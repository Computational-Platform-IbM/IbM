function integTime(simulation_file, directory)

    %% load preset file
    load(simulation_file, 'grid', 'bac', 'constants', 'init_params', 'settings')
    constants.debug.plotConvergence = false;
    constants.debug.plotDiffRegion = false;
    
    %% Overall settings
    if settings.parallelized
        nChunks_dir = ceil(sqrt(feature('numcores')));

        if isempty(gcp('nocreate')) % check if parpool is already started (development ease-of-work)
            pc = parcluster('local');
            %             pc.JobStorageLocation = './.tmp/';
            parpool(pc, feature('numcores'));
        end

        fprintf('Parallelisation enabled for %d cores\n', feature('numcores'));
    end

    %% determine initial simulation parameters
    backup_file = sprintf('%s/backup.mat', directory);
    profiling_file = sprintf('%s/profilingResults.mat', directory);

    if isfile(backup_file) && isfile(profiling_file)
        % load from files
        load(backup_file, 'bac', 'bulk_concs', 'invHRT', 'conc', 'reaction_matrix', 'pH');
        load(profiling_file, 'profiling', 'maxErrors', 'normOverTime', 'nDiffIters', 'bulk_history', 'Time');
    else
        % initiate from preset values
        [conc, bulk_concs, invHRT, reaction_matrix, pH, bac] = ...
            initTime(grid, bac, init_params, constants, settings);

        % initiate time and profiling information/storage from preset
        Time = struct;
        Time.current = 0;
        Time.steadystate = Time.current + (constants.nDiffusion_per_SScheck - 1) * constants.dT;
        Time.save = constants.dT_save;
        Time.backup = constants.dT_backup;
        Time.changed_dT = 0;
        Time.changed_dT_bac = 0;
        Time.dT = constants.dT;
        Time.dT_bac = constants.dT_bac;
        Time.bac = Time.dT_bac; % include dT_bac and dT_divide in one variable

        if settings.dynamicDT
            Time.maxDT = constants.dynamicDT.maxDT;
            Time.minDT = constants.dynamicDT.minDT;
            Time.maxDT_bac = constants.dynamicDT.maxDT_bac;
            Time.minDT_bac = constants.dynamicDT.minDT_bac;

            Time.history = zeros(ceil(constants.simulation_end / Time.minDT_bac), 1, 'single'); % save time at each Steady-state cycle
            profiling = zeros(ceil(constants.simulation_end / Time.minDT_bac), 11, 'single');
            maxErrors = zeros(ceil(constants.simulation_end / Time.minDT_bac), 1, 'single'); % store max error per dT_bac
            normOverTime = zeros(ceil(constants.simulation_end / Time.minDT_bac), 1, 'single'); % store norm of concentration differance per dT_bac
            nDiffIters = zeros(ceil(constants.simulation_end / Time.minDT_bac), 1, 'uint16'); % store number of diffusion iterations per steady state
            bulk_history = zeros(size(bulk_concs, 1), ceil(constants.simulation_end / Time.minDT_bac), 'single');
            maxInitRES = zeros(ceil(constants.simulation_end / Time.minDT_bac), 1, 'single');
        else
            Time.history = zeros(ceil(constants.simulation_end / constants.dT_bac), 1, 'single'); % save time at each Steady-state cycle
            Time.minDT = Time.dT;
            profiling = zeros(ceil(constants.simulation_end / constants.dT_bac), 11, 'single');
            maxErrors = zeros(ceil(constants.simulation_end / constants.dT_bac), 1, 'single'); % store max error per dT_bac
            normOverTime = zeros(ceil(constants.simulation_end / constants.dT_bac), 1, 'single'); % store norm of concentration differance per dT_bac
            nDiffIters = zeros(ceil(constants.simulation_end / constants.dT_bac), 1, 'uint16'); % store number of diffusion iterations per steady state
            bulk_history = zeros(size(bulk_concs, 1), ceil(constants.simulation_end / constants.dT_bac), 'single');
            maxInitRES = zeros(ceil(constants.simulation_end / constants.dT_bac), 1, 'single');
        end
        bulk_history(:, 1) = bulk_concs; % is added after += 1 of iProf, thus give first value already

        % initialise saving file
        save_slice(bac, conc, bulk_concs, pH, invHRT, 0, grid, constants, directory);
    end


    RESvalues = zeros(length(constants.compoundNames), 200); % reserve for n steady state checks beforehand (can be more)
    norm_diff = zeros(200, 1);
    res_bacsim = zeros(200, 2);

    iProf = find(profiling == 0, 1, 'first'); % keep track of index of profiling (every simulated dT_bac == +1 index)
    iDiffusion = 1; % keep track of index of diffusion (per 1 dT_bac: iDiffusion == cycles of diffusion)
    iRES = 0;

    % make bacterial-grid matrix
    [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);

    % determine diffusion layer and calculate ranges for focus mask
    [diffusion_region, focus_region] = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid);
    xRange = focus_region.x0:focus_region.x1;
    yRange = focus_region.y0:focus_region.y1;

    if settings.parallelized
        % create chunks
        chunks = create_chunks(nChunks_dir, focus_region);

        % sort bacteria
        bac = sort_bacteria_into_chunks(bac, grid, chunks, focus_region, nChunks_dir);

        % recalculate the grid2bac matrix
        [grid2bac, ~] = determine_where_bacteria_in_grid(grid, bac);
    else
        % set dummy values for chunk variables
        chunks = 0;
        nChunks_dir = 0;
    end

    %% time advancements (dT / dT_steadystate)
    prev_conc = conc;

    while Time.current < constants.simulation_end
        
        if mod(iDiffusion, constants.dynamicDT.nItersCycle) == 0
            fprintf('Currently at diffusion iteration %d (max error: %g)\n', iDiffusion, max(RESvalues(:, iRES)))
        end
        
        % diffuse (MG)
        tic;
        try
            conc(xRange, yRange, :) = diffusionMG(conc(xRange, yRange, :), reaction_matrix(xRange, yRange, :), ...
                bulk_concs, diffusion_region(xRange, yRange), grid, constants, Time);
        catch e
            switch e.identifier
                case 'Diffusion:NegativeConcentration'
                    Time = decrease_dT_diffusion(Time, e.message, grid.dx, constants);
%                     NegConc = 1;
                otherwise
                    rethrow(e)
            end
        end
        profiling(iProf, 1) = profiling(iProf, 1) + toc;

        % set bulk layer concentrations (in theory, not needed anymore with
        % correct diffusion model
        conc = set_concentrations(conc, bulk_concs, ~diffusion_region);

        % updata bacterial mass
        tic;
        bac = update_bacterial_mass(bac, Time.dT);
        profiling(iProf, 2) = profiling(iProf, 2) + toc;

        % calculate reaction matrix
        tic;
        [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
            grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange), chunks, nChunks_dir, settings);
        profiling(iProf, 3) = profiling(iProf, 3) + toc;

%         if mod(iDiffusion, constants.dynamicDT.nItersCycle) == 0 && ...
%                 non_convergent(iRES, RESvalues, constants)
%             
%             [bulk_concs, invHRT] = calculate_bulk_concentrations(bac, constants, bulk_concs, invHRT, reaction_matrix, Time.dT_bac - Time.current, settings);
%             conc = set_concentrations(conc, bulk_concs, ~diffusion_region);
% 
%             % ----- DEBUG -----
%             if constants.debug.plotConvergence
%                 plotConvergence(RESvalues, iRES, constants, Time.current)
%                 plotNormDiff(norm_diff, iRES, constants, Time.current)
%                 plotBacSimError(res_bacsim, iRES, constants, Time.current)
%                 drawnow()
%             end
%             fprintf('Bulk concentrations recalculated,\nmax RES value at the moment: %.2f %%\n\n', max(RESvalues(:,iRES))*100)
%             % ----- END DEBUG -----
%             
%             % recompute reaction matrix
%             tic;
%             [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
%                 grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange), chunks, nChunks_dir, settings);
%             profiling(iProf, 3) = profiling(iProf, 3) + toc;
%         end


%         if settings.dynamicDT && mod(iDiffusion, constants.dynamicDT.nItersCycle) == 0
%             if mod(iDiffusion, 2*constants.dynamicDT.nItersCycle) == 0
%                 Time = increase_dT_diffusion(Time, 'Diffusion takes a long time', grid.dx, constants);
%             else
%                 % if increasing does not work, try decreasing...
%                 Time = decrease_dT_diffusion(Time, 'Diffusion takes a long time', grid.dx, constants);
%             end
%         end
        
        
        % if T>T_ss: calculate residual
        if Time.current >= Time.steadystate

            % increase counter of RES checks performed
            iRES = iRES + 1;

            % perform check for steady state
            tic;
            [ssReached, RESvalues(:, iRES)] = steadystate_is_reached(conc(xRange, yRange, :), reaction_matrix(xRange, yRange, :), grid.dx, bulk_concs, diffusion_region(xRange, yRange), constants);
            norm_diff(iRES) = sqrt(sum((prev_conc - conc).^2, 'all'));
            res_bacsim(iRES, 1) = max(abs((prev_conc - conc) ./ Time.dT), [], 'all');
            res_bacsim(iRES, 2) = norm_diff(iRES) / Time.dT;
            profiling(iProf, 4) = profiling(iProf, 4) + toc;

            prev_conc = conc;

            % perform dynamic dT for diffusion
            if settings.dynamicDT
                if slow_convergence(iRES, RESvalues, constants) && Time.dT < Time.maxDT/2 %&& NegConc == 0 %-- Testing --%
                    Time = increase_dT_diffusion(Time, 'Slow convergence', grid.dx, constants);
%                     Time = decrease_dT_diffusion(Time, 'Diffusion takes a long time', grid.dx, constants);
                end
                if upward_trend(iRES, RESvalues)
                    Time = decrease_dT_diffusion(Time, 'Upward trend in RES values detected', grid.dx, constants);
                    Time = decrease_dT_diffusion(Time, 'Upward trend in RES values detected', grid.dx, constants);
                elseif non_convergent(iRES, RESvalues, constants.dynamicDT.tolerance_no_convergence)
                    Time = decrease_dT_diffusion(Time, 'Convergence is stuck', grid.dx, constants);
                end
                if iDiffusion > 5000
                    ssReached = true;
                end
            else
                if iDiffusion > 5000 && non_convergent_diffusion(iDiffusion, iRES, RESvalues, Time, constants)
                    % without dynamic timestep & negative concentrations
                    % due to too large step size, accept SS under
                    % non-convergent conditions
                    ssReached = true;
                end
                if iDiffusion > 10000
                    ssReached = true;
                end
            end

            if ssReached
%                 NegConc = 0;
                fprintf('Steady state reached after %d diffusion iterations\n', iDiffusion)
                fprintf('\twith at most %.4g mol/L/h off of steady state (norm = %e)\n', max(RESvalues(:, iRES)), norm_diff(iRES))

                % set time to next bacterial activity time
                previousTime = Time.current;
                Time.current = Time.bac;

                if Time.current > constants.simulation_end
                    Time.current = constants.simulation_end - Time.dT;
                end

                % perform dynamic dT for diffusion
                if settings.dynamicDT && multiple_high_iters(iDiffusion, iProf, nDiffIters, Time, constants)
                    Time = increase_dT_diffusion(Time, sprintf('Multiple steady states reached with more than %d diffusion iterations', constants.dynamicDT.iterThresholdIncrease), grid.dx, constants);
                end

                % calculate actual dT for integration of bacterial mass
                dT_actual = Time.current - previousTime;

                %% time advancements (dT_bac)
                if Time.current >= Time.bac

                    % ------ DEBUG ------
                    % plot convergence
                    if constants.debug.plotConvergence
                        plotConvergence(RESvalues, iRES, constants, Time.current)
                        plotNormDiff(norm_diff, iRES, constants, Time.current)
                        plotBacSimError(res_bacsim, iRES, constants, Time.current)
                        drawnow()

                        %                       plotBulkConcOverTime(bulk_history(:, iProf-10:iProf), Time.history(iProf-10:iProf), constants);
                        if iProf > 20
                            figure(19);
                            plot(Time.history(1:iProf), bulk_history(2, 1:iProf), 'LineWidth', 2);
                            legend('NO2 bulk_conc', 'Location', 'northwest')

                            figure(20);
                            plot(Time.history(iProf - 10:iProf - 1), diff(bulk_history(:, iProf - 10:iProf), 1, 2) ./ (bulk_history(:, iProf - 10:iProf - 1) + 1e-20), 'LineWidth', 2);
                            legend(constants.compoundNames, 'Location', 'northwest');
                        end

                        drawnow();
                    end
                    % ------- END DEBUG -------

                    maxErrors(iProf) = max(RESvalues(:, iRES));
                    normOverTime(iProf) = norm_diff(iRES);
                    nDiffIters(iProf) = iDiffusion;
                    maxInitRES(iProf) = max(RESvalues(:, 1));
                    iDiffusion = 1;
                    iRES = 0;
                    RESvalues = zeros(length(constants.compoundNames), 100); % reset RES value array

                    % reaction_matrix & mu & pH are already calculated (steady state, so still valid)

                    % update bacteria: mass
                    tic;
                    bac = update_bacterial_mass(bac, dT_actual);
                    profiling(iProf, 2) = profiling(iProf, 2) + toc;

                    %% time advancements (dT_divide)
                    % determine radius bacteria from mass
                    tic;
                    bac = update_bacterial_radius(bac, constants);
                    profiling(iProf, 2) = profiling(iProf, 2) + toc;

                    % bacteria: inactivate or die
                    tic;
                    if constants.inactivationEnabled
                        bac = bacteria_inactivate(bac, constants);
                    else
                        bac = bacteria_die(bac, constants);
                    end

                    % bacteria: divide
                    [bac, nDivCycles] = bacteria_divide(bac, constants);
                    if nDivCycles > 1
                        Time = decrease_dT_bac(Time, 'Bacteria are dividing too fast');
                    end
                    profiling(iProf, 5) = profiling(iProf, 5) + toc;

                    % shove bacteria
                    tic;
                    bac = bacteria_shove(bac, grid, constants);
                    bac = bacteria_shove(bac, grid, constants); % perform second shove call in case of overcrowding
                    profiling(iProf, 6) = profiling(iProf, 6) + toc;

                    % bacteria: detachment {for now only rough detachment is implemented}
                    tic;
                    bac = bacteria_detachment(bac, grid, constants, settings, Time.dT_bac);
                    profiling(iProf, 7) = profiling(iProf, 7) + toc;

                    % display number of bacteria in system
                    fprintf('current number of bacteria: %d (%d active)\n', length(bac.x), sum(bac.active))

                    % auto detect when to switch to parallel
                    % computation of rMatrix
                    % TODO: determine exact cutoff value (will be
                    % around 15000 approx)
                    if settings.parallelized == false && length(bac.x) > 15000
                        settings.parallelized = true;
                        nChunks_dir = ceil(sqrt(feature('numcores')));

                        if isempty(gcp('nocreate')) % check if parpool is already started (development ease-of-work)
                            pc = parcluster('local');
                            %                             pc.JobStorageLocation = './.tmp/';
                            parpool(pc, feature('numcores'));
                        end

                        fprintf('Parallelisation enabled for %d cores\n', feature('numcores'));
                    end

                    % update/re-determine where bacs
                    tic;
                    [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
                    profiling(iProf, 8) = profiling(iProf, 8) + toc;

                    % update diffusion region
                    tic;
                    [diffusion_region, focus_region] = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid);
                    xRange = focus_region.x0:focus_region.x1;
                    yRange = focus_region.y0:focus_region.y1;
                    profiling(iProf, 9) = profiling(iProf, 9) + toc;

                    if settings.parallelized
                        tic;
                        % create chunks
                        chunks = create_chunks(nChunks_dir, focus_region);

                        % sort bacteria
                        bac = sort_bacteria_into_chunks(bac, grid, chunks, focus_region, nChunks_dir);
                        profiling(iProf, 11) = profiling(iProf, 11) + toc;

                        % recalculate the grid2bac matrix
                        [grid2bac, ~] = determine_where_bacteria_in_grid(grid, bac);
                    end

                    %% apply dynamic dT_bac
                    if settings.dynamicDT && multiple_low_initRES(iProf, maxInitRES, Time, constants)
                        Time = increase_dT_bac(Time, sprintf('Multiple steady state cycles with initRES value of <%.1f %%', constants.dynamicDT.initRESThresholdIncrease * 100));
                    end

                    %% prepare for next steady-state cycle
                    fprintf('\n ================ \n')
                    fprintf('Current simulation time: %.2f h\n\n', Time.current)

                    % recompute reaction matrix for next cycle
                    tic;
                    [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
                        grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange), chunks, nChunks_dir, settings);
                    profiling(iProf, 3) = profiling(iProf, 3) + toc;
                    
                    % calculate and set bulk concentrations
                    tic;
                    while true % should be a "do while" loop, but Matlab doesn't have that functionality ...
                        [new_bulk_concs, invHRT] = calculate_bulk_concentrations(bac, constants, bulk_concs, invHRT, reaction_matrix, Time.dT_bac, settings);
                        if ~settings.dynamicDT || bulk_conc_diff_within_limit(new_bulk_concs, bulk_concs, constants)
                            break
                        end

                        if Time.dT_bac <= Time.minDT_bac
                            fprintf(2, 'Smallest dT_bac reached, cannot decrease smaller than %g h\n', Time.minDT_bac);
                            break
                        end

                        Time = decrease_dT_bac(Time, 'Too large bulk concentration jump detected');
                    end

                    bulk_concs = new_bulk_concs;
                    conc = set_concentrations(conc, bulk_concs, ~diffusion_region);
                    profiling(iProf, 10) = profiling(iProf, 10) + toc;

                    iProf = iProf + 1;
                    bulk_history(:, iProf) = bulk_concs;
                    Time.history(iProf) = Time.current;

                    % set next bacterial time
                    Time.bac = Time.bac + Time.dT_bac;
                end

                %% time advancements (dT_save)
                if Time.current >= Time.save
                    % set next save time
                    Time.save = Time.save + constants.dT_save;

                    % save all important variables
                    save_slice(bac, conc, bulk_concs, pH, invHRT, Time.current, grid, constants, directory);
%                     save_profile(bac, conc, bulk_concs, pH, invHRT, Time.current, grid, constants, directory); % entire plane of simulation

                    if Time.current >= Time.backup
                        % set next backup time
                        Time.backup = Time.backup + constants.dT_backup;

                        % save all important variables for continuing
                        % simulation from this point
                        save_backup(bac, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
                        save_profiling(profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time, directory)
                    end
                    
                    if strcmp(settings.detachment, 'SBR')
                        bac.molarMass(bac.active) = max(bac.molarMass(bac.active) - 0.1*constants.min_bac_mass_grams/constants.bac_MW, eps);
                    end
                end
            end

            % set next steadystate time
            Time.steadystate = Time.current + constants.nDiffusion_per_SScheck * Time.dT;
        end

        %% post-dT updates
        % advance current simulation time
        Time.current = Time.current + Time.dT;
        iDiffusion = iDiffusion + 1;
    end

    % save all important variables one last time?
    save_slice(bac, conc, bulk_concs, pH, invHRT, Time.current, grid, constants, directory);
    % save_profile(bac, conc, bulk_concs, pH, invHRT, Time.current, grid, constants, directory); % entire plane of simulation
    save_backup(bac, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
    save_profiling(profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time, directory)
end

%{
Figure reservations:
1) profiling
2) Bacteria
3) Diffusion region
4) Convergence to steady state
5) Max error over time
6) Bulk concentrations over time
7) concentration profiles (2D)
8) Norm(conc1 - conc0) over time
9) Norm(conc1 - conc0) to steady state
10) Reaction profiles (2D)
11) RES value (or residual of diffusion) profile (2D)
12) BacSim error (approximation of dC/dt)
13) Detach Time
%}

% IMPORTANT FOR PLOTTING!!!
% watch out with indexing in matrices: matrix(i, j) is at grid cell with:
%   xi = i & yi = j. THIS DOES NOT CORRELATE WITH VISUAL OF MATRIX:
%   increasing the first index increases the y coordinate, not the x coord.
%   Instead, visualise the transpose of the matrix!
