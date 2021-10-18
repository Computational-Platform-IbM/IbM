function integTime(simulation_file, directory)

% load simulation file
% check if there is a results file already
%   if not: initTime()
%           initProfiling()
%   else: load from backup file
% perform time looping


    %% load preset file
    load(simulation_file, 'grid', 'bac', 'constants', 'init_params', 'settings')
    
    %% Overall settings
    if settings.parallelized
        nChunks_dir = ceil(sqrt(feature('numcores')));
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
            initTime(grid, bac, init_params, constants);

        % initiate time and profiling information/storage from preset
        Time = struct;
        Time.current = 0;
        Time.steadystate = Time.current + (constants.nDiffusion_per_SScheck - 1)*constants.dT;
        Time.bac = constants.dT_bac;
        Time.divide = constants.dT_divide;
        Time.save = constants.dT_save;
        Time.backup = constants.dT_backup;

        profiling = zeros(ceil(constants.simulation_end / constants.dT_bac)+1, 11);
        maxErrors = zeros(ceil(constants.simulation_end / constants.dT_bac), 1); % store max error per dT_bac
        normOverTime = zeros(ceil(constants.simulation_end / constants.dT_bac), 1); % store norm of concentration differance per dT_bac
        nDiffIters = zeros(ceil(constants.simulation_end / constants.dT_bac), 1); % store number of diffusion iterations per steady state
        bulk_history = zeros(size(bulk_concs, 1), ceil(constants.simulation_end / constants.dT_bac)+1);        
        bulk_history(:,1) = bulk_concs; % is added after += 1 of iProf, thus give first value already
       
        % initialise saving file
        save_slice(bac, conc, bulk_concs, pH, invHRT, 0, grid, constants, directory);
    end


    RESvalues = zeros(sum(constants.isLiquid), 200); % reserve for n steady state checks beforehand (can be more)
    norm_diff = zeros(200,1);

    iProf = find(profiling == 0, 1, 'first');        % keep track of index of profiling (every simulated hour == +1 index)
    iDiffusion = 1;   % keep track of index of diffusion (per 1 dT_bac: iDiffusion == cycles of diffusion)
    iRES = 1;

    
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
    end
    
    
    %% time advancements (dT / dT_steadystate)
    prev_conc = conc;
    dT = constants.dT;
    
    while Time.current < constants.simulation_end
        % diffuse (MG)
        tic;
        conc(xRange, yRange, :) = diffusion(conc(xRange, yRange, :), reaction_matrix(xRange, yRange, :), ...
            bulk_concs, diffusion_region(xRange, yRange), grid, constants, dT);
        profiling(iProf, 1) = profiling(iProf, 1) + toc;
        
        % set bulk layer concentrations (in theory, not needed anymore with
        % correct diffusion model
        conc = set_concentrations(conc, bulk_concs, ~diffusion_region);
    
        % updata bacterial mass
        tic;
        bac = update_bacterial_mass(bac, dT);        
        profiling(iProf, 2) = profiling(iProf, 2) + toc;
        
        % calculate reaction matrix
        tic;
        if settings.parallelized
            [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = par_calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
                grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange), chunks, nChunks_dir);
        else
            [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
                grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange));
        end
        profiling(iProf, 3) = profiling(iProf, 3) + toc;


        % if T>T_ss: calculate residual
        if Time.current >= Time.steadystate
            
            % increase counter of RES checks performed
            iRES = iRES + 1; 
            
            % perform check for steady state
            tic;
            [ssReached, RESvalues(:,iRES)] = steadystate_is_reached(conc(xRange, yRange, :), reaction_matrix(xRange, yRange, :), grid.dx, bulk_concs, diffusion_region(xRange, yRange), constants);
            norm_diff(iRES) = sqrt(sum((prev_conc - conc).^2, 'all'));
            profiling(iProf, 4) = profiling(iProf, 4) + toc;
            
            prev_conc = conc;
            
            % increase dT if convergence is slow
            if settings.dynamicDT && Time.current > 5 && mod(iRES, ceil(100 / constants.nDiffusion_per_SScheck)) == 0 && norm_diff(iRES)/norm_diff(iRES-1) > 0.99
                dT = dT*1.2;
                fprintf('dT increased to %g at t=%g after %d diffusion iterations\n', dT, Time.current, iDiffusion)
            end
            
            
            % check if system is converging towards steady state
%             noLongerConverging = abs(max(RESvalues(:, iRES)) - max(RESvalues(:, iRES - 1))) < constants.convergence_accuracy;
            noLongerConverging = false;
            
            if ssReached || noLongerConverging || iDiffusion > 1500
                fprintf('Steady state reached after %d diffusion iterations\n', iDiffusion)
                if noLongerConverging
                    fprintf('\nNo longer converging (delta-RES = %g), steady state accepted with at most %.4f %% off of steady state\n\n', abs(max(RESvalues(:, iRES)) - max(RESvalues(:, iRES - 1))), max(RESvalues(:, iRES))*100)
                elseif iDiffusion > 1500
                    fprintf('\nNo longer converging (>1500 diffusion iterations), steady state accepted with at most %.4f %% off of steady state (norm = %e)\n\n', max(RESvalues(:, iRES))*100, norm_diff(iRES))
                else
                    fprintf('\twith at most %.4f %% off of steady state (norm = %e)\n', max(RESvalues(:, iRES))*100, norm_diff(iRES))
                end
                                
                % set time to next bacterial activity time
                previousTime = Time.current;
                Time.current = Time.bac;
                if Time.current > constants.simulation_end
                    Time.current = constants.simulation_end - constants.dT;
                end
                
                % (re)set dT to constant value each dT_bac
                dT = constants.dT;
                
                % calculate actual dT for integration of bacterial mass
                dT_actual = Time.current - previousTime;

                %% time advancements (dT_bac)
                if Time.current >= Time.bac
                    % set next bacterial time
                    Time.bac = Time.bac + constants.dT_bac;
                    
                    if constants.debug.plotConvergence
                        plotConvergence(RESvalues, iRES, constants, Time.current)
                        plotNormDiff(norm_diff, iRES, constants, Time.current)
                    end
                    
                    maxErrors(iProf) = max(RESvalues(:,iRES));
                    normOverTime(iProf) = norm_diff(iRES);
                    nDiffIters(iProf) = iDiffusion;
                    iDiffusion = 1;
                    iRES = 1;
                    RESvalues = zeros(sum(constants.isLiquid), 100); % reset RES value array
                    
                    % reaction_matrix & mu & pH are already calculated (steady state, so still valid)

                    % update bacteria: mass
                    tic;
                    bac = update_bacterial_mass(bac, dT_actual);
                    profiling(iProf, 2) = profiling(iProf, 2) + toc;


                        %% time advancements (dT_divide)
                        if Time.current >= Time.divide
                            % set next division time
                            Time.divide = Time.divide + constants.dT_divide;
                            
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
                            

                            %{
                              <E: In the PREV model we are able to include some evolution-adaptation of bacteria (for example Ks, Ki, Yield etc) />
                              <E: For this reason, in R struct you can see a matrix of Ks and Ki for every single bacterium, because everyone is an entity with its own kinetic properties. />
                              <E: In this case, I'm thinking some options here: (1) To include the Ks_Ki matrix, 
                                                                                (2) To include a matrix with Ks and Ki multipliers for every bacterium, 
                                                                                (3) Two IbM versions: one w/ evolution-adaptation and another w/o. />
                            %}
                            
                            %{
                            <C: I'm leaning towards the 3rd option, for now
                            we don't seem to be using evolution. So doesn't
                            seem fit to be using redundant memory space...
                            would be easy to introduce when we need it
                            %}
                            
                            % bacteria: divide
                            bac = bacteria_divide(bac, constants);
                            profiling(iProf, 5) = profiling(iProf, 5) + toc;
                            
                            % shove bacteria
                            % <E: we could only call bacteria_shove if
                            % division is done. />
                            tic;
                            bac = bacteria_shove(bac, grid, constants);
                            bac = bacteria_shove(bac, grid, constants); % add second shove to make ensure no overlap
                            profiling(iProf, 6) = profiling(iProf, 6) + toc;
                            
                            %{
                              <E: Here we could add a constant to define if we want include detachment and how we simulate this detachment. />
                              <E: For now, only the rough detachment in included. />
                              <C: agree, for now this works
                            %}
                            % bacteria: detachment
                            tic;
                            bac = bacteria_detachment(bac, grid, constants);
                            profiling(iProf, 7) = profiling(iProf, 7) + toc;
                            
                            % display number of bacteria in system
                            fprintf('current number of bacteria: %d \n', length(bac.x))
                            
                            if constants.debug.plotBacteria
                                plotBacs(grid, bac, constants)
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
                            
                            if constants.debug.plotDiffRegion
                                plotDiffRegion(grid, bac, diffusion_region, true)
                            end
                            
                        end
                        
                    % calculate and set bulk concentrations
                    tic;
                    [bulk_concs, invHRT] = calculate_bulk_concentrations(constants, bulk_concs, invHRT, reaction_matrix, constants.dT_bac);
                    conc = set_concentrations(conc, bulk_concs, ~diffusion_region);
                    profiling(iProf, 10) = profiling(iProf, 10) + toc;
                    
                    % recompute reaction matrix for next cycle
                    tic;
                    if settings.parallelized
                        [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = par_calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
                            grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange), chunks, nChunks_dir);
                    else
                        [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
                            grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange));
                    end
                    profiling(iProf, 3) = profiling(iProf, 3) + toc;
                    
                    iProf = iProf + 1;
                    bulk_history(:,iProf) = bulk_concs;
                    fprintf('\n\n ================ \n')
                    fprintf('Current simulation time: %.1f h\n', Time.current)
                end
                
                %% time advancements (dT_save)
                if Time.current >= Time.save
                    % set next save time
                    Time.save = Time.save + constants.dT_save;
                    
                    % save all important variables
                    save_slice(bac, conc, bulk_concs, pH, invHRT, Time.current, grid, constants, directory);
%                     save_plane(bac, conc, pH, Time, grid, constants, directory); % entire plane of simulation

                    if Time.current >= Time.backup
                        % set next backup time
                        Time.backup = Time.backup + constants.dT_backup;
                        
                        % save all important variables for continuing
                        % simulation from this point
                        save_backup(bac, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
                        save_profiling(profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time, directory)
                    end
                end
                
            end
            
            % set next steadystate time
            Time.steadystate = Time.current + constants.nDiffusion_per_SScheck*dT;
        end
        
        %% post-dT updates
        % advance current simulation time
        Time.current = Time.current + dT;
        iDiffusion = iDiffusion + 1;
    end
    
    if constants.debug.plotMaxErrors
        plotMaxErrorOverTime(maxErrors, constants.dT_bac)
        plotNorm(normOverTime, constants.dT_bac)
    end
    
    if constants.debug.plotBulkConcsOverTime
        plotBulkConcOverTime(bulk_history, constants)
    end
    
    % save all important variables one last time?
    save_slice(bac, conc, bulk_concs, pH, invHRT, Time.current, grid, constants, directory);
    % save_profile(bac, conc, bulk_concs, pH, invHRT, Time.current, grid, constants, directory); % entire plane of simulation
    save_backup(bac, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
    save_profiling(profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time, directory)
    
    
%     plotConcs(conc, constants, Time.current);
%     plotBacs(grid, bac, constants);
%     plotDiffRegion(grid, bac, diffusion_region, false);
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
%}


% IMPORTANT FOR PLOTTING!!!
% watch out with indexing in matrices: matrix(i, j) is at grid cell with:
%   xi = i & yi = j. THIS DOES NOT CORRELATE WITH VISUAL OF MATRIX:
%   increasing the first index increases the y coordinate, not the x coord.
%   Instead, visualise the transpose of the matrix!



