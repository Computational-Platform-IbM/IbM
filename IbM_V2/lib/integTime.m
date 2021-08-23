function profiling = integTime(grid, bac, conc, directory, constants, init_params)
    %% initialisation
    % calculate boundary conditions
    [bulk_concs, invHRT] = calculate_bulk_concentrations(constants, init_params.init_bulk_conc, init_params.invHRT, 0, constants.dT);
    
    % make bacterial-grid matrix
    [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
    
    % determine diffusion layer
    diffusion_region = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid); 
    
    % initialise concentrations
    conc = set_concentrations(conc, init_params.init_concs, diffusion_region);
    
    % set bulk layer concentrations
    conc = set_concentrations(conc, bulk_concs, ~diffusion_region);

    % calculate reaction matrix
    [reaction_matrix, bac.mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, grid, conc, constants, constants.pHsetpoint);
    
    % initiate times
    Time = struct;
    Time.current = 0;
    Time.steadystate = Time.current + constants.dT;
    Time.bac = constants.dT_bac;
    Time.divide = constants.dT_divide;
    Time.save = constants.dT_save;
   
    profiling = zeros(ceil(constants.simulation_end / constants.dT_bac)+1, 10);
    iProf = 1;
    %% time advancements (dT / dT_steadystate)
    
    while Time.current < constants.simulation_end
        % diffuse (MG)
        tic;
        conc = diffusion(conc, reaction_matrix, bulk_concs, grid, constants);
        profiling(iProf, 1) = profiling(iProf, 1) + toc;
        
        % updata bacterial mass
        tic;
        bac = update_bacterial_mass(bac, constants.dT);        
        profiling(iProf, 2) = profiling(iProf, 2) + toc;
        
        % calculate reaction matrix
        tic;
        [reaction_matrix, bac.mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, grid, conc, constants, pH);
        profiling(iProf, 3) = profiling(iProf, 3) + toc;


        % if T>T_ss: calculate residual
        if Time.current >= Time.steadystate
            
            tic;
            ssReached = steadystate_is_reached(conc, reaction_matrix, grid.dx, bulk_concs, constants);
            profiling(iProf, 4) = profiling(iProf, 4) + toc;
    
            if ssReached
                % set time to next bacterial activity time
                previousTime = Time.current;
                Time.current = Time.bac;
                if Time.current > constants.simulation_end
                    Time.current = constants.simulation_end - constants.dT;
                end
                
                % calculate actual dT for integration of bacterial mass
                dT_actual = Time.current - previousTime;

                %% time advancements (dT_bac)
                if Time.current >= Time.bac
                    % set next bacterial time
                    Time.bac = Time.bac + constants.dT_bac;

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
                            tic;
                            bac = bacteria_shove(bac, grid, constants);
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
                            
                            % update/re-determine where bacs
                            tic;
                            [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
                            profiling(iProf, 8) = profiling(iProf, 8) + toc;
                            
                            % update diffusion region
                            tic;
                            diffusion_region = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid);
                            profiling(iProf, 9) = profiling(iProf, 9) + toc;
                        end
                        
                    % calculate and set bulk concentrations
                    tic;
                    [bulk_concs, invHRT] = calculate_bulk_concentrations(constants, bulk_concs, invHRT, reaction_matrix, constants.dT);
                    conc = set_concentrations(conc, bulk_concs, ~diffusion_region);
                    profiling(iProf, 10) = profiling(iProf, 10) + toc;
                    
                    % recompute reaction matrix for next cycle
                    tic;
                    [reaction_matrix, bac.mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, grid, conc, constants, pH);
                    profiling(iProf, 3) = profiling(iProf, 3) + toc;
                    
                    iProf = iProf + 1;
                    fprintf('current simulation time: %.1f h', Time.current)
                end
                
                % set next steadystate time
                Time.steadystate = Time.current + 2*constants.dT;
                
                
                %% time advancements (dT_save)
                if Time.current >= Time.save
                    % set next save time
                    Time.save = Time.save + constants.dT_save;
                    
                    % save all important variables
                    save_slice(bac, conc(:, ceil(grid.nX/2)), bulk_concs, pH(:, ceil(grid.nX/2)), Time, grid, constants, directory); % along central horizontal axis
%                     save_plane(bac, conc, pH, Time, grid, constants, directory); % entire plane of simulation
                end
            end
        end
        
        %% post-dT updates
        % advance current simulation time
        Time.current = Time.current + constants.dT;
        

    end
end


% IMPORTANT FOR PLOTTING!!!
% watch out with indexing in matrices: matrix(i, j) is at grid cell with:
%   xi = i & yi = j. THIS DOES NOT CORRELATE WITH VISUAL OF MATRIX:
%   increasing the first index increases the y coordinate, not the x coord.

%{
constants:
    Keq, chrM, StNames, pH (setpoint), Vr, Vg, Dir_k, ...
        isLiquid, influent_concentrations, NH3sp, ...
        constantN                                           ==> bulk_conc
    Keq, chrM, T, Vg, StNames, react_v, Ki, Ks, ...
        MatrixMet, MatrixDecay                              ==> react_m
    diffusion_rates, diffusion_accuracy, Tol_a, dT          ==> diffusion
    dT, dT_bac, dT_division, dT_save                        ==> integ

init_params:
    init_bulk_conc (== Sbc_dir)                                  ==> bulk_conc
    invHRT                                                  ==> bulk_conc

bac:
    x, y                                                    ==> where bac
    x, y                                                    ==> diff region
    species, molarMass                                      ==> react_m
    radius

grid:
    dx, dy, nX, nY                                          ==> where bac
    nX, nY, blayer_thickness                                ==> diff region
    nX, nY                                                  ==> react_m

%}
