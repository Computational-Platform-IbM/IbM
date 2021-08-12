function integTime_fromScratch(grid, bac, conc, directory, constants, init_params)
    %% initialisation
    % calculate boundary conditions
    bulk_concs = calculate_bulk_concentrations(constants, init_params.init_conc, invHRT, reactionMatrix, dT);
    
    % make bacterial-grid matrix
    [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
    
    % determine diffusion layer
    diffusion_region = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid, constants); 
    
    % initialise concentrations
    conc = set_concentrations(conc, init_concs, diffusion_region);
    
    % set bulk layer concentrations
    conc = set_concentrations(conc, bulk_concs, ~diffusion_region);

    % calculate reaction matrix
    [reaction_matrix, bac.mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, grid, conc, constants, constants.pHsetpoint);
    
    % set time_indices
    T_indices = struct;
    T_indices.bac = 0;
    T_indices.divide = 0;
    T_indices.save = 0;
    
       
    %% time advancements (dT / dT_steadystate)
    
    % diffuse (MG)
    conc_old = conc;
    conc = diffusion(conc_old, reaction_matrix, bulk_concs, grid, constants, dT);
    
    % calculate reaction matrix
    [reaction_matrix, bac.mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, grid, conc, constants, pH);
    
    % if T>T_ss: calculate residual
    if T > T_ss && steadystate_is_reached(conc_0, conc)
        
        
        %% time advancements (dT_bac)
        % update bacteria: mass (+ update reaction matrix?)

            %% time advancements (dT_divide)
            % grow radius bacteria from mass
            % divide bacteria
            % shove bacteria
            % update diffusion region

        %% time advancements (dT_save)
        % save all important variables in R.mat
    
        % update boundary conditions
    
    end
    %% post-dT updates
    % 
    

    
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
    diffusion_rates, diffusion_accuracy, Tol_a              ==> diffusion

init_params:
    init_conc (== Sbc_dir)                                  ==> bulk_conc

bac:
    x, y                                                    ==> where bac
    x, y                                                    ==> diff region
    species, molarMass                                      ==> react_m

grid:
    dx, dy, nX, nY                                          ==> where bac
    nX, nY, blayer_thickness                                ==> diff region
    nX, nY                                                  ==> react_m

%}