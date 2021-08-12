function conc = diffusion(conc_old, reaction_matrix, diffusion_region, bulk_concentrations, grid, constants, dT)
    % Solve diffusion for all molecules using the multigrid method.
    %
    % conc_old: concentration of each molecule in the grid 
    %   (ix, iy, compound)
    % reaction_matrix: matrix with per grid cell and per compound the
    %   change [h-1] due to bacterial activity
    % diffusion_region: logical matrix with per grid cell whether it is
    %   in the diffusion region (true) or not (false).
    % bulk_concentrations: vector with the bulk concentration per compound
    % grid: struct containing all information regarding the grid
    % constants: struct containing all simulation constants
    % dT: time over which to solve the diffusion equations
    %
    % -> conc: concentrations after solving the diffusion equations

    % set concentrations of bulk layer
    
    % create right hand side
    % - boundary conditions
    % - reaction matrix
    
    % solve using multigrid


end
