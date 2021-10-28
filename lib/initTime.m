function [conc, bulk_concs, invHRT, reaction_matrix, pH, bac] = initTime(grid, bac, init_params, constants, settings)
    %% initialisation
    % initiate values from preset parameters
    
    % calculate boundary conditions
    [bulk_concs, invHRT] = calculate_bulk_concentrations(constants, init_params.init_bulk_conc, init_params.invHRT, 0, constants.dT_bac, settings);
    
    % <C: why calculate bulk concentration with dT_bac? />
    
    % make bacterial-grid matrix
    [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
    
    % determine diffusion layer and calculate ranges for focus mask
    [diffusion_region, focus_region] = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid);
    xRange = focus_region.x0:focus_region.x1;
    yRange = focus_region.y0:focus_region.y1;

    if constants.debug.plotDiffRegion
        plotDiffRegion(grid, bac, diffusion_region, true)
    end
    
    % initialise concentrations and pH
    conc = zeros(grid.nX, grid.nY, size(constants.isLiquid, 1));
    conc = set_concentrations(conc, init_params.init_concs, diffusion_region);
    reaction_matrix = zeros(grid.nX, grid.nY, size(conc, 3));
    pH = ones(grid.nX, grid.nY) * constants.pHsetpoint;
    
    % set bulk layer concentrations
    conc = set_concentrations(conc, bulk_concs, ~diffusion_region);

    % calculate reaction matrix
    [reaction_matrix(xRange, yRange, :), bac.mu, pH(xRange, yRange)] = calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
        grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, pH(xRange, yRange), settings);
    
end