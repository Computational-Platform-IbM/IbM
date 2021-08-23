function bac = bacteria_detachment(bac, grid, constants)
    % Rough detachment: Remove all bacteria that are outside of maximum radius of granule
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    bac_distance_from_center = sqrt((bac.x - grid.dx * grid.nX/2).*(bac.x - grid.dx * grid.nX/2) + (bac.y - grid.dy * grid.nY/2).*(bac.y - grid.dy * grid.nY/2));
    bac_detach = bac_distance_from_center > constants.max_granule_radius;
    nCellsDetach = sum(bac_detach);

    if nCellsDetach % <TODO: check if these are all struct-elements />
        bac.x(bac_detach) = [];
        bac.y(bac_detach) = [];
        bac.radius(bac_detach) = [];
        bac.species(bac_detach) = [];
        bac.molarMass(bac_detach) = [];
        bac.active(bac_detach) = [];
    end
end
