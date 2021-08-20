function bac = bacteria_detachment(bac, constants)
    % Rough detachment: Remove all bacteria that are outside of maximum radius of granule
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    bac_norm = sqrt((bac.x - constant.Gcenter)*(bac.x - constant.Gcenter) + (bac.y - constant.Gcenter)*(bac.y - constants.Gcenter));
    bac_detach = bac_norm > constants.Grmax;
    nCellsDetach = sum(bac_detach);

    if nCellsDetach % <TODO: check if these are all struct-elements />
        bac.x(bac_detach) = [];
        bac.y(bac_detach) = [];
        bac.radius(bac_detach) = [];
        bac.species(bac_detach) = [];
        bac.molarMass(bac_detach) = [];
    end
end