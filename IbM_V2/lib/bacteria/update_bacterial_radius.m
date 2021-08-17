function bac = update_bacterial_radius(bac, constants)
    % Update the radius of each bacterium
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    bac.radius = ((bac.molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))).^(1/3);
end