function bac = granule_density(bac, constants)
    % Update the density of granule
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    bac_m = bac.molarMass * constants.bac_MW;
    bac.bac_rho_bio = sum(bac_m)/((max(bac.y)- min(bac.y))*(max(bac.x)- min(bac.x))*(1e-6));
end