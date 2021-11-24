function bac = bacteria_die(bac, constants)
    % Remove all bacteria that are smaller than the mass threshold
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    mask_tooSmall = bac.molarMass * constants.bac_MW < constants.min_bac_mass_grams;
    nCellsTooSmall = sum(mask_tooSmall);

    if nCellsTooSmall
        bac = killBacs(bac, mask_tooSmall);
    end
end