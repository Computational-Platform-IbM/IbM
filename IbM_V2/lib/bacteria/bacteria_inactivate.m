function bac = bacteria_inactivate(bac, constants)
    % Inactivate any bacteria that are below the minimum mass threshold
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above

    mask_tooSmall = bac.molarMass * constants.bac_MW < constants.min_bac_mass_grams;

    bac.active(mask_tooSmall) = 0;
    bac.active(~mask_tooSmall) = 1;
    
    %{ 
    <E: Bacteria can be active again if local conditions are favourable. />
    maskd_reactivate = (~bac.active) * bac.molarMass * constants.bac_MW > constants.min_bac_mass_grams;
    bac.active(maskd_reactivate) = 1;
    %}
end