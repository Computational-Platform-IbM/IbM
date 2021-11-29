function bac = bacteria_inactivate(bac, constants)
    % Inactivate any bacteria that are below the minimum mass threshold
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above

    mask_tooSmall = bac.molarMass * constants.bac_MW < constants.min_bac_mass_grams;
    mask_growthRate = bac.mu < 0;
    mask_combined = mask_tooSmall & mask_growthRate;
    
    bac.active(mask_combined) = 0;
    bac.active(~mask_combined) = 1;
    
    %{ 
    Reactivation is implicitely done already by inactive bacteria being able to grow.
    If at any point those bacteria become larger than the threshold, they
    will be set to active again.
    %}
end