function bac = bacteria_inactivate(bac, constants)
    % Inactivate any bacteria that are below the minimum mass threshold
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above

    mask_tooSmall = bac.molarMass * constants.bac_MW < constants.min_bac_mass_grams;
    mask_positiveGrowthRate = bac.mu > 0;
    mask_possible_reactivation = ~bac.active & mask_positiveGrowthRate;
%     mask_tiny = bac.molarMass * constants.bac_MW < 0.4*constants.min_bac_mass_grams;
%     mask_combined = (mask_tooSmall & mask_growthRate) | mask_tiny;
    
    % too small & ~active -> ~active (trivial case, i.e. same value)
    bac.active(mask_tooSmall & bac.active) = 0; % too small & active -> inactivate
    random_reactivation = rand(sum(mask_possible_reactivation),1) < 0.1;
    bac.active(mask_possible_reactivation) = random_reactivation; % ~active & positive growth rate -> 50% activate
    % large enough & active -> active (trivial case, i.e. same value)
    
    %{ 
    Reactivation is implicitely done already by inactive bacteria being able to grow.
    If at any point those bacteria become larger than the threshold, they
    will be set to active again.
    %}
    fprintf('%d individuals reactivated\n', sum(random_reactivation));
end