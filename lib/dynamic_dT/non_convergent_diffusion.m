function bool = non_convergent_diffusion(iDiffusion, iRES, RESvalues, Time, constants)
    % Detect whether the diffusion is no longer converging after a certain
    % number of diffusion iteration in the steady-state cycle
    
%     after_threshold = mod(iDiffusion, constants.dynamicDT.nItersCycle) > constants.dynamicDT.iterThresholdDecrease; % +1 to also include the last iDiff
    changed_this_ss = Time.current - Time.changed_dT < Time.dT_bac;
    no_convergence = non_convergent(iRES, RESvalues, constants);
    
    bool = ~changed_this_ss && no_convergence;
end

