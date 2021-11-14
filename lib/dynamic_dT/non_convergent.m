function bool = non_convergent(iDiffusion, iRES, RESvalues, Time, constants)
    % Detect whether the diffusion is no longer converging after a certain
    % number of diffusion iteration in the steady-state cycle
    
    after_threshold = mod(iDiffusion, constants.dynamicDT.nItersCycle) > constants.dynamicDT.iterThresholdDecrease;
    changed_this_ss = Time.current - Time.changed_dT < Time.dT_bac;
    if iRES > 5
        no_convergence = abs(max(RESvalues(:, iRES)) - max(RESvalues(:, iRES-1))) < constants.dynamicDT.tolerance_no_convergence;
    else
        no_convergence = false;
    end
    bool = after_threshold && ~changed_this_ss && no_convergence;
end

