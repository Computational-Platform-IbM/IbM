function bool = slow_convergence(iRES, RESvalues, constants)
    % Detect whether the diffusion is no longer converging after a certain
    % number of diffusion iteration in the steady-state cycle
    
    at_cycle_time = mod(iRES, ceil(constants.dynamicDT.nItersCycle/constants.nDiffusion_per_SScheck)) == 0; % only at nCycle iterations
    if at_cycle_time
        
        no_conv = non_convergent(iRES, RESvalues, constants.dynamicDT.tolerance_no_convergence); % better convergence than no-convergence
        direction_of_conv = max(RESvalues(:,iRES-1)) - max(RESvalues(:,iRES)) > 0 && max(RESvalues(:,iRES-2)) - max(RESvalues(:,iRES)) > 0; % has to be some convergence
        little_conv =  max(RESvalues(:,iRES-1)) - max(RESvalues(:,iRES)) < 1e-3 && max(RESvalues(:,iRES-2)) - max(RESvalues(:,iRES)) < 2e-3; % has to be some convergence
        bool = ~no_conv && direction_of_conv && little_conv;
    else
        bool = false;
    end
end

