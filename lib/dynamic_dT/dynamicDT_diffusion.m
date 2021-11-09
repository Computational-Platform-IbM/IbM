function Time = dynamicDT_diffusion(Time, iDiffusion, iRES, iProf, RESvalues, nDiffIters, constants, dx, ssReached)
    % Function to dynamically reduce or increase the dT for diffusion after
    % the initial simulation hours, based on:
    % 1) whether there is an upward trend in the RES values (decrease)
    % 2) non-convergence in the diffusion (decrease) -> reduce oscillations
    % 3) high iDiffusion numbers for multiple steady states (increase)
    %
    % Also sets the changed_dT value to the current time when changed to
    % limit the amount of consecutive changes.

    % only perform dynamic dT after the initial simulation phase
    if Time.current > 5

        upward_trend = iDiffusion > 50 && ...
                max(RESvalues(:, iRES)) - max(RESvalues(:, iRES - 10)) > 0.01;
        non_convergence = iDiffusion > constants.dynamicDT.iterThresholdDecrease;
        changed_this_ss = Time.current - Time.changed_dT < Time.dT_bac;
        changed_previous_ss = Time.current - Time.changed_dT > constants.dynamicDT.nIterThresholdIncrease*Time.dT_bac;
        multiple_ss_with_high_nIter = all([nDiffIters(iProf-constants.dynamicDT.nIterThresholdIncrease+1:iProf-1)', iDiffusion] > constants.dynamicDT.iterThresholdIncrease);
            
            
        if upward_trend 
            Time.dT = Time.dT * 0.9;
            Time.changed_dT = Time.current;
            fprintf(2, 'upward trend in RES values detected, \n\tthus dT decreased to %g\n', Time.dT);
            fprintf(2, 'Neumann value of stability: %g\n', max(constants.diffusion_rates * Time.dT / (dx ^ 2)))
        elseif non_convergence && ~changed_this_ss
            Time.dT = Time.dT * 0.9;
            Time.changed_dT = Time.current;
            fprintf(2, 'Diffusion takes longer than %d diffusion iterations, \n\tthus dT decreased to %g\n', constants.dynamicDT.iterThresholdDecrease, Time.dT);                
            fprintf(2, 'Neumann value of stability: %g\n', max(constants.diffusion_rates * Time.dT / (dx ^ 2)))
        elseif ssReached && ~changed_previous_ss && multiple_ss_with_high_nIter
            Time.dT = Time.dT * 1.1;
            Time.changed_dT = Time.current;
            fprintf(2, 'Multiple steady states reached with more than %d diffusion iterations, \n\tthus dT increased to %g\n', constants.dynamicDT.iterThresholdIncrease, Time.dT);
            fprintf(2, 'Neumann value of stability: %g\n', max(constants.diffusion_rates * Time.dT / (dx ^ 2)))
        end
    end
end