function bool = multiple_high_iters(iDiffusion, iProf, nDiffIters, Time, constants)
    % Detect whether the previous steady states have been reached with a
    % high number of diffusion iterations
    
    changed_previous_ss = Time.current - Time.changed_dT < constants.dynamicDT.nIterThresholdIncrease*Time.dT_bac;
    if iProf > constants.dynamicDT.nIterThresholdIncrease
        multiple_ss_with_high_nIter = all([nDiffIters(iProf-constants.dynamicDT.nIterThresholdIncrease+1:iProf-1)', iDiffusion] > constants.dynamicDT.iterThresholdIncrease);
    else
        multiple_ss_with_high_nIter = false;
    end
    bool = ~changed_previous_ss && multiple_ss_with_high_nIter;
end
