function bool = multiple_low_initRES(iProf, maxInitRES, Time, constants)
    % Determine whether the previous initial RES values were all below the
    % threshold
    
    recently_changed = Time.current - Time.changed_dT_bac < constants.dynamicDT.nIterThresholdIncrease*Time.dT_bac;
    if iProf >= constants.dynamicDT.nIterThresholdIncrease % can only check if iProf > n otherwise negative indexing
        multiple_ss_with_low_initRES = all(maxInitRES(iProf - constants.dynamicDT.nIterThresholdIncrease + 1 : iProf ) < constants.dynamicDT.initRESThresholdIncrease);
    else
        multiple_ss_with_low_initRES = false;
    end
    
    bool = ~recently_changed && multiple_ss_with_low_initRES;

end