function Time = dynamicDT_bac(Time, maxInitRES, iProf, constants)
    % Dynamically adjust the dT_bac to limit the relative/absolute amount
    % of change in the bulk_concentrations
    % based on:
    % 1) high initial RES values (decrease)
    % 2) multiple low initial RES values (increase)
    %
    % Also sets the last time the dT_bac was changed, to prevent multiple
    % consecutive dT_bac changes.
    
    
    
    
    if Time.current > 5
%         initRES_above_threshold = maxInitRES(iProf) > constants.dynamicDT.initRESThresholdDecrease;
%         upward_trend = all(diff(maxInitRES((iProf - constants.dynamicDT.nIterThresholdIncrease + 1):iProf)) >= 0);
        recently_changed = Time.current - Time.changed_dT_bac < constants.dynamicDT.nIterThresholdIncrease*Time.dT_bac;
        if iProf >= constants.dynamicDT.nIterThresholdIncrease % can only check if iProf > n otherwise negative indexing
            multiple_ss_with_low_initRES = all(maxInitRES(iProf - constants.dynamicDT.nIterThresholdIncrease + 1 : iProf ) < constants.dynamicDT.initRESThresholdIncrease);
        else
            multiple_ss_with_low_initRES = false;
        end
%         if  initRES_above_threshold && upward_trend 
%             Time.dT_bac = Time.dT_bac * 0.7;
%             Time.changed_dT_bac = Time.current;
%             fprintf(2, 'large initial RES value detected, dT_bac decreased to %g\n', Time.dT_bac);
%         else
        if ~recently_changed && multiple_ss_with_low_initRES
            Time.dT_bac = min(Time.dT_bac / 0.8, Time.maxDT_bac);
            Time.changed_dT_bac = Time.current;
            fprintf(2, 'multiple cycles with initRES value of <%.1f %%, dT_bac increased to %g\n', constants.dynamicDT.initRESThresholdIncrease*100, Time.dT_bac);                        
        end
    end
end