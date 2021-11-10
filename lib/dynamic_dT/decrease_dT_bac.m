function Time = decrease_dT_bac(Time, msg)
    % Utility function for decreasing the dT for bacterial stepping
    
    Time.dT_bac = min(Time.dT_bac / 0.8, Time.maxDT_bac);
    Time.changed_dT_bac = Time.current;
    fprintf(2, '%s, \n\tthus dT_bac decreased to %g', msg, Time.dT_bac)
end