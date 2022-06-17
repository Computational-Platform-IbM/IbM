function Time = decrease_dT_diffusion(Time, msg, dx, constants)
    % Utility function for decreasing the dT for diffusion
    
    prev = Time.dT;
    
    Time.dT = max(Time.dT * 0.9, Time.minDT);
    if Time.dT ~= prev % only set changed dT time if actually changed
        Time.changed_dT = Time.current;
    end
    if Time.dT ~= constants.dynamicDT.minDT
        fprintf(2, '%s, \n\tthus dT decreased to %g\n', msg, Time.dT)
        fprintf(2, 'Neumann value of stability: %g\n', max(constants.diffusion_rates * Time.dT / (dx ^ 2)))
    end
end

