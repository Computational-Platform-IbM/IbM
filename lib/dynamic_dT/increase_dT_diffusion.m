function Time = increase_dT_diffusion(Time, msg, dx, constants)
    % Utility function for increasing the dT for diffusion
    
    Time.dT = min(Time.dT / 0.9, Time.maxDT);
    Time.changed_dT = Time.current;
    fprintf(2, '%s, \n\tthus dT increased to %g\n', msg, Time.dT)
    fprintf(2, 'Neumann value of stability: %g\n', max(constants.diffusion_rates * Time.dT / (dx ^ 2)))
end