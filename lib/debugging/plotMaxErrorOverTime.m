function plotMaxErrorOverTime(maxErrors, dT)
    % debug plotting of maximum error from steady state over the entire
    % simulation time
    %
    % maxErrors: maximum error per dT_bac over the entire simulation domain
    
    figure(5);
    plot((1:length(maxErrors))*dT, maxErrors*100, "LineWidth", 2);
    xlabel('Simulation time [h]')
    ylabel('Maximum error [%]')
    title('Maximum error over simulation time')
    yl = ylim();
    ylim([0, yl(2)])
end