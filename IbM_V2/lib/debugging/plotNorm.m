function plotNorm(Norm, dT)
    % debug plotting of norm of difference between conc1 and conc0. Representative of steady state over the entire
    % simulation time
    %
    % Norm: norm of diff(conc - conc_prev) per dT_bac over the entire simulation domain
    
    figure('Name', 'Difference (norm) between consecutive concentration profiles');
    plot((1:length(Norm))*dT, Norm, "LineWidth", 2);
    set(gca, 'YScale', 'log')
    xlabel('Simulation time [h]')
    ylabel('norm [mol/L]')
    title('Norm over simulation time')
%     yl = ylim();
%     ylim([0, yl(2)])
end

