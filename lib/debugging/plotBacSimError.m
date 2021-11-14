function plotBacSimError(res_bacsim, iRES, constants, Time)
    % Debug function to plot the residual as computed by BacSim (max diff
    % of concentrations divided by dT) and similarly the norm diff of concs
    % divided by dT
    
    figure(12); clf;
    plot((1:iRES-1)*constants.nDiffusion_per_SScheck, res_bacsim(2:iRES, 1), 'LineWidth', 2); hold on;
    yline(1, '--', 'Threshold', 'HandleVisibility','off'); 
    ylabel('bacsim residual');
    ylim([0, 10]);
    yyaxis right
    plot((1:iRES-1)*constants.nDiffusion_per_SScheck, res_bacsim(2:iRES, 2), 'LineWidth', 2); hold off;
    ylabel('norm / dT');
    ylim([0, 100])
    legend('bacsim res', 'norm/dT');
    xlabel('diffusion iteration');
    xlim([0, inf]);
    
%     set(gca, 'YScale', 'log')
    title(sprintf('[mol/h] residuals at t=%.1f', Time))
end