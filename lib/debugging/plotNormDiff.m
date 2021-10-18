function plotNormDiff(norm_diff, iRES, constants, Time)
    % Debug function to plot the difference in concentration as the
    % norm(conc - conc_prev).
    
    figure(9); clf;
    plot((1:iRES-1)*constants.nDiffusion_per_SScheck, norm_diff(2:iRES), 'LineWidth', 2);
    set(gca, 'YScale', 'log')
    title(sprintf('norm(conc_{prev} - conc) at t=%.1f', Time))
end