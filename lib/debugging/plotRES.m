function plotRES(RES, compoundName, diffRegion)
    
    % plot RES over entire domain
    figure(11); clf;
    ax1 = axes;
    imagesc(ax1, RES');
    max_val = max(abs(RES(diffRegion)));
    caxis([-max_val, max_val])
    title(compoundName)
    axis square;
    
    % plot diffRegion on top (with opacity)
    hold on;

    ax2 = axes;
    h = imagesc(ax2, diffRegion');
    h.AlphaData = 0.5;
    axis square;

    colorbar(ax1);

    linkaxes([ax1, ax2]);
    set(ax2, 'Position', ax1.Position)
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    colormap(ax1, redblue());
    colormap(ax2, [0 0 0; 1 1 1]);
    hold off;
    
    set(ax1, 'YDir', 'normal')
    set(ax2, 'YDir', 'normal')
end