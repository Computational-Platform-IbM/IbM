function plotDiffRegion(g, bac, diffRegion, plotBoundaryThickness) %, checked_cell, active_cell)
    % debugging function to visualise the determination of the diffusion region
    x = bac.x;
    y = bac.y;
    
    diffRegion = double(diffRegion);
%     diffRegion(checked_cell(1), checked_cell(2)) = 0.4;
%     diffRegion(active_cell(1), active_cell(2)) = 0.6;
    diffRegion = diffRegion';
    
    figure('Name', 'Diffusion region'); clf;
    
    % plot diff region
    im = imagesc([g.dx/2, g.nX*g.dx - g.dx/2], [g.dy/2, g.nY*g.dy - g.dy/2], diffRegion); hold on;
    im.AlphaData = 0.5;
    set(gca, 'YDir', 'normal');
    colormap([0 0 0; 1 0.1 0.1; .1 1 .1; 1 1 1]);
%     active_bac = floor(x / g.dx) + 1 == active_cell(1) & floor(y / g.dy) + 1 == active_cell(2);
    c = ones(size(x)) * [0 .7 .7];
%     c(active_bac, :) = ones(nnz(active_bac), 1) * [0 .55 .55];
    
    
    
    % plot bacs
    scatter(x, y, [], c, 'filled', 'MarkerEdgeColor', [0 .5 .5], 'LineWidth', 1.5); hold on;
    
    % plot circles around each bacterium with boundary layer thickness
    if plotBoundaryThickness
        r = g.blayer_thickness;
        for i = 1:length(x)
            rectangle('Curvature', [1 1], 'Position', [x(i) - r, y(i) - r, 2 * r, 2 * r], 'LineWidth', 1, 'EdgeColor', [0, 0.5, 0.5 1]);
        end
    end
    
    axis equal;

    xlim([0, g.nX*g.dx]);
    ylim([0, g.nY*g.dy]);
    xt = linspace(0, g.nX*g.dx, g.nX+1);
    yt = linspace(0, g.nY*g.dy, g.nY+1);
    xticks(xt);
    yticks(yt);    
    grid on;
end
