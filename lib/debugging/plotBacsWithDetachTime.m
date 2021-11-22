function plotBacsWithDetachTime(bac, g, T)
    % visualise inserted bacteria
        % testing function to visualise bacteria
    x = bac.x;
    y = bac.y;
    r = bac.radius;
    
    figure(2); clf; % same figure as plotting bacteria
    col = [0 0 0 0.2];
    
    if ~isscalar(T)
        % center of gridcells
        gx = linspace(g.dx / 2, g.dx*g.nX - g.dx/2, g.nX);
        gy = linspace(g.dy / 2, g.dy*g.nY - g.dy/2, g.nY);

        [X, Y] = meshgrid(gx,gy);

        contourf(X,Y,T',0:10); hold on;
        colormap(viridis());
        colorbar();
    end
    
    
    for i = 1:length(x)
        rectangle('Curvature', [1 1], 'Position', [x(i) - r(i), y(i) - r(i), 2 * r(i), 2 * r(i)], 'LineWidth', 0.1, 'EdgeColor', col*1.5, 'FaceColor', col);
    end
    
    axis equal;

    xlim([0, g.nX*g.dx]);
    ylim([0, g.nY*g.dy]);
    xt = linspace(0, g.nX*g.dx, g.nX+1);
    yt = linspace(0, g.nY*g.dy, g.nY+1);
    xticks(xt(1:10:end));
    yticks(yt(1:10:end));
    xticklabels(xt(1:10:end));
    yticklabels(yt(1:10:end));
    
    ax = gca;
    ax.XRuler.MinorTick = 'on'; %or 'off'
    ax.XRuler.MinorTickValues = xt; %just like major ticks
    ax.XRuler.MinorTickValuesMode = 'manual'; %or 'manual'
    
    ax.YRuler.MinorTick = 'on'; %or 'off'
    ax.YRuler.MinorTickValues = yt; %just like major ticks
    ax.YRuler.MinorTickValuesMode = 'manual'; %or 'manual'
    
    grid on;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    
    ylabel('y-position [μm]');
    xlabel('x-position [μm]');
    

end