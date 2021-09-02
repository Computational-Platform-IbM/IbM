function plotBacs(g, bac)
    % testing function to visualise bacteria
    x = bac.x;
    y = bac.y;
    r = bac.radius;
    s = bac.species;
    
    coloring = {'#D81B60', '#1E88E5', '#FFC107', '#004D40'}; % colorblind-friendly colours
    
    f = figure(2); clf;
%     f.Position = [-1800, 65, 1200, 900]; % desktop with two screens
    f.Position = [60 60 1000 700];
    
    for i = 1:size(x, 1)
        rectangle('Curvature', [1 1], 'Position', [x(i) - r(i), y(i) - r(i), 2 * r(i), 2 * r(i)], 'LineWidth', 1, 'EdgeColor', [0, 0, 0], 'FaceColor', coloring{s(i)});
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