%% testing
clear ans bac constants diffusion_region grid grid2bac grid2nBacs;

% init for determine_where_bacteria_in_grid
grid = struct;
grid.nX = 20;
grid.nY = 20;
grid.dx = .5;
grid.dy = .5;

% init for determine_diffusion_region
grid.blayer_thickness = 5/3 * grid.dx;

% bac init in circle
n = 4;
radius = 2;
bac = struct;
[bac.x, bac.y] = rand_circle(n, grid.nX/2*grid.dx, grid.nY/2*grid.dy, radius);
bac.r = ones(n, 1)*0.2;
clear n radius

% check manually for correct indices of bacteria
[grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
plotBacs(grid, bac);

% check manually for correct diffusion region
diffusion_region = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid);
plotDiffRegion(grid, bac, diffusion_region);


%%%% testing for reaction matrix:
bac.molarMass = ones(4, 1);
bac.species = 1:4;

constants = struct;
constants.T = 293;
constants.Keq = R.pTh.Keq;
constants.chrM = R.pTh.chrM;
constants.Vg = grid.dx ^ 3;
constants.StNames = R.St.StNames(1:8);
constants.react_v = R.pTh.react_v;
constants.Ks = R.pTh.Ks(:, 1:3);
constants.Ki = R.pTh.Ks(:, 4:6);
constants.MatrixMet = R.rm.MatrixMet;
constants.MatrixDecay = R.rm.MatrixDecay_mod;

conc = ones(grid.nX, grid.nY, 8);
for ci = 1:length(R.Sxy.Sbc_Dir)
    conc(:,:, ci) = R.Sxy.Sbc_Dir(ci);
end


[reaction_matrix, mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, grid, conc, constants, 7);





%% helper functions
function [X, Y] = rand_circle(N, x_center, y_center, r)
    Ns = round(4/pi * N + 2.5*sqrt(N) + 100);
    X = rand(Ns,1)*(2*r) - r;
    Y = rand(Ns,1)*(2*r) - r;
    I = find(sqrt(X.^2 + Y.^2)<=r);
    X = X(I(1:N)) + x_center;
    Y = Y(I(1:N)) + y_center;
end

function plotBacs(g, bac)
    % testing function to visualise bacteria
    x = bac.x;
    y = bac.y;
    
    dx_text = 0.1;
    dy_text = 0.1;
    
    labels = num2str([1:size(x, 1)]');
    
    figure(1); clf;
    
%     for i = 1:size(x, 1)
%         rectangle('Curvature', [1 1], 'Position', [x(i) - r(i), y(i) - r(i), 2 * r(i), 2 * r(i)], 'LineWidth', 1, 'EdgeColor', [0, 0, 0]);
%     end
    scatter(x, y, 'filled', 'MarkerFaceColor', [0, 0.7, 0.7], 'MarkerEdgeColor', [0 .5 .5], 'LineWidth', 1.5); hold on;
    text(x+dx_text, y+dy_text, labels);
    
    axis equal;

    xlim([0, g.nX*g.dx]);
    xticks(0:g.nX);
    ylim([0, g.nY*g.dy]);
    yticks(0:g.nY);
    grid on;
end

function plotDiffRegion(g, bac, diffRegion)
    % testing function to visualise bacteria and the diffusion region
    x = bac.x;
    y = bac.y;

    diffRegion = diffRegion';
    figure(2); clf;
    
    % plot diff region
    im = imagesc([g.dx/2, g.nX*g.dx - g.dx/2], [g.dy/2, g.nY*g.dy - g.dy/2], diffRegion); hold on;
    im.AlphaData = 0.5;
    set(gca, 'YDir', 'normal');
    colormap([0 0 0; 1 1 1]);
    
    % plot bacs
    scatter(x, y, 'filled', 'MarkerFaceColor', [0, 0.7, 0.7], 'MarkerEdgeColor', [0 .5 .5], 'LineWidth', 1.5); hold on;
    
    % plot circles around each bacterium with boundary layer thickness
    r = g.blayer_thickness;
    for i = 1:size(x, 1)
        rectangle('Curvature', [1 1], 'Position', [x(i) - r, y(i) - r, 2 * r, 2 * r], 'LineWidth', 1, 'EdgeColor', [0, 0.5, 0.5 0.5]);
    end
    
    axis equal;

    xlim([0, g.nX*g.dx]);
    xticks(0:g.nX*g.dx);
    ylim([0, g.nY*g.dy]);
    yticks(0:g.nY*g.dy);
    grid on;
end