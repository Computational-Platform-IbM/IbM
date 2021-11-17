clear all

% variable declaration
kDet = 20; % 1/(micrometer h)



if isfile('testScenario.mat')
    load('testScenario.mat', 'grid', 'bac')
    pBac(bac, grid);
else
    [grid, bac] = createTestScenario();
end

% % center of gridcells
% gx = linspace(grid.dx / 2, grid.dx*grid.nX - grid.dx/2, grid.nX);
% gy = linspace(grid.dy / 2, grid.dy*grid.nY - grid.dy/2, grid.nY);
% 
% [X, Y] = meshgrid(gx,gy);
% Z = Y ./ max(gy);

% figure();
% contourf(X,Y,Z,10);
% xlim([0, grid.dx*grid.nX])
% ylim([0, grid.dy*grid.nY])
% colorbar();

%% init detachment
% get grid2nBacs
[grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
biofilm = grid2nBacs > 0;

% plot biofilm
plotLogicalGrid(grid, biofilm, 'Biofilm'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}



% create three matrices with all grid cells
Visited = zeros(grid.nX, grid.nY, 'logical');
Narrow_band = zeros(grid.nX, grid.nY, 'logical');
Far = zeros(grid.nX, grid.nY, 'logical');
T = zeros(grid.nX, grid.nY);

T(~biofilm) = 0; % redundant, because initialised with 0
Visited = ~biofilm;
plotLogicalGrid(grid, Visited, 'Visited'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}


% find narrow band
kernel = zeros(3,3);
kernel([1, 3],2) = -1/4;
kernel(2,[1, 3]) = -1/4;
kernel(2,2) = 1;
temp = zeros(size(biofilm)+2, 'logical');
temp(2:end-1, 2:end-1) = biofilm;
temp(:,1) = 1;
temp(1, :) = temp(end-1, :);
temp(end, :) = temp(2, :);
Narrow_band = convn(temp, kernel, 'valid') > 0;
clear temp

plotLogicalGrid(grid, Narrow_band, 'Narrow band'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}


Far = biofilm & ~Narrow_band;
plotLogicalGrid(grid, Far, 'Far'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}

% set far away point to T infinity
T(Far) = Inf;

[i, j] = find(Narrow_band);
for k = 1:length(i)
    ii = i(k);
    jj = j(k);
    Fdetach = calculateLocalDetachmentRate(ii, jj, kDet, grid);
    
    % --------- IMPORTANT ----------
    % What is the impact of the number of free neighbours?
    % In idynomics it is calculated with the number of non-biomass
    % gridcells, but mathematically it does not make too much sense...?
    % --------END IMPORTANT --------
    
    nFreeNb = 4;
    
    T(ii, jj) = grid.dx / (Fdetach * nFreeNb);
end

plotDetachTime(grid, T);












%% create testScenario
function [grid, bac] = createTestScenario()

    % create test scene/scenario
    javaaddpath([pwd '/lib/bacteria/shovingQuadTreekDist.jar']);

    % create rectangular grid (1000 micrometer by 500 micrometer)
    grid.dx = 4;
    grid.dy = 4;

    grid.nX = 32;
    grid.nY = 16;

    % make new bac struct
    bac = struct();
    bac.x = zeros(20, 1);
    bac.y = zeros(20, 1);
    bac.radius = zeros(20, 1);
    bac.offset = 0.001;

    max_radius = 1;
    kDist = 1.5;

    rRange = [0.3 * max_radius, 0.8 * max_radius];
    bac = fill_rect(bac, [1, 127], [1, 20], rRange);
    bac = fill_rect(bac, [10, 30], [20, 40], rRange);
    bac = fill_rect(bac, [50, 70], [20, 40], rRange);
    bac = fill_rect(bac, [90, 110], [20, 40], rRange);
    qt = shoving.BiomassQuadtree(0, grid.dx*grid.nX, 0, grid.dy*grid.nY);
    r = qt.pushing2D(length(bac.x), bac.x, bac.y, bac.radius, 0.1, max_radius * 2, kDist); % .jar file in /lib is from granule, thus does not have kDist implemented yet
    bac.x = r.bac_x;
    bac.y = r.bac_y;
    
    % remove bacs outside of boundary
    rmv_mask = bac.x < 0 | bac.x > grid.dx*grid.nX | bac.y < 0;
    
    
    bac.x(rmv_mask) = [];
    bac.y(rmv_mask) = [];
    bac.radius(rmv_mask) = [];
    fprintf('removed %d bacteria outside domain\n', sum(rmv_mask));
    
    nBacs = nnz(bac.x);
    bac.x(nBacs+1:end) = [];
    bac.y(nBacs+1:end) = [];
    bac.radius(nBacs+1:end) = [];
    
    fprintf('%d bacteria remained in the domain\n', length(bac.x))
    


    pBac(bac, grid);
    save('testScenario.mat', 'bac', 'grid')
end

%% helper functions
function bac = fill_rect(bac, xRange, yRange, rRange)
    % fill a rectangle with bacteria of random size at blue noise positions
    % (rejection sampling)
    
    nBacs = nnz(bac.x);
    range = [xRange; yRange; rRange];
    cycles_without_addition = 0;
    cycles = 0;
    
    while cycles_without_addition < 1000 && cycles < 100000
        temp_xyr = rand(3, 1) .* (range(:,2) - range(:,1)) + range(:,1); % (x,y,r)
        if ~any(abs(bac.x - temp_xyr(1)) < bac.radius + temp_xyr(3) + bac.offset & abs(bac.y - temp_xyr(2)) < bac.radius + temp_xyr(3) + bac.offset)
            bac.x(nBacs + 1) = temp_xyr(1);
            bac.y(nBacs + 1) = temp_xyr(2);
            bac.radius(nBacs + 1) = temp_xyr(3);
            nBacs = nBacs + 1;
            cycles_without_addition = 0;
        else
            cycles_without_addition = cycles_without_addition + 1;
        end
        cycles = cycles + 1;
    end
    
    fprintf('ended because: ')
    if cycles_without_addition >= 1000
        fprintf('too many cycles without addition\n')
    elseif cycles >= 1000
        fprintf('too many cycles\n')
    else
        fprintf('error?\n')
    end
end

function Fdet = calculateLocalDetachmentRate(~, j, kDet, grid)
    % Calculate the local detachment speed: Fdet = kDet * d^2.
    % Using height in biofilm as measure for now.
    
    % calculate height/distance
    d = (j - 0.5)*grid.dy;
    
    % calculate Fdet
    Fdet = kDet * d^2;
end



%% plotting functions
function pBac(bac, g)
    % visualise inserted bacteria
        % testing function to visualise bacteria
    x = bac.x;
    y = bac.y;
    r = bac.radius;
    
    figure(1); clf;
    col = [0.5, 0.5, 0.5];
    
    for i = 1:length(x)
        rectangle('Curvature', [1 1], 'Position', [x(i) - r(i), y(i) - r(i), 2 * r(i), 2 * r(i)], 'LineWidth', 0.1, 'EdgeColor', col, 'FaceColor', col);
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

function plotLogicalGrid(grid, L, title_string)
    switch title_string
        case {'Biofilm', 'biofilm'}
            figure(2);
        case {'Visited', 'visited'}
            figure(3);
        case {'Narrow_band', 'narrow_band', 'Narrow band', 'narrow band'}
            figure(4);
        case {'Far', 'far'}
            figure(5);
        otherwise
            fprintf(2, 'no figure yet for %s, new figure created...\n', title_string);
            figure();
    end
            
    im = imagesc([grid.dx/2, grid.nX*grid.dx - grid.dx/2], [grid.dy/2, grid.nY*grid.dy - grid.dy/2], L'); hold on;
%     im.AlphaData = 0.5;
    colormap([0 0 0; 1 1 1]);
    axis equal;
    xlim([0, grid.nX*grid.dx]);
    ylim([0, grid.nY*grid.dy]);
    set(gca, 'YDir', 'normal'); 
    title(title_string);
end

function plotDetachTime(grid, T)
    figure(6); clf;
    imagesc([grid.dx/2, grid.nX*grid.dx - grid.dx/2], [grid.dy/2, grid.nY*grid.dy - grid.dy/2], T'); hold on;
    colormap(viridis());
    colorbar();
    caxis([0, 1e-3])
    axis equal;
    xlim([0, grid.nX*grid.dx]);
    ylim([0, grid.nY*grid.dy]);
    set(gca, 'YDir', 'normal'); 
    title('Detachment times [h]');
end




