clear all

% variable declaration
kDet = 200e-6; % 1/(micrometer h)



if isfile('testScenario.mat')
    load('testScenario.mat', 'grid', 'bac')
%     pBac(bac, grid, 0);
else
    [grid, bac] = createTestScenario();
%     pBac(bac, grid, 0);
end



%% init detachment
% grid.dx = 4;
% grid.dy = 4;
% grid.nX = 128;
% grid.nY = 128;

% get grid2nBacs
[grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
timer = tic;

% ----- IMPORTANT ------
% use diffregion here for the real version
biofilm = convn(grid2nBacs, ones(3), 'same') > 0;
% ----END IMPORTANT ----



% plot biofilm
% plotLogicalGrid(grid, biofilm, 'Biofilm'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}



% create three matrices with all grid cells
T = zeros(grid.nX, grid.nY);

% T(~biofilm) = 0; % redundant, because initialised with 0
Visited = ~biofilm;
% plotLogicalGrid(grid, Visited, 'Visited'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}


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

% plotLogicalGrid(grid, Narrow_band, 'Narrow band'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}


Far = biofilm & ~Narrow_band;
% plotLogicalGrid(grid, Far, 'Far'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}

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
    
    nFreeNb = getFreeNeighbourCount(ii, jj, Visited, grid);
    
    T(ii, jj) = grid.dx / (Fdetach * nFreeNb);
end

% plotDetachTime(grid, T, kDet);


% ------- IMPORTANT ------
% use of heapstack for better performance in narrow_band
% ----- END IMPORTANT ----


% ------------ Fast Marching --------------
while sum(Narrow_band, 'all')
    % take narrow-band value with the lowest T value
    [T_list, I] = sort(T(Narrow_band));
    [i, j] = find(Narrow_band);
    i = i(I);
    j = j(I);

    i_point = i(1);
    j_point = j(1); 

    % add point to visited
    Narrow_band(i_point, j_point) = 0;
    Visited(i_point, j_point) = 1;

    % for all neighbours of the point (from narrow-band or far):
    [Narrow_band, Far, T] = updateNeighbour(mod(i_point+1-1, grid.nX)+1, j_point, Narrow_band, Far, Visited, T, grid, kDet);
    [Narrow_band, Far, T] = updateNeighbour(mod(i_point-1-1, grid.nX)+1, j_point, Narrow_band, Far, Visited, T, grid, kDet);

    [Narrow_band, Far, T] = updateNeighbour(i_point, mod(j_point+1-1, grid.nY)+1, Narrow_band, Far, Visited, T, grid, kDet);
    if j_point ~= 1
        [Narrow_band, Far, T] = updateNeighbour(i_point, mod(j_point-1-1, grid.nY)+1, Narrow_band, Far, Visited, T, grid, kDet);
    end

    
% ------- DEBUG PLOTS -------
%     plotLogicalGrid(grid, Visited, 'Visited'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
%     plotLogicalGrid(grid, Narrow_band, 'Narrow band'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
%     plotLogicalGrid(grid, Far, 'Far'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
%     plotDetachTime(grid, T, kDet);
%     drawnow()
% ----- END DEBUG PLOTS -----


end
toc(timer)

% ------- DEBUG PLOTS -------
%     plotLogicalGrid(grid, Visited, 'Visited'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
%     plotLogicalGrid(grid, Narrow_band, 'Narrow band'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
%     plotLogicalGrid(grid, Far, 'Far'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
%     plotDetachTime(grid, T, kDet);
%     drawnow()
% ----- END DEBUG PLOTS -----



% ----- ANOTHER DEBUG PLOT --------
pBac(bac, grid, T);

figure(8); clf;
% center of gridcells
gx = linspace(grid.dx / 2, grid.dx*grid.nX - grid.dx/2, grid.nX);
gy = linspace(grid.dy / 2, grid.dy*grid.nY - grid.dy/2, grid.nY);

[X, Y] = meshgrid(gx,gy);

% smooth the graph a bit...
T_smooth = convn(T, (1/9)*ones(3), 'same');
T_smoother = convn(T_smooth, (1/9)*ones(3), 'same');

contourf(X,Y,T_smoother',[0,1,2,3,4,5,6,7,8,9,10]); hold on;
colormap(viridis());
colorbar();
xlim([0, grid.nX*grid.dx]);
ylim([0, grid.nY*grid.dy]);
axis equal
% ------ END DEBUG PLOT ------






%% create testScenario
function [grid, bac] = createTestScenario()

    % create test scene/scenario
    javaaddpath([pwd '/lib/bacteria/shovingQuadTreekDist.jar']);

    % create rectangular grid (1000 micrometer by 500 micrometer)
    grid.dx = 8;
    grid.dy = 8;

    grid.nX = 64;
    grid.nY = 64;

    % make new bac struct
    bac = struct();
    bac.x = zeros(20, 1);
    bac.y = zeros(20, 1);
    bac.radius = zeros(20, 1);
    bac.offset = 0.001;

    max_radius = 2;
    kDist = 1.3;

    rRange = [0.5 * max_radius, max_radius];
    bac = fill_rect(bac, [1, 512], [1, 100], rRange);
    bac = fill_rect(bac, [100, 300], [100, 400], rRange);
%     bac = fill_rect(bac, [50, 70], [100, 400], rRange);
%     bac = fill_rect(bac, [90, 110], [100, 400], rRange);
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
    


    pBac(bac, grid, 0);
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

function [Narrow_band, Far, T] = updateNeighbour(i, j, Narrow_band, Far, Visited, T, grid, kDet)
    % update status of grid cell with regards to detachment and recalculate
    % the Time-of-crossing using the quadratic approximation of the
    % gradient.
    
    if Visited(i, j)
        return
    elseif Far(i, j) % remove from far and add to narrow band
        Far(i, j) = 0;
        Narrow_band(i, j) = 1;
    end % if none of the above, then already in narrow band
    
    % recalculate T value
    T(i, j) = recalculateT(T, i, j, kDet, grid, Visited);
end

function T_new = recalculateT(T, i, j, kDet, grid, Visited)
    % Recalculate the T value at gridcell (i, j) using the quadratic 
    % approximation of the gradient.
    % using method of https://github.com/kreft/iDynoMiCS/blob/master/src/simulator/detachment/LevelSet.java
    
    right_visited = Visited(mod(i+1-1, grid.nX)+1, j); % x + 1 (right)
    left_visited = Visited(mod(i-1-1, grid.nX)+1, j); % x - 1 (left)
    if  left_visited
        if right_visited
            Tx = min(T(mod(i+1-1, grid.nX)+1, j), T(mod(i-1-1, grid.nX)+1, j));
        else
            Tx = T(mod(i-1-1, grid.nX)+1, j);
        end
    elseif right_visited
        Tx = T(mod(i+1-1, grid.nX)+1, j);
    else
        Tx = Inf;
    end
    
    top_visited = Visited(i, j+1); % j + 1 (top)
    bottom_visited = j ~= 1 && Visited(i, j-1); % j - 1 (bottom)
    
    if top_visited
        if bottom_visited
            Ty = min(T(i, j+1), T(i, j-1));
        else
            Ty = T(i, j+1);
        end
    elseif bottom_visited
        Ty = T(i, j-1);
    else
        Ty = Inf;
    end
    
    if isinf(Tx) && isinf(Ty)
        error('all neighbours have infinite time of crossing')
    end
    
    Fdet = calculateLocalDetachmentRate(i, j, kDet, grid);
    
    if Fdet == 0
        error('Detachment speed equals 0, thus infinite time of crossing')
    end
    
    T_new = computeRoot(Tx, Ty, Fdet, grid.dx);
end

function root = computeRoot(Tx, Ty, Fdet, dx)
    % Solve the quadratic equation for the gradient (approximation)
    
    % first init all parameters
    a = 0;
    b = 0;
    c = -(dx/Fdet)^2;
    
    % if Tx is finite, then add those respective terms to parameters
    if isfinite(Tx)
        a = a + 1;
        b = b - 2*Tx;
        c = c + Tx^2;
    end
    
    % same for Ty
    if isfinite(Ty)
        a = a + 1;
        b = b - 2*Ty;
        c = c + Ty^2;
    end
    
    % if all are infinite, then return inf
    if a == 0
        root = Inf;
        return
    end
    
    % now get the 2 solutions
    D = sqrt(b^2 - 4*a*c);
    
    if D < 0
        error('ValueError:Should always be above 0...')
    end
    
    % positive solution is only valid
    root = (-b+D) / (2*a);
end

function nFreeNb = getFreeNeighbourCount(i, j, Visited, grid)
    % count the number of non-biomass neighbouring gridcells are there for
    % gridcell (i, j)
    
    nFreeNb = 0;
    if Visited(mod(i+1-1, grid.nX)+1, j)
        nFreeNb = nFreeNb + 1;
    end
    if Visited(mod(i-1-1, grid.nX)+1, j)
        nFreeNb = nFreeNb + 1;
    end
    if Visited(i, j+1)
        nFreeNb = nFreeNb + 1;
    end
    if j>1 && Visited(i, j-1)
        nFreeNb = nFreeNb + 1;
    end
end

%% plotting functions
function pBac(bac, g, T)
    % visualise inserted bacteria
        % testing function to visualise bacteria
    x = bac.x;
    y = bac.y;
    r = bac.radius;
    
    figure(1); clf;
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

function plotLogicalGrid(grid, L, title_string)
    switch title_string
        case {'Biofilm', 'biofilm'}
            f = figure(2);
        case {'Visited', 'visited'}
            f = figure(3);
        case {'Narrow_band', 'narrow_band', 'Narrow band', 'narrow band'}
            f = figure(4);
        case {'Far', 'far'}
            f = figure(5);
        otherwise
            fprintf(2, 'no figure yet for %s, new figure created...\n', title_string);
            f = figure();
    end
            
    try
        f.Children.Children(1).CData = L';
        drawnow();
    catch e
        switch e.identifier
            case 'MATLAB:noSuchMethodOrField' % has not been made before
                im = imagesc([grid.dx/2, grid.nX*grid.dx - grid.dx/2], [grid.dy/2, grid.nY*grid.dy - grid.dy/2], L'); hold on;
        %           im.AlphaData = 0.5;
                colormap([0 0 0; 1 1 1]);
                axis equal;
                xlim([0, grid.nX*grid.dx]);
                ylim([0, grid.nY*grid.dy]);
                set(gca, 'YDir', 'normal'); 
                title(title_string);
            otherwise
                rethrow(e)
        end
                
    end
    
end

function plotDetachTime(grid, T, kDet)
    f = figure(6);
    
    try
        ax = gca;
        ax.Children(1).CData = T';
        caxis([0, 12]);
        drawnow();
    catch e
        switch e.identifier
            case {'MATLAB:noSuchMethodOrField', 'MATLAB:noPublicFieldForClass'} % has not been made before
                imagesc([grid.dx/2, grid.nX*grid.dx - grid.dx/2], [grid.dy/2, grid.nY*grid.dy - grid.dy/2], T'); hold on;
                colormap(viridis());
                colorbar();
                caxis([0, 12])
                axis equal;
                xlim([0, grid.nX*grid.dx]);
                ylim([0, grid.nY*grid.dy]);
                set(gca, 'YDir', 'normal'); 
                title('Detachment times [h]');
            otherwise
                rethrow(e)
        end
    end
end




