function [diffusion_region, extraction_region] = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid)
    % Determine which grid cells lay in the diffusion layer, i.e. the
    % region composed of all grid cells within the boundary distance from a
    % bacterium.
    % 
    % grid2bac: matrix with per grid cell which bacteria reside there
    % grid2nBacs: matrix with how many bacteria per grid cell
    % bac: struct containing all information regarding the bacteria
    % grid: struct containing all information regarding the grid
    % constants: struct containing all simulation constants
    %
    % -> diffusion_region: logical matrix with per grid cell whether it is
    %   in the diffusion region (true) or not (false).
    
    diffusion_region = zeros(grid.nX, grid.nY, 'logical');
    
    % get diffusion region
    [diffNodesX, diffNodesY] = getDiffusionNodes(bac, grid);    
    
    % in the diffusion region, apply convolution to find boundary of grid
    % cells with bacteria
    kernel = ones(3,3)/-8;
    kernel(2,2) = 1;
    hasBac = grid2nBacs(diffNodesX, diffNodesY) > 0;
    isBacBoundary = convn(hasBac, kernel, 'same') > 0;
    
    % for the boundary of bacterial grid cells, compute for neighbouring
    % grid cells whether they are in diffusion region
    maxOffsetX = ceil(grid.blayer_thickness / grid.dx);
    dx = find(diffNodesX == 1, 1, 'first') - 1;
    maxOffsetY = ceil(grid.blayer_thickness / grid.dy);
    dy = find(diffNodesY == 1, 1, 'first') - 1;
    
    [i, j] = find(isBacBoundary); % find indices in restricted matrix
    isDiffRegion = hasBac; % working matrix: 1 == diffregion, 0 == bulk, restricted from the entire matrix
    
    for n = 1:length(i)
        for di = -maxOffsetX:maxOffsetX
            for dj = -maxOffsetY:maxOffsetY
                if isDiffRegion(i(n)+di, j(n)+dj) ~= 0
                    continue % skip cell if already in diff region
                else
                    % perform actual check between gridcell and bacteria
                    bacs = nonzeros(grid2bac(i(n)+dx, j(n)+dy,:));
                    gridcell_center = [grid.dx*(dx+i(n)+di) - grid.dx/2 ...
                        grid.dy*(dy+j(n)+dj) - grid.dy/2];

% % <DEBUGGING>
% dr = diffusion_region;
% dr(diffNodesX, diffNodesY) = isDiffRegion;
% plotDiffRegion(grid, bac, dr, [dx+i(n)+di, dy+j(n)+dj], [dx+i(n), dy+j(n)]);
% pause(0.01);
% % </DEBUGGING>
                    for iBac = 1:length(bacs)
                        if isWithinBoundaryLayer(bac.x(bacs(iBac)), bac.y(bacs(iBac)), gridcell_center, grid.dx, grid.dy, di, dj, grid.blayer_thickness)
%                         if isWithinBoundaryLayer_v2(bac.x(bacs(iBac)), bac.y(bacs(iBac)), gridcell_center, grid.blayer_thickness)
                            isDiffRegion(i(n)+di, j(n)+dj) = 1;
                            break
                        end
                    end
                end
            end
        end
    end
    
    % put restricted matrix back in original
    diffusion_region(diffNodesX, diffNodesY) = isDiffRegion;
    
    
    % determine extraction region, which:
    % - includes all diffusion nodes
    % - has an odd number of nodes in x & y direction
    % - has at least one bulk layer node at each side
    
    extraction_region = determine_extraction_region(diffusion_region);
end



function extraction_region = determine_extraction_region(diffRegion)
    extraction_region = struct;
    extraction_region.mask = zeros(size(diffRegion), 'logical');
    
    x_diffRegion = sum(diffRegion, 1);
    first_x = find(x_diffRegion, 1, 'first') - 1;
    last_x = find(x_diffRegion, 1, 'last') + 1;
    
    dx = last_x - first_x + 1; % number of diffusion cells
    if mod(dx, 2) == 0
        first_x = first_x - 1;
    end
    
    y_diffRegion = sum(diffRegion, 2);
    first_y = find(y_diffRegion, 1, 'first') - 1;
    last_y = find(y_diffRegion, 1, 'last') + 1;
    
    dy = last_y - first_y + 1; % number of diffusion cells
    if mod(dy, 2) == 0
        first_y = first_y - 1;
    end
    if first_x < 1 || last_x > size(diffRegion, 1) || first_y < 1 || last_y > size(diffRegion, 2)
        warning('DEBUG:actionRequired', 'debug: not enough bulk liquid present around granule.');
    end

    extraction_region.x0 = first_x;
    extraction_region.x1 = last_x;
    extraction_region.y0 = first_y;
    extraction_region.y1 = last_y;
    extraction_region.mask(first_x:last_x, first_y:last_y) = true;
end

function [diffusionNodesX, diffusionNodesY] = getDiffusionNodes(bac, grid)
    % Determine in both x and y direction which of the nodes are
    % potentially in the diffusion region.
    %
    % bac: struct containing all information regarding the bacteria
    % grid: struct containing all information regarding the grid
    %
    % -> diffusionNodesX, diffusionNodesY: logical arrays that signify
    %   per direction which nodes are potentially in the diffusion region.
    
    x_min = min(bac.x);
    x_max = max(bac.x);
    y_min = min(bac.y);
    y_max = max(bac.y);
    dx = grid.dx; dy = grid.dy;

    % define an offset, such that the boundary layer always falls in the right node
    offsetX = grid.blayer_thickness + dx;
    offsetY = grid.blayer_thickness + dy;

    % create an array (columnvector) with all end coordinates of the nodes
    nodeEndCoordinatesX = (1:grid.nX)' * dx;
    nodeEndCoordinatesY = (1:grid.nY)' * dy;

    % check which nodes have either the boundary layer or biofilm in them
    diffusionNodesX = nodeEndCoordinatesX > x_min - offsetX ...
        & nodeEndCoordinatesX - dx < x_max + offsetX;
    diffusionNodesY = nodeEndCoordinatesY > y_min - offsetY ...
        & nodeEndCoordinatesY - dy < y_max + offsetY;
end

function isBLayer = isWithinBoundaryLayer(bac_x, bac_y, gridcell_center, dx, dy, di, dj, blayer_thickness)
    % Determine whether a gridcell is within the boundary layer thickness 
    % of a specific bacterium. 8 distinct scenarios, one for each
    % direction. For each scenario a different measure is applied.
    %
    % bac_x, bac_y: bacterial coordinates
    % gridcell_center: center coordinates (x, y) of the gridcell
    % dx, dy: width and height of each gridcell
    % di, dj: respective grid position of checked cell w.r.t. bacterial
    %   containing cell
    % blayer_thickness: size of the boundary layer
    %
    % -> isBLayer: boolean with whether the checked gridcell has any
    % boundary layer in it.
    
    if di == 0
        if dj > 0
            % check bottom y coordinate
            grid_edge_y = gridcell_center(2) - dy/2;
            isBLayer = blayer_thickness > grid_edge_y - bac_y;
            
        else % dj < 0
            % check top y coordinate
            grid_edge_y = gridcell_center(2) + dy/2;
            isBLayer = blayer_thickness > bac_y - grid_edge_y;

        end
    elseif dj == 0
        if di < 0
            % check right x coordinate
            grid_edge_x = gridcell_center(1) + dx/2;
            isBLayer = blayer_thickness > bac_x - grid_edge_x;
            
        else % di > 0
            % check left x coordinate
            grid_edge_x = gridcell_center(1) - dx/2;
            isBLayer = blayer_thickness > grid_edge_x - bac_x;
            
        end
    elseif di < 0 && dj < 0
        % check top right corner
        grid_corner = gridcell_center + [dx dy]/2;
        isBLayer = blayer_thickness^2 > sum(([bac_x bac_y] - grid_corner).^2);
        
    elseif di > 0 && dj < 0
        % check top left corner
        grid_corner = gridcell_center + [-dx dy]/2;
        isBLayer = blayer_thickness^2 > sum(([bac_x bac_y] - grid_corner).^2);
        
    elseif di < 0 && dj > 0
        % check bottom right corner
        grid_corner = gridcell_center + [dx -dy]/2;
        isBLayer = blayer_thickness^2 > sum(([bac_x bac_y] - grid_corner).^2);
        
    else % di > 0 && dj > 0
        % check bottom left corner
        grid_corner = gridcell_center - [dx dy]/2;
        isBLayer = blayer_thickness^2 > sum(([bac_x bac_y] - grid_corner).^2);
        
    end
end

function isBLayer = isWithinBoundaryLayer_v2(bac_x, bac_y, gridcell_center, blayer_thickness)
    % Determine whether a gridcell is within the boundary layer thickness 
    % of a specific bacterium. Only look at whether the center of a
    % gridcell is in the boundary layer.
    %
    % bac_x, bac_y: bacterial coordinates
    % gridcell_center: center coordinates (x, y) of the gridcell
    % blayer_thickness: size of the boundary layer
    %
    % -> isBLayer: boolean with whether the checked gridcell is within the
    %   boundary layer.
    
    isBLayer = blayer_thickness^2 > sum(([bac_x bac_y] - gridcell_center).^2);
end

function plotDiffRegion(g, bac, diffRegion, checked_cell, active_cell)
    % debugging function to visualise the determination of the diffusion region
    x = bac.x;
    y = bac.y;
    
    diffRegion = double(diffRegion);
    diffRegion(checked_cell(1), checked_cell(2)) = 0.4;
    diffRegion(active_cell(1), active_cell(2)) = 0.6;
    diffRegion = diffRegion';
    
    figure(2); clf;
    
    % plot diff region
    im = imagesc([g.dx/2, g.nX*g.dx - g.dx/2], [g.dy/2, g.nY*g.dy - g.dy/2], diffRegion); hold on;
    im.AlphaData = 0.5;
    set(gca, 'YDir', 'normal');
    colormap([0 0 0; 1 0.1 0.1; .1 1 .1; 1 1 1]);
    active_bac = floor(x / g.dx) + 1 == active_cell(1) & floor(y / g.dy) + 1 == active_cell(2);
    c = ones(size(x)) * [0 .7 .7];
    c(active_bac, :) = ones(nnz(active_bac), 1) * [0 .55 .55];
    
    
    
    % plot bacs
    scatter(x, y, [], c, 'filled', 'MarkerEdgeColor', [0 .5 .5], 'LineWidth', 1.5); hold on;
    
    % plot circles around each bacterium with boundary layer thickness
    r = g.blayer_thickness;
    si = sum(active_bac);
    ii = find(active_bac);
    for j = 1:si
        i = ii(j);
        rectangle('Curvature', [1 1], 'Position', [x(i) - r, y(i) - r, 2 * r, 2 * r], 'LineWidth', 1, 'EdgeColor', [0, 0.5, 0.5 1]);
    end
    
    axis equal;

    xlim([0, g.nX*g.dx]);
    xticks(0:g.nX);
    ylim([0, g.nY*g.dy]);
    yticks(0:g.nY);
    grid on;
end
