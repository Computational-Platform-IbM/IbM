function concentrations = get_concentrations_for_children(ci, grid, t)
    %%% returns the concentrations in the 4 child nodes of a grid cell
    % +--------+
    % | 1    3 |
    % |   ci   |
    % | 2    4 |
    % +--------+
    
    % get difference of new x, y positions of centers w.r.t. current center
    dxy = [-1,-1;-1,1;1,-1;1,1]'.*[grid.dx; grid.dy]./(2^(grid.cells.lvl(ci)+2)); % [dx; dy]
    
    if t ~= 0
        dFluxes = slope_limited_gradient(ci, grid);
        concentrations = grid.cells.conc(ci) + dFluxes * dxy;
    else
        center_x = (grid.nodes.x(grid.cells.nod(ci, 1)) + grid.nodes.x(grid.cells.nod(ci, 4)))/2;
        center_y = (grid.nodes.y(grid.cells.nod(ci, 1)) + grid.nodes.y(grid.cells.nod(ci, 4)))/2;
        new_centers = [center_x; center_y] + dxy;
        concentrations = GaussianDiffusionExact(new_centers(1,:), new_centers(2,:), 0);
    end
end

function F = slope_limited_gradient(ci, grid)
    % determine the gradient with slope limiter
    % +-- 1--+    fluxes determined at the cell boundary
    % |      |    returns the decomposed gradient at the center of the cell:
    % 4  ci  2           
    % |      |       F = [Fx, Fy]
    % +-- 3--+           
    if grid.dx ~= grid.dy
        error('dx and dy are not the same')
    end
    
    dx = grid.dx / 2^(grid.cells.lvl(ci) + 1);
    self_conc = grid.cells.conc(ci);
    
    % find all neighbours of the given cell (ci)
    neighbours = zeros(8, 1); % store all orthogonal neighbours
    neighbourIndices = [2, 3; 1, 4; 1, 4; 2, 3];
    nb_dir = [2,5; 6,8; 4,7; 1,3]; % in the neighbour array: which neighbours correspond per direction (1=north, then clockwise)
    sub_nb_dir = [2,4; 1,2; 1,3; 3,4]; % if neighbour is one level lower, which subconc is required per neighbour direction

    for dir = 1:4 % loop over all 4 directions for nodes
        ni = grid.cells.nod(ci, dir); % find node from this cell in given direction
        neighbours([1, 2] + (dir - 1) * 2) = grid.nodes.cel(ni, neighbourIndices(dir, :));
    end
    
    flux = zeros(4,1);
    for dir = 1:4 % loop over each neighbour direction (north, east, south, west)
         % find neighbours using both nodes that touch north face
        % -> neighbour indices from nb_dir
        direct_neighbours = neighbours(nb_dir(dir,:));
        % if both are 0, then next to boundary of simulation, thus set
        % concentration to boundary concentration
        if all(direct_neighbours == 0)
            conc_neighbour = 1; % TODO: this is only through if neumann boundary with gamma=1;

        % level of neighbour is lower than level(ci)
        elseif any(direct_neighbours == 0)
            temp = get_concentrations_for_children(direct_neighbours(direct_neighbours ~= 0), grid, t);
%             if direct_neighbours(1) ~= 0 % neighbour 2 ~= 0 -> get subconcentration of 2
            conc_neighbour = temp(sub_nb_dir(dir, direct_neighbours ~= 0));
%             else % neighbour 5 ~= 0 -> get subconcentration of 5
%                 conc_neighbour = temp(sub_nb_dir(dir, 2));
%             end

        % level of neighbour is higher than level(ci)
        elseif direct_neighbours(1) ~= direct_neighbours(2)
            siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
            have_same_lvl = all(~grid.cells.children(siblings, 1));
            if ~have_same_lvl
                error('cannot simply take the average over neighbours parents children as these have different levels again');
            else
                conc_neighbour = mean(grid.cells.conc(siblings));
            end

        % level of neighbour is the same as level(ci)
        else
            conc_neighbour = grid.cells.conc(direct_neighbours(1));
        end
        
        if dir == 1 || dir == 4
            flux(dir) = (self_conc - conc_neighbour)/dx;
        else
            flux(dir) = (conc_neighbour - self_conc)/dx;
        end
    end
    
%     %----- flux 1 (north)
%     % find neighbours using both nodes that touch north face
%     % -> neighbour indices 2, 5
%     direct_neighbours = neighbours([2,5]);
%     % if both are 0, then next to boundary of simulation, thus set
%     % concentration to boundary concentration
%     if all(direct_neighbours == 0)
%         conc_neighbour = 1; % TODO: this is only through if neumann boundary with gamma=1;
%         
%     % level of neighbour is lower than level(ci)
%     elseif any(direct_neighbours == 0)
% %         error('neighbour is lower level')
%         temp = get_concentrations_for_children(direct_neighbours(direct_neighbours ~= 0), grid, t);
%         if direct_neighbours(1) ~= 0 % neighbour 2 ~= 0 -> get subconcentration of 2
%             conc_neighbour = temp(2);
%         else % neighbour 5 ~= 0 -> get subconcentration of 5
%             conc_neighbour = temp(4);
%         end
%         
%     % level of neighbour is higher than level(ci)
%     elseif direct_neighbours(1) ~= direct_neighbours(2)
%         siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
%         have_same_lvl = all(grid.cells.lvl(siblings) == grid.cells.lvl(direct_neighbours(1)));
%         if ~have_same_lvl
%             error('cannot simply take the average over neighbours parents children as these have different levels again');
%         else
%             conc_neighbour = mean(grid.cells.conc(siblings));
%         end
%     
%     % level of neighbour is the same as level(ci)
%     else
%         conc_neighbour = grid.cells.conc(direct_neighbours(1));
%     end
%     flux1 = (self_conc - conc_neighbour)/dy;
%     
%     
%     %----- flux 2 (east)
%     % find neighbours using both nodes that touch north face
%     % -> neighbour indices 6, 8
%     direct_neighbours = neighbours([6,8]);
%     % if both are 0, then next to boundary of simulation, thus set
%     % concentration to boundary concentration
%     if all(direct_neighbours == 0)
%         conc_neighbour = 1;
%         
%     % level of neighbour is lower than level(ci)
%     elseif any(direct_neighbours == 0)
% %         error('neighbour is lower level')
% 
%         temp = get_concentrations_for_children(direct_neighbours(direct_neighbours ~= 0), grid, t);
%         if direct_neighbours(1) ~= 0 % neighbour 6 ~= 0 -> get subconcentration of 6
%             conc_neighbour = temp(1);
%         else % neighbour 8 ~= 0 -> get subconcentration of 8
%             conc_neighbour = temp(2);
%         end
%         
%     % level of neighbour is higher than level(ci)
%     elseif direct_neighbours(1) ~= direct_neighbours(2)
%         siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
%         have_same_lvl = all(grid.cells.lvl(siblings) == grid.cells.lvl(direct_neighbours(1)));
%         if ~have_same_lvl
%             error('cannot simply take the average over neighbours parents children as these have different levels again');
%         else
%             conc_neighbour = mean(grid.cells.conc(siblings));
%         end
%     
%     % level of neighbour is the same as level(ci)
%     else
%         conc_neighbour = grid.cells.conc(direct_neighbours(1));
%     end
%     flux2 = (conc_neighbour - self_conc)/dx;
%     
%     
%     %----- flux 3 (south)
%     % find neighbours using both nodes that touch north face
%     % -> neighbour indices 4, 7
%     direct_neighbours = neighbours([4,7]);
%     % if both are 0, then next to boundary of simulation, thus set
%     % concentration to boundary concentration
%     if all(direct_neighbours == 0)
%         conc_neighbour = 1;
%         
%     % level of neighbour is lower than level(ci)
%     elseif any(direct_neighbours == 0)
% %         error('neighbour is lower level')
%         temp = get_concentrations_for_children(direct_neighbours(direct_neighbours ~= 0), grid, t);
%         if direct_neighbours(1) ~= 0 % neighbour 4 ~= 0 -> get subconcentration of 4
%             conc_neighbour = temp(1);
%         else % neighbour 7 ~= 0 -> get subconcentration of 7
%             conc_neighbour = temp(3);
%         end
%         
%     % level of neighbour is higher than level(ci)
%     elseif direct_neighbours(1) ~= direct_neighbours(2)
%         siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
%         have_same_lvl = all(grid.cells.lvl(siblings) == grid.cells.lvl(direct_neighbours(1)));
%         if ~have_same_lvl
%             error('cannot simply take the average over neighbours parents children as these have different levels again');
%         else
%             conc_neighbour = mean(grid.cells.conc(siblings));
%         end
%     
%     % level of neighbour is the same as level(ci)
%     else
%         conc_neighbour = grid.cells.conc(direct_neighbours(1));
%     end
%     flux3 = (conc_neighbour - self_conc)/dy;
%     
%     
%     %----- flux 4 (west)
%     % find neighbours using both nodes that touch north face
%     % -> neighbour indices 1, 3
%     direct_neighbours = neighbours([1,3]);
%     % if both are 0, then next to boundary of simulation, thus set
%     % concentration to boundary concentration
%     if all(direct_neighbours == 0)
%         conc_neighbour = 1;
%         
%     % level of neighbour is lower than level(ci)
%     elseif any(direct_neighbours == 0) 
% %         error('neighbour is lower level')
%         temp = get_concentrations_for_children(direct_neighbours(direct_neighbours ~= 0), grid, t);
%         if direct_neighbours(1) ~= 0 % neighbour 1 ~= 0 -> get subconcentration of 1
%             conc_neighbour = temp(3);
%         else % neighbour 3 ~= 0 -> get subconcentration of 3
%             conc_neighbour = temp(4);
%         end
%         
%     % level of neighbour is higher than level(ci)
%     elseif direct_neighbours(1) ~= direct_neighbours(2)
%         siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
%         have_same_lvl = all(grid.cells.lvl(siblings) == grid.cells.lvl(direct_neighbours(1)));
%         if ~have_same_lvl
%             error('cannot simply take the average over neighbours parents children as these have different levels again');
%         else
%             conc_neighbour = mean(grid.cells.conc(siblings));
%         end
%     
%     % level of neighbour is the same as level(ci)
%     else
%         conc_neighbour = grid.cells.conc(direct_neighbours(1));
%     end
%     flux4 = (self_conc - conc_neighbour)/dx;   
    
    
    Fx = minmod(flux(2), flux(4));
    Fy = minmod(flux(1), flux(3));
%     Fx = mean([flux(2), flux(4)]);
%     Fy = mean([flux(1), flux(3)]);
    F = [Fx, Fy];
end

function m = minmod(a, b)
    m = (sign(a)+sign(b))/2.*min(abs(a), abs(b));
end