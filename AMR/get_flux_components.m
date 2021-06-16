function [x_ind, x_frac, y_ind, y_frac] = get_flux_components(ci, grid)
    %%% Function that returns the indices of cells and their respective
    %%% fractions of which the flux in cell ci are made up, using linear
    %%% interpolation (excludes the one time self concentration)
    % initialise arrays of indices and fractions for x and y direction
    x_ind = zeros(9, 1); % if all neighbours are one level higher, maximum of 9 indices (including 1 self)
    x_frac = zeros(9, 1);
    y_ind = zeros(9,1);
    y_frac = zeros(9,1);
    
    cur_index = 1;
    dx2 = (grid.dx / 2^(grid.cells.lvl(ci) + 1))^2;
    
    % find all neighbours of the given cell (ci)
    neighbours = zeros(8, 1); % store all orthogonal neighbours
    neighbourIndices = [2, 3; 1, 4; 1, 4; 2, 3];
    nb_dir = [2,5; 6,8; 4,7; 1,3]; % in the neighbour array: which neighbours correspond per direction (1=north, then clockwise)
%     sub_nb_dir = [2,4; 1,2; 1,3; 3,4]; % if neighbour is one level lower, which subconc is required per neighbour direction

    for dir = 1:4 % loop over all 4 directions for nodes
        ni = grid.cells.nod(ci, dir); % find node from this cell in given direction
        neighbours([1, 2] + (dir - 1) * 2) = grid.nodes.cel(ni, neighbourIndices(dir, :));
    end
    
    for dir = 2:2:4 % loop over each neighbour direction (east, west)
         % find neighbours using both nodes that touch north face
        % -> neighbour indices from nb_dir
        direct_neighbours = neighbours(nb_dir(dir,:));
        % if both are 0, then next to boundary of simulation, thus set
        % concentration to boundary concentration
        if all(direct_neighbours == 0)
            error('all neighbours are 0, to be implemented still');
            
        % level of neighbour is lower than level(ci)
        elseif any(direct_neighbours == 0)
            error('neighbour is level lower than current cell, to be implemented still');

        % level of neighbour is higher than level(ci)
        elseif direct_neighbours(1) ~= direct_neighbours(2)
            siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
            have_same_lvl = all(~grid.cells.children(siblings, 1));
            if ~have_same_lvl
                error('cannot simply take the average over neighbours parents children as these have different levels again');
            else
                x_ind(cur_index:cur_index+3) = siblings;
                x_frac(cur_index:cur_index+3) = 0.25;
                cur_index = cur_index + 4;
            end

        % level of neighbour is the same as level(ci)
        else
            x_ind(cur_index) = direct_neighbours(1);
            x_frac(cur_index) = 1;
            cur_index = cur_index + 1;
        end
    end
    
    % add ci fractions as well
    x_ind(cur_index) = ci;
    x_frac(cur_index) = -2;

    cur_index = 1;
    for dir = 1:2:3 % loop over each neighbour direction (north, south)
         % find neighbours using both nodes that touch north face
        % -> neighbour indices from nb_dir
        direct_neighbours = neighbours(nb_dir(dir,:));
        % if both are 0, then next to boundary of simulation, thus set
        % concentration to boundary concentration
        if all(direct_neighbours == 0)
            error('all neighbours are 0, to be implemented still');
            
        % level of neighbour is lower than level(ci)
        elseif any(direct_neighbours == 0)
            error('neighbour is level lower than current cell, to be implemented still');

        % level of neighbour is higher than level(ci)
        elseif direct_neighbours(1) ~= direct_neighbours(2)
            siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
            have_same_lvl = all(~grid.cells.children(siblings, 1));
            if ~have_same_lvl
                error('cannot simply take the average over neighbours parents children as these have different levels again');
            else
                y_ind(cur_index:cur_index+3) = siblings;
                y_frac(cur_index:cur_index+3) = 0.25;
                cur_index = cur_index + 4;
            end

        % level of neighbour is the same as level(ci)
        else
            y_ind(cur_index) = direct_neighbours(1);
            y_frac(cur_index) = 1;
            cur_index = cur_index + 1;
        end
    end
    
    % add ci fractions as well
    y_ind(cur_index) = ci;
    y_frac(cur_index) = -2;
    
    
    % correct for distance
    x_frac = x_frac / dx2;
    y_frac = y_frac / dx2;

end