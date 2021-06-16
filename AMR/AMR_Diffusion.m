function conc_new = AMR_Diffusion(grid, nCells, dt)

    neighbourIndices = [2, 3; 1, 4; 1, 4; 2, 3]; % indices of directions from a node per direction for orthogonal neighbours
    nb_dir = [2,5; 6,8; 4,7; 1,3]; % in the neighbour array: which neighbours correspond per direction (1=north, then clockwise)
    sub_nb_dir = [2,4; 1,2; 1,3; 3,4]; % if neighbour is one level lower, which subconc is required per neighbour direction
    sub_cell_directions = [-1,-1;-1,1;1,-1;1,1]; % for each of the sub_directions, which direction is the center located
    
    cells_in_use = 1:nCells;
%     Lxy = sparse(nCells, nCells); % create Lxy
    ghost_cell_conc = zeros(4*nCells, 1); % create array to store ghost cells, assume that maximally 1 ghost cell per normal cell
    ghost_ref = zeros(nCells, 1); % reference of parent cell to first ghost cell in the ghost array
    gi = 1; % ghost index
    
    lxy_constr = zeros(5*nCells, 3); % array for constructing sparse Lxy matrix: 1=x_index, 2=y_index, 3=value
    index = 1;

%     Lxy = cell(grid.maxLevel+1, 1); % cell array to store all constructed Laplacian matrices
%     Lxy_constr = cell(grid.maxLevel+1, 1); % cell array to store all construction matrices for the Laplacian matrices
%     cells_per_level = cell(grid.maxLevel+1, 1); % save which cells per level

    
    
    % go over each level and create Lxy
    for lvl = grid.maxLevel:-1:0
        % find cells in this level
        cell_in_lvl = find(grid.cells.lvl(cells_in_use) == lvl & ~grid.cells.children(cells_in_use,1));
        fprintf('%d cells in level %d\n', length(cell_in_lvl), lvl);
        
%         lxy_constr = zeros(17*length(cell_in_lvl), 3); % array for constructing sparse Lxy matrix: 1=x_index, 2=y_index, 3=value
%         index = 1;
        for cii = 1:length(cell_in_lvl)
            ci = cell_in_lvl(cii);
            self_val = 0;
            
            neighbours = zeros(8,1);
            % collect all neighbours
            for dir = 1:4 % loop over all 4 directions for nodes
                ni = grid.cells.nod(ci, dir); % find node from this cell in given direction
                neighbours([1, 2] + (dir - 1) * 2) = grid.nodes.cel(ni, neighbourIndices(dir, :));            
            end
                
            for dir = 1:4 % north, east, south, west
                direct_neighbours = neighbours(nb_dir(dir,:));
                if all(direct_neighbours == 0)
                    % if both neighbours are 0, then either continue?
                    % (neumann boundary) or find neighbours at other side
                    % of the simulation (periodic boundary)
                    if ci <= grid.nX^2 % cell is in level 0
                        % add periodic boundary condition
                        
                        % find row and column of ci
                        [row, column] = ind2sub([grid.nX, grid.nY], ci);
                        if dir == 1 % looking for north neighbour, loop around to most south cell in same column
                            nb = sub2ind([grid.nX, grid.nY], grid.nY, column);
                        elseif dir == 2 % looking for east neighbour, loop around to most west cell in same row
                            nb = sub2ind([grid.nX, grid.nY], row, 1);
                        elseif dir == 3 % looking for south neighbour, loop around to most north cell in same column
                            nb = sub2ind([grid.nX, grid.nY], 1, column);
                        else % looking for west neighbour, loop arond to most east cell in same row
                            nb = sub2ind([grid.nX, grid.nY], row, grid.nX);
                        end
                        
                        if grid.cells.children(nb, 1)
                            error('Periodic cell at the other side has a different level, to be implemented')
                        end
                        
                        self_val = self_val - 1;
                        lxy_constr(index, :) = [ci, nb, 1];
                        index = index + 1;

                    else
                        error('both neighbours are not found, add periodic boundary if needed?')
                    end
                    
                elseif any(direct_neighbours == 0)
                    % if one of the neighbours is zero, then neighbour is
                    % one level lower, thus create ghost cells of non-zero
                    % neighbour. nb = [nCells+gi+sub_nb_dir-1];
                    
                    % which of the ghost subconcentration is required
%                     sub_index = sub_nb_dir(dir, direct_neighbours ~= 0); 
                    offset = sub_nb_dir(dir, direct_neighbours ~= 0) - 1; % 1 subtracted, because first subconc starts at gi, not gi+1
                    
                    % which cell is parent
                    parent = direct_neighbours(direct_neighbours ~= 0);
                    
                    % check if ghost already made of this cell
                    % if not made, make ghost cells
                    if ~ghost_ref(parent)
                        ghost_cell_conc(gi+[0,1,2,3]) = get_concentrations_for_children(parent, grid, 1);
                        ghost_ref(parent) = gi;
                        nb = nCells + gi + offset;
                        
                        % for the ghost cell, add linear interpolation
                        % formula (based on the flux at the center of the
                        % parent cell):
                        
                        %    sub_conc = parent_conc + dFluxes * dxy
                        [x_ind, x_frac, y_ind, y_frac] = get_flux_components(parent, grid);
                        dx = grid.dx / 2^(lvl + 1); dy = dx;

                        for jj = 1:nnz(x_ind) % add x-flux components * dx
                            val = x_frac(jj) * sub_cell_directions(:, 1) * dx;
                            lxy_constr(index:index+3, :) = [nCells + gi+[0,1,2,3]', x_ind(jj)*ones(4,1), val];
                            index = index + 4;
                        end

                        for jj = 1:nnz(y_ind) % add y-flux components * dx
                            val = y_frac(jj) * sub_cell_directions(:, 2) * dy;
                            lxy_constr(index:index+3, :) = [nCells + gi+[0,1,2,3]', y_ind(jj)*ones(4,1), val];
                            index = index + 4;
                        end

                        % add concentration of parent
                        lxy_constr(index:index+3, :) = [nCells + gi+[0,1,2,3]', parent*ones(4,1), ones(4,1)];
                        index = index + 4;

                        gi = gi + 4;

                    else
                        nb = nCells + ghost_ref(parent) + offset;
                    end
                    
                    lxy_constr(index, :) = [ci, nb, 1];
                    index = index + 1;

                    self_val = self_val - 1;


                elseif direct_neighbours(1) ~= direct_neighbours(2)
                    % if neighbours are different, then neighbour is one level
                    % higher, thus take the average (i.e. 1/4 per sibling of
                    % neighbour)
                    siblings = grid.cells.children(grid.cells.parent(direct_neighbours(1)),:);
                    
                    self_val = self_val - 1;
                    for i = 1:4
                        lxy_constr(index, :) = [ci, siblings(i), 0.25];
                        index = index + 1;
                    end

%                     % create ghost cell with mean concentration of siblings
%                     ghost_cell_conc(gi) = mean(grid.cells.conc(siblings));
%                     nb = nCells + gi;
%                     gi = gi + 1;
% 
%                     lxy_constr(index, :) = [ci, nb, 1];
%                     index = index + 1;
                    
                    
                    
                else
                    % if both neighbours are the same, then neighbour is of the
                    % same level, thus add 1 to neighbour, -1 from self
                    
                    self_val = self_val - 1;
                    lxy_constr(index, :) = [ci, direct_neighbours(1), 1];
                    index = index + 1;                
                end

            end
            
            if self_val ~= -4
                disp(self_val);
                error('something went wrong?')
            end
            lxy_constr(index, :) = [ci, ci, self_val];
            index = index + 1;
            
        end
        
%         % create Lxy for this level if there where any cells
%         if ~isempty(cell_in_lvl)
%             Lxy_constr{lvl+1} = lxy_constr;
%             cells_per_level{lvl+1} = cell_in_lvl;
%         end

    end % end init
    
    
%     % create Laplacian matrices for each level
%     for i = 1:grid.maxLevel
%         % how many elements are in the constructor?
%         nEl = nnz(Lxy_constr{i}(:,1));
%         Lxy{i} = sparse(Lxy_constr{i}(1:nEl,1), Lxy_constr{i}(1:nEl,2), Lxy_constr{i}(1:nEl,3), nCells+nnz(ghost_cell_conc), nCells+nnz(ghost_cell_conc));
%         figure(6); clf;
%         spy(Lxy{i});
%     end

    nGhost = gi - 1;
    nEl = index - 1;
    
    % create matrix for use in: 
    % (1 - alpha*Lxy)*phi_new = (1 + alpha*Lxy)*phi_old
    Lxy = sparse(lxy_constr(1:nEl,1), lxy_constr(1:nEl,2), lxy_constr(1:nEl,3), nCells+nGhost, nCells+nGhost);
    figure(6);
    clf;
    spy(Lxy);
    
    
    
    
    

    % diffuse
    conc_new = diffuse(Lxy, grid, ghost_cell_conc, nCells, nGhost, dt, cells_per_level);
    
    % update ghost cell concentrations
    
    
    % check solution
    
   
    
    

end

function conc_new = diffuse(Lxy, grid, ghost_cell_conc, nCells, nGhost, dt, cells_per_level)
    conc_new = zeros(nCells, 1);
    
    % solve diffusion
    % solve diffusion itself (naive method)
    D = 1;
    
    temp_conc = [grid.cells.conc(1:nCells); ghost_cell_conc(1:nGhost)];
    I = speye(nCells + nGhost);
    
    %             alpha = dt * D / (2 * grid.dx^2);
    for lvl = 0:length(Lxy)-1
        if ~isempty(Lxy{lvl+1})
            alpha = (dt * D * 2^lvl) / (grid.dx^2);

            rhs = (I + alpha * Lxy{lvl+1})*temp_conc;
            temp_conc = (I - alpha * Lxy{lvl+1}) \ rhs;
            conc_new(cells_per_level{lvl+1}) = temp_conc(cells_per_level{lvl+1});
        end
    end
end
