function [grid, nCells] = refine_grid(grid, nCells, nNodes, t)
    % refresher on data structure:
    % 4 nodes connected to each cell (ci = cell_index)
    %    |    |
    % -- 1 -- 3 --
    %    | ci |
    % -- 2 -- 4 --
    %    |    |

    % 4 Cells connected per node (ni = node_index)
    % 1  ||  3
    % --|ni|--
    % 2  ||  4

    toRefine = zeros(grid.maxCells, 1);
    toCoarsen = zeros(grid.maxCells, 1);

    minValueForRefinement = 0.5; % minimum concentration difference between grid cells for refinement
    maxValueForCoarsen = 0.1; % maximum concentration difference between grid cells for coarsening
    grid.maxLevel = 3;


    % neighbourIndices => conversion matrix/lookup table for neighbouring cells based on
    % direction of found node of specific cell. E.g. if top-left node of a cell is
    % checked to find neighbours, those are found in grid.nodes.cel(ni,
    % [2,3])
    % direction of ni -> neighbour index (non-diagonal)
    %               1 -> 2,3
    %               2 -> 1,4
    %               3 -> 1,4
    %               4 -> 2,3
    neighbourIndices = [2, 3; 1, 4; 1, 4; 2, 3];
    neighbours = zeros(grid.maxCells, 8); % store all orthogonal neighbours per cell
    
    iter = 0;
    while sum(toRefine) + sum(toCoarsen) > 0 || iter < 2
        %% mark cells for coarsening/refinement based on concentration gradient
        
        for ci = 1:nCells

            if grid.cells.children(ci, 1) || toRefine(ci) % if cell is parent or removed or already marked, then continue (only interested in leaf-cells of tree)
                continue
            end

            % check all cells around for concentration difference
            for dir = 1:4 % loop over all 4 directions for nodes
                ni = grid.cells.nod(ci, dir); % find node from this cell in given direction
                neighbours(ci, [1, 2] + (dir - 1) * 2) = grid.nodes.cel(ni, neighbourIndices(dir, :));
            end

            % what are the concentrations at neighbours
            nz_neighbours = nonzeros(neighbours(ci,:)); % non-zero neighbours
            neighbour_concentrations = grid.cells.conc(nz_neighbours);
            self_concentration = grid.cells.conc(ci);

            % determine if refinement or coarsening is required
%             val = abs(neighbour_concentrations - self_concentration) ./ (min(neighbour_concentrations, self_concentration) + 1e-5);
%             val = abs(neighbour_concentrations - self_concentration) ./ (self_concentration + 1e-5);
%             val = abs(neighbour_concentrations - self_concentration);
            val = (neighbour_concentrations) ./ (self_concentration - 1e-2);
%             d1 = grid.dx ./ 2.^(grid.cells.lvl(nonzeros(neighbours(ci,:))) + 1);
%             d2 = grid.dx ./ 2.^(grid.cells.lvl(nonzeros(neighbours(ci,:))) + 1);
%             val = 25 * abs(neighbour_concentrations - self_concentration) ./ (d1 + d2); % backward finite difference -> gradient
            
            if grid.cells.lvl(ci) < grid.maxLevel
                toRefine(ci) = any(val > minValueForRefinement);
                if toRefine(ci)
                    % check if neighbours need to be refined as well
                    % if level of neighbours < this cell's level - 1, then also refine
                    % those neighbours
                    mask = grid.cells.lvl(nz_neighbours) < grid.cells.lvl(ci);
                    toRefine(nz_neighbours(mask)) = 1;
                end
            end
                
            if iter == 0 && grid.cells.lvl(ci) % only check for coarsening if level > 0 and on first iteration
                toCoarsen(ci) = all(val < maxValueForCoarsen);
            end

        end

        %% visualise marked cells
        visualise_grid(grid, nCells, toCoarsen, toRefine);

        %% correct marked cells
        % remove twice marked
        % before coarsening: check level of neighbours (if not sufficient:
        %       cancel coarsening
        % before refining: check level of neighbours (if not sufficient: mark
        %       for refinement as well

        disp(find(toCoarsen(toRefine == 1)));

        %% refine cells
        % from lowest level to highest level
        cells_in_use = 1:nCells;
        
        for lvl = 0:grid.maxLevel
            % isolate all cells that are of this level
            cell_in_lvl = find(grid.cells.lvl(cells_in_use) == lvl);
            fprintf('%d cells in level %d\n', length(cell_in_lvl), lvl);
            
            
            for i = 1:length(cell_in_lvl)
                ci = cell_in_lvl(i);
                if ~toRefine(ci)
                    continue
                end

                % refine cell
                [grid, nCells, nNodes] = splitcell(grid, ci, nCells, nNodes, t);

                % remove mark for refined cell
                toRefine(ci) = 0;
            end
        end

        %% possible duplicate coarsening and refining?
        disp(find(toCoarsen(toRefine == 1)));

        %% coarsen cells
        if iter == 0
            for ci = 1:nCells

                if ~toCoarsen(ci)
                    continue
                end

                % check whether all children of parent are marked for coarsening
                parent = grid.cells.parent(ci);
                siblings = nonzeros(grid.cells.children(parent, :));

                if ~all(toCoarsen(siblings)) % if not all siblings are marked, remove those that are and continue
                    toCoarsen(siblings) = 0;
                    continue
                end

                % check level of neighbours and remove if level(ci) <
                % level(neighbours). If so, don't coarsen
                if any(grid.cells.lvl(ci) < grid.cells.lvl(nonzeros(neighbours(ci,:))))
                    toCoarsen(ci) = 0;
                    continue
                end

                % coarsen parent of cell
                [grid] = mergecell(grid, parent);

                % remove cell and siblings from marked cells
                toCoarsen(siblings) = 0;
            end
        end

        % reset coarsening marked cells
        % toCoarsen = zeros(grid.maxCells, 1);
        iter = iter + 1;
        disp(iter);
    end
end

function visualise_grid(grid, nCells, toCoarsen, toRefine)
    f = figure(4);
    clf;
    f.Position = [-469 57 425 400]; % fig1: [-1773 164 700 600]
    alpha = 0.5;

    % find dx and dy for base level
    dy = grid.dy;
    dx = grid.dx;

    % for all cells, draw rectangle with color based on whether
    % refinement/coarsening is required
    for ci = 1:nCells

        if ~grid.cells.children(ci,1)
            % find bottom left node
            ni = grid.cells.nod(ci, 1);
            x = grid.nodes.x(ni);
            y = grid.nodes.y(ni);

            if toRefine(ci)
                c = [1, 0, 0, alpha];
            elseif toCoarsen(ci)
                c = [0, 1, 0, alpha];
            else
                c = [0, 0, 0, alpha];
            end

            % draw rectangle for cell
            rectangle('Position', [x, y, dx / (2^grid.cells.lvl(ci)), dy / (2^grid.cells.lvl(ci))], 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', c);
        end

    end

    axis equal; axis off;
    xlim([grid.xmin, grid.xmax]);
    ylim([grid.ymin, grid.ymax]);
end
