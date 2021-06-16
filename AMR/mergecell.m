function [grid] = mergecell(grid, parent)
    % function to merge all children of a parent cell into one cell

    %----- mark cells for rey_centercling (index in grid.nodes arrays can be reused)
    grid.cells.nRec = grid.cells.nRec + 1;
    grid.cells.rec(grid.cells.nRec) = grid.cells.children(parent, 1); % mark first cell of children (implicit = 3 other children recycled too)

    %----- delete children reference and connected nodes
    for child_index = 1:4 % loop over all children

        for ni = 1:4
            % remove grid.nodes.cel reference to children
            if ni == child_index % retain references for the 4 nodes of parent cell
                grid.nodes.cel(grid.cells.nod(parent, ni), 5 - ni) = parent;
            else
                % remove reference to child cel
                grid.nodes.cel(grid.cells.nod(grid.cells.children(parent, child_index), ni), 5 - ni) = 0;

                % rey_centercle the node if there is no more references to any
                % cell left
                if all(grid.nodes.cel(grid.cells.nod(grid.cells.children(parent, child_index), ni), [1, 2, 3, 4]) == 0)
                    grid.nodes.nRec = grid.nodes.nRec + 1;
                    grid.nodes.rec(grid.nodes.nRec) = grid.cells.nod(grid.cells.children(parent, child_index), ni);
                end
            end
        end

        % set children value of removed children to -1, so that they are skipped in loops over all
        % cells
        grid.cells.children(grid.cells.children(parent, child_index), 1) = -1;
    end

    %----- update parent cell
    grid.cells.conc(parent) = mean(grid.cells.conc(grid.cells.children(parent,:)));
    grid.cells.children(parent, 1) = 0;

end