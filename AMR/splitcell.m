function [grid, nCells, nNodes] = splitcell(grid, parent, nCells, nNodes, t)
    % function to split a given cell into 4 children
    lvl = grid.cells.lvl(parent) + 1;
    x_center = (grid.nodes.x(grid.cells.nod(parent, 1)) + grid.nodes.x(grid.cells.nod(parent, 4))) / 2;
    y_center = (grid.nodes.y(grid.cells.nod(parent, 1)) + grid.nodes.y(grid.cells.nod(parent, 4))) / 2;

    %----- check if there are cells recycled and remember old nCells
    cells_recycled = grid.cells.nRec > 0;

    if cells_recycled
        cell_index = grid.cells.rec(grid.cells.nRec) - 1;
        grid.cells.nRec = grid.cells.nRec - 1;
    else
        cell_index = nCells;
    end

    %----- create new children cells
    % +------+
    % | 1  3 |
    % | 2  4 |
    % +------+
    i = [1, 2, 3, 4];

    % reference to children cells
    grid.cells.children(parent, :) = cell_index + i;
    % level of children
    grid.cells.lvl(cell_index + i) = lvl;
    % set no children for children cells
    grid.cells.children(cell_index + i, 1) = zeros(4, 1);
    % set parent of children
    grid.cells.parent(cell_index + i) = parent;
    
    % set concentration in children
    grid.cells.conc(cell_index + i) = get_concentrations_for_children(parent, grid, t);
    

    %----- create new nodes
    % central node
    if ~grid.nodes.nRec
        nNodes = nNodes + 1;
        node_index = nNodes;
    else
        node_index = grid.nodes.rec(grid.nodes.nRec);
        grid.nodes.nRec = grid.nodes.nRec - 1;
    end

    grid.nodes.x(node_index) = x_center; % x coordinate
    grid.nodes.y(node_index) = y_center; % y coordinate

    % Add  parent and central nodes
    for i = 1:4 % looping over children/directions/nodes
        % change references to parent node to reference child node
        grid.cells.nod(cell_index + i, i) = grid.cells.nod(parent, i);
        grid.nodes.cel(grid.cells.nod(cell_index + i, i), 5 - i) = cell_index + i;

        % create references to/from the central node
        grid.cells.nod(cell_index + i, 5 - i) = node_index;
        grid.nodes.cel(node_index, i) = cell_index + i;
    end

    % other nodes ----> TODO: optimize to not do it so manually.... OOP?
    % Cell 1--------------------------------------
    % 4 nodes connected to the cell
    % x    3
    %   ci
    % 2    x
    % Node 2

    % check if node already exists
    if grid.nodes.cel(grid.cells.nod(parent, 1), 2) > 0 && ...
            grid.cells.lvl(grid.nodes.cel(grid.cells.nod(parent, 1), 2)) == lvl
        grid.cells.nod(cell_index + 1, 2) = grid.cells.nod(grid.nodes.cel(grid.cells.nod(parent, 1), 2), 4);
    else

        if ~grid.nodes.nRec
            nNodes = nNodes + 1;
            node_index = nNodes;
        else
            node_index = grid.nodes.rec(grid.nodes.nRec);
            grid.nodes.nRec = grid.nodes.nRec - 1;
        end

        grid.nodes.x(node_index) = grid.nodes.x(grid.cells.nod(parent, 1)); % x coordinate
        grid.nodes.y(node_index) = y_center; % y coordinate
        grid.cells.nod(cell_index + 1, 2) = node_index;
    end

    grid.nodes.cel(grid.cells.nod(cell_index + 1, 2), 3) = cell_index + 1;
    % Node 3
    if (grid.nodes.cel(grid.cells.nod(parent, 1), 3) > 0 && grid.cells.lvl(grid.nodes.cel(grid.cells.nod(parent, 1), 3)) == lvl)
        grid.cells.nod(cell_index + 1, 3) = grid.cells.nod(grid.nodes.cel(grid.cells.nod(parent, 1), 3), 4);
    else

        if ~grid.nodes.nRec
            nNodes = nNodes + 1;
            node_index = nNodes;
        else
            node_index = grid.nodes.rec(grid.nodes.nRec);
            grid.nodes.nRec = grid.nodes.nRec - 1;
        end

        grid.nodes.x(node_index) = x_center; % x coordinate
        grid.nodes.y(node_index) = grid.nodes.y(grid.cells.nod(parent, 1)); % y coordinate
        grid.cells.nod(cell_index + 1, 3) = node_index;
    end

    grid.nodes.cel(grid.cells.nod(cell_index + 1, 3), 2) = cell_index + 1;

    % Cell 2--------------------------------------
    % 4 nodes connected to the cell
    % 1    x
    %   ci
    % x    4
    % Node 1
    grid.cells.nod(cell_index + 2, 1) = grid.cells.nod(cell_index + 1, 2);
    grid.nodes.cel(grid.cells.nod(cell_index + 2, 1), 4) = cell_index + 2;
    % Node 4
    if (grid.nodes.cel(grid.cells.nod(parent, 2), 4) > 0 && grid.cells.lvl(grid.nodes.cel(grid.cells.nod(parent, 2), 4)) == lvl)
        grid.cells.nod(cell_index + 2, 4) = grid.cells.nod(grid.nodes.cel(grid.cells.nod(parent, 2), 4), 3);
    else

        if ~grid.nodes.nRec
            nNodes = nNodes + 1;
            node_index = nNodes;
        else
            node_index = grid.nodes.rec(grid.nodes.nRec);
            grid.nodes.nRec = grid.nodes.nRec - 1;
        end

        grid.nodes.x(node_index) = x_center; % x coordinate
        grid.nodes.y(node_index) = grid.nodes.y(grid.cells.nod(parent, 2)); % y coordinate
        grid.cells.nod(cell_index + 2, 4) = node_index;
    end

    grid.nodes.cel(grid.cells.nod(cell_index + 2, 4), 1) = cell_index + 2;

    % Cell 3--------------------------------------
    % 4 nodes connected to the cell
    % 1    x
    %   ci
    % x    4
    % Node 1
    grid.cells.nod(cell_index + 3, 1) = grid.cells.nod(cell_index + 1, 3);
    grid.nodes.cel(grid.cells.nod(cell_index + 3, 1), 4) = cell_index + 3;
    % Node 4
    if (grid.nodes.cel(grid.cells.nod(parent, 3), 4) > 0 && grid.cells.lvl(grid.nodes.cel(grid.cells.nod(parent, 3), 4)) == lvl)
        grid.cells.nod(cell_index + 3, 4) = grid.cells.nod(grid.nodes.cel(grid.cells.nod(parent, 3), 4), 2);
    else

        if ~grid.nodes.nRec
            nNodes = nNodes + 1;
            node_index = nNodes;
        else
            node_index = grid.nodes.rec(grid.nodes.nRec);
            grid.nodes.nRec = grid.nodes.nRec - 1;
        end

        grid.nodes.x(node_index) = grid.nodes.x(grid.cells.nod(parent, 3)); % x coordinate
        grid.nodes.y(node_index) = y_center; % y coordinate
        grid.cells.nod(cell_index + 3, 4) = node_index;
    end

    grid.nodes.cel(grid.cells.nod(cell_index + 3, 4), 1) = cell_index + 3;

    % Cell 4--------------------------------------
    % 4 nodes connected to the cell
    % x    3
    %   ci
    % 2    x
    % Node 2
    grid.cells.nod(cell_index + 4, 2) = grid.cells.nod(cell_index + 2, 4);
    grid.nodes.cel(grid.cells.nod(cell_index + 4, 2), 3) = cell_index + 4;
    % Node 3
    grid.cells.nod(cell_index + 4, 3) = grid.cells.nod(cell_index + 3, 4);
    grid.nodes.cel(grid.cells.nod(cell_index + 4, 3), 2) = cell_index + 4;

    % Add number of cells
    if ~cells_recycled
        nCells = nCells + 4;
    end

end