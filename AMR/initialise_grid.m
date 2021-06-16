function [grid, nCells, nNodes] = initialise_grid(grid)
    % every grid vector is linearised as per columnbased for level 0, afterwards assigned by order of encountering:
    %
    %   1 | 4 | 7
    %  ---+---+---
    %   2 | 5 | 8
    %  ---+---+---
    %   3 | 6 | 9

    % nX & nY ==> number of CELLS per direction (nodes = grid.nX + 1)
    grid.dx = (grid.xmax - grid.xmin) / (grid.nX); 
    grid.dy = (grid.ymax - grid.ymin) / (grid.nY);

    % Coordinate vectors
    x=grid.xmin:grid.dx:grid.xmax;
    y=grid.ymin:grid.dy:grid.ymax;

    % Coordinates of cell-centers
    xc=grid.xmin+grid.dx/2:grid.dx:grid.xmax-grid.dx/2;
    yc=grid.ymin+grid.dy/2:grid.dy:grid.ymax-grid.dy/2;

    % Max amount of nodes
    grid.maxNodes = (grid.nY+1)*(grid.nX+1)*100; % maximum based on maxLevel: (2^(m+1) - 1)^2 * (grid.nX + 1) * (grid.nY + 1);
    % Max amount of cells
    grid.maxCells = grid.nY*grid.nX*100; % maximum based on maxLevel: (2^(m+1) - 1)^2 * grid.nX * grid.nY;

    % Nodal arrays
    grid.nodes.x=zeros(grid.maxNodes,1); % x coordinates of nodes
    grid.nodes.y=zeros(grid.maxNodes,1); % y coordinates of nodes
    grid.nodes.cel=zeros(grid.maxNodes,4); % 4 cells connected to each node
    grid.nodes.rec=zeros(grid.maxNodes,1); % listing of recycled nodes
    grid.nodes.nRec=0; % Number of nodes for recycling

    % Cell arrays
    grid.cells.nod=zeros(grid.maxCells,4); % 4 nodes surrounding each cell
    grid.cells.lvl=zeros(grid.maxCells,1); % resolution level for each cell
    grid.cells.parent=zeros(grid.maxCells,1); % parent cell
    grid.cells.children=zeros(grid.maxCells,4); % 4 daughter cells
    grid.cells.rec=zeros(grid.maxCells,1); % listing of recycled cells
    grid.cells.nRec=0; % Number of cells for recycling
    grid.cells.conc=zeros(grid.maxCells,1);

    %%%%% Level 0 init
    % Initial number of nodes and cells
    nNodes=(grid.nX + 1) * (grid.nY + 1);
    nCells=grid.nX*grid.nY;

    % Node coordinates
    grid.nodes.x(1:nNodes) = kron(x, ones(1, grid.nY+1))';
    grid.nodes.y(1:nNodes) = repmat(y, 1, grid.nX+1)';

    % Node-Cell connections
    % 4 nodes connected to each cell (ci = cell_index)
    % 1 -- 3
    % | ci |
    % 2 -- 4
    basis = (1:grid.nY)' + [0, 1, grid.nY+1, grid.nY+2]; % basic element
    grid.cells.nod(1:nCells, :) = kron((0:grid.nX-1)', ones(grid.nY, 4)) * (grid.nY+1) + repmat(basis, grid.nX, 1);

    % 4 Cells connected per node (ni = node_index)
    % 1  ||  3
    % --|ni|--
    % 2  ||  4
    basis = (0:grid.nY)';
    basis2 = repmat(basis, grid.nX, 1) + kron(0:grid.nX-1, [0, ones(1, grid.nY)])' * (grid.nY);
    grid.nodes.cel(grid.nY+2:nNodes, 1) = basis2;
    grid.nodes.cel(grid.nY+1:nNodes-1, 2) = basis2;
    grid.nodes.cel(1:nNodes-grid.nY-1, 3) = basis2;
    grid.nodes.cel(1:nNodes-grid.nY-2, 4) = basis2(2:end);

    % set initial concentration in grid cells
    [X,Y] = meshgrid(xc, yc);
    grid.cells.conc(1:nCells) = reshape(GaussianDiffusionExact(X, Y, 0), nCells, 1);

end