function [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac)
    % Create a matrix (nX * nY * ?) with reference to which bacteria reside
    % in each grid-cell
    %
    % grid: struct containing all static information about the grid
    % bac: struct containing all information regarding the bacteria
    %
    % -> grid2bac: matrix (nx*ny) with per grid cell the bacterial indices
    % that that grid cell contains
    % -> grid2nBacs: matrix (nx*ny). Number of bacteria in each node
    % <E: Am I correct with grid2nBacs? />

    maxBacPerGrid = 4; % estimate of how many bacteria can be in one grid cell at most
    % <TODO: Could we actually estimate the maxBacPerGrid using a typical density
    % value of granules. maxBacPerGrid = f(dx, GranuleDensity, BacterialDensity...)/>
    % <E: Idk if it is necessary to apply a formula here, but I guess that
    % we should check this value at least one time. />
   
    %% determine grid cell per bacteria
    bac.x;
    bac.y;
    ix = floor(bac.x / grid.dx) + 1; % +1 because of MATLAB indexing
    iy = floor(bac.y / grid.dy) + 1;
    bac_grid = [ix, iy]; 
    
    
    %% invert to be which bacs per grid cell
    grid2nBacs = zeros(grid.nX, grid.nY, 'uint8'); % uint8 for reduced storage requirements (assumes no more than 255 bacs per grid cell)
    grid2bac = zeros(grid.nX, grid.nY, maxBacPerGrid, 'uint32'); % uint32 for reduced storage requirements -> lower storage requirements than sparse matrix (nBacs*nGridcells)
    
    for i = 1:size(bac_grid, 1)
        ix = bac_grid(i, 1);
        iy = bac_grid(i, 2);
        grid2bac(ix, iy, grid2nBacs(ix, iy) + 1) = i;
        grid2nBacs(ix, iy) = grid2nBacs(ix, iy) + 1;
    end
end