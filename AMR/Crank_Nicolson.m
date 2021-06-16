function sol = Crank_Nicolson(grid, dt, nCells)
    conc = grid.cells.conc(1:nCells);
    % check if dx == dy
    if grid.dx ~= grid.dy || grid.nX ~= grid.nY
        error('dx and dy are not the same... or nX and nY are not the same');
    end
    
    D = 1; % diffusion coefficient
    
    % create Lxy matrix
    alpha = dt * D / (2* grid.dx^2);
    
    nx = grid.nX; % size of Lx matrix
    e = ones(nx, 1); % diagonal unit
    
    L = spdiags([e -2*e e], -1:1, nx, nx);
    bc = sparse([1, nx], [nx, 1], 1, nx, nx);
    L = L + bc;
    
    I = speye(nx);
    Lxy = kron(I, L) + kron(L, I);
    
    rhs = (speye(nx^2) + alpha*Lxy)*conc;
    
    sol = (speye(nx^2) - alpha*Lxy) \ rhs;

end