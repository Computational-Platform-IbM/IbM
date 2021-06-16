function conc = ADI(grid, dt, nCells)
    % check if dx == dy
    if grid.dx ~= grid.dy || grid.nX ~= grid.nY
        error('dx and dy are not the same... or nX and nY are not the same');
    end
    conc = grid.cells.conc(1:nCells);
    
    D = 1; % diffusion coefficient
    
    % create B matrices
    alpha = dt * D / (2 * grid.dx^2);
        
    nx = grid.nX; % size of Lx matrix
    e = ones(nx, 1);
    
    B = spdiags([1*e -2*e 1*e], -1:1, nx, nx);
    bc = sparse([1, nx], [nx, 1], 1, nx, nx);
    B = B + bc;
    B_orthogonal = kron(B, speye(nx));
    
    % keep y-direction fixed -> concentration linearized columnwise
    temp = zeros(nCells, 1);
    I = speye(nx);
    I_orthogonal = speye(nx^2);
    for ix = 1:nx
        ind = (ix-1)*nx+1:(ix*nx);
        rhs = (I_orthogonal(ind,:) + alpha*B_orthogonal(ind,:))*conc;
        temp(ind) = (I - alpha*B) \ rhs;
    end
        
    % keep x-direction fixed -> remake the concentration vector column-wise
    temp = reshape(reshape(temp, nx, nx)', nx^2, 1); % recast into row-wise
    for iy = 1:nx
        ind = (iy-1)*nx+1:(iy*nx);
        rhs = (I_orthogonal(ind,:) + alpha*B_orthogonal(ind,:))*temp;
        conc(ind) = (I - alpha*B) \ rhs;
    end
    
    conc = reshape(reshape(conc, nx, nx)', nx^2, 1); % recast back into column-wise vector

end