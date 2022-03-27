function T_new = recalculateT(T, i, j, kDet, grid, Visited, xcenter, ycenter)
    % Recalculate the T value at gridcell (i, j) using a quadratic 
    % approximation of the gradient.
    % 
    % T: matrix (nX*nY) with the current times of detachment
    % i, j: gridcell indices for the gridcell that needs a recalculation of
    %   the T value
    % kDet: Detachment constant determining how fast detachment takes place
    % grid: struct containing all information regarding spacial
    %   discretization
    % Visited: logical matrix (nX*nY) with per gridcell whether the correct
    %   T value has already been calculated
    % xcenter, ycenter: x and y coordinates of the center of the granule
    %
    % -> T_new: new time of detachment (T) for the respective gridcell
    
    right_visited = Visited(i+1, j); % x + 1 (right)
    left_visited = Visited(i-1, j); % x - 1 (left)
    if left_visited
        if right_visited
            Tx = min(T(i+1, j), T(i-1, j));
        else
            Tx = T(i-1, j);
        end
    elseif right_visited
        Tx = T(i+1, j);
    else
        Tx = Inf;
    end
    
    top_visited = Visited(i, j+1); % j + 1 (top)
    bottom_visited = j ~= 1 && Visited(i, j-1); % j - 1 (bottom)
    
    if top_visited
        if bottom_visited
            Ty = min(T(i, j+1), T(i, j-1));
        else
            Ty = T(i, j+1);
        end
    elseif bottom_visited
        Ty = T(i, j-1);
    else
        Ty = Inf;
    end
    
    if isinf(Tx) && isinf(Ty)
        error('all neighbours have infinite time of crossing')
    end
    
    Fdet = calculateLocalDetachmentRate(i, j, kDet, grid, xcenter, ycenter);
    
    if Fdet == 0
        warning('Detachment speed equals 0, thus infinite time of crossing')
        T_new = Inf;
        return
    end
    
    T_new = computeRoot(Tx, Ty, Fdet, grid.dx);
end
