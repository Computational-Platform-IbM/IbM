function root = computeRoot(Tx, Ty, Fdet, dx)
    % Solve the quadratic equation for the gradient (approximation) of the
    % detachment front
    %
    % Tx, Ty: Time of detachment value of the neighbouring gridcell in the
    %   x and y direction respectively
    % Fdet: speed of detachment at the specific gridcell
    % dx: discretization resolution
    %
    % -> root: solution of the quadratic equation
    
    % first init all parameters
    a = 0;
    b = 0;
    c = -(dx/Fdet)^2;
    
    % if Tx is finite, then add those respective terms to parameters
    if isfinite(Tx)
        a = a + 1;
        b = b - 2*Tx;
        c = c + Tx^2;
    end
    
    % same for Ty
    if isfinite(Ty)
        a = a + 1;
        b = b - 2*Ty;
        c = c + Ty^2;
    end
    
    % if all are infinite, then return inf
    if a == 0
        root = Inf;
        return
    end
    
    % now get the 2 solutions
    D = sqrt(b^2 - 4*a*c);
    
    if D < 0
        error('ValueError:Should always be above 0...')
    end
    
    % positive solution is only valid
    root = (-b+D) / (2*a);
end
