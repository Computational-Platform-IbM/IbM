function conc = multigrid_diffusion(grid, nCells, dt);


end

%% V_cycle function

function phi = V_Cycle(phi,f,h)
	% Recursive V-Cycle Multigrid for solving the Diffusion equation on a uniform grid of spacing h

	% Pre-Smoothing
	phi = smoothing(phi,f,h);
	
	% Compute Residual Errors
	r = residual(phi,f,h);
	
	% Restriction
	rhs = restriction(r);

	eps = zeros(size(rhs));

	% stop recursion at smallest grid size, otherwise continue recursion
	if smallest_grid_size_is_achieved
        	eps = smoothing(eps,rhs,2*h);
	else        
        	eps = V_Cycle(eps,rhs,2*h);        
	end
	
	% Prolongation and Correction
	phi = phi + prolongation(eps);
	
	% Post-Smoothing
	phi = smoothing(phi,f,h);    
end

%% smoothing function

function x = gauss_seidel(A, b, x, iters)
    % perform one iteration of smoothing using the Gauss-Seidel method on
    % the equation: Ax=b
    for i = 1:iters
        for j = 1:size(A,1)
            x(j) = (b(j) - sum(A(j,:).*x) + A(j,j)*x(j))/A(j,j);
        end
    end
end
    
    
    
    
end