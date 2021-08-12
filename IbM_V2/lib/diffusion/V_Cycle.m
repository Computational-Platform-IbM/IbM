function phi = V_Cycle(phi, f, L_0, L_restriction, L_prolongation, maxDepth, depth, iter_pre, iter_post, iter_final)
	% Recursive V-Cycle Multigrid for solving the Diffusion equation on a uniform grid
    % 
    % phi: matrix x in Ax=b
    % f: matrix b in Ax=b
    % L_0: basis of the laplacian stencil
    % L_restriction: stencil for restriction
    % L_prolongation: stencil to prolongate error correction
    % maxDepth: minimum number of gridcells in the highest recursion of the
    %   V-cycle, i.e. in the deepest level the grid is at least 
    %   (maxDepth, maxDepth)
    % depth: current depth (starting at 0)
    % iter_pre: number of smoothing steps before recursion
    % iter_post: number of smoothing steps after error correction
    % iter_final: number of smoothing steps on the deepest level
    %
    % -> phi: advanced matrix x towards the solution in Ax=b

    % Create correct left-hand side stencil
    L_lhs = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*depth))*L_0;
    
	% Pre-Smoothing
    for i = 1:iter_pre
        phi = smoothing(phi,f,L_lhs);
    end
	
	% Compute Residual Errors
    r = residual(phi,f,L_lhs);
    
	% Restriction
	rhs = restriction(r, L_restriction);

	eps = zeros(size(rhs));

	% stop recursion at smallest grid size, otherwise continue recursion
    if ceil(sqrt(numel(phi))) <= maxDepth
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*(depth+1)))*L_0;
        for i = 1:iter_final
            eps = smoothing(eps,rhs,L_lhs_deeper);
        end
    else
        eps = V_Cycle(eps, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, iter_pre, iter_post, iter_final);        
    end
    
	% Prolongation and Correction
	phi = phi + prolongation(eps, L_prolongation, size(phi));
    
	% Post-Smoothing
    for i = 1:iter_post
        phi = smoothing(phi,f,L_lhs);
    end
end