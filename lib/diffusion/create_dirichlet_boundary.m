function x_dirichlet = create_dirichlet_boundary(x, value)
    % Construct dirichlet boundary around the given matrix x
    %
    % x: matrix of (:,:) around which to construct boundary
    % value: boundary value for this compound
    %
    % x_dirichlet: padded matrix x with the boundary values
    
    x_dirichlet = zeros(size(x) + 2);
    x_dirichlet(2:end-1, 2:end-1) = x;

    x_dirichlet([1, end], :) = value;
    x_dirichlet(:, [1, end]) = value;
end
