function rhs = restriction(r, L_restriction)
    % creates new rhs matrix for a coarser mesh, based on the residuals.
    rhs = convn(r, L_restriction, 'same');
    rhs = rhs(1:2:end, 1:2:end); % *4?
end