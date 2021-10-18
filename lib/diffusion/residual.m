function r = residual(phi, rhs, L_lhs)
    % calculate the residual, given a concentration matrix and a lhs
    % stencil.
    r = rhs - convn(phi, L_lhs, 'same');
end