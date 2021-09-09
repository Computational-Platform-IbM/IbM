function r = residual(phi, rhs, L_lhs, diffRegion)
    % calculate the residual, given a concentration matrix and a lhs
    % stencil. Only returns values in the diffusion region.
    r = rhs - convn(phi, L_lhs, 'same');
    r = r(diffRegion);
end