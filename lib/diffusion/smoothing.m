function phi = smoothing(phi, rhs, L_lhs)
    % apply a Jacobi iteration on the system using convolution
    L_sm = L_lhs;
    L_sm(2, 2) = 0;
    phi = (rhs - convn(phi, L_sm, 'same'))/L_lhs(2,2);
end