function phi_fine = prolongation(phi_coarse, L_prolongation, sz)
    % linearly interpolates the coarse values to the fine values
    phi_fine = zeros(sz);
    phi_fine(1:2:end, 1:2:end) = phi_coarse;
    phi_fine = convn(phi_fine, L_prolongation, 'same');
end