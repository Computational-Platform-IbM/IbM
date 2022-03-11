function bac = update_bacterial_mass(bac, dT)
    % Integrate bacterial growth during dT and update the respective mass
    % of each bacterium.
    %
    % bac: struct containing all information regarding the bacteria
    % dT: time over which to integrate bacterial growth
    %
    % -> bac: see above

    % mu >= 0: grow active cells
    bac.molarMass(bac.mu >= 0 & bac.active) = bac.molarMass(bac.mu >= 0 & bac.active) + dT * bac.mu(bac.mu >= 0 & bac.active) .* bac.molarMass(bac.mu >= 0 & bac.active);

    % mu < 0: shrink only active
    bac.molarMass(bac.mu < 0 & bac.active) = bac.molarMass(bac.mu < 0 & bac.active) + dT * bac.mu(bac.mu < 0 & bac.active) .* bac.molarMass(bac.mu < 0 & bac.active);
end
