function bac = update_bacterial_mass(bac, dT)
    % Integrate bacterial growth during dT and update the respective mass
    % of each bacterium.
    %
    % bac: struct containing all information regarding the bacteria
    % dT: time over which to integrate bacterial growth
    %
    % -> bac: see above

    bac.molarMass = bac.molarMass + dT * bac.mu .* bac.molarMass;
end
