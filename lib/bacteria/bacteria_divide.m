function [bac, cycle] = bacteria_divide(bac, constants)
    % Divide bacteria that are above the mass threshold
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    cycle = 0;
    
    while sum(bac.molarMass * constants.bac_MW > constants.max_bac_mass_grams) > 1
        cycle = cycle + 1;
        
        mask_tooBig = bac.molarMass * constants.bac_MW > constants.max_bac_mass_grams;
        nCellsTooBig = sum(mask_tooBig);

        %% copy x, y, and active from parents
        fi = rand(nCellsTooBig, 1) * 2 * pi;
        new_x = bac.x(mask_tooBig) + bac.radius(mask_tooBig) .* cos(fi);
        new_y = bac.y(mask_tooBig) + bac.radius(mask_tooBig) .* sin(fi);
        new_species = bac.species(mask_tooBig);
        new_mu = bac.mu(mask_tooBig);
        new_active = ones(nCellsTooBig, 1, 'logical');

        %% split mass over parent and child
        % mass of child and parent
        new_molarMass = bac.molarMass(mask_tooBig) .* (0.45 + 0.1 * rand(nCellsTooBig, 1));
        bac.molarMass(mask_tooBig) = bac.molarMass(mask_tooBig) - new_molarMass;
        % radius of child and parent
        new_radius = ((new_molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))).^(1/3);
        bac.radius(mask_tooBig) = ((bac.molarMass(mask_tooBig) * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))).^(1/3);

        %% update variables
        bac.x = [bac.x; new_x];
        bac.y = [bac.y; new_y];
        bac.species = [bac.species; new_species];
        bac.molarMass = [bac.molarMass; new_molarMass];
        bac.radius = [bac.radius; new_radius];
        bac.mu = [bac.mu; new_mu];
        bac.active = [bac.active; new_active];
        
    end
end