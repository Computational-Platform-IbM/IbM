function bac = killBacs(bac, indices)
    % Remove cells at specific indices from the simulation domain
    
    bac.x(indices) = [];
    bac.y(indices) = [];
    bac.radius(indices) = [];
    bac.species(indices) = [];
    bac.molarMass(indices) = [];
    bac.active(indices) = [];
    bac.mu(indices) = [];
end