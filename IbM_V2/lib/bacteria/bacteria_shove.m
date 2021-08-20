function bac = bacteria_shove(bac, constants)
    % Shove bacteria so that no overlap is present in the bio-aggregate
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    r = quadtree.pushing2D(length(bac.x), bac.x, bac.y, bac.radius, 0.1, constants.bac_max_radius * 2, constants.kDist);
    bac.x = r.bac_x;
    bac.y = r.bac_y;
end