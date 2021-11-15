function bac = bacteria_shove(bac, grid, constants)
    % Shove bacteria so that no overlap is present in the bio-aggregate
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    qt = shoving.BiomassQuadtree(0, grid.dx*grid.nX, 0, grid.dy*grid.nY);
    r = qt.pushing2D(length(bac.x), bac.x, bac.y, bac.radius, 0.1, constants.bac_max_radius * 2, constants.kDist); % .jar file in /lib is from granule, thus does not have kDist implemented yet
    bac.x = r.bac_x;
    bac.y = r.bac_y;
end