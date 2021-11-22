function bac = bacteria_detachment(bac, grid, constants, settings, timestep)
    % Rough detachment: Remove all bacteria that are outside of maximum radius of granule
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    switch settings.detachment
        case 'naive'
            bac_distance_from_center = sqrt((bac.x - grid.dx * grid.nX/2).*(bac.x - grid.dx * grid.nX/2) + (bac.y - grid.dy * grid.nY/2).*(bac.y - grid.dy * grid.nY/2));
            bac_detach = bac_distance_from_center > constants.max_granule_radius;
            nCellsDetach = sum(bac_detach);

            if nCellsDetach % <TODO: check if these are all struct-elements />
                bac.x(bac_detach) = [];
                bac.y(bac_detach) = [];
                bac.radius(bac_detach) = [];
                bac.species(bac_detach) = [];
                bac.molarMass(bac_detach) = [];
                bac.active(bac_detach) = [];
                bac.mu(bac_detach) = [];
            end
            
        case 'mechanistic'
            % update/generate grid2nBacs struct (after division previous is
            % not valid anymore)
            [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
            T = calcTimeOfDetach(bac, grid, grid2nBacs, constants);
            
            ratio = timestep ./ T;
            
            % decrease mass of bacteria (radius will be updated when
            % needed)
            [i, j] = find(ratio > 0 & ratio < 1);
            for k=1:length(i)
                ii = i(k);
                jj = j(k);
                r = ratio(ii,jj);
                iBacs = nonzeros(grid2bac(ii, jj, :));
                for n = 1:length(iBacs)
                    iBac = iBacs(n);        
                    bac.molarMass(iBac) = bac.molarMass(iBac) * (1 - r);
                end    
            end

            % remove bacteria with T < timestep
            [i, j] = find(ratio >= 1 & ratio < Inf);
            bac_detach = zeros(length(i)*size(grid2bac, 3), 1, 'uint32');
            nDetach = 0;
            for k = 1:length(i)
                ii = i(k);
                jj = j(k);
                iBacs = nonzeros(grid2bac(ii, jj, :));
                n_temp = length(iBacs);
                bac_detach(nDetach + 1:nDetach + n_temp) = iBacs;
                nDetach = nDetach + n_temp;
            end
            bac_detach = bac_detach(1:nDetach);
            nCellsDetach = sum(bac_detach);



            bac.x(bac_detach) = [];
            bac.y(bac_detach) = [];
            bac.radius(bac_detach) = [];
            bac.species(bac_detach) = [];
            bac.molarMass(bac_detach) = [];
            bac.active(bac_detach) = [];
            bac.mu(bac_detach) = [];            
            
            
            
        otherwise
            error('Detachment of %s is unknown.', settings.detachment)
    end
    
    fprintf('%d cells detached from the granule\n', nCellsDetach)

            
end
