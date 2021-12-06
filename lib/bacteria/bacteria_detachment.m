function bac = bacteria_detachment(bac, grid, constants, settings, timestep)
    % Rough detachment: Remove all bacteria that are outside of maximum radius of granule
    %
    % bac: struct containing all information regarding the bacteria
    % constants: struct containing all simulation constants
    %
    % -> bac: see above
    
    switch settings.detachment
        case 'none'
            return
        case 'naive'
            bac_distance_from_center = sqrt((bac.x - grid.dx * grid.nX/2).*(bac.x - grid.dx * grid.nX/2) + (bac.y - grid.dy * grid.nY/2).*(bac.y - grid.dy * grid.nY/2));
            bac_detach = bac_distance_from_center > constants.max_granule_radius;
            nCellsDetach = sum(bac_detach);

            if nCellsDetach % <TODO: check if these are all struct-elements />
                bac = killBacs(bac, bac_detach);
            end
            
        case 'mechanistic'
            % update/generate grid2nBacs struct (after division previous is
            % not valid anymore)
            [grid2bac, grid2nBacs] = determine_where_bacteria_in_grid(grid, bac);
            T = calcTimeOfDetach(bac, grid, grid2nBacs, constants);
            
            ratio = timestep ./ T;
            
            % decrease mass of bacteria (radius will be updated when
            % needed) -> Erosion
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

            % remove bacteria with T < timestep -> Detachment
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

            if nCellsDetach
                bac = killBacs(bac, bac_detach);
            end
            
            % remove bacteria that are way too small (only on outside, due to erosion)
            % factor 2 smaller than inactive bacteria should only be reached with erosion,
            % cells on the inside of the granule that are too small, are
            % not removed
            mask_tooSmall = bac.molarMass * constants.bac_MW < constants.min_bac_mass_grams / 2;
            xc = mean(bac.x(bac.active));
            yc = mean(bac.y(bac.active));
            dist = sqrt((bac.x - xc).^2 + (bac.y - yc).^2);
            mask_outside = dist > max(dist(bac.active & ~mask_tooSmall)) - grid.blayer_thickness;
            mask_outsideCellRemoval = mask_tooSmall & mask_outside; 
            nCellsRemoved = sum(mask_outsideCellRemoval);

            if nCellsRemoved
                bac = killBacs(bac, mask_outsideCellRemoval);
                nCellsDetach = nCellsDetach + nCellsRemoved;
            end

            
        otherwise
            error('Detachment of %s is unknown.', settings.detachment)
    end
    
    fprintf('%d cells detached from the granule\n', nCellsDetach)

            
end
