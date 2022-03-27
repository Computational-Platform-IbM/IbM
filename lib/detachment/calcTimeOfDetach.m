function T = calcTimeOfDetach(bac, grid, grid2bac, grid2nBacs, constants)
    % Calculate the time of detachment for each gridcell in the simulation
    % domain
    %
    % bac: struct containing all information regarding the bacteria
    % grid: struct containing all information regarding spacial
    %   discretizations
    % grid2nBacs: matrix with per gridcell the amount of bacteria in there
    %
    % -> T: matrix with per gridcell the time of detachment

    % extract variables
    kDet = constants.kDet;
    
    % calculate center of the granule
    xcenter = mean(bac.x);
    ycenter = mean(bac.y);

    % where are bacteria located in the grid? apply some smoothing of
    % logical grid with bacteria (grid2nBacs) for the case of high
    % resolution grid
    detachment_grid = grid;
    detachment_grid.blayer_thickness = constants.kDist*constants.bac_max_radius*2;
    [aggregate, ~] = determine_diffusion_region(grid2bac, grid2nBacs, bac, detachment_grid);
    biofilm = aggregate > 0;

    % create three matrices with all grid cells
    T = zeros(size(grid2nBacs));

    % T(~biofilm) = 0; % redundant, because initialised with 0
    Visited = ~biofilm;
    % plotLogicalGrid(grid, Visited, 'Visited'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}


    % find narrow band
    kernel = zeros(3,3);
    kernel([1, 3],2) = -1/4;
    kernel(2,[1, 3]) = -1/4;
    kernel(2,2) = 1;
    Narrow_band = convn(biofilm, kernel, 'same') > 0;

    Far = biofilm & ~Narrow_band;

    % set far away point to T infinity
    T(Far) = Inf;

    [i, j] = find(Narrow_band);
    for k = 1:length(i)
        ii = i(k);
        jj = j(k);
        Fdetach = calculateLocalDetachmentRate(ii, jj, kDet, grid, xcenter, ycenter);

        % --------- IMPORTANT ----------
        % What is the impact of the number of free neighbours?
        % In idynomics it is calculated with the number of non-biomass
        % gridcells, but mathematically it does not make too much sense...?
        % --------END IMPORTANT --------

        nFreeNb = getFreeNeighbourCount(ii, jj, Visited);

        T(ii, jj) = grid.dx / (Fdetach * nFreeNb);
    end


    % initialise stacks (not a real minheap, but we are not going to get closer
    % than this with MATLAB... <sadface>
    nStack = sum(Narrow_band, 'all'); % number of points in T_stack
    T_stack = zeros(nStack,1);
    index_stack = zeros(nStack, 2);
    T_stack(1:nStack,1) = T(Narrow_band);
    [i, j] = find(Narrow_band);
    index_stack(1:nStack, :) = [i, j];

    % ------------ Fast Marching --------------
    while nStack
        % take narrow-band value with the lowest T value
        [T_stack(1:nStack), I] = sort(T_stack(1:nStack));
        index_stack = index_stack(I, :);

        i_point = index_stack(1,1);
        j_point = index_stack(1,2); 

        % remove from Narrow_band, T_stack and index_stack
        Narrow_band(i_point, j_point) = 0;
        T_stack(1) = [];
        index_stack(1,:) = [];
        nStack = nStack - 1;

        % add point to visited
        Visited(i_point, j_point) = 1;

        offSet = [0, 0, +1, -1;  % i offset
                  +1, -1, 0, 0]; % j offset

        % update neighbours with new values and correct classifications
        for nb = 1:length(offSet)
            i_nb = i_point + offSet(1, nb);
            j_nb = j_point + offSet(2, nb);

            if Visited(i_nb, j_nb) % if neighbour was already visited, then continue to next neighbour
                continue
            end

            % for all neighbours of the point (from narrow-band or far)
            if Narrow_band(i_nb, j_nb)
                % get index in stack
                old_index = find(T_stack == T(i_nb, j_nb));
            else
                old_index = 0;
            end

            % recalculate T value
            T_val = recalculateT(T, i_nb, j_nb, kDet, grid, Visited, xcenter, ycenter);


            % if it was already in stack, then update value, otherwise add it
            % to the stack        
            if Far(i_nb, j_nb) % remove from far and add to narrow band and add to stack
                Far(i_nb, j_nb) = 0;
                Narrow_band(i_nb, j_nb) = 1;
                T_stack(nStack + 1) = T_val;
                index_stack(nStack+1, :) = [i_nb, j_nb];
                nStack = nStack + 1;
            else % update stack
                T_stack(old_index) = T_val;
            end
            T(i_nb, j_nb) = T_val;
        end


    % ------- DEBUG PLOTS FOR ANIMATION -------
    %     plotLogicalGrid(grid, Visited, 'Visited'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
    %     plotLogicalGrid(grid, Narrow_band, 'Narrow band'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
    %     plotLogicalGrid(grid, Far, 'Far'); % {'Biofilm', 'Visited', 'Narrow band', 'Far'}
    %     plotDetachTime(grid, T, kDet);
    %     drawnow()
    % ------------ END DEBUG PLOTS ------------


    end
end



%% helper functions
function nFreeNb = getFreeNeighbourCount(i, j, Visited)
    % count the number of non-biomass neighbouring gridcells are there for
    % gridcell (i, j)
    
    nFreeNb = 0;
    if Visited(i+1, j)
        nFreeNb = nFreeNb + 1;
    end
    if Visited(i-1, j)
        nFreeNb = nFreeNb + 1;
    end
    if Visited(i, j+1)
        nFreeNb = nFreeNb + 1;
    end
    if Visited(i, j-1)
        nFreeNb = nFreeNb + 1;
    end
end





