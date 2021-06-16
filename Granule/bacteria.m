% %%%% bacteria.m -Generates the IbM: Position of the cells and division
% %%%% R.bac is the structure that keeps all the information related with this
% subscript

function [R, StVLiq, BacterialChange] = bacteria(L, G, X, R)
    detaching = 1; % set detachment to default naive way
    St = R.St;
    Sxy = R.Sxy;
    bac = R.bac;
    numStVLiq = St.numStVLiq;
    numStVLiq2 = St.numStVLiq2;

    U = zeros(St.numSt, 1);

    for k = 1:St.numSt

        if k <= numStVLiq2
            U(k) = Sxy.Sbc_Dir(k);
        elseif k <= numStVLiq
            U(k) = G(Sxy.nT * ((k - numStVLiq2) - 1) + 1);
        elseif k > numStVLiq
            U(k) = sum((X .* (bac.atrib(:, 5) == (k - numStVLiq))) / Sxy.Vg) / Sxy.nT;
        end

    end

    % %%%% Update of the variables
    St.StV = U;
    St.StVLiq = St.StV(1:numStVLiq); % Liquid and gas variables
    St.StVLiq2 = St.StV(1:numStVLiq2); % Liquid variables
    St.StVX = St.StV(numStVLiq + 1:end); % Biomass
    StVLiq = [L; G];

    % Sxy.StVLiq = [L; G];
    % Sxy.StVLiq2 = L;
    % Sxy.StVGas = G;
    R.bac.atrib(:, 3) = X;
    R.bac.atrib(:, 6) = ((X * bac.bac_MW / bac.bac_rho) * (3 / (4 * pi))).^(1/3);

    R.Sxy = Sxy;
    R.St = St;

    % check if the number of any bacteria have to divide or die
    bac_m = R.bac.atrib(:, 3) * bac.bac_MW; % Using grams here
    BacChange = (sum(bac_m > bac.bac_mmax) >= 1) + (sum(bac_m < bac.bac_mmin) >= 1);

    if BacChange
        % let those bacteria divide and die
        R.bac = bac_division(R.bac);
    end

    % shove the cells (always, even if no bacteria divide/die as growth can
    % cause overlap as well)
    qt = shoving.BiomassQuadtree(0, Sxy.maxxSys, 0, Sxy.maxySys);
    R.bac = shove(R.bac, qt);

    % calculate derived bacterial variable(s)
    R.bac.bac_rho_bio = sum(R.bac.atrib(:, 3) * R.bac.bac_MW) / ((max(R.bac.atrib(:, 2)) - min(R.bac.atrib(:, 2))) * (max(R.bac.atrib(:, 1)) - min(R.bac.atrib(:, 1))) * (1e-6));

    if detaching
        % perform detachment
        [R.bac, detachChange] = detachCells(R.bac);
    else
        detachChange = 0;
    end
    
    % update number of each cell type
    for i = 1:length(R.bac.bac_ns)
        R.bac.bac_ns(i) = sum(R.bac.atrib(:, 5) == i);
    end

    BacterialChange = BacChange || detachChange;
end

function bac = bac_division(bac)
    bac_n = bac.bac_n;
    bac_x = bac.atrib(:, 1);
    bac_y = bac.atrib(:, 2);
    bac_m = bac.atrib(:, 3) * bac.bac_MW;
    bac_a = bac.atrib(:, 4);
    bac_s = bac.atrib(:, 5);
    bac_r = bac.atrib(:, 6);
    bac_yield = bac.atrib(:, 7);
    bac_Ks = bac.bac_Ks;
    z = 1;

    while sum(bac_m > bac.bac_mmax) >= 1 || sum(bac_m < bac.bac_mmin) >= 1

        % if there are cells that are too small, die
        mask_tooSmall = bac_m < bac.bac_mmin;
        nCellsTooSmall = sum(mask_tooSmall);

        if nCellsTooSmall
            bac_x(mask_tooSmall) = [];
            bac_y(mask_tooSmall) = [];
            bac_a(mask_tooSmall) = [];
            bac_s(mask_tooSmall) = [];
            bac_r(mask_tooSmall) = [];
            bac_m(mask_tooSmall) = [];
            bac_Ks(mask_tooSmall, :) = [];
            bac_yield(mask_tooSmall) = [];
        end

        % if there are cells that are too big, divide
        mask_tooBig = bac_m > bac.bac_mmax;
        nCellsTooBig = sum(mask_tooBig);

        if nCellsTooBig
            fi = (1 + (-1 - 1) .* rand(nCellsTooBig, 1)) * 2 * pi;
            new_x = bac_x(mask_tooBig) + bac_r(mask_tooBig) .* cos(fi);
            new_y = bac_y(mask_tooBig) + bac_r(mask_tooBig) .* sin(fi);
            new_a = bac_a(mask_tooBig);
            new_s = bac_s(mask_tooBig);
            new_yield = bac_yield(mask_tooBig);
            new_Ks = bac_Ks(mask_tooBig, :);
            % mass of child and parent
            new_m = bac_m(mask_tooBig) .* (0.45 + 0.1 * rand(nCellsTooBig, 1));
            bac_m(mask_tooBig) = bac_m(mask_tooBig) - new_m;
            % radius of child and parent
            new_r = ((new_m / bac.bac_rho) * (3 / (4 * pi))).^(1/3);
            bac_r(mask_tooBig) = ((bac_m(mask_tooBig) / bac.bac_rho) * (3 / (4 * pi))).^(1/3);

            % update variables
            bac_x = [bac_x; new_x];
            bac_y = [bac_y; new_y];
            bac_a = [bac_a; new_a];
            bac_s = [bac_s; new_s];
            bac_m = [bac_m; new_m];
            bac_yield = [bac_yield; new_yield];
            bac_Ks = [bac_Ks; new_Ks];
            bac_r = [bac_r; new_r];
        end

        % update number of bacteria
        bac_n = bac_n + nCellsTooBig - nCellsTooSmall;

        % save variables
        bac.bac_n = bac_n;
        bac.bac_Ks = bac_Ks;
        bac.atrib = [bac_x, bac_y, bac_m / bac.bac_MW, bac_a, bac_s, bac_r, bac_yield];

        z = z + 1;
    end

end

function bac = shove(bac, quadtree)
    bac_x = bac.atrib(:, 1);
    bac_y = bac.atrib(:, 2);
    bac_r = bac.atrib(:, 6);

    tic
    r = quadtree.pushing2D(length(bac_x), bac_x, bac_y, bac_r, 0.1, bac.s_dist);
    bac_x_final = r.bac_x;
    bac_y_final = r.bac_y;
    toc

    bac.atrib(:, 1) = bac_x_final;
    bac.atrib(:, 2) = bac_y_final;
end

function [bac, nCellsDetached] = detachCells(bac)
    bac_n = bac.bac_n;
    bac_x = bac.atrib(:, 1);
    bac_y = bac.atrib(:, 2);
    bac_m = bac.atrib(:, 3); % back to using moles here
    bac_a = bac.atrib(:, 4);
    bac_s = bac.atrib(:, 5);
    bac_r = bac.atrib(:, 6);
    bac_yield = bac.atrib(:, 7);
    bac_Ks = bac.bac_Ks;
    detach_method = 'naive'; % set to default detachment

    % Remove bacteria outside of max-granule size (naive detachment)
    mask_detach = bacterialDetachment(detach_method, bac_x, bac_y, bac);

    nCellsDetached = sum(mask_detach);
    fprintf("%d cells detached. ", nCellsDetached);

%     if nCellsDetached
%         figure(1);
%         clf;
%         colors = hsv(2);
%         %                         colors = [1,1,1 ; 1,1,1; 0.8, 0.1, 0.1];
%         for i = 1:2
%             max_mu(i) = max(bac_a(bac_s == i));
%             min_mu(i) = min(bac_a(bac_s == i)) - 1e-20;
%         end
% 
%         maxOut = 1;
%         minOut = 0.3;
% 
%         for i = 1:bac.bac_n
%             mapping = (bac_a(i) - min_mu(bac_s(i))) / (max_mu(bac_s(i)) - min_mu(bac_s(i))) * (maxOut - minOut) + minOut;
% 
%             if mask_detach(i)
%                 c = [1, 1, 1];
%             else
%                 c = colors(bac_s(i), :) * mapping;
%             end
% 
%             rectangle('Curvature', [1 1], 'Position', [bac_x(i) - bac_r(i), bac_y(i) - bac_r(i), 2 * bac_r(i), 2 * bac_r(i)], 'LineWidth', 1, 'FaceColor', c, 'EdgeColor', [0, 0, 0, 0.5]);
%         end
% 
%         axis equal;
%     end

    % delete detached bacteria
    bac_m(mask_detach) = [];
    bac_r(mask_detach) = [];
    bac_x(mask_detach) = [];
    bac_y(mask_detach) = [];
    bac_s(mask_detach) = [];
    bac_a(mask_detach) = [];
    bac_yield(mask_detach) = [];
    bac_Ks(mask_detach, :) = [];

    % update number of bacteria
    bac_n = bac_n - nCellsDetached;

    % save variables
    bac.bac_n = bac_n;
    bac.bac_Ks = bac_Ks;
    bac.atrib = [bac_x, bac_y, bac_m, bac_a, bac_s, bac_r, bac_yield];
end

function mask = bacterialDetachment(method, x, y, bac)
    %%% create boolean mask for which cells are to be detached
    center = bac.bac_c;
    dx = x - center;
    dy = y - center;
    normd = sqrt(dx .* dx + dy .* dy);
    mask = normd > (bac.bac_ymax / 2); % default: remove bacteria outside of max-granule size

    if strcmp(method, 'linear')
        mask_outside = findOutsideBacs(dx, dy, normd, bac.bac_rmax * 2);
        r1 = bac.bac_ystart / 2;
        r2 = bac.bac_ymax / 2;

        detach_chance = min(max((normd - r1) ./ (r2 - r1), 0), 1); % clamped between 0 and 1;
        mask = mask_outside & (detach_chance > rand(length(normd), 1));

    elseif strcmp(method, 'quadratic')
        mask_outside = findOutsideBacs(dx, dy, normd, bac.bac_rmax * 2);
        r1 = bac.bac_ystart / 2;
        r2 = bac.bac_ymax / 2;

        detach_chance = min(max((normd - r1).^2 ./ ((r2 - r1)^2), 0), 1); % clamped between 0 and 1;
        mask = mask_outside & (detach_chance > rand(length(normd), 1));
    end

    function bac_mask_outsideBacs = findOutsideBacs(dx, dy, normd, resolution)
        % section into different regions based on the angle
        d_theta = (resolution / max(normd)); % divide into sections with this delta_angle
        nSections = 2 * pi / d_theta + 1; % number of sections that need to be created. Also add one as endpoint in linspace to ensure even spacing
        sectionStarts = linspace(-pi, pi, nSections); % create list of section starts (also includes pi itself, but is not used)

        % where is each bacterium located (angle)
        bacterialAngles = atan2(dy, dx);

        bac_mask_outsideBacs = zeros(length(normd), 1);

        for i = 1:length(sectionStarts) - 1
            % find all bacteria per section
            ind_theta = (bacterialAngles > sectionStarts(i)) .* (bacterialAngles <= sectionStarts(i + 1));
            [~, index_bac_outside] = max(normd .* ind_theta);
            bac_mask_outsideBacs(index_bac_outside) = true;
        end

    end

end
