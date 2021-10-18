clear all;

%{
Main Takeaway: overhead of parfor is too large for efficient
parallelisation... (for a 32x32 system with 500 bacteria), at a certain
point this does become better with more bacteria (and more chunks);
%}
% load('testingRMatrix.mat')
load('testingRMatrix_medium.mat')
% load('testingRMatrix_large.mat')
% diffRegion = diffusion_region;

% grid.dx = 4e-6;
% grid.dy = grid.dx;
% grid.nX = grid.;%size(diffRegion, 1); %32;
% grid.nY = 255;%size(diffRegion, 2); %32;
% grid.blayer_thickness = 5e-6;


%% determine diffRegion and focus region
[diffRegion, focus_region] = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid); 


%% create chunks
% reorganise bacs into chunks (based on diffRegion)
% 3-by-3 chunks
nChunks_dir = 2;
nChunks = nChunks_dir^2;
dx = focus_region.x1 - focus_region.x0 + 1;
% dx = size(diffRegion, 1);
dx_chunk = ceil(dx/nChunks_dir);
dy = focus_region.y1 - focus_region.y0 + 1;
% dy = size(diffRegion, 2);
dy_chunk = ceil(dy/nChunks_dir);

chunk_ind_x = zeros(nChunks_dir, 2);

temp_x = focus_region.x0;
for ixChunk = 1:nChunks_dir
    chunk_ind_x(ixChunk,[1,2]) = [temp_x, temp_x + dx_chunk - 1];
    temp_x = temp_x + dx_chunk;
end
chunk_ind_x(end) = focus_region.x1;
% chunk_ind_x(end) = size(diffRegion, 1);

chunk_ind_y = zeros(nChunks_dir, 2);

temp_y = focus_region.y0;
for iyChunk = 1:nChunks_dir
    chunk_ind_y(iyChunk,[1,2]) = [temp_y, temp_y + dy_chunk - 1];
    temp_y = temp_y + dy_chunk;
end
chunk_ind_y(end) = focus_region.y1;
% chunk_ind_y(end) = size(diffRegion, 2);


%% reorganize indices in bacterial struct based on respective chunk
chunk2nBacs = zeros(nChunks, 1);
for iChunk = 1:nChunks
    ix = floor((iChunk - 1) / nChunks_dir) + 1;
    iy = mod(iChunk - 1, nChunks_dir) + 1;
    xRange = chunk_ind_x(ix,1):chunk_ind_x(ix,2);
    yRange = chunk_ind_y(iy,1):chunk_ind_y(iy,2);
    chunk2nBacs(iChunk) = sum(grid2nBacs(xRange,yRange), 'all');
end

ix = floor(bac.x / grid.dx) + 1; % +1 because of MATLAB indexing
iy = floor(bac.y / grid.dy) + 1;
bac_grid_ixiy = [ix, iy];

ixChunk = floor((ix - focus_region.x0) / dx_chunk) + 1;
iyChunk = floor((iy - focus_region.y0) / dy_chunk) + 1;

bac_chunk = (iyChunk - 1) + nChunks_dir * (ixChunk - 1) + 1; % why is Matlab indexing this stupid...?!
[~, sortChunkIndex] = sort(bac_chunk);

bacOffset = zeros(nChunks + 1, 1);
bacOffset(2:end) = cumsum(chunk2nBacs);


% reorganise bac struct
bac.x = bac.x(sortChunkIndex);
bac.y = bac.y(sortChunkIndex);
bac.species = bac.species(sortChunkIndex);
bac.molarMass = bac.molarMass(sortChunkIndex);
bac.radius = bac.radius(sortChunkIndex);
bac.active = bac.active(sortChunkIndex);
% bac.mu

%% redo grid2bac with new indices
[grid2bac, ~] = determine_where_bacteria_in_grid(grid, bac);

tStart = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate_reaction_matrix_chunky
    Keq = constants.Keq;
    chrM = constants.chrM;
    Vg = constants.Vg;
    compounds = constants.StNames;
    reactive_form = constants.react_v;
    Ks = constants.Ks;
    Ki = constants.Ki;
    mMetabolism = constants.MatrixMet;
    mDecay = constants.MatrixDecay;
    
    iNH3 = find(strcmp(compounds, 'NH3'));
    iNO2 = find(strcmp(compounds, 'NO2'));
    iO2 = find(strcmp(compounds, 'O2'));
    
    reaction_matrix = zeros(size(conc));
    mu = zeros(size(bac.x));
    
    % pre-compute for bulk-liquid (at 1,1 there should never be diffusion
    % layer)
    Sh_bulk = 10^-pH(1, 1);
    [spcM, Sh_bulk] = solve_pH(Sh_bulk, [reshape(conc(1,1,:), [], 1, 1); 1; 0], Keq, chrM, constants.constantpH, constants.pHtolerance); % <C: why [...; 1; 0]? />
    pH_bulk = -log10(Sh_bulk);
    
    % check if pH is solved correctly
    if any(spcM < 0)
        warning('DEBUG:actionRequired', 'debug: negative concentration encountered after pH calculation...');
    end

    % group constants for easy passing to multiple cores
    constantValues = [pH_bulk, constants.constantpH, constants.pHtolerance, constants.T, iNH3, iNO2, iO2];
    grouped_bac = [bac.species, bac.molarMass];
    
    % create storage variables for output per chunk
    nCompounds = length(compounds);
    xlens = kron(ones(nChunks_dir,1), diff(chunk_ind_x, [], 2)+1);
    ylens = kron(diff(chunk_ind_y, [], 2)+1, ones(nChunks_dir, 1));
    
    chunk_rMatrix = arrayfun(@(xlen, ylen) zeros(xlen, ylen, nCompounds), xlens, ylens, 'UniformOutput', false);
    chunk_pH = arrayfun(@(xlen, ylen) zeros(xlen, ylen), xlens, ylens, 'UniformOutput', false);
    chunk_mu = arrayfun(@(x) zeros(x, 1), chunk2nBacs, 'UniformOutput', false);
    ix = @(iChunk) floor((iChunk - 1) / nChunks_dir) + 1;
    iy = @(iChunk) mod(iChunk - 1, nChunks_dir) + 1;
    chunk_xRange = arrayfun(@(iChunk) chunk_ind_x(ix(iChunk), 1):chunk_ind_x(ix(iChunk), 2), 1:nChunks, 'UniformOutput', false);
    chunk_yRange = arrayfun(@(iChunk) chunk_ind_y(iy(iChunk), 1):chunk_ind_y(iy(iChunk), 2), 1:nChunks, 'UniformOutput', false);
    chunk_pH_OG = cellfun(@(xR, yR) pH(xR, yR), chunk_xRange, chunk_yRange, 'UniformOutput', false);
    chunk_conc = cellfun(@(xR, yR) conc(xR, yR, :), chunk_xRange, chunk_yRange, 'UniformOutput', false);
    chunk_grid2bac = cellfun(@(xR, yR) grid2bac(xR, yR, :), chunk_xRange, chunk_yRange, 'UniformOutput', false);
    chunk_grid2nBacs = cellfun(@(xR, yR) grid2nBacs(xR, yR), chunk_xRange, chunk_yRange, 'UniformOutput', false);
    chunk_diffRegion = cellfun(@(xR, yR) diffRegion(xR, yR), chunk_xRange, chunk_yRange, 'UniformOutput', false);
    chunk_bacRange = arrayfun(@(iChunk) bacOffset(iChunk)+1:bacOffset(iChunk+1), 1:nChunks, 'UniformOutput', false);
    chunk_grouped_bac = cellfun(@(bacRange) grouped_bac(bacRange, :), chunk_bacRange, 'UniformOutput', false);
    
    % do all stuff in chunks
    for i = 1:10
    parfor iChunk = 1:nChunks
        [chunk_rMatrix{iChunk}, chunk_mu{iChunk}, chunk_pH{iChunk}] = ...
            rMatrix_chunk(chunk_pH_OG{iChunk}, chunk_conc{iChunk}, chunk_grid2bac{iChunk}, chunk_grid2nBacs{iChunk}, chunk_diffRegion{iChunk}, ...
            chunk_grouped_bac{iChunk}, chunk2nBacs(iChunk), bacOffset(iChunk), ...
            reactive_form, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constantValues);
    end
    end
    
    % put everything back into correct matrix/vector
    pH = ones(size(diffRegion))*pH_bulk; % just to shut up Matlab...
    for iChunk = 1:nChunks
        ix = floor((iChunk - 1) / nChunks_dir) + 1;
        iy = mod(iChunk - 1, nChunks_dir) + 1;
        xRange = chunk_ind_x(ix,1):chunk_ind_x(ix,2);
        yRange = chunk_ind_y(iy,1):chunk_ind_y(iy,2);
        bacRange = bacOffset(iChunk)+1:bacOffset(iChunk+1);

        reaction_matrix(xRange, yRange, :) = chunk_rMatrix{iChunk};
        pH(xRange, yRange) = chunk_pH{iChunk};
        mu(bacRange) = chunk_mu{iChunk};
    end
    
    
    % do final unit correction/calculation
    reaction_matrix = reaction_matrix / Vg;                                                                             % [molS/L/h]

tTaken = toc(tStart);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parallelized code (with focus) for %d-by-%d system with %d bacs: %.4f\n', grid.nX, grid.nY, length(bac.x), tTaken);


function [reaction_matrix, mu, pH_new] = rMatrix_chunk(pH, conc, grid2bac, grid2nBacs, diffRegion, ...
    grouped_bac, nBacs, bacOffset, ...
    reactive_form, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constants)

    pH_bulk = constants(1);
    constantpH = constants(2);
    pHtolerance = constants(3);
    T = constants(4);
    iNH3 = constants(5);
    iNO2 = constants(6);
    iO2 = constants(7);

    bac_species = grouped_bac(:,1);
    bac_molarMass = grouped_bac(:,2);

    mu = zeros(nBacs, 1);
    reaction_matrix = zeros(size(conc));
    pH_new = zeros(size(pH));

    % for each gridcell
    for ix = 1:size(conc, 1) % parfor?
        for iy = 1:size(conc, 2)
            if ~diffRegion(ix, iy) % bulk layer
                pH_new(ix, iy) = pH_bulk;
            
            else % in diffusion layer, thus pH calculation needs to be performed
                % calculate pH & speciation
                Sh_old = 10^-pH(ix, iy);
                [spcM, Sh] = solve_pH(Sh_old, [reshape(conc(ix,iy,:), [], 1, 1); 1; 0], Keq, chrM, constantpH, pHtolerance); % <C: why [...; 1; 0]? />
                pH_new(ix, iy) = -log10(Sh);

                % check if pH is solved correctly
                if any(spcM < 0)
                    warning('DEBUG:actionRequired', 'debug: negative concentration encountered after pH calculation...');
                end

                if grid2nBacs(ix, iy) % if also cells are found here, then update reaction matrix
                    % get bacteria in this grid cell
                    iBacs = reshape(grid2bac(ix, iy, 1:grid2nBacs(ix, iy)), [], 1); % n-by-1 vector of bacterial indices

                    % correct for chunk indexing...
                    iBacs = iBacs - bacOffset;
                    
                    speciesGrid = bac_species(iBacs);
                    unique_species = unique(speciesGrid); % which species
                    for i = 1:length(unique_species)
                        curr_species = unique_species(i);

                        % update mu_max per species based on pH
                        [mu_max, maint] = determine_max_growth_rate_and_maint(curr_species, T, Sh);

                        % get reactive concentrations for soluble components
                        reactive_conc = [spcM(iNH3, reactive_form(iNH3)), ...
                            spcM(iNO2, reactive_form(iNO2)), ...
                            spcM(iO2, reactive_form(iO2))];

                        % set mu for all bacteria of same species in that gridcell
                        M = calculate_monod(Ks(curr_species,:), Ki(curr_species, :), reactive_conc);                        % [DN]
                        mu_noMaintenance = mu_max * M;                                                                      % [1/h]
                        mu_withMaintenance = mu_noMaintenance - maint;                                                      % [1/h]
                        mu(iBacs(speciesGrid == curr_species)) = mu_withMaintenance;                                        % [1/h]

                        % update reaction_matrix element for this grid cell
                        concentrationChange = mMetabolism(:, curr_species) * mu_noMaintenance;                              % [molS/molX/h]
                        if mu_withMaintenance < 0
                            concentrationChange = concentrationChange - mDecay(:, curr_species) * mu_withMaintenance;       % [molS/molX/h]
                        end
                        concentrationChange = concentrationChange * sum(bac_molarMass(iBacs(speciesGrid == curr_species))); % [molS/h]

                        reaction_matrix(ix, iy, :) = reaction_matrix(ix, iy, :) + reshape(concentrationChange, 1,1,[]);     % [molS/h]
                    end
                end
            end
        end
    end
end
    

% end

% function [spcM, Sh] = solve_pH(Sh_ini, StV, Keq, chrM, keepConstantpH, Tol)
%     % Solve the pH and speciation per grid cell
%     %
%     % <>
%     % original as copied from latest version, to be improved?
%     % </>
%     %
%     % -> spcM: species matrix for the respective pH
%     % -> Sh: proton concentration
%     % -> flagpH: flagpH = 1 -> Fix pH; flagpH = 0 -> Variable pH 
%     
%     w = 1;
% 
%     if keepConstantpH
%         spcM = zeros(size(chrM));
%         Sh = Sh_ini;
%         Denm = (1 + Keq(:, 1) / w) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);
% 
%         spcM(:, 1) = ((Keq(:, 1) / w) .* StV * Sh^3) ./ Denm;
%         spcM(:, 2) = (StV * Sh^3) ./ Denm;
%         spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
%         spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
%         spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;
% 
%     else
%         % Newton-Raphson method
%         Sh = Sh_ini;
%         % Counter of convergences
%         ipH = 1; 
%         maxIter = 20;
%         err = 1; % initial error
%         F = 1; % initial electroneutrality error
% 
%         % Inicialization of matrix of species
%         spcM = zeros(size(chrM)); dspcM = spcM;
% 
%         while (abs(err) > Tol) && (abs(F) > Tol) && ipH <= maxIter
%             Denm = (1 + Keq(:, 1)) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);
%             spcM(:, 1) = ((Keq(:, 1)) .* StV * Sh^3) ./ Denm;
%             spcM(:, 2) = (StV .* Sh^3) ./ Denm;
%             spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
%             spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
%             spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;
% 
%             % Evaluation of the charge balance for the current Sh value, F(Sh)
%             F = Sh + sum(spcM .* chrM, 'all');
% 
%             % Calculation of all derivated functions
%             dDenm = Denm.^2;
%             aux = 3 * Sh^2 * (Keq(:, 1) + 1) + 2 * Sh * Keq(:, 2) + Keq(:, 2) .* Keq(:, 3);
% 
%             dspcM(:, 1) = (3 * Sh^2 * Keq(:, 1) .* StV) ./ (Denm) - ((Keq(:, 1) .* StV * Sh^3) .* aux) ./ (dDenm);
%             dspcM(:, 2) = (3 * Sh^2 * StV) ./ Denm - (StV * Sh^3 .* aux) ./ dDenm;
%             dspcM(:, 3) = (2 * Sh * Keq(:, 2) .* StV) ./ Denm - ((Keq(:, 2) .* StV * Sh^2) .* aux) ./ dDenm;
%             dspcM(:, 4) = (Keq(:, 2) .* Keq(:, 3) .* StV) ./ Denm - ((Keq(:, 2) .* Keq(:, 3) .* StV * Sh) .* aux) ./ dDenm;
%             dspcM(:, 5) = -(Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4) .* StV .* aux) ./ dDenm;
% 
%             % Evaluation of the charge balance for the current Sh value, dF(Sh)
%             dF = 1 + sum(dspcM .* chrM, 'all');
%             %Error
%             err = F / dF;
%             % Newton-Raphson algorithm
%             Sh = max(Sh - err, 1e-10);
% 
%             ipH = ipH + 1;
%         end
%         if any(spcM < 0)
%             warning('DEBUG:actionRequired', 'debug: negative concentration encountered after pH calculation...');
%         end
%         if (Sh < 1e-14) || (Sh > 1)
%             warning('DEBUG:actionRequired', 'debug: pH found outside of 1-14 range...');
%         end
%     end
% end

function [mu_max, maint] = determine_max_growth_rate_and_maint(species, T, Sh)
    % Determine the maximum growth rate and maintenance [h-1] for a
    % specific species under certain conditions
    %
    % species: species of bacterium
    % T: temperature (in Kelvin)
    % Sh: concentration of protons
    %
    % -> mu_max: max growth rate
    % -> maint: maintenance requirement

    switch species
        case 1 % AOB
            mu_max = ((1.28*10^(12)*exp(-8183/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
            maint = (1.651*10^(11)*exp(-8183/T))/24;

        case 2 % NOB; Nitrobacter
            mu_max = ((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
            maint = (8.626*10^(6)*exp(-5295/T))/24;

        case 3 % NOB; Nitrospira
            mu_max = 0.63 * ((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
            maint = 0.63 * (8.626*10^(6)*exp(-5295/T))/24;

        case 4 % AMX; Brocadia spp [Puyol 2014]
            mu_max = 1.89*10^(8)*exp(-7330/T);
            maint = 0.05 * mu_max;
            
        otherwise
            error('Bacterial species not implemented: %d', species);
    end
end

function M = calculate_monod(Ks, Ki, conc)
    % Calculate the Monod coefficient given the Ks's, Ki's and the
    % corresponding concentrations
    %
    % Ks: vector with all Ks values (1-by-n)
    % Ki: vector with all Ki values (1-by-n)
    % conc: vector with all concentrations corresponding to the Ks and Ki
    %   values (1-by-n)
    %
    % -> M: Monod block
    
    % apply Ks
    M = prod((conc + 1e-25) ./ (conc + Ks + 1e-25)); % + 1e-25 to prevent NaN when conc == 0 and Ks == 0
    
    % apply Ki
    for i = 1:length(Ki)
        if Ki(i) ~= 0
            M = M * Ki(i) / (Ki(i) + conc(i));
        end
    end
end
