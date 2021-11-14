function [reaction_matrix, mu, pH] = par_calculate_reaction_matrix(grid2bac, grid2nBacs, bac, diffRegion, conc, constants, pH, chunks, nChunks_dir, settings)
    % Calculate how much compounds are consumed/produced per grid cell due
    % to bacterial activity. Also updates the growth rate of the respective
    % bacteria.
    % 
    % grid2bac: matrix with per grid cell which bacteria reside there
    % grid2nBacs: matrix with how many bacteria per grid cell
    % bac: struct containing all information regarding the bacteria
    % grid: struct containing all information regarding the grid
    % conc: matrix containing all concentrations per grid cell as of
    %   (ix, iy, compound)
    % constants: struct containing all simulation constants
    % pH_old: matrix with per grid cell the pH
    % chunks: struct containing the info regarding the chunks for
    %   parallelisation
    % nChunks_dir: the number of chunks per direction
    % settings: struct containing simulation wide settings
    %
    % -> reaction_matrix: matrix with per grid cell and per compound the
    %   change [h-1] due to bacterial activity
    % -> mu: vector with updated growth rates per bacterium
    % -> pH: matrix with per grid cell the pH
    
    structure_model = settings.structure_model;
    pHincluded = settings.pHincluded;

    nChunks = nChunks_dir^2;
    
    % convert init pH to pH matrix
    if isscalar(pH)
        pH = ones(size(grid2nBacs))*pH;
    end
    
    Keq = constants.Keq;
    chrM = constants.chrM;
    Vg = constants.Vg;
    compounds = constants.StNames;
    reactive_form = constants.react_v;
    Ks = constants.Ks;
    Ki = constants.Ki;
    mMetabolism = constants.MatrixMet;
    mDecay = constants.MatrixDecay;
    
    if structure_model
        iA = find(strcmp(compounds, 'A'));
        iB = find(strcmp(compounds, 'B'));
        iC = find(strcmp(compounds, 'C'));
        iO2 = find(strcmp(compounds, 'O2'));        
    else
        iNH3 = find(strcmp(compounds, 'NH3'));
        iNO2 = find(strcmp(compounds, 'NO2'));
        iO2 = find(strcmp(compounds, 'O2'));
    end
    
    reaction_matrix = zeros(size(conc));
    mu = zeros(size(bac.x));
    
    % pre-compute for bulk-liquid 
    % -> at 1,1 there should never be diffusion layer
    Sh_bulk = 10^-pH(1, 1);
    
    if pHincluded
        [~, Sh_bulk] = solve_pH(Sh_bulk, [reshape(conc(1,1,:), [], 1, 1); 1; 0], Keq, chrM, constants.constantpH, constants.pHtolerance); % <C: why [...; 1; 0]? />
        pH_bulk = -log10(Sh_bulk);
    else
        pH_bulk = pH(1, 1);
    end

    % group constants for easy passing to multiple cores
    if structure_model
        constantValues = [pH_bulk, constants.constantpH, constants.pHtolerance, constants.T, iA, iB, iC, iO2];
    else
        constantValues = [pH_bulk, constants.constantpH, constants.pHtolerance, constants.T, iNH3, iNO2, iO2];
    end
    grouped_bac = [bac.species, bac.molarMass];

    if settings.parallelized
        % ================ PARALLEL CALCULATION ====================
        
        % create cell arrays with input variables per chunk
        ix = @(iChunk) floor((iChunk - 1) / nChunks_dir) + 1;
        iy = @(iChunk) mod(iChunk - 1, nChunks_dir) + 1;
        chunk_xRange = arrayfun(@(iChunk) chunks.indices_x(ix(iChunk), 1):chunks.indices_x(ix(iChunk), 2), 1:nChunks, 'UniformOutput', false);
        chunk_yRange = arrayfun(@(iChunk) chunks.indices_y(iy(iChunk), 1):chunks.indices_y(iy(iChunk), 2), 1:nChunks, 'UniformOutput', false);
        chunk_pH_OG = cellfun(@(xR, yR) pH(xR, yR), chunk_xRange, chunk_yRange, 'UniformOutput', false);
        chunk_conc = cellfun(@(xR, yR) conc(xR, yR, :), chunk_xRange, chunk_yRange, 'UniformOutput', false);
        chunk_grid2bac = cellfun(@(xR, yR) grid2bac(xR, yR, :), chunk_xRange, chunk_yRange, 'UniformOutput', false);
        chunk_grid2nBacs = cellfun(@(xR, yR) grid2nBacs(xR, yR), chunk_xRange, chunk_yRange, 'UniformOutput', false);
        chunk_diffRegion = cellfun(@(xR, yR) diffRegion(xR, yR), chunk_xRange, chunk_yRange, 'UniformOutput', false);

        % chunk2nBacs: number of bacteria per chunk
        chunk2nBacs = zeros(nChunks, 1);
        for iChunk = 1:nChunks
            xRange = chunk_xRange{iChunk};
            yRange = chunk_yRange{iChunk};
            chunk2nBacs(iChunk) = sum(grid2nBacs(xRange,yRange), 'all');
        end

        % bacOffset
        bacOffset = zeros(nChunks + 1, 1);
        bacOffset(2:end) = cumsum(chunk2nBacs);

        chunk_bacRange = arrayfun(@(iChunk) bacOffset(iChunk)+1:bacOffset(iChunk+1), 1:nChunks, 'UniformOutput', false);
        chunk_grouped_bac = cellfun(@(bacRange) grouped_bac(bacRange, :), chunk_bacRange, 'UniformOutput', false);

        % create storage variables for output per chunk
        nCompounds = length(compounds);
        xlens = kron(ones(nChunks_dir,1), diff(chunks.indices_x, [], 2)+1);
        ylens = kron(diff(chunks.indices_y, [], 2)+1, ones(nChunks_dir, 1));

        chunk_rMatrix = arrayfun(@(xlen, ylen) zeros(xlen, ylen, nCompounds), xlens, ylens, 'UniformOutput', false);
        chunk_pH = arrayfun(@(xlen, ylen) zeros(xlen, ylen), xlens, ylens, 'UniformOutput', false);
        chunk_mu = arrayfun(@(x) zeros(x, 1), chunk2nBacs, 'UniformOutput', false);


        % calculate reaction matrix per chunk in parallel
        parfor iChunk = 1:nChunks
            [chunk_rMatrix{iChunk}, chunk_mu{iChunk}, chunk_pH{iChunk}] = ...
                rMatrix_chunk(chunk_pH_OG{iChunk}, chunk_conc{iChunk}, chunk_grid2bac{iChunk}, chunk_grid2nBacs{iChunk}, chunk_diffRegion{iChunk}, ...
                chunk_grouped_bac{iChunk}, chunk2nBacs(iChunk), bacOffset(iChunk), ...
                reactive_form, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constantValues, structure_model, pHincluded);
        end

        % put everything back into correct matrix/vector
        pH = ones(size(diffRegion))*pH_bulk; % everything outside of diffRegion
        for iChunk = 1:nChunks
            xRange = chunk_xRange{iChunk};
            yRange = chunk_yRange{iChunk};
            bacRange = chunk_bacRange{iChunk};

            reaction_matrix(xRange, yRange, :) = chunk_rMatrix{iChunk};
            pH(xRange, yRange) = chunk_pH{iChunk};
            mu(bacRange) = chunk_mu{iChunk};
        end
    
    else
        % ================ SEQUENTIAL CALCULATION ====================
        bacOffset = 0;
        % calculate by using just one chunk
        [reaction_matrix, mu, pH] = ...
            rMatrix_chunk(pH, conc, grid2bac, grid2nBacs, diffRegion, ...
            grouped_bac, length(bac.x), bacOffset, ...
            reactive_form, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constantValues, structure_model, pHincluded);
    end
    
    
    % do final unit correction/calculation
    reaction_matrix = reaction_matrix / Vg;
end


