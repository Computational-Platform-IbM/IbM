function [reaction_matrix, mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, diffRegion, conc, constants, pH, settings)
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
    %
    % -> reaction_matrix: matrix with per grid cell and per compound the
    %   change [h-1] due to bacterial activity
    % -> mu: vector with updated growth rates per bacterium
    % -> pH: matrix with per grid cell the pH

    structure_model = settings.structure_model;
    pHincluded = settings.pHincluded;
    
    % convert init pH to pH matrix
    if isscalar(pH)
        pH = ones(size(grid2nBacs))*pH;
    end
    
    % init
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
    
    % pre-compute for bulk-liquid (at 1,1 there should never be diffusion
    % layer)
    Sh_bulk = 10^-pH(1, 1);
    
    if pHincluded
        [~, Sh_bulk] = solve_pH(Sh_bulk, [reshape(conc(1,1,:), [], 1, 1); 1; 0], Keq, chrM, constants.constantpH, constants.pHtolerance); % <C: why [...; 1; 0]? /> => H2O & H+
        pH_bulk = -log10(Sh_bulk);
    else
        pH_bulk = pH(1, 1);
    end
    
    % prep for call to rMatrix_chunk
    grouped_bac = [bac.species, bac.molarMass];
    bacOffset = 0;
    
    if structure_model
        constantValues = [pH_bulk, constants.constantpH, constants.pHtolerance, constants.T, iA, iB, iC, iO2];
    else
        constantValues = [pH_bulk, constants.constantpH, constants.pHtolerance, constants.T, iNH3, iNO2, iO2];
    end
    
    % calculate by using just one chunk
    [reaction_matrix, mu, pH] = ...
            rMatrix_chunk(pH, conc, grid2bac, grid2nBacs, diffRegion, ...
            grouped_bac, length(bac.x), bacOffset, ...
            reactive_form, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constantValues, structure_model, pHincluded);
    
    reaction_matrix = reaction_matrix / Vg;                                                                                 % [molS/L/h]

end
