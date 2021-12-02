filename = './planning/Excels/Templates/AOBNOBAMXCMX_template.xlsx';
% function [grid, bac, constants, settings, init_params] = loadPresetFile(filename)
    % Generate the structs for running the IbM model from the provided
    % Excel file.
    %
    % filename: filename (with correct path) to the preset Excel file
    %
    % -> grid: struct containing all information regarding spacial
    %   discretisation
    % -> bac: struct containing all information regarding the bacteria
    % -> constants: struct containing all constants
    % -> settings: struct with settings for switching between different
    %   methods within the IbM model
    % -> init_params: struct containing the initial values from which the
    %   model will run
    
    %% initialisation of structs
    grid = struct;
    bac = struct;
    constants = struct;
    settings = struct;
    init_params = struct;

    
    %% grid
    [vals, names] = xlsread(filename, 'Discretization');
    grid.dx = vals(strcmp(names, 'dx'));                                    % [m]
    grid.dy = vals(strcmp(names, 'dy'));                                    % [m]
    grid.nX = vals(strcmp(names, 'nx'));                                    % [m]
    grid.nY = vals(strcmp(names, 'ny'));                                    % [m]
    grid.blayer_thickness = vals(strcmp(names, 'Boundary layer thickness'));% [m]
    constants.Vg = (grid.dx ^ 3) * 1000;                                    % [L]
    constants.max_granule_radius = (grid.dx * (grid.nX - 4)) / 2;           % [m]
    
    
    %% constants (Time)
    constants.simulation_end = vals(strcmp(names, 'Simulation end'));       % [h]
    constants.dT = vals(strcmp(names, 'Initial dT diffusion'));             % [h]
    constants.dT_bac = vals(strcmp(names, 'Initial dT bacteria'));          % [h]
    constants.dT_save = vals(strcmp(names, 'dT save'));                     % [h]
    constants.dT_backup = vals(strcmp(names, 'dT backup'));                 % [h]
    
    settings.dynamicDT = logical(vals(strcmp(names, 'Dynamic dT')));
    
    if settings.dynamicDT
        constants.dynamicDT.nIterThresholdIncrease = vals(strcmp(names, 'nIterThreshold'));
        constants.dynamicDT.iterThresholdDecrease = vals(strcmp(names, 'iterThresholdDecrease'));
        constants.dynamicDT.iterThresholdIncrease = vals(strcmp(names, 'iterThresholdIncrease'));
        constants.dynamicDT.initRESThresholdIncrease = vals(strcmp(names, 'initial RES threshold increase'));
        constants.dynamicDT.nItersCycle = vals(strcmp(names, 'nIters per cycle'));
        constants.dynamicDT.tolerance_no_convergence = vals(strcmp(names, 'tolerance no-convergence')); 
        constants.dynamicDT.maxRelDiffBulkConc = vals(strcmp(names, 'maximum relative bulk conc change'));
        
        constants.dynamicDT.maxDT = vals(strcmp(names, 'Maximum dT diffusion'));    % [h]
        constants.dynamicDT.minDT = vals(strcmp(names, 'Minimum dT diffusion'));    % [h]
        constants.dynamicDT.maxDT_bac = vals(strcmp(names, 'Maximum dT bacteria')); % [h]
        constants.dynamicDT.minDT_bac = vals(strcmp(names, 'Minimum dT bacteria')); % [h]
    end

    
    %% constants (Operational parameters)
    [vals, names] = xlsread(filename, 'Parameters');
    constants.pHsetpoint = vals(strcmp(names, 'pH setpoint'));              % [-]
    constants.T = vals(strcmp(names, 'Temperature')) + 273.15;              % [K]
    constants.Vr = vals(strcmp(names, 'Representative volume')) * 1000;     % [L]

    settings.variableHRT = vals(strcmp(names, 'Variable HRT')); % -> changed from constants.constantN
    init_params.invHRT = 1/vals(strcmp(names, 'HRT'));                      % [1/h]

    if settings.variableHRT
        constants.bulk_setpoint = vals(strcmp(names, 'Setpoint'));          % [mol/L]    % -> changed from constants.pOp.NH3sp
        compound_name = names{strcmp(names, 'Compound setpoint'), 2};
%         constants.setpoint_compound_index = find(strcmp(%%%%%%%%%%%%%%%%%)
    end
    
    % load some constant values for use later in preset file generator
    
    
    %% constants (Diffusion)
    [vals, names] = xlsread(filename, 'Diffusion');
    constants.compoundNames = names(:,1);
    constants.diffusion_rates = vals;                                       % [m2/h]
    
    
    %% constants (Bacteria)
    [vals, names] = xlsread(filename, 'Bacteria');
    constants.bac_MW = vals(strcmp(names, 'Molecular weight bacterium'));                   % [g/mol]
    constants.bac_rho = vals(strcmp(names, 'Density bacterium'));                           % [g/m3]
    constants.max_nBac = vals(strcmp(names, 'Maximum nBacteria')); % <TODO: calculate with better method?>  % [-]
    constants.inactivationEnabled = logical(vals(strcmp(names, 'Inactivation enabled')));   % [-]
    constants.min_bac_mass_grams = vals(strcmp(names, 'Minimum mass bacterium'));           % [g]
    constants.max_bac_mass_grams = vals(strcmp(names, 'Maximum mass bacterium'));           % [g]
    constants.bac_max_radius = vals(strcmp(names, 'Maximum radius bacterium'));             % [m]
    constants.kDist = vals(strcmp(names, 'kDist'));                                         % [-]
    settings.detachment = names{strcmp(names, 'Detachment method'), 2};                     % {mechanistic, naive, none}
    constants.kDet = vals(strcmp(names, 'Detachment constant'));                            % [1/m2.h]
    constants.max_granule_radius = min(constants.max_granule_radius, vals(strcmp(names, 'Maximum granule radius'))); % [m]
    
    
    %% constants (Solver)
    [vals, names] = xlsread(filename, 'Solver');
    
    constants.diffusion_accuracy = vals(strcmp(names, 'Diffusion tolerance')); 
    constants.steadystate_tolerance = vals(strcmp(names, 'Steady state RES threshold'));
    tol_abs = vals(strcmp(names, 'Concentration tolerance'));
    constants.correction_concentration_steady_state = tol_abs / constants.steadystate_tolerance; % [mol/L]
    constants.RESmethod = names{strcmp(names, 'RES determination method'), 2};

    constants.nDiffusion_per_SScheck = vals(strcmp(names, 'nIters diffusion per SS check'));
    
    settings.pHbulkCorrection = logical(vals(strcmp(names, 'pH bulk concentration corrected')));
    settings.pHincluded = logical(vals(strcmp(names, 'pH solving included')));
    if settings.pHincluded
        constants.pHtolerance = vals(strcmp(names, 'pH solver tolerance'));
    end
    settings.speciation = logical(vals(strcmp(names, 'Speciation included')));

    settings.structure_model = logical(vals(strcmp(names, 'Structure model')));
    if settings.structure_model
        settings.type = names{strcmp(names, 'Structure model type'), 2};
    end
    
    settings.parallelized = logical(vals(strcmp(names, 'Parallelisation')));

    
    %% Constants (Influent & initial condition)
    % influent concentrations
    [vals, names] = xlsread(filename, 'Influent');
    
    assert(areEqual(constants.compoundNames, names(:,1)), 'ERROR: Compounds do not have the same name, or are not in the same order.')
    
    constants.Dir_k = strcmp(names(:,4), 'D');
    constants.influent_concentrations = vals;                               % [mol/L]

    % initial conditions
    [vals, names] = xlsread(filename, 'Initial condition');
    
    assert(areEqual(constants.compoundNames, names(:,1)), 'ERROR: Compounds do not have the same name, or are not in the same order.')

    init_params.init_concs = vals;                                          % [mol/L]
    init_params.init_bulk_conc = init_params.init_concs;                    % [mol/L]
    init_params.init_bulk_conc(constants.Dir_k) = constants.influent_concentrations(constants.Dir_k);
    
    
    %% constants (Equilibrium constants & charge matrix)
    [vals, names] = xlsread(filename, 'ThermoParam');
    nColumns = size(vals, 2) - 1; % one column required for preferred subspecies
    
    section_starts = find(strcmp(constants.compoundNames{1}, names(:,1)));
    section_ends = find(strcmp(constants.compoundNames{end}, names(:,1)));
    H2O_index = find(strcmp(names(:,1), 'H2O'));
    H_index = find(strcmp(names(:,1), 'H'));
    
    dG_rows = section_starts(1):section_ends(1);
    assert(areEqual(constants.compoundNames, names(dG_rows,1)), 'ERROR: Compounds do not have the same name, or are not in the same order.')
    dG_rows = [dG_rows, H2O_index(1), H_index(1)];
    
    % extract preferred state per compound
    constants.preferred_state = vals(dG_rows-1, end);

    % extract equilibrium constants
    Keq_rows = section_starts(2):section_ends(2);
    Keq_rows = [Keq_rows, H2O_index(2), H_index(2)];
    
    constants.Keq = vals(Keq_rows-1, 1:nColumns-1); % dG_rows - 1, because first row in names is for header; nColumns - 1, because always 1 less Keq than subspecies
    
    
    
    if settings.pHincluded && settings.speciation
        % add all dG values
    else 
        % only add dG values for preferred state
    end
    
    
    % extract charge matrix
    charge_rows = section_starts(3):section_ends(3);
    assert(areEqual(constants.compoundNames, names(charge_rows,1)), 'ERROR: Compounds do not have the same name, or are not in the same order.')

    charge_rows = [charge_rows, H2O_index(3), H_index(3)];
    constants.chrM = vals(charge_rows-1, 1:nColumns);
    constants.chrM(isnan(constants.chrM)) = 0;
    
    
    %% constants (Ks & Ki)
    % create Ks matrix (species * compounds)
    [vals, names] = xlsread(filename, 'Ks');
    
    constants.speciesNames = names(2:end, 1);
    constants.Ks = zeros(length(constants.speciesNames), length(constants.compoundNames));
    for i = 1:length(constants.compoundNames)
        compoundName = constants.compoundNames{i};
        iColumn = find(strcmp(names(1,:), compoundName));
        if iColumn
            constants.Ks(:, i) = vals(:, iColumn-1); % -1, because first column is header
        end
    end
    
    % create Ks matrix (species * compounds)
    [vals, names] = xlsread(filename, 'Ki');
    
    assert(areEqual(constants.speciesNames, names(2:end, 1)), 'ERROR: Bacterial species do not have the same name, or are not in the same order.');
    constants.Ki = zeros(length(constants.speciesNames), length(constants.compoundNames));
    for i = 1:length(constants.compoundNames)
        compoundName = constants.compoundNames{i};
        iColumn = find(strcmp(names(1,:), compoundName));
        if iColumn
            constants.Ki(:, i) = vals(:, iColumn-1); % -1, because first column is header
        end
    end
    
    
    %% constants (Yield, eDonor, mu_max, maint)
    [vals, names] = xlsread(filename, 'Yield');
    
    assert(areEqual(constants.speciesNames, names(2:end,1)), 'ERROR: Bacterial species do not have the same name, or are not in the same order.');
    yield = vals(:, find(strcmp(names(1,:), 'Yield C/N')) - 1);
    eD = names(2:end, strcmp(names(1,:), 'eD'));
    constants.maintenance = vals(:, find(strcmp(names(1,:), 'Maintenance')) - 1);
    constants.mu_max = vals(:, find(strcmp(names(1,:), 'Max growth rate')) - 1);
    
    
    
    %% constants (ReactionMatrix)
    [vals, names] = xlsread(filename, 'ReactionMatrix');
    assert(areEqual(constants.speciesNames, names(1, 2:3:end)'), 'ERROR: Bacterial species do not have the same name, or are not in the same order.');
    compound_rows = find(strcmp(constants.compoundNames{1}, names(:, 1))):find(strcmp(constants.compoundNames{end}, names(:, 1)));
    assert(areEqual(constants.compoundNames, names(compound_rows,1)), 'ERROR: Compounds do not have the same name, or are not in the same order.')

    % find catabolism, anabolism & decay parts for each bacterium
    compound_rows = compound_rows - 2; % correct for two lines of headers
    iCat = find(strcmp(names(2,:), 'Cat')) - 1; % correct for one column of header
    iAnab = find(strcmp(names(2,:), 'Anab')) - 1; % idem
    iDecay = find(strcmp(names(2,:), 'Decay')) - 1; % idem
    
    % initialise matrix for metabolisms & decay
    constants.MatrixMet = zeros(length(constants.compoundNames), length(constants.speciesNames));
    constants.MatrixDecay = zeros(length(constants.compoundNames), length(constants.speciesNames));
    
    % loop over species
    for species = 1:length(constants.speciesNames)
        % get catabolism & anabolism
        cat = vals(compound_rows, iCat(species));
        ana = vals(compound_rows, iAnab(species));
        
        % get eDonor & yield
        eD_species = eD{species};
        Y = yield(species);
        
        % check whether eD has value -1 in cat
        eD_index = strcmp(constants.compoundNames, eD_species);
        assert(cat(eD_index) == -1, 'ERROR, catabolism is not normalised towards the electron donor.');
        
        % calculate metabolism matrix entry
        constants.MatrixMet(:, species) = cat/Y + ana;
        
        % set decay matrix entry
        constants.MatrixDecay(:, species) = vals(compound_rows, iDecay(species));
    end
    
    
    %% initialise bacteria
    % what initialisation method?
    % how big initial granule? / how many individuals?
    
    
    
    
    
    
    
    
    
settings.structure_model = false;
if settings.structure_model
    settings.type = 'Neut'; % {'Neut': Neutralism, 'Comp': Competition, 'Comm': Commensalism, 'Copr': Co-protection}
end


% constants.isLiquid = strcmp(R.St.Phase(1:8), 'L') | strcmp(R.St.Phase(1:8), 'P'); % should no longer be needed...
% constants.StNames = R.St.StNames(1:8); --> changed to constants.compoundNames
% check all settings in model if they are applied correctly...
% remove N2 from all bulk_conc, reaction_matrix, etc. in the model
% how to implement the Ks and Ki generally speaking with matrix instead of
% values?
% check pH algorithms with new structs... is Kd now working?
% constants.speciesNames = {'AOB', 'NOB (Nitrob.)', 'NOB (Nitrosp.)', 'AMX'};
% constants.react_v = R.pTh.react_v; % --> renamed to constants.prefered state
% think about what to do with mu_max?
% think about what to do with maintenance?
constants.constantpH = false; % check what this is used for? and put into excel then?


    



% end

function b = areEqual(cellarray1, cellarray2)
    % compare two cell arrays with strings and return whether they are identical
    b = all(cellfun(@strcmp, cellarray1, cellarray2));
end