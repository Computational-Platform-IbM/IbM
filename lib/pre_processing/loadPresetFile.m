function [grid, bac_init, constants, settings, init_params] = loadPresetFile(filename)
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
    constants = struct;
    settings = struct;
    init_params = struct;
    bac_init = struct;

    
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

    
    %% constants (Diffusion)
    [vals, names] = xlsread(filename, 'Diffusion');
    constants.compoundNames = names(:,1);
    nCompounds = length(constants.compoundNames);
    constants.diffusion_rates = vals;                                       % [m2/h]
    
    
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
        constants.setpoint_index = find(strcmp(constants.compoundNames, compound_name));
    end
        
    
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
    else
        constants.pHtolerance = nan;
    end
    
    % if pH is variable, speciation is always set to true. Otherwise, read
    % from excel
    settings.speciation = logical(vals(strcmp(names, 'Speciation included'))) || settings.pHincluded;

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
%     init_params.init_bulk_conc = constants.influent_concentrations; % TODO: keep original?
    
    
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
    preferred_state = vals(dG_rows-1, end);

    % extract equilibrium constants
    Keq_rows = section_starts(2):section_ends(2);
    Keq_rows = [Keq_rows, H2O_index(2), H_index(2)];
    
    constants.Keq = vals(Keq_rows-1, 1:nColumns-1); % dG_rows - 1, because first row in names is for header; nColumns - 1, because always 1 less Keq than subspecies
    
    % extract charge matrix
    charge_rows = section_starts(3):section_ends(3);
    assert(areEqual(constants.compoundNames, names(charge_rows,1)), 'ERROR: Compounds do not have the same name, or are not in the same order.')

    charge_rows = [charge_rows, H2O_index(3), H_index(3)];
    constants.chrM = vals(charge_rows-1, 1:nColumns);
    constants.chrM(isnan(constants.chrM)) = 0;
    
    
    %% constants (Ks & Ki)
    % create Ks matrix (species * compounds)
    [vals_ks, names_ks] = xlsread(filename, 'Ks');
    [vals_ki, names_ki] = xlsread(filename, 'Ki');

    constants.speciesNames = names_ks(2:end, 1);
    assert(areEqual(constants.speciesNames, names_ki(2:end, 1)), 'ERROR: Bacterial species do not have the same name, or are not in the same order.');

    % check which unique compounds in Ks & Ki compounds
    comp_names = {names_ks(1,2:end-1), names_ki(1, 2:end-1)};
    uniq_compounds = unique(cat(2, comp_names{:}));
    
    % create Ks and Ki matrix
    constants.Ks = zeros(length(constants.speciesNames), length(uniq_compounds));
    constants.Ki = zeros(length(constants.speciesNames), length(uniq_compounds));
    
    for i = 1:length(uniq_compounds)
        compoundName = uniq_compounds{i};
        iColumn = find(strcmp(names_ks(1,:), compoundName));
        if iColumn
            constants.Ks(:, i) = vals_ks(:, iColumn-1); % -1, because first column is header
        end
        iColumn = find(strcmp(names_ki(1,:), compoundName));
        if iColumn
            constants.Ki(:, i) = vals_ki(:, iColumn-1); % -1, because first column is header
        end
    end
    
    % create reactive_indices
    constants.reactive_indices = zeros(length(uniq_compounds), 1);

    if settings.speciation
        sz = [nCompounds+2, nColumns]; % spcM also includes rows for H2O and H

        for i = 1:length(uniq_compounds)
            I = find(strcmp(constants.compoundNames, uniq_compounds{i}));
            J = preferred_state(I);
            constants.reactive_indices(i) = sub2ind(sz, I, J);
        end
    else
        for i = 1:length(uniq_compounds)
            constants.reactive_indices(i) = find(strcmp(constants.compoundNames, uniq_compounds{i}));
        end
    end

    
    %% constants (Yield, eDonor, mu_max, maint)
    [vals, names] = xlsread(filename, 'Yield');
    
    assert(areEqual(constants.speciesNames, names(2:end,1)), 'ERROR: Bacterial species do not have the same name, or are not in the same order.');
    yield = vals(:, find(strcmp(names(1,:), 'Yield C/N')) - 1);
    eD = names(2:end, strcmp(names(1,:), 'eD'));
    try
        constants.maintenance = vals(:, find(strcmp(names(1,:), 'Maintenance')) - 1);
        constants.mu_max = vals(:, find(strcmp(names(1,:), 'Max growth rate')) - 1);
    catch e
        switch e.identifier
            case 'MATLAB:badsubscript'
                constants.maintenance = nan;
                constants.mu_max = nan;
                fprintf('Maintenance and maximum growth rate are not set, thus calculating dynamically. \nPlease make sure the equations and species match up in the code.\n')
            otherwise
                rethrow(e)
        end
    end
    
    
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
    constants.MatrixMet = zeros(nCompounds, length(constants.speciesNames));
    constants.MatrixDecay = zeros(nCompounds, length(constants.speciesNames));
    
    % loop over species
    for species = 1:length(constants.speciesNames)
        % get catabolism & anabolism
        cata = vals(compound_rows, iCat(species));
        ana = vals(compound_rows, iAnab(species));
        
        % get eDonor & yield
        eD_species = eD{species};
        Y = yield(species);
        
        % check whether eD has value -1 in cata
        eD_index = strcmp(constants.compoundNames, eD_species);
        assert(cata(eD_index) == -1, 'ERROR, catabolism is not normalised towards the electron donor.');
        
        % calculate metabolism matrix entry
        constants.MatrixMet(:, species) = cata/Y + ana;
        
        % set decay matrix entry
        constants.MatrixDecay(:, species) = vals(compound_rows, iDecay(species));
    end
    
    
    %% initialisation of bacteria
    [vals, names] = xlsread(filename, 'Bacteria');
    bac_init.method = names{strcmp(names, 'Initialisation method'), 2};
    switch (bac_init.method)
        case 'granule'
            bac_init.granule_radius = vals(strcmp(names, 'Starting granule radius'));
            bac_init.start_nBac = vals(strcmp(names, 'Starting number of bacteria (granule)'));
        case 'suspension'
            bac_init.start_nBac = vals(strcmp(names, 'Starting number of bacteria (suspension)'));
        otherwise
            error(['Initialisation method <', bac_init.method,'> is not a valid method.'])
    end
        
end
    
    
    
    
    
    
    
    

% check all settings in model if they are applied correctly...
% check pH algorithms with new structs... is Kd now working?
% think about what to do with mu_max?
% think about what to do with maintenance?


    



% end

function b = areEqual(cellarray1, cellarray2)
    % compare two cell arrays with strings and return whether they are identical
    b = all(cellfun(@strcmp, cellarray1, cellarray2));
end