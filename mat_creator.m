function mat_creator(num_replicates)
    % Creates .mat files & registrate them in "simulation_log.md"
    % > num_replicates: number of replicates for each simulation
    
    %% argument validation
    arguments
        num_replicates (1,1) {mustBeInteger, mustBeInRange(num_replicates, 1, 9999)}
    end
    
    addpath(genpath('lib')); % make every subfolder with functions accessible to the code
    
    %% IbM version
    version = 'v2.3.0';
    
    %% Last simulation number
    fileparts(which('readwrite_md.py'));
    matID = double(py.readwrite_md.readmd());
    
    %% Creation & registration of simulations 
    excels = dir('planning\Excels\*.xlsx');
    if size(excels,1) == 0
        disp('>> Error 404: Excel files not found.')
    else
        for excel = excels'
            for r = 1:num_replicates
                fprintf('\n')
                disp(['>> Excel: ' excel.name ' (Run ' num2str(r) ')'])
                matID = matID + 1;
                [simname, siminfo, simgoal] = mc(excel.name, matID);
                clear R
                py.readwrite_md.writemd(matID, simname, siminfo, simgoal, version, r);
            end
        end
    end
end

function [simname, siminfo, simgoal] = mc(filename, sim_number)
    %{
    Mimicing the integTime from before to see if everything is working.
    Manually setting the structs for bac, grid, constants
    %}

    mat_name = sprintf('sim_%04d.mat', sim_number);
    folder_name = 'planning';
    path = ['planning\Excels\' filename];

    R = loadModelXlsx(path);

    warning('on', 'DEBUG:noActionRequired');
    warning('on', 'DEBUG:actionRequired');

    simname = R.info.simname;
    siminfo = R.info.siminfo;
    simgoal = R.info.simgoal;
    
    settings = struct;
    settings.parallelized = false;
    settings.structure_model = R.settings.structure_model;
    if settings.structure_model
        settings.type = R.settings.type; % {'Neut': Neutralism, 'Comp': Competition, 'Comm': Commensalism, 'Copr': Co-protection}
    end
    settings.pHincluded = R.settings.pHincluded;

    %% Manually creating the structs required for integTime (i.e. loadExcel mimic)
    rng(2021);

    grid = struct;
    grid.dx = R.Sxy.dx;
    grid.dy = R.Sxy.dy;
    grid.nX = R.Sxy.nx; % for now let's not worry about adding 1 gridcell at the beginning and end of the simulation domain...
    grid.nY = R.Sxy.ny;
    grid.blayer_thickness = R.Sxy.T_blayer;

    bac = struct;
    bac.x = R.bac.atrib(:,1);
    bac.y = R.bac.atrib(:,2);
    n = length(bac.x);
    % n = 5000; radius = R.Sxy.dx * 19; % arbitrary %500 bacs => 4, 5000 bacs => 19
    % [bac.x, bac.y] = rand_circle(n, grid.nX/2*grid.dx, grid.nY/2*grid.dy, radius);
    bac.species = randi(3, size(bac.x)); % random for now
    % bac.species = R.bac.atrib(:,5);
    bac.molarMass = R.bac.atrib(1,3) * ones(n, 1);
    bac.radius = R.bac.atrib(1,6) * ones(n, 1);
    bac.active = ones(size(bac.x), 'logical');
    clear n radius;

    % % retain only a handful of bacteria for faster debugging
    % retain = rand(size(bac.x)) < 0.2;
    % bac.x = bac.x(retain);
    % bac.y = bac.y(retain);
    % bac.species = bac.species(retain);
    % bac.molarMass = bac.molarMass(retain);
    % bac.radius = bac.radius(retain);
    % bac.active = bac.active(retain);

    constants = struct;
    constants.constantpH = R.settings.constantpH;
    constants.pHsetpoint = R.pOp.pH;
    constants.T = R.pOp.T;
    constants.isLiquid = strcmp(R.St.Phase(1:5), 'L') | strcmp(R.St.Phase(1:5), 'P');
    constants.Keq = R.pTh.Keq;
    constants.chrM = R.pTh.chrM;
    constants.Dir_k = logical(R.Inf.Dir_k);
    constants.Vr = R.pOp.Vr;                    % still don't know what this variable means, or more: why it is max_nBacs volume? seems arbitrary
    constants.Vg = (grid.dx ^ 3) * 1000; % L
    constants.StNames = R.St.StNames(1:R.St.numStVLiq2);
    constants.speciesNames = R.rm.rNamesX'; %{'B1', 'B2', 'B3'};
    constants.Ks = R.pTh.Ks; %R.pTh.Ks(:, 1:4);
    constants.Ki = R.pTh.Ki; %R.pTh.Ks(:, 5:7); % Ki > 0 only when co-protection is assumed
    constants.react_v = R.pTh.react_v;
    constants.MatrixMet = R.rm.MatrixMet;
    constants.MatrixDecay = R.rm.MatrixDecay_mod;
    constants.influent_concentrations = R.Inf.St;
    constants.pOp.NH3sp = R.pOp.SP; % Setpoint of substrate S in reactor
    constants.constantN = R.flagN;
    constants.kDist = 1.0;                         % Extra distance between bacteria, when kDist > 1.
    constants.max_granule_radius = R.bac.bac_ymax; % R.bac.bac_ymax;  % C: see excel          % [m] Maximum radius of granule. To compute the detachment of bacteria when bac_norm > r_max
    constants.dT = R.Sxy.dT;
    constants.dT_bac = R.Sxy.dT_bac;
    constants.dT_divide = R.Sxy.dT_Div;
    constants.dT_save = R.Sxy.dT_Print;
    constants.simulation_end = R.Sxy.maxT;
    constants.diffusion_rates = R.kTr.Diffn;
    constants.diffusion_accuracy = 1e-8; % to be tweaked still
    constants.Tol_a = 1e-12; % in [mol/m3], not [mol/L]! (prev. R.kTr.Tolabs)
    constants.pHtolerance = 1e-15;
    constants.steadystate_tolerance = 0.005; % [0, 1] -> relative/absolute tolerance of steady state
    % constants.correction_concentration_steady_state = 1e-4; % [mol/L]
    constants.correction_concentration_steady_state = R.kTr.Tolabs / constants.steadystate_tolerance; % [mol/L]
    constants.RESmethod = 'max'; % {'mean', 'max', 'norm'}
    constants.bac_MW = R.bac.bac_MW;
    constants.bac_rho = R.bac.bac_rho;
    constants.max_nBac = R.bac.bac_nmax; % <TODO: calculate with better method?>
    constants.inactivationEnabled = R.settings.inactivationBac;
    constants.min_bac_mass_grams = R.bac.bac_mmin;
    constants.max_bac_mass_grams = R.bac.bac_mmax;
    constants.bac_max_radius = R.bac.bac_rmax;
    constants.convergence_accuracy = 1e-6; % should be around 1e-6... the maximum difference between absolute maximum RES values of two diffusion cycles before it is said to no longer converge [RES% / 100]
    constants.nDiffusion_per_SScheck = 2;

    Neumann = 0.2;
    constants.dT = min(grid.dx^2./constants.diffusion_rates * Neumann);

    constants.debug.plotBacteria = false;
    constants.debug.plotConvergence = false; %
    constants.debug.plotMaxErrors = false; %
    constants.debug.plotDiffRegion = false;
    constants.debug.plotBulkConcsOverTime = false; %
    constants.debug.plotProfiling = false; %

    init_params = struct;
    init_params.init_bulk_conc = R.Sxy.Sbc_Dir;
    init_params.init_concs = R.St.StVLiq;
    init_params.invHRT = R.pOp.invHRT;

    %% actual call to integTime
    % directory = 'Testing';
    % constants.simulation_end = 24*7*3; % 3 weeks

    constants.dT_backup = 7*24;

    %{
    bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
    bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
    bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
    bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
    %}
    if constants.debug.plotBacteria
        plotBacs(grid, bac, constants, 0)
    end

    % temp settings of variables: pls remove
    constants.dynamicDT.iterThreshholdDecrease = 500;
    constants.dynamicDT.iterThreshholdIncrease = 50;
    constants.dynamicDT.nIterThreshholdIncrease = 3;
    settings.dynamicDT = true;
    
    save(fullfile(folder_name, mat_name), 'grid', 'bac', 'constants', 'init_params', 'settings', '-v7.3')
    
end