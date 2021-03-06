%{
Mimicing the integTime from before to see if everything is working.
Manually setting the structs for bac, grid, constants
%}

clearvars -except R;
clc;
addpath(genpath('lib'));
javaaddpath([pwd '\lib\bacteria\shovingQuadTreekDist.jar']);

if ~exist('R', 'var')
    R = loadModelXlsx('./planning/Excels/Templates/AOBNOBAMX_template.xlsx');
%     R = loadModelXlsx('Testing.xlsx');
    clc;
end

warning('on', 'DEBUG:noActionRequired');
warning('on', 'DEBUG:actionRequired');

%% Manually creating the structs required for integTime (i.e. loadExcel mimic)
grid = struct;
% grid.dx = R.Sxy.dx;
% grid.dy = R.Sxy.dy;
% grid.nX = R.Sxy.nx; % for now let's not worry about adding 1 gridcell at the beginning and end of the simulation domain...
% grid.nY = R.Sxy.ny;
grid.dx = R.Sxy.dx;
grid.dy = R.Sxy.dy;
grid.nX = 1025; % for now let's not worry about adding 1 gridcell at the beginning and end of the simulation domain...
grid.nY = 1025;
grid.blayer_thickness = R.Sxy.T_blayer;

bac = struct;

% ----- load from R-struct
% bac.x = R.bac.atrib(:,1);
% bac.y = R.bac.atrib(:,2);
% n = length(bac.x);

% ----- initiate circle with bacteria (random)
% n = 5000; radius = R.Sxy.dx * 19; % arbitrary %500 bacs => 4, 5000 bacs => 19
% [bac.x, bac.y] = rand_circle(n, grid.nX/2*grid.dx, grid.nY/2*grid.dy, radius);

% ----- initiate suspension (semi-random)
n = 500; margin = 0.1*grid.dx*grid.nX;
xRange = [margin, grid.dx*grid.nX - margin];
yRange = xRange; % assume square domain
[bac.x, bac.y] = distribute_rectangle(n, xRange, yRange);



bac.species = randi(4, size(bac.x)); % random for now
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
constants.pHsetpoint = R.pOp.pH;
constants.T = 293;
constants.isLiquid = strcmp(R.St.Phase(1:8), 'L') | strcmp(R.St.Phase(1:8), 'P');
constants.Keq = R.pTh.Keq;
constants.chrM = R.pTh.chrM;
constants.Dir_k = logical(R.Inf.Dir_k);
constants.Vr = R.pOp.Vr;                    % still don't know what this variable means, or more: why it is max_nBacs volume? seems arbitrary
constants.Vg = (grid.dx ^ 3) * 1000; % L
constants.StNames = R.St.StNames(1:8);
constants.speciesNames = {'AOB', 'NOB (Nitrob.)', 'NOB (Nitrosp.)', 'AMX'};
constants.react_v = R.pTh.react_v;
constants.Ks = R.pTh.Ks;
constants.Ki = R.pTh.Ki;
constants.MatrixMet = R.rm.MatrixMet;
constants.MatrixDecay = R.rm.MatrixDecay_mod;
constants.influent_concentrations = R.Inf.St;
constants.pOp.NH3sp = R.pOp.SP; % Setpoint of NH3 in reactor
constants.constantN = R.flagN;
constants.kDist = 1.5;                                    % Extra distance between bacteria, when kDist > 1.
constants.max_granule_radius = 1000*10^(-6);  % C: see excel          % [m] Maximum radius of granule. To compute the detachment of bacteria when bac_norm > r_max
constants.dT = R.Sxy.dT; % AOB/NOB/AMX -> 1e-6
constants.dT_bac = R.Sxy.dT_bac;
constants.dT_divide = R.Sxy.dT_Div;
constants.dT_save = 24;%R.Sxy.dT_Print;
constants.constantpH = false;
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
constants.inactivationEnabled = true;
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

settings = struct;
settings.parallelized = false;
settings.structure_model = false;
if settings.structure_model
    settings.type = 'Neut'; % {'Neut': Neutralism, 'Comp': Competition, 'Comm': Commensalism, 'Copr': Co-protection}
end
settings.pHincluded = true; % true -> pH resolution included; false -> pH resolution not included

init_params = struct;
init_params.init_bulk_conc = R.Sxy.Sbc_Dir;
init_params.init_concs = R.St.StVLiq;
init_params.invHRT = R.pOp.invHRT;

%% actual call to integTime
% directory = 'Testing';
% constants.simulation_end = 24*7*3; % 3 weeks

constants.dT_backup = 7*24;

bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...
bac = bacteria_shove(bac, grid, constants); % otherwise bacteria might overlap at start...

% if constants.debug.plotBacteria
    plotBacs(grid, bac, constants, 0)
% end

% overall settings dynamic dT
settings.dynamicDT = true;
constants.dynamicDT.nIterThresholdIncrease = 3;

% dynamic dT diffusion
constants.dynamicDT.iterThresholdDecrease = 200;
constants.dynamicDT.iterThresholdIncrease = 25;

% dynamic dT bac
constants.dynamicDT.initRESThresholdIncrease = 20/100;

constants.dynamicDT.nItersCycle = 500; % in diffusion to steady state, after how many diffusion iterations should bulk conc be recalculated
constants.dynamicDT.tolerance_no_convergence = 1e-6; % maximum difference between RES values between diffusion iterations to be considered not converging
constants.dynamicDT.maxRelDiffBulkConc = 0.02; % maximum relative difference between bulk concentration values

constants.simulation_end = 5000;
constants.dT_bac = 1;
settings.detachment = 'none'; % {'mechanistic', 'naive', 'none'}
% constants.kDet = 1; % [1/(m.h)]


% save('sim_xxxx.mat', 'grid', 'bac', 'constants', 'init_params', 'settings', '-v7.3')
% clear ans
%% 
function [X, Y] = rand_circle(N, x_center, y_center, r)
    Ns = round(4/pi * N + 2.5*sqrt(N) + 100);
    X = rand(Ns,1)*(2*r) - r;
    Y = rand(Ns,1)*(2*r) - r;
    I = find(sqrt(X.^2 + Y.^2)<=r);
    X = X(I(1:N)) + x_center;
    Y = Y(I(1:N)) + y_center;
end

function [x, y] = distribute_rectangle(n, xRange, yRange)
    nSections = ceil(sqrt(n)*1.2);
    xlist = linspace(xRange(1), xRange(2), nSections)';
    ylist = linspace(yRange(1), yRange(2), nSections)';
    
    space_margin = xlist(2) - xlist(1);
    
    x = repmat(xlist, nSections, 1);
    y = kron(ylist, ones(nSections, 1));
    
    % add some random noise
    noise_x = rand(size(x)) * 0.6*space_margin - 0.3*space_margin;
    noise_y = rand(size(y)) * 0.6*space_margin - 0.3*space_margin;
    
    x = x + noise_x;
    y = y + noise_y;
    
    
    % reduce the number of bacteria
    i = randperm(length(x));
    i = i(1:n);
    
    x = x(i);
    y = y(i);
end





