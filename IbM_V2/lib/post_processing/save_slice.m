function save_slice(bac, conc, bulk_concentrations, pH, Time, grid, constants, directory)
    % Save important variables along the central axis of the bio-aggregate
    %
    % bac: struct containing all information regarding the bacteria
    % conc: matrix containing all concentrations per grid cell as of
    %   (ix, iy, compound)    % pH:
    % Time: simulation time
    % grid: struct containing all information regarding the grid
    % constants: struct containing all simulation constants
    % directory: directory where results are to be stored in
    
    %% initialise or load previous values
    
    results_file = [directory, '/results.mat'];
    if Time == 0
        Results = init_save(constants, grid);
    else
        Results = load(results_file);
    end
    
    %% set values
    iSave = ceil(Time / constants.dT_save) + 1;
    
    % bacterial variables
    nBacs = length(bac.x);
    Results.bac.x(iSave, 1:nBacs) = bac.x;
    Results.bac.y(iSave, 1:nBacs) = bac.y;
    Results.bac.radius(iSave, 1:nBacs) = bac.radius;
    Results.bac.species(iSave, 1:nBacs) = bac.species;
    Results.bac.active(iSave, 1:nBacs) = bac.active;
    
    % concentration variable
    Results.conc(iSave, :) = conc;
    
    % pH variable
    Results.pH(iSave, :) = pH;
    
    % reactor properties
    Results.reactor.bulk_concs(iSave, :) = bulk_concentrations;
    Results.reactor.HRT(iSave) = 1 / invHRT;
    Results.reactor.granule_density(iSave) = sum(bac.molarMass * constants.bac_MW) / ...
        ( (max(bac.y) - min(bac.y)) * (max(bac.x) - min(bac.x)) * 1e-6 );
    
    %% save struct to file
    save(results_file, 'Results', '-v7.3'); % <C: -v7.3 is default? v7.3 file is more than twice the size of v7 file... but still only 100 kB />
    
end

function Results = init_save(constants, grid)
    % initialise the results struct for the number of saves that are going
    % to happen. Utilizes datatypes within the required precision with the
    % lowest storage requirements.
    %
    % constants: struct containing all simulation constants
    % grid: struct containing all information regarding the grid
    % 
    % -> Results: struct containing all variables to be saved, initialised
    %   for the number of saves that are going to be made
    
    nSaves = ceil(constants.simulation_end / constants.dT_save) + 1;

    Results = struct;
    
    % bacterial variables
    Results.bac.x = zeros(nSaves, constants.max_nBac, 'single');
    Results.bac.y = zeros(nSaves, constants.max_nBac, 'single');
    Results.bac.radius = zeros(nSaves, constants.max_nBac, 'single');
    Results.bac.species = zeros(nSaves, constants.max_nBac, 'int8');
    Results.bac.active = zeros(nSaves, constants.max_nBac, 'logical');
    
    % concentration variable
    Results.conc = zeros(nSaves, grid.nX, 'single');
    
    % pH variable
    Results.pH = zeros(nSaves, grid.nX, 'single');
    
    % reactor properties
    Results.reactor.bulk_concs = zeros(nSaves, sum(constants.isLiquid), 'single');
    Results.reactor.HRT = zeros(nSaves, 1, 'single');
    Results.reactor.granule_density = zeros(nSaves, 1, 'single');
end