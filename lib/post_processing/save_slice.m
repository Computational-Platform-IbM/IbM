function save_slice(bac, conc, bulk_concentrations, pH, invHRT, Time, grid, constants, directory)
    % Save important variables along the central axis of the bio-aggregate
    %
    % bac: struct containing all information regarding the bacteria
    % conc: matrix containing all concentrations per grid cell as of
    %   (ix, iy, compound)
    % bulk_concentrations: vector of the bulk liquid concentration of all
    %   compounds
    % pH: matrix containing the pH value per grid cell as (ix, iy)
    % Time: simulation time
    % grid: struct containing all information regarding the grid
    % constants: struct containing all simulation constants
    % directory: directory where results are to be stored in
    
    %% initialise or load previous values
    
    results_file = [directory, '/results1D.mat'];
    if Time == 0
        [bac_saved, conc_saved, pH_saved, reactor_saved] = init_save(constants, grid);
    else
        load(results_file, 'bac_saved', 'conc_saved', 'pH_saved', 'reactor_saved');
    end
    
    %% set values
    iSave = ceil((Time+0.01) / constants.dT_save);
    
    % bacterial variables
    nBacs = length(bac.x);
    bac_saved.nBacs(iSave) = nBacs;
    bac_saved.x(iSave, 1:nBacs) = bac.x;
    bac_saved.y(iSave, 1:nBacs) = bac.y;
    bac_saved.radius(iSave, 1:nBacs) = bac.radius;
    bac_saved.species(iSave, 1:nBacs) = bac.species;
    bac_saved.active(iSave, 1:nBacs) = bac.active;
    bac_saved.mu(iSave, 1:nBacs) = bac.mu;
    
    % concentration variable
    conc_saved(iSave, :, :) = conc(:, ceil(grid.nY / 2), :); % save horizontal slice through center of granule
    
    % pH variable
    pH_saved(iSave, :) = pH(:, ceil(grid.nY / 2));
    
    % reactor properties
    reactor_saved.bulk_concs(iSave, :) = bulk_concentrations;
    reactor_saved.HRT(iSave) = 1 / invHRT;
    reactor_saved.granule_density(iSave) = sum(bac.molarMass * constants.bac_MW) / ...
        ( (max(bac.y) - min(bac.y)) * (max(bac.x) - min(bac.x)) * 1e-6 );
    
    %% save struct to file
    save(results_file, 'bac_saved', 'conc_saved', 'pH_saved', 'reactor_saved', '-v7.3'); % <C: -v7.3 is default? v7.3 file is more than twice the size of v7 file... but still only 100 kB />
    
end

function [bac_saved, conc_saved, pH_saved, reactor_saved] = init_save(constants, grid)
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

    % bacterial variables
    bac_saved = struct;
    bac_saved.nBacs = zeros(nSaves, 1, 'uint32');
    bac_saved.x = zeros(nSaves, constants.max_nBac, 'single');
    bac_saved.y = zeros(nSaves, constants.max_nBac, 'single');
    bac_saved.radius = zeros(nSaves, constants.max_nBac, 'single');
    bac_saved.species = zeros(nSaves, constants.max_nBac, 'uint8');
    bac_saved.active = zeros(nSaves, constants.max_nBac, 'logical');
    bac_saved.mu = zeros(nSaves, constants.max_nBac, 'single');
    
    % concentration variable
    nCompounds = length(constants.compoundNames);
    conc_saved = zeros(nSaves, grid.nX, nCompounds, 'single');
    
    % pH variable
    pH_saved = zeros(nSaves, grid.nX, 'single');
    
    % reactor properties
    reactor_saved.bulk_concs = zeros(nSaves, nCompounds, 'single');
    reactor_saved.HRT = zeros(nSaves, 1, 'single');
    reactor_saved.granule_density = zeros(nSaves, 1, 'single');
end
