function save_backup(bac, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
    % Save important variables required for restart at this point in time
    % Will overwrite the last backup in order to always have the latest
    % file
    %
    % bac: struct containing all information regarding the bacteria
    % conc: matrix containing all concentrations per grid cell as of
    %   (ix, iy, compound)
    % bulk_concs: vector of the bulk liquid concentration of all
    %   compounds
    % pH: matrix containing the pH value per grid cell as (ix, iy)
    % Time: struct containing all time parameters
    % invHRT: current 1/HRT value [1/h]
    % reaction_matrix: matrix containing all reaction rates per grid cell
    %   as (ix, iy, compound) [mol/L/h]
    % directory: directory where results are to be stored in
    
    %% initialise or load previous values
    
    results_file = [directory, '/backup.mat'];

    %% save struct to file
    save(results_file, 'bac', 'bulk_concs', 'invHRT', 'conc', 'reaction_matrix', 'pH');
end