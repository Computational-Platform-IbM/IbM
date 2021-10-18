function save_profiling(profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time, directory)
    % Save important profiling arrays with information on the performance
    % of the simulation. Will overwrite the last backup in order to always 
    % have the latest file.
    %
    % profiling: matrix with per dT_bac the time per functionality of the
    %   model
    % maxErrors: vector with per dT_bac the maximum RES value
    % normOverTime: vector with per dT_bac the norm of delta-concentrations
    % nDiffIters: vector with per dT_bac the number of diffusion iterations
    % bulk_history: vector with per dT_bac the bulk concentration of each
    %   compound
    % Time: struct containing all information regarding timings
    % directory: directory where results are to be stored in
    
    %% initialise or load previous values
    
    results_file = [directory, '/profilingResults.mat'];

    %% save struct to file
    save(results_file, 'profiling', 'maxErrors', 'normOverTime', 'nDiffIters', 'bulk_history', 'Time');
end