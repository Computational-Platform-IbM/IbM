function plot2d_S(simulation_number, substrates_number, Time)
    % Function to create 2D substrates profiles for a given simulation number
    % Time (days)
        
    % check if there is a results folder for this simulation
    output_dir = sprintf('./Results/%04d', simulation_number);
    if ~isfolder(output_dir)
        error('Simulation %04d is not present', simulation_number)
    end
    
    % check if simulation file exists (meaning the simulation has finished)
    simulation_file = sprintf('%s/sim_%04d.mat', output_dir, simulation_number);
    if ~isfile(simulation_file)
        error('The simulation file %s has not finished running', simulation_file)
    end
    
    % check if simulation result file exists
    simulation_result = sprintf('%s/results1D.mat', output_dir);
    if ~isfile(simulation_result)
        error('The results file %s does not exist', simulation_result)
    end
    
    load(simulation_result, 'conc_saved');                 % result1D.mat
    load(simulation_file, 'grid', 'constants');             % sim_####.mat
    
    plotConcs2D(grid, conc_saved, constants, Time, substrates_number)
    
end