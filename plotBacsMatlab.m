function plotBacsMatlab(simulation_number, Time)
    % Function to create bacteria animation for a given simulation number
    % Time (days)

    addpath(genpath('lib')); % make every subfolder with functions accessible to the code
    
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
    
    load(simulation_result, 'bac_saved');               % result1D.mat
    load(simulation_file, 'grid', 'constants');         % sim_####.mat 
    
    maxT = constants.simulation_end;
    dt_save = constants.dT_save;
    maxi = maxT/dt_save;
    i = ceil((Time*24)/24) + 2;
    if (i > maxi), i = maxi;  end
    
    nBac = bac_saved.nBacs(i);
    bac_plot = struct;
    bac_plot.x = bac_saved.x(i,1:nBac);
    bac_plot.y = bac_saved.y(i,1:nBac); 
    bac_plot.radius = bac_saved.radius(i,1:nBac);
    bac_plot.species = bac_saved.species(i,1:nBac);
    bac_plot.active = bac_saved.active(i,1:nBac);
    
    plotBacs(grid, bac_plot, constants, Time * 24)
    
end