function IbM(simulation_number)
%%%%% Main function to run the IbM model for a given preset file (naming
% convention: 'sim_xxxx.mat', where xxxx is the simulation number)
%
% -> Creates a folder in the Results directory for the output of the
%   simulation.
% -> Run the simulation based on the simulation file (save results in the
%   corresponding results folder)
% -> Save profiling results in the corresponding results folder
% -> Move the preset file in the corresponding results folder (signifying
%   that the simulation has been run) after running the model

    %% argument validation
    arguments
        simulation_number (1,1) {mustBeInteger, mustBeInRange(simulation_number, 1, 9999)}
    end
        
    % check if simulation file exists
    simulation_file = sprintf('sim_%04d.mat', simulation_number);
    if ~isfile(simulation_file)
        error('The simulation file %s does not exist', simulation_file)
    end
    
    % create output directory for results
    output_dir = sprintf('./Results/%04d', simulation_number);
    if ~isfolder(output_dir)
        mkdir(output_dir)
    end
    
    
    %% java import for shoving
%     javaaddpath([pwd '/lib/shovingQuadTree.jar']);
    javaaddpath([pwd '/lib/bacteria/shovingQuadTreekDist.jar']);
    addpath(genpath('lib')); % make every subfolder with functions accessible to the code
    
    %% enable/disable debug disp/warning
    warning('off', 'DEBUG:noActionRequired');
    warning('on', 'DEBUG:actionRequired');
    
    
    %% ========== Time advancements ==========
    fprintf('> SIMULATION RUNNING >>>>>\n');
        
    tTime = tic;
    integTime(simulation_file, output_dir);
    totalTime = toc(tTime);
    
    fprintf('> SIMULATION FINISHED >>>>>\n');
    

    %% save profiling information in output directory
    load(simulation_file, 'constants');
    fprintf('\n\nTotal time for simulation of %.2f hours:\n\t%.2f seconds\n', constants.simulation_end, totalTime)

    
    %% cleanup of root directory
    movefile(simulation_file, output_dir);
end


