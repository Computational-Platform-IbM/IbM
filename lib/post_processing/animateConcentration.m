function animateConcentration(simulation_number, finished)
    % Function to create animations for a given simulation number

    %% argument validation
    arguments
        simulation_number (1,1) {mustBeInteger, mustBeInRange(simulation_number, 1, 9999)}
        finished
    end
    
    
    %% load correct files
    % check if there is a results folder for this simulation
    output_dir = sprintf('./Results/%04d', simulation_number);
    if ~isfolder(output_dir)
        error('Simulation %04d is not present')
    end

    
    % check if simulation file exists (meaning the simulation has finished)
    if finished
        simulation_file = sprintf('%s/sim_%04d.mat', output_dir, simulation_number);
    else
        simulation_file = sprintf('sim_%04d.mat', simulation_number);
    end
    if ~isfile(simulation_file)
        error('The simulation file %s has not finished running', simulation_file)
    end
    
    % check if simulation result file exists
    simulation_result = sprintf('%s/results1D.mat', output_dir);
    if ~isfile(simulation_result)
        error('The results file %s does not exist', simulation_result)
    end
    
    % check if simulation profiling file exists
    profiling_result = sprintf('%s/profilingResults.mat', output_dir);
    if ~isfile(profiling_result)
        error('The profiling results file %s does not exist', profiling_result)
    end

    
    %% add results to path
    addpath(genpath('Results')); % add Results folder and all subfolders to path
    addpath(genpath('lib'));
    
    
    %% load results
    load(simulation_result, 'bac_saved', 'conc_saved', 'pH_saved', 'reactor_saved');
%     load(profiling_result, 'profiling', 'maxErrors', 'normOverTime', 'nDiffIters', 'bulk_history', 'Time');
    load(simulation_file, 'constants', 'grid');
    
    
    %% pH over time animation
    fprintf('\n===== granule concentration ANIMATION =====\n')
%     pH_animation = sprintf('%s/pH.avi', output_dir);
    timer = tic;

    % determine size of figure frame for plotting
    last_nonzero = find(bac_saved.nBacs ~= 0, 1, 'last');
    xlim_conc = [min(bac_saved.x(last_nonzero, 1:bac_saved.nBacs(last_nonzero))) - 5*grid.dx, ...
        max(bac_saved.x(last_nonzero, 1:bac_saved.nBacs(last_nonzero))) + 5*grid.dx] * 1e6;
%     ylim_pH = [min(pH_saved(pH_saved ~= 0), [], 'all') - 1, ...
%         max(pH_saved(pH_saved ~= 0), [], 'all') + 1];
%     ylims = max(conc_saved, [], [1,2])*1.1;
    

    % start animator object
%     writerObjB = VideoWriter(pH_animation);
%     writerObjB.FrameRate = 1;
%     open(writerObjB);

    for i = 1:last_nonzero
        T = constants.dT_save * (i - 1);

        f = figure(10);
        clf;
        for ii = 1:4
            subplot(2,2,ii)
            plot(((0:grid.nX-1)*grid.dx - .5*grid.dx)*1e6, conc_saved(i, :, ii), 'LineWidth', 2);
            xlabel('Position along central axis in granule [μm]');
            ylabel('Concentration [mol/L]');
            title(sprintf('%s', constants.compoundNames{ii}))
            xlim(xlim_conc);
            ylim([0, max(conc_saved(i, :, ii))]);
        end
        drawnow();
%         frame = getframe(f);
%         writeVideo(writerObjB,frame);
    end

    % close animation object
%     close(writerObjB);

    timeTaken = toc(timer);

    % report back to console
    fprintf('Composed in %.2f seconds\n', timeTaken);
end

