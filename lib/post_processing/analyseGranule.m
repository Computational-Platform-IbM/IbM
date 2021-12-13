function analyseGranule(simulation_number, finished)
    % create detachment analysis for granule simulations
    
    
    %% load correct files
    output_dir = sprintf('./Results/%04d', simulation_number);
    if finished
        simulation_file = sprintf('%s/sim_%04d.mat', output_dir, simulation_number);
    else
        simulation_file = sprintf('sim_%04d.mat', simulation_number);
    end
    simulation_result = sprintf('%s/results1D.mat', output_dir);
    profiling_result = sprintf('%s/profilingResults.mat', output_dir);

    load(simulation_file);
    load(simulation_result);
    load(profiling_result);
    
    
    %% plot number of bacteria over time
    figure(1); clf;
    last_nonzero = find(bac_saved.nBacs ~= 0, 1, 'last');
    save_times = (1:last_nonzero) * constants.dT_save;
    plot(save_times, bac_saved.nBacs(1:last_nonzero), 'LineWidth', 2); hold on
    plot(save_times, sum(bac_saved.active(1:last_nonzero, 1:bac_saved.nBacs(last_nonzero)), 2));
    title('Bacteria over time');
    xlabel('Simulation time [h]');
    ylabel('Number of bacteria');
    
    
    %% plot bacteria
    
    %% plot granule size over time
    figure(3); clf;
    granule_radius = zeros(last_nonzero, 1);
    
    
    for i = 1:last_nonzero
        nBac = bac_saved.nBacs(i);
        x = bac_saved.x(i, bac_saved.active(i, 1:nBac));
        y = bac_saved.y(i, bac_saved.active(i, 1:nBac));
        center_x = mean(x);
        center_y = mean(y);
        granule_radius(i) = max(sqrt((x - center_x).^2 + (y - center_y).^2));
    end
    plot(save_times, granule_radius*1e6, 'LineWidth', 2);
    title('Granule size over time')
    xlabel('Simulation time [h]')
    ylabel('Granule radius [um]')
    
    
    
    %% plot stratification
    % only at last point in time to check if it works
    figure(4); clf;
    subplot(2,1,1);
    
    % create bins for radial distance
    nBins = 50;
    bin_size = granule_radius(end) / nBins;
    
    center_x = grid.dx*grid.nX / 2;
    center_y = grid.dy*grid.nY / 2;
    nBac = bac_saved.nBacs(i);
    x = bac_saved.x(i, 1:nBac);
    y = bac_saved.y(i, 1:nBac);
    radial_dist = sqrt((x - center_x).^2 + (y - center_y).^2);
    
    bin = ceil(radial_dist ./ bin_size);
    
    
    % create data table (nBins * nSpecies)
    nSpecies = length(constants.speciesNames);
    strat_data = zeros(nBins, nSpecies);
    for b = 1:nBins
        for s = 1:nSpecies
            strat_data(b, s) = sum(bac_saved.species(i, 1:nBac) == s & bin == b & bac_saved.active(i, 1:nBac));
        end
    end
    
    % calculate per bin the relative abundance of the bacteria
    
    % plot in stacked bar chart
    coloring = {'#D81B60', '#1E88E5', '#FFC107', '#004D40'};
    coloring = {'#E69F00', '#56B4E9','#009E73','#F0E442','#0072B2','#D55E00'};

    ar = area((1:nBins) * bin_size * 1e6, strat_data);
    legend(constants.speciesNames);
    for s = 1:nSpecies
        ar(s).FaceColor = coloring{s};
    end
    title('Absolute bacterial numbers per radial segment')
    xlabel('Distance from center of granule [um]')
    ylabel('Number of bacteria')
    
    subplot(2,1,2);
    strat_data_relative = strat_data ./ sum(strat_data, 2);
    ar = area((1:nBins) * bin_size * 1e6, strat_data_relative);
    legend(constants.speciesNames);
    for s = 1:nSpecies
        ar(s).FaceColor = coloring{s};
    end
    title('Relative bacterial numbers per radial segment')
    xlabel('Distance from center of granule [um]')
    ylabel('Relative portion of the number of bacteria')
    
    
    
    %% plot HRT
    figure(11);
    clf;
    plot(Time.history(1:last_nonzero), reactor_saved.HRT(1:last_nonzero), 'LineWidth', 2);


end