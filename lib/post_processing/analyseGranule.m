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
%     profiling_result = sprintf('%s/profilingResults.mat', output_dir);

    load(simulation_file, 'constants', 'grid', 'settings');    
    load(simulation_result, 'bac_saved', 'conc_saved', 'pH_saved', 'reactor_saved');
%     load(profiling_result);
    
    
    %% plot number of bacteria over time
    figure(1); clf;
    last_nonzero = find(bac_saved.nBacs ~= 0, 1, 'last');
    save_times = (1:last_nonzero) * constants.dT_save;
    plot(save_times, bac_saved.nBacs(1:last_nonzero), 'LineWidth', 2); hold on
    plot(save_times, sum(bac_saved.active(1:last_nonzero, 1:bac_saved.nBacs(last_nonzero)), 2));
    title('Bacteria over time');
    xlabel('Simulation time [h]');
    ylabel('Number of bacteria');
    
        
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
    colors_species_raw = {'#E69F00','#56B4E9','#33b190','#F0E442','#0072B2','#D55E00'};
    %                      orange   light blue  green    yellow   dark blue    red
    species_per_color = { 'An-NRMX', 'CMX',    'NOB',    'AOB',    'NRMX',    'AMX'};
    
    species_index = zeros(nSpecies, 1);
    for s = 1:nSpecies
        species_index(s) = find(strcmp(species_per_color, constants.speciesNames{s}));
    end   
    coloring = colors_species_raw(species_index);
    
    
    ar = area((1:nBins) * bin_size * 1e6, strat_data);
    legend(constants.speciesNames, 'Location', 'northwest');
    for s = 1:nSpecies
        ar(s).FaceColor = coloring{s};
    end
    
    xlim([0, granule_radius(end)*1e6])
    ax = gca;
    ax.FontSize = 11;
%     title('Number of individuals per radial segment')
    xlabel('Distance from center of granule [μm]', 'FontSize', 14)
    ylabel({'Number of', 'individuals'}, 'FontSize', 14)
    
    subplot(2,1,2);
    strat_data_relative = strat_data ./ sum(strat_data, 2) * 100;
    ar = area((1:nBins) * bin_size * 1e6, strat_data_relative);
    legend(constants.speciesNames, 'Location', 'northwest');
    for s = 1:nSpecies
        ar(s).FaceColor = coloring{s};
    end
    
    xlim([0, granule_radius(end)*1e6])
    ax = gca;
    ax.FontSize = 11;
%     title('Relative abundance per radial segment')
    xlabel('Distance from center of granule [μm]', 'FontSize', 14)
    ylabel({'Relative','abundance [%]'}, 'FontSize', 14)
    
    
    
    %% plot HRT
    figure(11);
    clf;
    plot(save_times, reactor_saved.HRT(1:last_nonzero), 'LineWidth', 2);
    title('HRT over time')
    xlabel('Simulation time [h]')
    ylabel('Hydrolic Retention Time (HRT) [h]')
    
    
    %% plot bacterial density in reactor
    cumul_mass = sum(4/3 * pi * bac_saved.radius.^3 * constants.bac_rho, 2);
    cumul_mass_active = zeros(size(cumul_mass));
    for t=1:last_nonzero
        cumul_mass_active(t) = sum(4/3 * pi * bac_saved.radius(t, bac_saved.active(t, :)).^3 * constants.bac_rho, 2);
    end
    if strcmp(settings.model_type, 'suspension')
        f = 1;
    else
        f = 4 * granule_radius ./ (3 * (constants.bac_max_radius * 2));
    end
    figure(12);
    clf;
    plot(save_times, cumul_mass(1:last_nonzero) .* f / constants.Vr, 'LineWidth', 2); hold on;
    plot(save_times, cumul_mass_active(1:last_nonzero) .* f/ constants.Vr, 'LineWidth', 2);
    title('Bacterial density in reactor')
    xlabel('Simulation time [h]')
    ylabel('density in reactor [g/L]')
    
    
    %% plot bacterial population over time
    figure(13);
    clf;
    for s = 1:nSpecies
        temp = sum(bac_saved.species == s & bac_saved.active, 2);
        p = plot(save_times, temp(1:last_nonzero), 'LineWidth', 2);
        p.Color = coloring{s};
        hold on;
    end
    hold off;
    legend(constants.speciesNames, 'Location', 'northwest')
    xlabel('Simulation time [h]')
    ylabel('Number of individuals')

    
    figure(14);
    clf;
    final_mass = zeros(nSpecies, 1);
    for s = 1:nSpecies
        mass = zeros(last_nonzero, 1);
        for t=1:last_nonzero
            mask = bac_saved.species(t,:) == s & bac_saved.active(t,:);
            temp = (constants.bac_rho * 4/3 * pi) * bac_saved.radius(t,:).^3;
            mass(t) = sum(temp(mask));
        end
        final_mass(s) = mass(end);
        p = plot(save_times, mass, 'LineWidth', 2);
        p.Color = coloring{s};
        hold on;
    end
    hold off;
    legend(constants.speciesNames, 'Location', 'northwest')
    xlabel('Simulation time [h]')
    ylabel('Cumulative active mass [g]')

    ix_AMX = strcmp(constants.speciesNames, 'AMX');
    fprintf('AMX : Nitrospira = %f : 1\n', final_mass(ix_AMX)/sum(final_mass(~ix_AMX)));
    
    
    %% plot consumption per compound
    nCompounds = numel(constants.compoundNames);
    R_cumulative = zeros(last_nonzero, nCompounds, nSpecies);
    for s = 1:nSpecies
        maint = constants.maintenance(s);
%         mask = bac_saved.species(t,:) == s & bac_saved.active(t,:);
        [t, bix] = find(bac_saved.species == s & bac_saved.active);
        
        for n=1:length(bix)
            mu_withMaintenance = bac_saved.mu(t(n), bix(n));
            mu_noMaintenance = mu_withMaintenance + maint;

            concentrationChange = constants.MatrixMet(:, s) * mu_noMaintenance;                              % [molS/molX/h]
            if mu_withMaintenance < 0
                concentrationChange = concentrationChange - constants.MatrixDecay(:, s) * mu_withMaintenance;       % [molS/molX/h]
            end
            concentrationChange = concentrationChange * (constants.bac_rho / constants.bac_MW * 4/3 * pi) * bac_saved.radius(t(n), bix(n))^3;
            R_cumulative(t(n), :, s) = R_cumulative(t(n), :, s) + reshape(concentrationChange, 1, []);
        end
    end
    R_cumulative = R_cumulative ./ constants.Vr; % mol/L/h
    
    R_sum = sum(R_cumulative(:, 1:4, :), 3);
    R_sum_o2 = R_sum(:,4); % mol/L/h
    infl_o2 = reactor_saved.HRT(1:last_nonzero) .* -R_sum_o2 + reactor_saved.bulk_concs(1:last_nonzero, 4);
    Q = 50e-3; % L/h
    Vm = 24; % L/mol
    O2_content_air = 0.2;
    Q_o2 = Q * infl_o2;
    F_air = Q_o2 * Vm / O2_content_air; % [L/h]
    F_gas_reactor = 0.6; % L/h -> 10 mL/min = 600 mL/h = 0.6 L/h
    
    
    figure(23); clf;
    plot(save_times, F_air ./ F_gas_reactor * 100, 'LineWidth', 2);
    xlabel('Simulation time [h]')
    ylabel('Percentage of air in gas inflow')
    title('air-% in gas inflow, from mass balance')
    
    
    H_o2 = 1.39e-3; % mol/L/atm
    
%     figure(22); clf;
%     plot(save_times, sum(R_cumulative(:,1:4,:), 3), 'LineWidth', 2); 
%     legend(constants.compoundNames{1:4})
%     xlabel('Simulation time [h]')
%     ylabel('Influent consumption/production [mol/L]')

    P_o2 = @(kLa) (-R_sum_o2./kLa + reactor_saved.bulk_concs(1:last_nonzero, 4))./H_o2;
    colors = {'#a6cee3','#1f78b4','#b2df8a','#33a02c'};

    figure(22); clf; hold on;
    kLa = [10 20 30];
    for i = 1:length(kLa)
        plot(save_times, P_o2(kLa(i))*100 / O2_content_air, 'LineWidth',2, 'Color', colors{i})
    end
    xlabel('Simulation time [h]')
    ylabel('Percentage of air in gas inflow')
    title('air-% in gas inflow, from oxygen transfer to liquid (kLa)')
    legend('kLa = 10', 'kLa = 20', 'kLa = 30')
    

end