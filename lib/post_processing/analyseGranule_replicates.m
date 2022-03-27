function analyseGranule_replicates(simulation_number, nReplicates, finished)
    % create detachment analysis for replicated simulations
    
    if isscalar(finished)
        finished = ones(nReplicates, 1) * finished;
    else
        assert(numel(finished) == nReplicates, 'Array of finished status should match the number of replicates');
    end
    
    %% create all figures:
    % 1- number of bacteria in the simulation domain over time (total &
    % active)
    % 2- active mass per species over time
    % 3- number of active bacteria per species over time
    % 4- density of the reactor over time
    % 5- oxygen uptake rate in reactor
    % 6- bulk concentrations over time
    for n = 1:6
        f = figure(n); hold off;
        clf(n); cla(n);
        close(f);
        clear f
    end
    
    %% set global parameters
    % colors
    global colors_species
    global colors_qualitative
    global linestyles

    global last_nonzero
    global save_times
    global nSpecies
    global nCompounds
    
%     colors_species_raw = {'#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00'};
    colors_species_raw = {'#E69F00','#56B4E9','#33b190','#F0E442','#0072B2','#D55E00'};
    %                      orange   light blue  green    yellow   dark blue    red
    species_per_color = { 'An-NRMX', 'CMX',    'NOB',    'AOB',    'NRMX',    'AMX'};
    colors_qualitative = {'#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c'};
    linestyles = {'-','--',':','-.'};

    
    
    %% plot data for each of the replicates
    [constants, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = load_data(simulation_number, finished(1));
    nSpecies = length(constants.speciesNames);
    nCompounds = numel(constants.compoundNames);
    
    species_index = zeros(nSpecies, 1);
    for s = 1:nSpecies
        species_index(s) = find(strcmp(species_per_color, constants.speciesNames{s}));
    end   
    colors_species = colors_species_raw(species_index);
    final_mass = zeros(nSpecies, nReplicates);

    for n = 1:nReplicates
        % load data for each replicate
        [constants, grid, settings, bac_saved, conc_saved, pH_saved, reactor_saved, profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time] = load_data(simulation_number+n-1, finished(n));

        last_nonzero = find(bac_saved.nBacs ~= 0, 1, 'last');
        save_times = (1:last_nonzero) * constants.dT_save / (7 * 24); % [weeks]

        % add each replicate to all graphs
        plot_nBacs(bac_saved, n)
        final_mass(:, n) = plot_active_mass(bac_saved, constants, n);
        plot_active_nBacs(bac_saved, constants, n)
        plot_reactor_density(bac_saved, constants, settings, n)
        if ~isnan(constants.maintenance)
            plot_reactor_oxygen_uptake(bac_saved, constants, reactor_saved, n)
        end
        plot_bulk_concs(reactor_saved, constants, n)
    end
    
  
    %% show AMX:Nitrospira ratio
    ix_AMX = strcmp(constants.speciesNames, 'AMX');
    AMX_nitrospira_ratio = final_mass(ix_AMX, :)./sum(final_mass(~ix_AMX, :));
    fprintf('AMX : Nitrospira = %f (+- %f): 1\n', mean(AMX_nitrospira_ratio), 2*std(AMX_nitrospira_ratio));
end

%% helper functions
function [constants, grid, settings, bac_saved, conc_saved, pH_saved, reactor_saved, profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time] = load_data(simulation_number, finished)
    %% load correct files
    output_dir = sprintf('./Results/%04d', simulation_number);
    if finished
        simulation_file = sprintf('%s/sim_%04d.mat', output_dir, simulation_number);
    else
        simulation_file = sprintf('sim_%04d.mat', simulation_number);
    end
    simulation_result = sprintf('%s/results1D.mat', output_dir);
    profiling_result = sprintf('%s/profilingResults.mat', output_dir);

    load(simulation_file, 'constants', 'grid', 'settings');    
    load(simulation_result, 'bac_saved', 'conc_saved', 'pH_saved', 'reactor_saved');
    load(profiling_result, 'profiling', 'maxErrors', 'normOverTime', 'nDiffIters', 'bulk_history', 'Time');
end

%% plotting functions
function plot_nBacs(bac_saved, replicate)
    % plot the number of bacteria (active and total) over time
    
    global save_times
    global last_nonzero
    global linestyles
    global colors_qualitative
    
    figure(1); hold on;
    plot(save_times, bac_saved.nBacs(1:last_nonzero), 'LineWidth', 2, 'Color', colors_qualitative{2}, 'LineStyle', linestyles{replicate});
    plot(save_times, sum(bac_saved.active(1:last_nonzero, 1:bac_saved.nBacs(last_nonzero)), 2), 'LineWidth', 2, 'Color', colors_qualitative{1}, 'LineStyle', linestyles{replicate});
    if replicate == 1
        legend('Total', 'Active', 'Location', 'northwest','AutoUpdate','off')
        title('Number of individuals over time');
        xlabel('Simulation time [weeks]');
        ylabel('Number of individuals');
        xlim([0, save_times(end)+5])
        ylim([0 inf])
        xtickies = 0:50:save_times(end);
        xticks(xtickies)
    end
end

function plot_active_nBacs(bac_saved, constants, replicate)
    % plot the number of active individuals per species over time
    
    global nSpecies
    global save_times
    global linestyles
    global last_nonzero
    global colors_species
    
    figure(2); hold on;
    for s = 1:nSpecies
        temp = sum(bac_saved.species == s & bac_saved.active, 2);
        p = plot(save_times, temp(1:last_nonzero), 'LineWidth', 2, 'LineStyle', linestyles{replicate});
        p.Color = colors_species{s};
    end
    
    if replicate == 1
        legend(constants.speciesNames, 'Location', 'northwest','AutoUpdate','off')
        xlabel('Simulation time [weeks]')
        ylabel('Number of active individuals')
        xlim([0, save_times(end)+5])
        ylim([0 inf])
        xtickies = 0:50:save_times(end);
        xticks(xtickies)
    end
end

function final_mass = plot_active_mass(bac_saved, constants, replicate)
    % plot the cumulative active mass per species over time
    
    global nSpecies
    global save_times
    global linestyles
    global last_nonzero
    global colors_species
    
    figure(3); hold on;
    final_mass = zeros(nSpecies, 1);
    for s = 1:nSpecies
        mass = zeros(last_nonzero, 1);
        for t=1:last_nonzero
            mask = bac_saved.species(t,:) == s & bac_saved.active(t,:);
            temp = (constants.bac_rho * 4/3 * pi) * bac_saved.radius(t,:).^3;
            mass(t) = sum(temp(mask));
        end
        final_mass(s) = mass(end);
        p = plot(save_times, mass, 'LineWidth', 2, 'LineStyle', linestyles{replicate});
        p.Color = colors_species{s};
    end
    
    if replicate == 1
        legend(constants.speciesNames, 'Location', 'northwest','AutoUpdate','off')
        xlabel('Simulation time [weeks]')
        ylabel('Cumulative active mass [g]')
        xlim([0, save_times(end)+5])
        xtickies = 0:50:save_times(end);
        xticks(xtickies)
        ylim([0 inf])
    end
end

function plot_reactor_density(bac_saved, constants, settings, replicate)
    % plot the reactor density (total and active) over time
    
    global last_nonzero
    global save_times
    global colors_qualitative
    global linestyles
    
    granule_radius = zeros(last_nonzero, 1);
        
    for i = 1:last_nonzero
        nBac = bac_saved.nBacs(i);
        x = bac_saved.x(i, bac_saved.active(i, 1:nBac));
        y = bac_saved.y(i, bac_saved.active(i, 1:nBac));
        center_x = mean(x);
        center_y = mean(y);
        granule_radius(i) = max(sqrt((x - center_x).^2 + (y - center_y).^2));
    end
    
    if strcmp(settings.model_type, 'suspension')
        f = 1;
    else
        f = 4 * granule_radius ./ (3 * (constants.bac_max_radius * 2));
    end
    
    figure(4); hold on;
    cumul_mass = sum(4/3 * pi * bac_saved.radius.^3 * constants.bac_rho, 2);
    cumul_mass_active = zeros(size(cumul_mass));
    for t=1:last_nonzero
        cumul_mass_active(t) = sum(4/3 * pi * bac_saved.radius(t, bac_saved.active(t, :)).^3 * constants.bac_rho, 2);
    end
    
    plot(save_times, cumul_mass(1:last_nonzero) .* f / constants.Vr, 'LineWidth', 2, 'LineStyle', linestyles{replicate}, 'Color', colors_qualitative{1}); 
    plot(save_times, cumul_mass_active(1:last_nonzero) .* f / constants.Vr, 'LineWidth', 2, 'LineStyle', linestyles{replicate}, 'Color', colors_qualitative{2});
    
    if replicate == 1
        legend('Total', 'Active','AutoUpdate','off', 'Location', 'northwest')
        title('Microbial density in reactor')
        xlabel('Simulation time [weeks]')
        ylabel('density in reactor [g/L]')
        xlim([0, save_times(end)+5])
        xtickies = 0:50:save_times(end);
        xticks(xtickies)
        ylim([0 inf])
    end
end

function plot_reactor_oxygen_uptake(bac_saved, constants, reactor_saved, replicate)
    % plot the uptake rate of oxygen in the reactor over time

    global last_nonzero
    global save_times
    global nCompounds
    global nSpecies
    global colors_qualitative
    global linestyles
    H_o2 = 1.39e-3; % mol/L/atm
    O2_content_air = 0.2;

    figure(5); hold on;


    R_cumulative = zeros(last_nonzero, nCompounds, nSpecies);
    for s = 1:nSpecies
        maint = constants.maintenance(s);
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
    R_sum_o2 = R_sum(:,4); % mol/L/h oxygen
    
    P_o2 = @(kLa) (-R_sum_o2./kLa + reactor_saved.bulk_concs(1:last_nonzero, 4))./H_o2;

    kLa = [10 20 30];
    for i = 1:length(kLa)
        plot(save_times, P_o2(kLa(i))*100 / O2_content_air, 'LineWidth',2, 'Color', colors_qualitative{i}, 'LineStyle', linestyles{replicate})
    end
    
    if replicate == 1
        xlabel('Simulation time [weeks]')
        ylabel('Percentage of air in gas inflow')
        title('air-% in gas inflow, from oxygen transfer to liquid (kLa)')
        legend('kLa = 10', 'kLa = 20', 'kLa = 30','AutoUpdate','off')
        xlim([0, save_times(end)+5])
        xtickies = 0:50:save_times(end);
        xticks(xtickies)
        ylim([0 inf])
    end
end

function plot_bulk_concs(reactor_saved, constants, replicate)
    % plot the concentration of all compounds in the bulk liquid over time

    global save_times
    global last_nonzero
    global colors_qualitative
    global linestyles
    global nCompounds
    
    f = figure(6); hold on;
    f.Position = [50, 50, 1000, 700];
    
    nPlots = [4,2]; % 3 plots in x direction, 2 in y direction (6 total)

    for subix = 1:nPlots(1)
        for subiy = 1:nPlots(2)
            plot_index = (subix - 1)*nPlots(2)+subiy;
            if plot_index > nCompounds
                return
            end
            subplot(nPlots(1), nPlots(2), plot_index); hold on;
            plot(save_times, reactor_saved.bulk_concs(1:last_nonzero, plot_index), 'LineWidth', 2, 'Color', colors_qualitative{2}, 'LineStyle', linestyles{replicate});
            
%             if plot_index == 1 || plot_index == 2
%                 ylim([0 reactor_saved.bulk_concs(last_nonzero, plot_index)*1.05])
%             end
            if replicate == 1
                title(sprintf('Bulk concentration [%s] over time', constants.compoundNames{plot_index}));
                xlabel('Simulation time [weeks]')
                ylabel('Concentration [mol/L]')
                xlim([0, save_times(end)+5])
                xtickies = 0:50:save_times(end);
                xticks(xtickies)
                ylim([0 inf])
            end
        end
    end    
end
