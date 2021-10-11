function plotProfiling(profiling, simulation_end)
    % debug plotting of profiling values
    %
    % profiling: time [seconds] taken per 1 dT_bac of simulation time
    
    f = figure('Name', 'Profiling of simulation'); clf;
    % f.Position = [-1600 100 1400 800]; % desktop with 2 screens
    f.Position = [50, 50, 1000, 700];

    names = {'Diffusion', 'update bacterial properties', 'reaction matrix', 'check if steady state', 'bacs: divide/die/inactivate', 'shoving', 'detachment', 'create bacterial position matrices', 'determine diffusion region', 'calculate bulk concentration', 'Creating chunks & reorganising bacs'};

    [~, sortIndex] = sort(profiling(min(10, ceil(size(profiling, 1)/6)), :), 'ascend');
    names_cat = categorical(names(sortIndex));
    names_cat = reordercats(names_cat, names(sortIndex));

    nPlots = [3,2]; % 3 plots in x direction, 2 in y direction (6 total)
    dt_plot = simulation_end / prod(nPlots);
    plot_time = ceil((1:prod(nPlots)) * dt_plot);

    for subix = 1:nPlots(1)
        for subiy = 1:nPlots(2)
            iPlot = (subix - 1)*nPlots(2)+subiy;
            subplot(nPlots(1), nPlots(2), iPlot)
            barh(names_cat, profiling(plot_time(iPlot), sortIndex), 0.7);
            title(sprintf('Profiling at t = %.1f', plot_time(iPlot)));
            xlabel('time taken [s]')
        end
    end
end

