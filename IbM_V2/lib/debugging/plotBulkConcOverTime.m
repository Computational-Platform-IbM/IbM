function plotBulkConcOverTime(bulk_history, constants)
    % plot the concentration of all compounds in the bulk liquid over time
    %
    % bulk_history: matrix (nCompounds-by-n) with bulk concentration per
    %   compound at each dT_bac timepoint
    f = figure(6);
    f.Position = [50,50, 1000, 700];
    
    nPlots = [4,2]; % 3 plots in x direction, 2 in y direction (6 total)

    for subix = 1:nPlots(1)
        for subiy = 1:nPlots(2)
            plot_index = (subix - 1)*nPlots(2)+subiy;
            subplot(nPlots(1), nPlots(2), plot_index)
            plot((0:length(bulk_history(1,:))-1)*constants.dT_bac, bulk_history(plot_index,:), 'LineWidth', 2);
            title(sprintf('Bulk concentration [%s] over time', constants.StNames{plot_index}));
            xlabel('Simulation time [h]')
            ylabel('Concentration [mol/L]')
        end
    end
    
    
end