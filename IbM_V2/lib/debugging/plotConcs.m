function plotConcs(conc, constants, T)
    f = figure('Name', 'Concentration profile');
    f.Position = [50, 80, 1450 680];
    clf;
    
    tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 

    for iPlot = 1:5
        nexttile
        imagesc(conc(:,:,iPlot)); 
        colormap(viridis()); 
        colorbar(); 
        title(sprintf('Concentration profile of %s', constants.StNames{iPlot}));
        cval = caxis();
        if cval(1) < 0
            warning('concentration below zero detected in %s', constants.StNames{iPlot})
        end
        caxis([0, cval(2)]);
        axis square;
    end
    
    sgtitle(sprintf('Concentration profiles at t = %.2f', T));
end

