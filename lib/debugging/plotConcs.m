function plotConcs(conc, constants, T)
    f = figure(7);
    f.Position = [50, 80, 1450 680];
    clf;
    
    tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 

    for iPlot = 1:5
        nexttile
        imagesc(conc(:,:,iPlot)'); 
        colormap(viridis()); 
        colorbar(); 
        title(sprintf('Concentration profile of %s', constants.compoundNames{iPlot}));
        cval = caxis();
        if cval(1) < 0
            warning('concentration below zero detected in %s', constants.compoundNames{iPlot})
        end
        caxis([0, cval(2)]);
        axis square;
        set(gca,'YDir','normal')
    end
    
    sgtitle(sprintf('Concentration profiles at t = %.2f', T));
end

