function plotRMatrix(rMatrix, constants, T)
    f = figure(10);
    f.Position = [50, 80, 1450 680];
    clf;
    
    tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 

    for iPlot = 1:5
        nexttile
        imagesc(rMatrix(:,:,iPlot)'); 
        colormap(redblue()); 
        colorbar(); 
        title(sprintf('Reaction profile of %s', constants.StNames{iPlot}));
        cval = caxis();
        cval_maxabs = max(abs(cval));

        caxis([-cval_maxabs, cval_maxabs]);
        axis square;
        set(gca, 'YDir', 'normal')
    end
    
    sgtitle(sprintf('Reaction profiles at t = %.2f', T));
end

