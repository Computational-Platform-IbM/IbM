function plotConcs2D(g, conc, constants, T, substrates_number)
    f = figure('Name', sprintf('Concentration profiles at t = %.1fd', T));
    f.Position = [60 60 750 600];
    clf;
    
    coloring = {'#CC66FF', '#00B04F', '#FFA200', '#010101', '#010101'}; % Old colours
%     coloring = {'#D81B60', '#1E88E5', '#FFC107', '#004D40'}; % colorblind-friendly colours
    
    xt = linspace(0, g.nX*g.dx, g.nX) * 1e6;    %um
%     yt = linspace(0, g.nY*g.dy, g.nY+1);
    maxY = max(max(max(conc * 1000)));
    if substrates_number > 1 && substrates_number < 4 % Neutralism, Commensalism and co-protection
        for iSubs = 1:substrates_number

            plot(xt, conc(T, :, iSubs) * 1000, 'Color', coloring{iSubs}, 'LineWidth', 2.5)      %mM vs um
            hold on
        end
        hold off
    elseif substrates_number == 4
        for iSubs = 1:substrates_number+1
            if iSubs == 4
                continue;
            end
            plot(xt, conc(T, :, iSubs) * 1000, 'Color', coloring{iSubs}, 'LineWidth', 2.5)      %mM vs um
            hold on
        end
        hold off
    else
%         plot(xt, conc(T, :, 1) * 1000 + 0.07, 'Color', coloring{1}, 'LineWidth', 2.5)                                   %mM vs um
        hold on
        plot(xt, conc(T, :, 1) * 1000 + 0.00, 'Color', coloring{1}, 'LineWidth', 2.5)                                   %mM vs um
        hold on
%         plot(xt, conc(T, :, 1) * 1000 - 0.06, 'Color', coloring{3}, 'LineWidth', 2.5)                                   %mM vs um
        hold off
    end
    xlim([200 330])
    ylim([0 1.1*maxY]) %
    box on
    set(gca,'TickDir','out');
    set(gcf,'color','w');
    set(gca,'FontSize',20)
    set(gca,'linewidth',1.5)
    xlabel('x (\mum)')
    ylabel('[S] (mM)')
    if substrates_number < 4
        legend(constants.compoundNames{1:substrates_number})
    else
        legend(constants.compoundNames{1:3},constants.compoundNames{5:substrates_number+1})
    end
end

