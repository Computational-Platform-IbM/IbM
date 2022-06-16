function f = plotBacs(g, bac, constants, Time)
    % testing function to visualise bacteria
    x = bac.x;
    y = bac.y;
    r = bac.radius;
    s = bac.species;
    act = bac.active;
    
    %% Colour sets
    % Nitrospira's metabolism niches
%     colors_species_raw = {'#E69F00','#56B4E9','#33b190','#F0E442','#0072B2','#D55E00'};
%     %                      orange   light blue  green    yellow   dark blue    red
%     colors_species_rgb = [230 159 0; 86 180 233; 51 177 144; 240 228 66; 0 114 178; 213 94 0]./255;
%     species_per_color = { 'An-NRMX', 'CMX',    'NOB',    'AOB',    'NRMX',    'AMX'};
    % NOB genera diversity
    colors_species_raw = {'#CC00FF','#00B050','#FFC000','#6292FE','#FF017A'};
    %                      purple     green     yellow    blue      pink
    colors_species_rgb = [204 0 255; 0 176 80; 255 192 0; 98 146 254; 255 1 122]./255;
    species_per_color = { 'AOB', 'Nitrobacter', 'Nitrospira', 'Nitrotoga', 'AMX'};    
    
    %% pLoting
    nSpecies = length(constants.speciesNames);
    species_index = zeros(nSpecies, 1);
    for si = 1:nSpecies
        species_index(si) = find(strcmp(constants.speciesNames{si}, species_per_color));
    end   
    coloring = colors_species_raw(species_index);
    coloring_rgb = colors_species_rgb(species_index, :);
    
    
    f = figure(); clf;
%     f.Position = [-1800, 65, 1200, 900]; % desktop with two screens
    f.Position = [60 60 1000 700];
    
    for i = 1:length(x)
        if act(i)
%             col = coloring{s(i)};	% HEX
            col = [coloring_rgb(s(i), :) 1];
            ir = r(i);
        else
            col = [0 0 0];          % HEX
%             col = [coloring_rgb(s(i), :) 0.3];
            ir = 1.0 * r(i);
        end
        rectangle('Curvature', [1 1], 'Position', [x(i) - ir, y(i) - ir, 2 * ir, 2 * ir], 'LineWidth', 0.1, 'EdgeColor', col, 'FaceColor', col);
    end
    
    axis equal;

    xlim([0, g.nX*g.dx]);
    ylim([0, g.nY*g.dy]);
    xt = linspace(0, g.nX*g.dx, g.nX+1);
    yt = linspace(0, g.nY*g.dy, g.nY+1);
    xticks(xt(1:10:end));
    yticks(yt(1:10:end));
    xticklabels(xt(1:10:end)*1e6);
    yticklabels(yt(1:10:end)*1e6);
    
    ax = gca;
    ax.XRuler.MinorTick = 'on'; %or 'off'
    ax.XRuler.MinorTickValues = xt; %just like major ticks
    ax.XRuler.MinorTickValuesMode = 'manual'; %or 'manual'
    
    ax.YRuler.MinorTick = 'on'; %or 'off'
    ax.YRuler.MinorTickValues = yt; %just like major ticks
    ax.YRuler.MinorTickValuesMode = 'manual'; %or 'manual'
    
    grid off;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    
    % add dummy scatters to make legend 
    hold on;
    uniq_s = unique(s);
    for i = 1:length(uniq_s)
        if uniq_s(i) == 0
            continue
        end
        scatter(NaN, NaN, 'Marker', 'o', 'MarkerFaceColor', coloring{uniq_s(i)}, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
    end
    
    % increase marker size in legend - Unsupported behaviour in Matlab, but does still work...
    [~, icons] = legend(constants.speciesNames{uniq_s}, 'FontSize', 12, 'Location', 'bestoutside');
    icons = findobj(icons, 'type', 'patch');
    set(icons, 'MarkerSize', 12);
    set(icons, 'LineWidth', 1);
    hold off;
    
    ylabel('y-position [μm]');
    xlabel('x-position [μm]');
    title(sprintf('Bacteria at t=%.1f', Time))

end