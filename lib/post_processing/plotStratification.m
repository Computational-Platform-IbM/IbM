function plotStratification(simulation_number, finished, TimeInput)
% Function to plot radial distribution (layered stratification) 
% and angular distribution (columned stratification) of microbial colonies.
    
    BoxChartPr = 0; % Draw BoxChart of relative perimeter size of microcolonies
    MATLABdraw_layered = 1; % Draw plot -> 1 | No draw plot -> 0
    ATLABdraw_columned = 1; % Draw plot -> 1 | No draw plot -> 0
    
    %% load correct files
    output_dir = sprintf('./Results/%04d', simulation_number);
    if finished
        simulation_file = sprintf('%s/sim_%04d.mat', output_dir, simulation_number);
    else
        simulation_file = sprintf('sim_%04d.mat', simulation_number);
    end
    simulation_result = sprintf('%s/results1D.mat', output_dir);
    profiling_result = sprintf('%s/profilingResults.mat', output_dir);
    
    folder_route = './Results/';

    load(simulation_file);      %#ok<LOAD>
    load(simulation_result);    %#ok<LOAD>
    load(profiling_result);     %#ok<LOAD>
    
    %% plot number of bacteria over time
    last_nonzero = find(bac_saved.nBacs ~= 0, 1, 'last');
%     save_times = (1:last_nonzero) * constants.dT_save;
%     figure(1); clf;
%     plot(save_times, bac_saved.nBacs(1:last_nonzero), 'LineWidth', 2); hold on
%     plot(save_times, sum(bac_saved.active(1:last_nonzero, 1:bac_saved.nBacs(last_nonzero)), 2));
%     title('Bacteria over time');
%     xlabel('Simulation time [h]');
%     ylabel('Number of bacteria');
    

    
    %% plot granule size over time
    granule_radius = zeros(last_nonzero, 1);
    
    for i = 1:last_nonzero
        nBac = bac_saved.nBacs(i);
        x = bac_saved.x(i, bac_saved.active(i, 1:nBac));
        y = bac_saved.y(i, bac_saved.active(i, 1:nBac));
        center_x = mean(x);
        center_y = mean(y);
        granule_radius(i) = max(sqrt((x - center_x).^2 + (y - center_y).^2));
    end
%     figure(3); clf;
%     plot(save_times, granule_radius*1e6, 'LineWidth', 2);
%     title('Granule size over time')
%     xlabel('Simulation time [h]')
%     ylabel('Granule radius [um]')
    
    i = TimeInput + 1; %days
    
    %% plot Layered Stratification
    
    % create bins for radial distance
    nBins = 50;
    bin_size = granule_radius(i-1) / nBins;
    
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
    strat_data_relative = strat_data ./ sum(strat_data, 2);
    
    layeredStratRelBac = figure('Name', 'Layered stratification (Relative bacteria)', 'Color', 'w', 'Position', [100, 100, 1200, 350], 'Visible', MATLABdraw_layered); clf;   
        coloring = {'#CC66FF', '#00B04F', '#FFA200', '#FF1482'}; %Old Colors
%         ar = area((1:nBins) * bin_size * 1e6, strat_data_relative * 100); % Granule radius
        rectangle('Position', [0 0 1 100], 'FaceColor', '#333333'); hold on
        ar = area(((1:nBins) * bin_size) ./ (nBins * bin_size), strat_data_relative * 100); % Dimensionless radius
        for s = 1:nSpecies
            ar(s).FaceColor = coloring{s};
        end
        axis([0 1 0 100])
        xticks(0:0.05:1)
        yticks(0:10:100)
        xtickformat('%.2f')
        ytickformat('percentage')
        set(gca, 'TickDir', 'out', 'Box', 'off')
        xlabel('R_{L}')
        ylabel('Relative abundance')
%         legend(constants.speciesNames);
        saveas(layeredStratRelBac,[folder_route, sprintf('%04d/', simulation_number), '/LayeredStrat_RelBac'],'tiffn')

        
    %% plot Columned Stratification
    
    % create bins for angular distance
    nBins = 180;
    nBins_column = nBins;
    bin_size = (360 / nBins); %º (d_theta)
    bin_size_column = bin_size;
    
    theta_rad = atan2((y - center_y), (x - center_x));                                                              %rad
    theta_deg = (theta_rad >= 0) .* theta_rad * (180 / pi) + (theta_rad < 0) .* (theta_rad + 2 * pi) * (180 / pi);	%º
    
    bin = ceil(theta_deg ./ bin_size);
    
    % create data table (nBins * nSpecies)
    nSpecies = length(constants.speciesNames);
    strat_data = zeros(nBins, nSpecies);
    BinRadius = zeros(nBins, 1);
    for b = 1:nBins
        for s = 1:nSpecies
            strat_data(b, s) = sum(bac_saved.species(i, 1:nBac) == s & bin == b & bac_saved.active(i, 1:nBac));
        end
        BinRadius(b) = max(radial_dist(bin == b)) * 10^6; %um
    end
    
    % calculate per bin the relative abundance of the bacteria
    strat_data_relative = strat_data ./ sum(strat_data, 2);
    
    columnedStratRelBac = figure('Name', 'Columned stratification (Relative bacteria)', 'Color', 'w', 'Position', [100, 100, 1200, 350], 'Visible', ATLABdraw_columned); clf;   
        coloring = {'#CC66FF', '#00B04F', '#FFA200', '#FF1482'}; %Old Colors
%         ar = area((1:nBins) * bin_size * 1e6, strat_data_relative * 100); % Granule radius
        rectangle('Position', [0 0 360 100], 'FaceColor', '#333333'); hold on
        ar = area(((0:nBins) * bin_size), [strat_data_relative(end, :); strat_data_relative] * 100); % Theta (º)
        for s = 1:nSpecies
            ar(s).FaceColor = coloring{s};
        end
        axis([0 360 0 100])
        xticks(0:10:360)
        yticks(0:10:100)
%         xtickformat('%.2f')
        ytickformat('percentage')
        set(gca, 'TickDir', 'out', 'Box', 'off')
        xlabel('\theta (º)')
        ylabel('Relative abundance')
%         legend(constants.speciesNames);
        saveas(columnedStratRelBac,[folder_route, sprintf('%04d/', simulation_number), '/ColumnedStrat_RelBac'],'tiffn')
    
    if BoxChartPr
        prompt = 'Do you want continue with BoxChart of Microcolonies? [Y = 1/N = 0]:\n>> ';
        fprintf('')
        r = input(prompt); %days
        
        if r == 1
            %% BoxChart of microcolonies            
            file_name = './BoxChart_microcolonies.xlsx';
            minTol = 0.9; %Minimum relative abundance to assume that is a unique species section
            prev_spc = nSpecies + 10;
            pos = [];
            spc = zeros(nBins_column,1);
            indx = [];
            for i = 1:nBins_column
                for s = 1:nSpecies
                    if strat_data_relative(i,s) >= minTol
                        spc(i) = s;
                        continue;
                    end
                end
                %- DEBUGGING -%
%                 if i >= 27 && i <= 29
%                     disp(i)
%                     disp(prev_spc)
%                     disp(spc(i))
%                 end
%                 if i >= 100 && i <= 103
%                     disp(i)
%                     disp(prev_spc)
%                     disp(spc(i))
%                 end
                %-------------%
                if spc(1) ~= 0 && prev_spc > nSpecies
                    pos = 1;
                end
                if prev_spc ~= spc(i) && prev_spc <= nSpecies
                    if prev_spc ~= 0 && spc(i) == 0
                        pos_corr = i - 1;
                    elseif prev_spc ~= 0 && spc(i) ~= 0
                        pos_corr = [i - 1 i];
                    else
                        pos_corr = i;
                    end
                    pos = [pos pos_corr]; %#ok<AGROW>
                end
                prev_spc = spc(i);
                
                indx = [indx; i]; %#ok<AGROW>
            end
            pos = [pos nBins_column];
            %- DEBUGGING -%
%             mPr = [strat_data_relative spc indx BinRadius];
%             disp(pos)
%             disp(mPr)
%             pause()
            %-------------%

            raw_result = [];
            perimeter_total = 2 * pi() * granule_radius(TimeInput) * (360 / 360) * 10^6; %um
            for iP = 1:2:length(pos)-1
                %- DEBUGGING -%
%                 disp(spc(pos(iP):pos(iP+1)))
                %-------------%
                UnSpc = unique(spc(pos(iP):pos(iP+1)));
                AvRadius = mean(BinRadius(pos(iP):pos(iP+1)));
                theta = (pos(iP+1) - pos(iP) + 1) * bin_size_column;
                perimeter = 2 * pi() * AvRadius * (theta / 360); %um
                rel_perimeter = perimeter/perimeter_total;
                raw_result = [raw_result; UnSpc AvRadius theta perimeter rel_perimeter]; %#ok<AGROW>
            end
            %% Continuity check of microcolony
            %- DEBUGGING -%
%             disp(raw_result)
%             pause()
            %-------------%
            if raw_result(1,1) == raw_result(end,1)
                raw_result(1,3:5) = raw_result(1,3:5) + raw_result(end,3:5);
                raw_result = raw_result(1:end-1,:);
            end
            %- DEBUGGING -%
%             disp(raw_result)
%             disp(perimeter_total)
%             disp(median(raw_result(:,end)))
            %-------------%
            writematrix(raw_result(:,1), file_name, 'Sheet', 'Species', 'WriteMode', 'append')
            writematrix(1*ones(size(raw_result(:,1))), file_name, 'Sheet', 'Concentration', 'WriteMode', 'append')
            writematrix(raw_result(:,4), file_name, 'Sheet', 'Perimeter', 'WriteMode', 'append')
            writematrix(raw_result(:,5), file_name, 'Sheet', 'Relative perimeter', 'WriteMode', 'append')
        end
        fprintf('>> Done!\n')
    end
end