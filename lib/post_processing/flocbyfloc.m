function flocbyfloc(simulation_number, nReplicates, finished)
    % Function to do the analysis of metabolisms/species floc by floc
    
    if isscalar(finished)
        finished = ones(nReplicates, 1) * finished;
    else
        assert(numel(finished) == nReplicates, 'Array of finished status should match the number of replicates');
    end
    
    addpath(genpath('lib')); % make every subfolder with functions accessible to the code

    %% colour set
    colors_species_raw = {'#E69F00','#56B4E9','#33b190','#F0E442','#0072B2','#D55E00'};
    %                      orange   light blue  green    yellow   dark blue    red
    species_per_color = { 'An-NRMX', 'CMX',    'NOB',    'AOB',    'NRMX',    'AMX'};
%     colors_qualitative = {'#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c'};
%     linestyles = {'-','--',':','-.'};

    %% Settings
    writeExcel = 1;
    plotMass = 0; 
    
    %% plot data
    [constants, grid, ~, bac_saved, ~, ~, ~, ~, ~, ~, ~, ~, ~] = load_data(simulation_number, finished(1));
    nSpecies = length(constants.speciesNames);
    
    species_index = zeros(nSpecies, 1);
    for s = 1:nSpecies
        species_index(s) = find(strcmp(species_per_color, constants.speciesNames{s}));
    end   
    colors_species = colors_species_raw(species_index); %#ok<NASGU>
    
%     floc2nBacs = determine_where_bacteria_in_flocs(1, grid, bac_saved);
    [floc2nBacs, floc_centers] = determine_where_bacteria_in_flocs(1, grid, bac_saved, 0);
    nFlocs = max(floc2nBacs);
    [r,~] = size(bac_saved.x);
    mass_matrix = single(zeros(r, nSpecies, nFlocs * nReplicates));
    
    for n = 1:nReplicates
        % load data for each replicate
        [~, grid, ~, bac_saved, ~, ~, ~, ~, ~, ~, ~, ~, ~] = load_data(simulation_number+n-1, finished(n));
        
        last_nonzero = find(bac_saved.nBacs ~= 0, 1, 'last');
        
        for t = 1:last_nonzero
            % Determine in which floc are all bacteria
            [floc2nBacs, floc_centers] = determine_where_bacteria_in_flocs(t, grid, bac_saved, floc_centers);
            floc2nBacs = floc2nBacs';
            nBacs = bac_saved.nBacs(t);
            masS = (constants.bac_rho * 4/3 * pi) * bac_saved.radius(t,1:nBacs).^3;
            actS = bac_saved.active(t,1:nBacs);
            speS = bac_saved.species(t,1:nBacs);
                        
            % Active mass for each species
            for f = 1:nFlocs
                for s = 1:nSpecies
                    mask = speS == s & actS == 1 & floc2nBacs == f;
                    mass = sum(masS(mask));
                    mass_matrix(t,s,f+(n-1)*nFlocs) = mass;
                end
            end
        end
    end
    
    if writeExcel
        write_excel(mass_matrix)
    end
    if plotMass
        save_times = (1:last_nonzero) * constants.dT_save / (7 * 24); % [weeks]
        plot_mass_by_floc(save_times*constants.dT_save/(7*24), mass_matrix, nFlocs, nReplicates, nSpecies, colors_species)
    end
    fprintf('>>> Done! \n')
    
    %-- TESTING --%
%     x = 0;
%     main_mass_matrix = mass_matrix;
%     while x == 1
%         prompt = '> Some jumping relative abundances? [Y/N]: ';
%         x = input(prompt);
%         if x == 1
%             correction(floc1, floc2, [species], TimeS)
%             
%             plot_mass_by_floc(save_times*constants.dT_save/(7*24), mass_matrix, nFlocs, nReplicates, nSpecies, colors_species)
%         else
%             fprintf('>> Done! \n')
%         end
%     end
    %-------------%
end

%% helper functions
function write_excel(mass_matrix)
    % Writing of results in Excel file
    path = 'RawData_flocbyfloc.xlsx';
    
    for i = 1:size(mass_matrix,3)
        % FlocByFloc
        %{
        sheet_name = sprintf('Floc %d', i);
        writematrix(mass_matrix(:,:,i), path, 'Sheet', sheet_name)
        %}
        % t_0 and t_max directly
        sheet_name = 't_0';
            writematrix(mass_matrix(1,:,i), path, 'Sheet', sheet_name, 'WriteMode', 'append')
        sheet_name = 't_max';
            last_sum_nonzero = find(sum(mass_matrix(:,:,i),2), 1, 'last');
            writematrix(mass_matrix(last_sum_nonzero,:,i), path, 'Sheet', sheet_name, 'WriteMode', 'append')
    end
end

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

function [floc2nBacs, floc_centers] = determine_where_bacteria_in_flocs(TimeS, grid, bac_saved, floc_center)
    % Create a vector with same order than bac_saved.x with reference to 
    % establish in what floc each bacterium reside in 
    nBacs = bac_saved.nBacs(TimeS);
    
    bacAUX = struct;
    bacAUX.x = bac_saved.x(TimeS,1:nBacs)'; bacAUX.y = bac_saved.y(TimeS,1:nBacs)';
    bacAUX.s = bac_saved.species(TimeS,1:nBacs)';
    floc2nBacs = zeros(size(bacAUX.x));
    
    [~, grid2nBacs] = determine_where_bacteria_in_grid(grid, bacAUX);
    
    kED = -1/8 * ones(3); kED(2,2) = 1;
    bacED = convn(grid2nBacs, kED, 'same') < 0; %logical ('outside' nodes of floc)
    [Bx, By] = find(bacED);
    
    unBx = unique(Bx) * grid.dx;
    unBy = unique(By) * grid.dy;
    
    diff_unBx = diff(unBx);
    diff_unBy = diff(unBy);
    
    potX = sort([unBx(1); unBx(diff_unBx > (1.1*grid.dx)); unBx(find(diff_unBx > (1.1*grid.dx)) + 1); unBx(end)]);
    potY = sort([unBy(1); unBy(diff_unBy > (1.1*grid.dy)); unBy(find(diff_unBy > (1.1*grid.dy)) + 1); unBy(end)]);
    potX_r = reshape(potX, [2, size(potX,1)/2])';
    potY_r = reshape(potY, [2, size(potY,1)/2])';
    
    if floc_center == 0
        c = 0;
    end
    for i = 1:size(potX_r, 1)
        for j = 1:size(potY_r, 1)
            bac_c = (bacAUX.x > potX_r(i,1))&(bacAUX.x < potX_r(i,2))&(bacAUX.y > potY_r(j,1))&(bacAUX.y < potY_r(j,2));
            if any(bac_c)
                if floc_center == 0
                    c = c + 1;
                    floc_centers(c,:) = [c, min(potX_r(i,:))+diff(potX_r(i,:))/2, min(potY_r(j,:))+diff(potY_r(j,:))/2]; %#ok<AGROW>
                else
                    %-- TESTING --%
%                     floc_centers = floc_center;
%                     diff_x = potX_r(i,:) - floc_centers(:,2);
%                     diff_y = potY_r(j,:) - floc_centers(:,3);
%                     [~, I] = min(sqrt((diff_x).^2+diff_y.^2));
%                     c = I(1);
                    %-------------%
                    
                    % update of center every time (considering the closest center)
                    new_center_X = min(potX_r(i,:))+diff(potX_r(i,:))/2;
                    new_center_Y = min(potY_r(j,:))+diff(potY_r(j,:))/2;
                    diff_c_x = new_center_X - floc_center(:,2);
                    diff_c_y = new_center_Y - floc_center(:,3);
                    [~,I_c] = min(sqrt((diff_c_x).^2+diff_c_y.^2));
                    floc_centers(I_c, 2:3) = [new_center_X, new_center_Y]; %#ok<AGROW>
                    
                    % c -> nearest center of inoculum flocs!
                    diff_x = potX_r(i,2) - floc_centers(:,2);
                    diff_y = potY_r(j,2) - floc_centers(:,3);
                    [~, I] = min(sqrt((diff_x).^2+diff_y.^2));
                    c = I;
                end
                floc2nBacs(bac_c) = c;
            end
        end
    end
end

function plot_mass_by_floc(time_plot, mass_matrix, nFlocs, nReplicates, nSpecies, colors_species) %#ok<DEFNU>
    % plot acumulative mass of active individuals per floc over time
    close all;
    
    for f = 1:nFlocs*nReplicates
        figure(900 + f);
        totalMass = sum(mass_matrix(:,:,f),2);
        relAbund = mass_matrix(:,:,f)./totalMass;
        
        for s = 1:nSpecies
            p = plot(time_plot, relAbund(1:end-1,s)*100, 'LineWidth', 2);
            p.Color = colors_species{s};
            ytickformat('percentage')
            ylabel('Relative abundances (wt. %)')
            ylim([0 100]);
            xlabel('Time (weeks)')
            set(gca,'TickDir','in');
            set(gcf,'color','w');
            hold on;
        end
    end
end