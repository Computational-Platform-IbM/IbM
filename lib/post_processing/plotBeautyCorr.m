function plotBeautyCorr
    % Function to do Spearman's rho analysis between species and plotting
    % Remember to paste all .xlsx file on IbM folder!!!!
    % Excel: Row -> nReplicates; Columns -> nSpecies
    % Do you need reorganize some Excel columns??

    format
    
    DataPlotting = 0; % Do u want to plot data? [Y:1/N:0]
    CorrMethod = 'Kendall'; % {'Pearson','Kendall','Spearman'}
    
    %% Excel name and path
    file_name = '1.1.[NH3];[NO2]=1;0_[O2]=1uM.xlsx';
%     file_name = '1.2.[NH3];[NO2]=1;0_[O2]=1.5uM.xlsx';
%     file_name = '1.3.[NH3];[NO2]=1;0_[O2]=3uM.xlsx';
%     file_name = '1.4.[NH3];[NO2]=1;0_[O2]=31.25uM+93.75uM.xlsx';
%     file_name = '2.1.[NH3];[NO2]=1;1_[O2]=1uM.xlsx';
%     file_name = '2.2.[NH3];[NO2]=1;1_[O2]=1.5uM.xlsx';
%     file_name = '2.3.[NH3];[NO2]=1;1_[O2]=3uM.xlsx';
%     file_name = '2.4.[NH3];[NO2]=1;1_[O2]=31.25uM+93.75uM.xlsx';
%     file_name = '3.1.[NH3];[NO2]=1.33;1_[O2]=1uM.xlsx';
%     file_name = '3.2.[NH3];[NO2]=1.33;1_[O2]=1.5uM.xlsx';
%     file_name = '3.3.[NH3];[NO2]=1.33;1_[O2]=3uM.xlsx';
%     file_name = '3.4.[NH3];[NO2]=1.33;1_[O2]=31.25uM+93.75uM.xlsx';
%     file_name = '4.1.[NH3];[NO2]=5;1_[O2]=1uM.xlsx';
%     file_name = '4.2.[NH3];[NO2]=5;1_[O2]=1.5uM.xlsx';
%     file_name = '4.3.[NH3];[NO2]=5;1_[O2]=3uM.xlsx';
%     file_name = '4.4.[NH3];[NO2]=5;1_[O2]=31.25uM+93.75uM.xlsx';
    fprintf(['>>> Excel: ' file_name '.xlsx \n\n'])
    
    % .mat order
%     nSpecies_name = {'AO', 'NO', 'AMX', 'CMX', 'NRMX', 'An-NRMX'};
    % Desired order
    nSpecies_name = {'AO', 'NO', 'CMX', 'NRMX', 'An-NRMX', 'AMX'};
    
    %% Read data (mass or relative abundances??)
    % t_0
    t0 = readmatrix(file_name, 'Sheet', 't_0');         % mass [g]
    t0_TM = sum(t0,2);                                  % total mass [g]
    t0_RA = (t0./t0_TM)*100;                            % Rel. abund. [%]
    % Reorganize columns 
    % From (AO,NO,AMX,CMX,NRMX,An-NRMX) to (AO,NO,CMX,NRMX,An-NRMX,AMX)
    t0_RA_aux = t0_RA;
    t0_RA(:,3) = t0_RA_aux(:,4);
    t0_RA(:,4:5) = t0_RA_aux(:,5:6);
    t0_RA(:,6) = t0_RA_aux(:,3);
    
    % t_max
    tmax = readmatrix(file_name, 'Sheet', 't_max');     % mass [g]
    tmax_TM = sum(tmax,2);                              % total mass [g]
    tmax_RA = (tmax./tmax_TM)*100;                      % Rel. abund. [%]
    % Reorganize columns
    % From (AO,NO,AMX,CMX,NRMX,An-NRMX) to (AO,NO,CMX,NRMX,An-NRMX,AMX)
    tmax_RA_aux = tmax_RA;
    tmax_RA(:,3) = tmax_RA_aux(:,4);
    tmax_RA(:,4:5) = tmax_RA_aux(:,5:6);
    tmax_RA(:,end) = tmax_RA_aux(:,3);
    
    [r0, c0] = size(t0_RA); [rmax, cmax] = size(tmax_RA);
    if r0 == rmax && c0 == cmax        
        %% Matrices definition
        nSpecies = c0;
        rhoM_aux = zeros(nSpecies+1);
        pvalM_aux = zeros(nSpecies+1);
        nM_aux = zeros(nSpecies+1);
        
        if DataPlotting
            cFlippedPlot = zeros(size(tmax,1),2);
            iFP = 1;
        end
        
        %% Correlation analysis
        for s1 = 1:nSpecies
            for s2 = 1:nSpecies
                s1_inocula = t0_RA(:,s1);
                s2_inocula = t0_RA(:,s2);
                check_inocula = s1_inocula > 0 & s2_inocula > 0;
                s1_final = tmax_RA(:,s1); s1_final = s1_final(check_inocula);
                s2_final = tmax_RA(:,s2); s2_final = s2_final(check_inocula);

                [rho,pval] = corr(s1_final, s2_final, 'Type', CorrMethod); 
                
                rhoM_aux(s1,s2) = rho;
                pvalM_aux(s1,s2) = pval;
                nM_aux(s1,s2) = length(s1_final);
                
                %% Data plotting
                if DataPlotting && s1 ~= s2 && ne(isnan(rhoM_aux(s1,s2)),1)
                    cFP = ismember(cFlippedPlot(1:iFP,:), [s1 s2], 'rows') | ismember(cFlippedPlot(1:iFP,:), [s2 s1], 'rows');
                    if nnz(cFP) == 0
                        cFlippedPlot(iFP,:) = [s1 s2];
                        iFP = iFP + 1;
                        if strcmp(CorrMethod,'Kendall')
                            TitleDataPlotting = sprintf('%s vs %s; 𝜏 = %0.4f; p-value = %0.4f; n = %d)',nSpecies_name{s1},nSpecies_name{s2},rhoM_aux(s1,s2),pvalM_aux(s1,s2),nM_aux(s1,s2));
                        else
                            TitleDataPlotting = sprintf('%s vs %s; ρ = %0.4f; p-value = %0.4f; n = %d)',nSpecies_name{s1},nSpecies_name{s2},rhoM_aux(s1,s2),pvalM_aux(s1,s2),nM_aux(s1,s2));
                        end
                        figure()
                        plot(s1_final, s2_final, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
                        title(TitleDataPlotting)
                        xlabel(sprintf('Active relative abundance of %s (%%)',nSpecies_name{s1}))
                        ylabel(sprintf('Active relative abundance of %s (%%)',nSpecies_name{s2}))
                        set(gca,'TickDir','out');
                        set(gcf,'color','w');
                    end
                end
            end
        end
        %% Cleaning of correlation data
        rhoM = tril(rhoM_aux, -1);
        rhoM(isnan(rhoM)) = 0;
        rhoM(logical(eye(size(rhoM)))) = 23;
        pvalM = tril(pvalM_aux, -1);
        pvalM(isnan(pvalM)) = 23.23;
        pvalM(logical(eye(size(pvalM)))) = 23;
        nM = tril(nM_aux, -1);
        nM(isnan(nM)) = 0;
        nM(logical(eye(size(nM)))) = 23.23;
        gM = 23*eye(size(pvalM));
                
        %% Plotting settings
        vec = [       -1;     -0.75;      -0.5;     -0.25;         0;      0.25;       0.5;      0.75;         1];
        hex = ['#DF8239'; '#E9A06B'; '#F7BA8E'; '#F7D7C0'; '#FFFFFF'; '#F4FBF4'; '#D1E9D2'; '#B3D5B4'; '#8FC891'];
        raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
        map = interp1(vec, raw, linspace(-1,1,1000),'pchip');
        map = [map; 235/255 235/255 235/255];
        
        [X,Y] = meshgrid(1:nSpecies+1);
        
        %% Correlation test plotting
        figure(); hold on;
        % Surface plot
        s = pcolor(X,flip(Y),rhoM);
        s.EdgeColor = 'none';
        colormap(map)
        caxis([-1 1])

        % Plot Borders
        plot(1:nSpecies+1,ones(1,nSpecies+1),'k')
        plot(ones(1,nSpecies+1),1:nSpecies+1,'k')
        plot(1:nSpecies+1,(nSpecies+1)*ones(1,nSpecies+1),'k')
        plot((nSpecies+1)*ones(1,nSpecies+1),1:nSpecies+1,'k')
        
        % "Diagonal" lines
        for i = 1:nSpecies
            %xLines
            plot([i i+1], ((nSpecies+1)-(i-1))*ones(1,2), 'k')
            %yLines
            plot(((nSpecies+1)-(i-1))*ones(1,2), [i i+1], 'k')
        end
        
        % p-value labels     
        pVal = pvalM(1:nSpecies,1:nSpecies);
        plabels = zeros(size(pVal)); plabels = string(plabels);
        plabels(pVal > 0.05) = "";
        plabels(pVal < 0.05 & pVal > 0.01) = "▪";
        plabels(pVal < 0.01 & pVal > 0.001) = "▪▪";
        plabels(pVal < 0.001) = "▪▪▪";
        plabels(pVal == 0) = "";
        plabels(pVal == 23.23) = "˟";
        
        for i = 1:nSpecies
            for j = 1:nSpecies
                if plabels(j,i) == "˟"
                    pY = 0.43;
                    FS = 19;
                else
                    pY = 0.51;
                    FS = 16;
                end
                text(i+0.5, nSpecies-(j-1)+pY, plabels(j,i), 'HorizontalAlignment', 'center', 'FontSize', FS)
            end
        end
        
        % plot settings
        axis([1 7 1 7])
        xticks(1.5:1:nSpecies+0.5)
        xticklabels(nSpecies_name);
        yticks(1.5:1:nSpecies+0.5)
        yticklabels(flip(nSpecies_name));
        set(gca,'TickDir','out');
        set(gcf,'color','w');
        
        %% Sample size plot
        figure(); hold on;
        % Surface plot
        s = pcolor(X,flip(Y),gM);
        s.EdgeColor = [235/255 235/255 235/255];
        s.LineWidth = 2;
        colormap(map)
        caxis([-1 1])
        
        % "Upper" white lines
        for i = 3:nSpecies
            %xLines
            plot(i:nSpecies+1, (nSpecies-(i-3))*ones(nSpecies-(i-3+1)), 'w', 'LineWidth', 2)
            %yLines
            plot((nSpecies-(i-3))*ones(nSpecies-(i-3+1)), i:nSpecies+1, 'w', 'LineWidth', 2)
        end
        
        % Plot Borders
        plot(1:nSpecies+1,ones(1,nSpecies+1),'k')
        plot(ones(1,nSpecies+1),1:nSpecies+1,'k')
        plot(1:nSpecies+1,(nSpecies+1)*ones(1,nSpecies+1),'k')
        plot((nSpecies+1)*ones(1,nSpecies+1),1:nSpecies+1,'k')
        
        % "Diagonal" lines
        for i = 1:nSpecies
            %xLines
            plot([i i+1], ((nSpecies+1)-(i-1))*ones(1,2), 'k')
            %yLines
            plot(((nSpecies+1)-(i-1))*ones(1,2), [i i+1], 'k')
        end
        
        % Sample size labels
        for i = 1:nSpecies
            for j = 1:nSpecies
                if nM(j,i) ~= 23.23 && nM(j,i) ~= 0
                    text(i+0.5, nSpecies-(j-1)+0.5, string(nM(j,i)), 'HorizontalAlignment', 'center', 'FontSize', 16)
                else
                    continue;
                end                
            end
        end
        
        % plot settings
        axis([1 7 1 7])
        xticks(1.5:1:nSpecies+0.5)
        xticklabels(nSpecies_name);
        yticks(1.5:1:nSpecies+0.5)
        yticklabels(flip(nSpecies_name));
        set(gca,'TickDir','out');
        set(gcf,'color','w');
        
        %% Display Command Window
        fprintf('>>> ρ values: (''23'' depicts diagonal correlation)\n')
        disp(rhoM(1:nSpecies,1:nSpecies))
        fprintf('>>> p-value: \n')
        disp(pvalM(1:nSpecies,1:nSpecies))
        fprintf('>>> Labels p-value: \n')
        disp(plabels)
        fprintf('>>> Sample sizes: \n')
        disp(nM(1:nSpecies,1:nSpecies))
    else
       fprintf('>>> Matrix dimensions between t_0 and t_max DO NOT match! \n') 
    end
    
    fprintf('>>> Done! \n\n')
end