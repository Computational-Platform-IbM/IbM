function draw
%%%%%% VISUALIZATION: Generates plots of cells and substrates
load('ResultsSim.mat')
load('R.mat')

fprintf('\n VISUALIZATION >>>> \n')

prompt = ['\n>> Graphics (default [0]):\n [1] - Bacteria\n', ...
                                        ' [2] - Histogram\n', ...
                                        ' [3] - Average of Substrates, Products and Biomass with time\n', ...
                                        ' [4] - Diffusion of Substrates or Products in the liquid phase\n', ... 
                                        ' [5] - pH\n', ... 
                                        ' [8] - Stratification proof \n', ...
                                        ' [9] - CC & DD matrices \n\n', ... 
                                        ' [0] - Nothing - End of the program\n\n>> Answer: '];
x = input(prompt);

if isempty(x)
    x = 0;
end

while x ~= 0
    if x == 1
        draw_bacteria(R, B)
    elseif x == 2
        draw_biomass(B,  All_StatesVar(:,1))
    elseif x == 3
        draw_substrate2(R,  All_StatesVar, B) %#ok<*NODEF>
    elseif x == 4
        draw_substrate1(R, SSxy, B, All_StatesVar(:,1))
    elseif x == 5
        draw_pH(R, pHs, B, All_StatesVar(:,1))
    elseif x == 8
        draw_stratification(R)
    elseif x == 9
        draw_CCDD(R, B)
    end
    x = input(prompt);
end
fprintf('\n>>>>>>>>> END of the program. See you later alligator !\n')
return

end

function draw_bacteria(R, B)

n_sp = R.St.numX;
CM = [0.6 0.2 1; 0 0.8 0;1 0.6 0;1 0 0.4]; %[AOB, Nitrob, Nitros, AMX];

 writerObjB = VideoWriter('bacteria.avi');
 writerObjB.FrameRate = 2.3;
 open(writerObjB);

for k=1:n_sp % Assignation of a colour to a microbial species
    fprintf('\n Microbial Species %d: %s  %d bacteria Activas: %d.>> COLOUR: ',k,char(R.rm.rNamesX(k)), R.bac.bac_ns(k),  sum((R.bac.atrib(:,5)==k).*(R.bac.atrib(:,4)> 0)))
    cprintf(CM(k,:),['GROUP - ',char(R.rm.rNamesX(k)),'\n']);
end
fprintf('\nTOTAL Number of cells: %d\n', R.bac.bac_n)

fprintf('\n>>> Bacterial growth with time ... ')

% Plot and video of bacteria
fB=figure('Name','bacterias','Position',get(0, 'Screensize'),'Color','w');
  
maxix = R.Sxy.maxxSys*1e6; %um
maxiy = R.Sxy.maxySys*1e6; %um
 minix = 0;
 miniy = 0;
%     for i = size(B,3)-1:1:size(B,3)
%     for i = 2:1:size(B,3)
    i = size(B,3); %no for loop
%     i = 120/4 + 2; %no for loop
        figure(fB)
        clf;
        Time = B(1,1,i)/24;
        bac_x = B(:,2,i)*1e6; %um
        bac_y = B(:,3,i)*1e6; %um
        bac_s = B(:,4,i);
        bac_r = B(:,5,i)*1e6; %um
        bac_a = B(:,6,i);
        act = B(:,7,i);
    
    for k=1:length(bac_x)          
        if bac_s(k)>0
            BM = max(bac_a(bac_s == bac_s(k)));
            if  BM <= 0
                BM = 1;
            end
            bacA = bac_a(k)/BM;
%             bacA = 1;
            if Time == 0
                bacA = 1;
            end
            if act(k) == 1
                col = CM(bac_s(k),:);
            else
                col = [0 0 0];
%                 col = CM(bac_s(k),:);
            end
            rectangle('Curvature',[1 1],'Position',[bac_x(k)-bac_r(k) bac_y(k)-bac_r(k) 2*bac_r(k) 2*bac_r(k)],'EdgeColor',col,'FaceColor',(bacA> 0)*bacA*col); %'LineWidth',1,
%             rectangle('Curvature',[1 1],'Position',[bac_x(k)-bac_r(k) bac_y(k)-bac_r(k) 2*bac_r(k) 2*bac_r(k)],'EdgeColor','k','LineWidth',1,'FaceColor',(bacA> 0)*bacA*col); %'LineWidth',1,
        end
    end

    axis equal; axis([minix maxix miniy maxiy]);
%     axis equal; axis([minix maxix miniy 50]); %um
%     numBac = sum(ne(bac_s,0));
    xlabel('x (\mum)','FontSize',20) 
    ax = gca;
    ax.FontSize = 20;
    ylabel('y (\mum)','FontSize',20)
    title([...
        '\color[rgb]{' sprintf('%1.1f,%1.1f,%1.1f', CM(1,:)) '} ' char(R.rm.rNamesX(1)), ...
        '\color[rgb]{' sprintf('%1.1f,%1.1f,%1.1f', CM(2,:)) '} ' char(R.rm.rNamesX(2)), ...
        '\color[rgb]{' sprintf('%1.1f,%1.1f,%1.1f', CM(3,:)) '} ' char(R.rm.rNamesX(3)), ...
        '\color[rgb]{' sprintf('%1.1f,%1.1f,%1.1f', CM(4,:)) '} ' char(R.rm.rNamesX(4))], ...
        'FontSize',24, 'FontWeight', 'bold', 'Interpreter', 'tex');
    subtitle({[num2str(Time),' days. ']},'FontSize',20)
     drawnow;
    frame = getframe(fB);
    writeVideo(writerObjB,frame);

%     end %for i = ...
    writeVideo(writerObjB,frame);
    writeVideo(writerObjB,frame);
    writeVideo(writerObjB,frame);
    writeVideo(writerObjB,frame);
    close(writerObjB);
    fprintf('\nTotal size of the space discretization x: %.3e y: %.3e', R.Sxy.maxxSys, R.Sxy.maxySys)
    fprintf('\nTotal size occupied by bacteria x: %.3e y: %e', maxix - minix, maxiy - miniy)
    fprintf('\nDiscretization size Dx: %.3e Dy: %.3e', R.Sxy.dx, R.Sxy.dx)
    fprintf('\nDiscretization steps in the biofilm Dx: %.3f Dy: %.3f\n', (maxix - minix)/R.Sxy.dx, (maxiy - miniy)/R.Sxy.dx)
    fprintf('END\n')

end
function draw_CCDD(R, B)        
    %- b_layer() -%
    Sxy = R.Sxy;
    
    xSys = 0:Sxy.dx:Sxy.maxxSys;                       
    ySys = 0:Sxy.dy:Sxy.maxySys;

    dx = Sxy.dx; dy = Sxy.dy;
    
    pos_minx = 1; pos_miny = 1;
    maxy = Sxy.maxy + dy + Sxy.T_blayer;
%     maxy = Sxy.maxy;
    pos_maxy = find((ySys >= maxy),1);
    if isempty(pos_maxy)
        pos_maxy = length(ySys);
    end

    x = xSys;               %biofilm + BDL space               
    y = ySys(1:pos_maxy);   %biofilm + BDL space
    nx = length(x)-1+2;
    ny = length(y)-1+2;

    %Check if it is necessary construct the Laplacian matrices again
    ax = zeros(Sxy.nxSys-2,1);
    ax(pos_minx:pos_minx+nx-3) = 1;
    ay = zeros(Sxy.nySys-2,1);
    ay(pos_miny:pos_miny+ny-3) = 1;
    Sxy.pos_xySys = sparse(kron(ax, ay));
    xp = x + dx/2;
    x_pn = kron(xp(1:end-1)', ones(ny-2,1));
    yp = y + dy/2;
    y_pn = kron(ones(nx-2,1),yp(1:end-1)');
    
    %- CC & DD construction + Figure -%    
    nT = R.Sxy.nT;
    nx = R.Sxy.nxSys - 2;
    ny = R.Sxy.nySys - 2;
    bac_n = R.bac.bac_n;
    
    aux_x = x_pn - R.Sxy.dx/2; aux_x1 = x_pn + R.Sxy.dx/2;
    aux_y = y_pn - R.Sxy.dy/2; aux_y1 = y_pn + R.Sxy.dy/2;
    DD = sparse(R.Sxy.nT, 1);
    CC = sparse(R.bac.bac_n, R.Sxy.nT);

    bac_x = R.bac.atrib(:,1);
    bac_y = R.bac.atrib(:,2);
%     bac_r = R.bac.atrib(:,6);
    T_blayer = R.Sxy.T_blayer;
    sM = 0; %Smoothness level
    
    Time = B(1,1,end)/24;
    fDD = figure('Name','DD','Position',get(0, 'Screensize'));
    figure(fDD)
    
    %Boundary layer construction
    My_aux = zeros(nx,1);
    x_aux = reshape(aux_x,[ny,nx]);     x_aux = x_aux(1,:);
    x1_aux = reshape(aux_x1,[ny,nx]);   x1_aux = x1_aux(1,:);
    
    for i = 1:nx
        c1 = ((bac_x) > x_aux(i)).*((bac_x) <= x1_aux(i));
        Mnorm_nTi = max(c1.*(bac_y)) + T_blayer;
        Mnorm_nTb = 0; Mnorm_nTf = 0;
        if ne(sM,0)
            Mnorm_nTb_aux = Mnorm_nTb; Mnorm_nTf_aux = Mnorm_nTf;
            for s = 1:sM
                if (i - s) > 0
                    c1b = ((bac_x) > x_aux(i-s)).*((bac_x) <= x1_aux(i-s));
                    Mnorm_nTb_aux = max(c1b.*(bac_y)) + T_blayer;
                end
                if (i + s) < nx
                    c1f = ((bac_x) > x_aux(i+s)).*((bac_x) <= x1_aux(i+s));
                    Mnorm_nTf_aux = max(c1f.*(bac_y)) + T_blayer;
                end
                Mnorm_nTb = max([Mnorm_nTb Mnorm_nTb_aux]);
                Mnorm_nTf = max([Mnorm_nTf Mnorm_nTf_aux]);
            end
        end
        %Only biofilm space is included
        Mnorm_nT = max([Mnorm_nTi Mnorm_nTb Mnorm_nTf]);
        My_aux(i,1) = Mnorm_nT;
        %All space is included
%         Mnorm_nT = R.Sxy.maxy;
    end
    My = kron(My_aux,ones(ny,1));
    
    %CC & DD construction
    for i = 1:R.Sxy.nT
        c1 = ((bac_x) > aux_x(i)).*((bac_x) <= aux_x1(i));
        c2 = ((bac_y) > aux_y(i)).*((bac_y) <= aux_y1(i));
        vv = (c1.*c2); vv(vv==1) = find(vv);
        ci = sparse(1: bac_n, i, vv, bac_n, nT); %sparce(position i, position j, value, rows sparse, columns sparse);
        di = sparse(i, 1, My(i,1) > aux_y(i), nT, 1); %value = 1 or 0 -> 1: inside biofilm + BDL | 0: outside biofilm + BDL
        DD = DD + di;
        CC = CC + ci;

        if ne(sum(vv),0)
            col = sum(full(di)).*[0 0.8 0];
        else
            col = sum(full(di)).*[1 1 1];
        end
        rectangle('Position',[aux_x(i) aux_y(i) R.Sxy.dx R.Sxy.dy],'FaceColor',col)
        title(['Numbac:',num2str(R.bac.bac_n),'.Time:',num2str(Time)])
        axis([0 R.Sxy.maxxSys 0 R.Sxy.maxySys]);
    end
    drawnow;
    saveas(fDD,[num2str(Time),'.jpg']);
    
    fprintf('\nDone!\n');
end
function draw_biomass(B, time)

figure('Name','Activity Prop')
bW = 1e-6;
[~, I]=max(max(B(:,3,:)));
[~,edges, bin] = histcounts(B(:,3,I), 'BinWidth', bW);
act1 = zeros(max(bin),1);
act2 = zeros(max(bin),1);
d_str = zeros(size(B,3));
bac_sTot = zeros(size(B,3));
bac_s1 = zeros(size(B,3));
bac_s2 = zeros(size(B,3));
for i = 1:1:size(B,3)
    bac_x = B(:,2,i);
    bac_y = B(:,3,i);
    bac_norm = sqrt(bac_x.*bac_x+bac_y.*bac_y);
    bac_s = B(:,4,i);
    bac_a = B(:,6,i).*(B(:,6,i)> 0);
    a1 = bac_a(bac_s == 1)/sum(bac_a(bac_s == 1));
    b1 = bac_norm(bac_s == 1);
    a2 = bac_a(bac_s == 2)/sum(bac_a(bac_s == 2));
    b2 = bac_norm(bac_s == 2);
    [N1,~, bin1] = histcounts(b1, 'BinWidth',bW);
    [~,~, bin2] = histcounts(b2, 'BinWidth',bW);
    if ne(N1,0)    
       act1(1:max(bin1)) = accumarray(bin1,a1);
       act2(1:max(bin2))  = accumarray(bin2,a2);
       k1 = 0;acc1 = 0;
       [sb1, I1] = sort(b1);
       sa1 = a1(I1);
       while acc1 < 0.5
           k1 = k1 +1;
           acc1 = sum(sa1(1:k1));
       end
       x = (edges(1)+bW/2):bW:edges(end);
       subplot(3,1,1)
       ax1 = gca;
       bar(x , act1,'r')
      hold on
      plot(sb1(k1)*ones(11,1), 0:ax1.YLim(2)/10:ax1.YLim(2), 'k')
        hold off
       subplot(3,1,2)
       k2 = 0;acc2 = 0;
       [sb2, I2] = sort(b2);
       sa2 = a2(I2);
       while acc2 < 0.5
           k2 = k2 +1;
           acc2 = sum(sa2(1:k2));
       end
       ax2 = gca;
       bar(x , act2,'c')
       hold on
      plot(sb2(k2)*ones(11,1), 0:ax2.YLim(2)/10:ax2.YLim(2), 'k')
        hold off
       subplot(3,1,3)
       plot(x , act1,'r')
       hold on
       plot(x , act2,'c')
       hold off
       drawnow; 
       d_str(i) = (sb1(k1)- sb2(k2))/(1e-6);
       bac_s1(i) = sum(bac_s == 1);
       bac_s2(i) = sum(bac_s == 2);
       bac_sTot(i) = bac_s1(i)+bac_s2(i);
    end
end
aux = ne(bac_s1,0) | ne(bac_s2,0); 
time = time(aux);
d_str = d_str(aux);
bac_s1 = bac_s1(aux);
bac_s2 = bac_s2(aux);
bac_sTot = bac_sTot(aux);
figure('Name','Stratification Distance')
plot(time,d_str)
ylabel('y \mum')
xlabel('time (h)')
figure('Name','Number bacteria')
plot(time,bac_s1, 'r')
hold on
plot(time,bac_s2, 'c')
hold on
plot(time,bac_sTot, 'k')
hold off
ylabel('y \mum')
xlabel('time (h)')
fprintf('\nPlease pulse any key to continue >>>>\n')
pause
figure('Name','Number of cells')
[~, I]=max(max(B(:,3,:)));
[~,edges, ~] = histcounts(B(:,3,I), 'BinWidth', bW);

for i = 1:size(B,3)
    bac_x = B(:,2,i);
    bac_y = B(:,3,i);
    bac_norm = sqrt(bac_x.*bac_x+bac_y.*bac_y);
    bac_s = B(:,4,i);
    bac_a = B(:,6,i);
    bac_y1 = bac_norm(bac_s == 1);
    bac_y2 = bac_norm(bac_s == 2);
    act1 = (bac_a(bac_s == 1)/max(bac_a(bac_s == 1)))> 1e-6;
    act2 = (bac_a(bac_s == 2)/max(bac_a(bac_s == 2)))> 1e-6;
    b1 =  bac_y1(act1);
    b2 =  bac_y2(act2);
    [N1,~, ~] = histcounts(b1, edges);
    [N2,~, ~] = histcounts(b2, edges);
    if ne(sum(N1),0)
       x = (edges(1)+bW/2):bW:edges(end);
       subplot(3,1,1)
       bar(x , N1,'r')
       subplot(3,1,2)
       bar(x , N2,'c')
       subplot(3,1,3)
       plot(x , 100*(N1/sum(bac_s == 1)),'r')
       hold on
       plot(x , 100*(N2/sum(bac_s == 2)),'c')
       hold off
       drawnow; 
    end
end

   
%    [M I]= sort(bac_y(bac_s == 1))
%    histogram(a(I))
%    pause


%     h1 = bac_y(bac_s == 1);
%     h2 = bac_y(bac_s == 2);
%     H1 = histogram(h1,10,'Normalization','probability','FaceColor','r');
%     hold on
%     H2 = histogram(h2,10,'Normalization','probability','FaceColor','b');
%     hold off
%     drawnow;
%     [~, I] = max(H1.Values);
%     Y1(i) = H1.BinEdges (I);
%     [~, I] = max(H2.Values);
%     Y2(i) = H2.BinEdges (I);

%     H1.bin
     


% plot(time, Y1,'r')
% hold on
% plot(time, Y2,'b')
% hold off
% legend('AOB', 'NOB')
end

function draw_substrate1(R, SSxy, B, time)

    fprintf('\n>>> Diffusion of substrates')
    nT = size(SSxy,3);
    SSxy = SSxy*(14*1000);
%     Si = reshape(SSxy(:,j,:),R.Sxy.nySys-2,nT);
    S1 = reshape(SSxy(:,1,:),R.Sxy.nySys-2,nT);
    S2 = reshape(SSxy(:,2,:),R.Sxy.nySys-2,nT);
    S3 = reshape(SSxy(:,3,:),R.Sxy.nySys-2,nT);
    S4 = reshape(SSxy(:,4,:)*(32/14),R.Sxy.nySys-2,nT);
    S5 = reshape(SSxy(:,5,:)*(44/14),R.Sxy.nySys-2,nT);
    S6 = reshape(SSxy(:,6,:)*(96/14),R.Sxy.nySys-2,nT);
    S7 = reshape(SSxy(:,7,:)*(70/14),R.Sxy.nySys-2,nT);
%     Smax = max(R.Sxy.Sbc_Dir(j),max(max(Si)));
     Smax1 = max(max(max([S1;S2;S3;S4])));
     Smax2 = max(max(max([S5;S6;S7])));
     Smax3 = max(max(max(S4)));
    Smax1 = Smax1 + Smax1*0.1;
    Smax2 = Smax2 + Smax2*0.1;
    SmaxO2 = Smax3 + Smax3*0.1;
    writerObjS = VideoWriter('diffusion.avi');
    open(writerObjS);
    miniy = 0;
    maxiy = R.Sxy.maxySys;
    
    ySys = 0:R.Sxy.dy:(R.Sxy.maxySys);
 
    yp = (ySys + R.Sxy.dy/2);
    y_pnSys = yp(1:end-1)';

    fS=figure('Name','Difusion of substrates and products','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 50 200 500 500 ],'Color','w');
    for i = 2:nT
        figure(fS)
        bac_s = B(:,4,i);
        numBac = sum(ne(bac_s,0));
        bac_y = B(1:numBac,3,i);
        minbac_y = min(bac_y);
        maxbac_y = max(bac_y);
        
        plot(y_pnSys*1e6,S1(:,i), 'b');       %plotting the field variable
        hold on
        plot(y_pnSys*1e6,S2(:,i),'r'); 
        hold on
        plot(y_pnSys*1e6,S3(:,i),'y'); 
        hold on
        plot(y_pnSys*1e6,S4(:,i),'k'); 
        hold on  
        plot(minbac_y*1e6,R.Sxy.Sbc_Dir(1)*(14*1000), 'k*')
        hold on
        plot(maxbac_y*1e6,R.Sxy.Sbc_Dir(1)*(14*1000), 'k*')
        hold off      
        axis ([miniy*1e6 maxiy*1e6 0 Smax1])
        title({'2-D Diffusion of substrates';['time (\itt) = ',num2str(time(i))]})
        xlabel('radius ({\mu}m)')
        ylabel('Concentration (mgN/L // mgO2/L)')
        legend('NH_{3}','NO_{2}^{-}','NO_{3}^{-}', 'O_{2}')
        drawnow; 
%           pause
        frame = getframe(fS);
        writeVideo(writerObjS,frame);
    end
    close(writerObjS);  
    
    writerObjS1 = VideoWriter('diffusion_O2.avi');
    open(writerObjS1);
    
    fS1=figure('Name','Difusion of substrates and products','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 50 200 500 500 ],'Color','w');

    for i = 2:nT
        figure(fS1)
        bac_s = B(:,4,i);
        numBac = sum(ne(bac_s,0));
        bac_y = B(1:numBac,3,i);
        minbac_y = min(bac_y);
        maxbac_y = max(bac_y);
        
        plot(y_pnSys*1e6,S4(:,i),'k'); 
        hold on
%         plot(y_pnSys*1e6,S6(:,i),'c'); 
%         hold on
%         plot(y_pnSys*1e6,S7(:,i),'g'); 
%         hold on
        plot(minbac_y*1e6,R.Sxy.Sbc_Dir(4)*(32*1000), 'k*')
        hold on
        plot(maxbac_y*1e6,R.Sxy.Sbc_Dir(4)*(32*1000), 'k*')
        hold off
        axis ([miniy*1e6 maxiy*1e6 0 SmaxO2])
        
        title({'2-D Diffusion of substrates';['time (\itt) = ',num2str(time(i))]})
        xlabel('radius ({\mu}m)')
        ylabel('Concentration (mgO2/L)')
        legend('O_{2}')
%         legend('CO_{2}','SO_{4}^{2-}', 'Na^{+}')
        drawnow; 
         frame = getframe(fS1);
        writeVideo(writerObjS1,frame);
    end
    close(writerObjS1);  
    
    fprintf('\nEND\n')
end
function draw_pH(R, pH, B, time)

    fprintf('\n>>> Diffusion of substrates')
    nT = size(pH,2);
   
    Smax = 14;
    writerObjS = VideoWriter('pH.avi');
    open(writerObjS);
    miniy = 0;
    maxiy = R.Sxy.maxySys;
    
    ySys = 0:R.Sxy.dy:R.Sxy.maxySys;   
   
    yp = (ySys + R.Sxy.dy/2);
    y_pnSys = yp(1:end-1);

    fS=figure('Name','Difusion of substrates and products','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 50 200 500 500 ],'Color','w');
    for i = 2:nT
        figure(fS)
        bac_s = B(:,4,i);
        numBac = sum(ne(bac_s,0));
        bac_y = B(1:numBac,3,i);
        minbac_y = min(bac_y);
        maxbac_y = max(bac_y);
        
        plot(y_pnSys*1e6,pH(:,i), 'k');       %plotting the field variable
        
        hold on
        plot(minbac_y*1e6,R.pOp.pH, 'k*')
        hold on
        plot(maxbac_y*1e6,R.pOp.pH, 'k*')
        hold off
        axis ([miniy*1e6 maxiy*1e6 0 Smax])
        
        title({'pH';['time (\itt) = ',num2str(time(i))]})
        xlabel('radius ({\mu}m)')
        ylabel('pH')
        drawnow; 
        frame = getframe(fS);
        writeVideo(writerObjS,frame);
    end
    close(writerObjS);  
    fprintf('\nEND\n')
end

function draw_substrate2(R,  All_StatesVar, B)
time = All_StatesVar(:,1);
    fprintf('\n>>> Evolution of Substrates and Products with time ... ')
    % Plot of liquid substrates
    fL = figure('Name','LIQUID Chemicals','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 50 200 500 500 ],'Color','w');
    figure(fL)
    n_Liq = sum(strcmp(R.St.Phase,'L'));
    l = R.St.StNames(1:n_Liq,1);
    MW = [14*1000 14*1000 14*1000 32*1000 45*1000 96*1000];
    for i = 1: n_Liq
        plot(time, All_StatesVar(:,i+1)*MW(i))
        hold on
    end
    legend(l);
%     savefig(fL,'SPLiquids.fig')

    % Plot of Biofilm density
    fDensity = figure('Name','Density of the biofilm','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 75 200 500 500 ],'Color','w');
    figure(fDensity)
    nt = {'bioflm density g/m3'};
    plot(time, All_StatesVar(:,end));
    legend(nt);

    % Plot of Ntotal
    fNtotal = figure('Name','Ntotal','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 75 200 500 500 ],'Color','w');
    figure(fNtotal)
    nt = {'Ntotal','NH3 Inf'};
    plot(time, All_StatesVar(:,end-1));
    hold on
    plot(time, R.Inf.St(1)*ones(size(time)), 'r');
    hold off
    legend(nt);
    
    
    % Plot of HRT
    HRT = figure('Name','HRT','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 75 200 500 500 ],'Color','w');
    figure(HRT)
    nhrt = 'HRT';
    plot(time, All_StatesVar(:,end-2));
    legend(nhrt);

    % Plot of liquid substrates
    fP = figure('Name','CHARGES Chemicals','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 75 200 500 500 ],'Color','w');
    figure(fP)
    n_Ph = sum(strcmp(R.St.Phase,'P'));
    p = R.St.StNames(n_Liq+1:n_Liq+n_Ph,1);
    for i = 1: n_Ph
        plot(time, All_StatesVar(:,n_Liq+i+1))
        hold on
    end
    legend(p);
%     savefig(fP,'SPCharges.fig')
   
    % Plot of gas substrates
    fG=figure('Name','GAS Chemicals','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 100 200 500 500 ],'Color','w');
    figure(fG)
    n_Gas = sum(strcmp(R.St.Phase,'G'));
    g = R.St.StNames(n_Liq+n_Ph+1:n_Liq+n_Ph+n_Gas,1);
    for i = 1: n_Gas
        plot(time, All_StatesVar(:,n_Liq+n_Ph+i+1))
        hold on
    end
    legend(g);
%     savefig(fG,'SPGas.fig')
    
    % Plot of biomass
    fX=figure('Name','Biomass overall','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 100 200 500 500 ],'Color','w');
    figure(fX)
    n_X = sum(strcmp(R.St.Phase,'S'));
    x = R.St.StNames(n_Liq+n_Ph+n_Gas+1:n_Liq+n_Ph+n_Gas+n_X,1);
    for i = 1: n_X
        plot(time, All_StatesVar(:,n_Liq+n_Ph+n_Gas+i+1))
        hold on
    end
    legend(x);
    
    % Plot of biomass
    fX=figure('Name','Biomass Number','PaperPosition',[0.05 0.0 3.5 3.5],'Position',[ 100 200 500 500 ],'Color','w');
    figure(fX)
    n_X = sum(strcmp(R.St.Phase,'S'));
    x = R.St.StNames(n_Liq+n_Ph+n_Gas+1:n_Liq+n_Ph+n_Gas+n_X,1);
    ns = zeros(size(B,3),n_X);
    for j = 1:size(B,3)
        bac_s = B(:,4,j);
        for i = 1: n_X
            ns(j,i) = sum(bac_s == i);
        end
    end
    for i = 1: n_X
        plot(time, ns(:,i))
        hold on
    end
    legend(x);
    

    fprintf('END\n') 
end


function draw_stratification(R)
figure('Name', 'Stratification?')

bac_y = R.bac.atrib(:,2);
bac_s = R.bac.atrib(:,5);

vv1 = bac_y(bac_s==1);
vv2 = bac_y(bac_s==2);
ySys = 0:R.Sxy.dy:R.Sxy.maxySys;

pp1 = zeros(size(ySys));
pp2 = zeros(size(ySys));
pp1a = zeros(size(ySys));
pp2a = zeros(size(ySys));
bac_a = R.bac.atrib(:,4);
for i = 2: length(ySys)
    pp1(i-1) = sum((vv1 > ySys(i-1)).*(vv1 < ySys(i)));
    pp2(i-1) = sum((vv2 > ySys(i-1)).*(vv2 < ySys(i)));
    pp1a(i-1) = sum(bac_a(((vv1 > ySys(i-1)).*(vv1 < ySys(i))>0)))/pp1(i-1);
    pp2a(i-1) =sum(bac_a(((vv2 > ySys(i-1)).*(vv2 < ySys(i))>0)))/pp2(i-1);
end
Y = 1e6*ySys(ne((pp1+pp2),0));
% i_maxAOB = find(pp1 == max(pp1));
% i_maxNOB = find(pp2 == max(pp2));

pp1n = pp1(ne((pp1+pp2),0));
pp2n = pp2(ne((pp1+pp2),0));
% pp1na = pp1a(ne((pp1+pp2),0));
% pp2na = pp2a(ne((pp1+pp2),0));
stem(Y, pp1n,'Color', 'r', 'LineWidth',0.2,'LineStyle','-','Marker','o','MarkerSize',3,'MarkerEdgeColor','r','MarkerFaceColor','r')
hold on
stem(Y, pp2n,'Color', 'b', 'LineWidth',0.2,'LineStyle','-','Marker','o','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b')
xlabel('y \mum')
ylabel('Number of individuals')
axis([ Y(1) Y(end) 0 max([pp1,pp2])])
legend('AOB','NOB')
% subplot(2,1,2);

figure('Name','Activity over y coordinate')
[Y1, I]=sort(bac_y);
bac_mu_max = R.bac.atrib(:,7);
mu_rel = bac_a./bac_mu_max;
mu_rel = mu_rel(I);
mu_rel1 = mu_rel(bac_s==1);
mu_rel2 = mu_rel(bac_s==2);
plot(Y1(bac_s==1), mu_rel1, 'r')
hold on
plot(Y1(bac_s==2), mu_rel2, 'b')
axis([Y1(1) Y1(end) 0 1])
xlabel('y \mum')
ylabel('\mu/\mu^{max}')
legend('AOB','NOB')


% stem(Y, pp1na,'Color', 'r', 'LineWidth',0.2,'LineStyle','-','Marker','o','MarkerSize',3,'MarkerEdgeColor','r','MarkerFaceColor','r')
% hold on
% stem(Y, pp2na,'Color', 'b', 'LineWidth',0.2,'LineStyle','-','Marker','o','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b')
% xlabel('y \mum')
% ylabel('\mu/\mu^{max}')
% axis([ Y(1) Y(end) 0 1])
% legend('AOB','NOB')

% %%% 2 Oxygen profile
% figure()
% SO2 = reshape(Sxy(:,4,end),R.Sxy.nySys-2,R.Sxy.nxSys-2,1);
% SO2_v = SO2(:,round(end/2))*32*1000;
% SO2_v = SO2_v(ne((pp1+pp2),0));
% plot (Y, SO2_v, 'k')
% hold on
% plot(Y(i_maxAOB(1)), SO2_v(i_maxAOB(1)),'r*','MarkerSize',14)
% hold on
% plot(Y(i_maxNOB(end)), SO2_v(i_maxNOB(end)),'b*','MarkerSize',14)
% hold on
% Yaob = 1e6*R.Sxy.ySys(ne((pp1),0));
% SO2_vaob = SO2_v(ne((pp1),0));
% plot(Yaob(end), SO2_vaob(end) ,'r')
% hold on
% Ynob = 1e6*R.Sxy.ySys(ne((pp2),0));
% SO2_vnob = SO2_v(ne((pp2),0));
% plot(Ynob(end), SO2_vnob(end),'b')
% hold off
% legend('[O_{2}]','Postion (y) of Maximum number of AOB','Postion (y) of Maximum number of NOB', 'Maximum Postion (y) for AOB in the biofilm', 'Maximum Postion (y) for NOB in the biofilm','Location','northwest')
% axis([Y(1) Y(end) 0 max(SO2_v)])
% xlabel('y \mum')
% ylabel('[O_{2}] mg/L')


figure()
area(Y, pp1n, 'FaceColor', 'r')
hold on
area(Y, pp2n, 'FaceColor', 'b')
hold off
legend('AOB','NOB')
xlabel('y \mum')
ylabel('Number of individuals')
% figure()
% y_d = pp2n.*(pp2n<=pp1n) + pp1n.*(pp1n<pp2n);
% area(Y, y_d, 'FaceColor', 'g')
% area_int  = 100*trapz(y_d)/(min([trapz(pp1n), trapz(pp2n)]));
% fprintf('\nSolapamiento %.2f \n', area_int)

end

function count = cprintf(style,format,varargin)

%%%%%%% This function only generates the necessary code to be able to write with
%%%%%%% colours (necessary for the visualization)



% CPRINTF displays styled formatted text in the Command Window
%
% Syntax:
%    count = cprintf(style,format,...)
%
% Description:
%    CPRINTF processes the specified text using the exact same FORMAT
%    arguments accepted by the built-in SPRINTF and FPRINTF functions.
%
%    CPRINTF then displays the text in the Command Window using the
%    specified STYLE argument. The accepted styles are those used for
%    Matlab's syntax highlighting (see: File / Preferences / Colors / 
%    M-file Syntax Highlighting Colors), and also user-defined colors.
%
%    The possible pre-defined STYLE names are:
%
%       'Text'                 - default: black
%       'Keywords'             - default: blue
%       'Comments'             - default: green
%       'Strings'              - default: purple
%       'UnterminatedStrings'  - default: dark red
%       'SystemCommands'       - default: orange
%       'Errors'               - default: light red
%       'Hyperlinks'           - default: underlined blue
%
%       'Black','Cyan','Magenta','Blue','Green','Red','Yellow','White'
%
%    STYLE beginning with '-' or '_' will be underlined. For example:
%          '-Blue' is underlined blue, like 'Hyperlinks';
%          '_Comments' is underlined green etc.
%
%    STYLE beginning with '*' will be bold (R2011b+ only). For example:
%          '*Blue' is bold blue;
%          '*Comments' is bold green etc.
%    Note: Matlab does not currently support both bold and underline,
%          only one of them can be used in a single cprintf command. But of
%          course bold and underline can be mixed by using separate commands.
%
%    STYLE also accepts a regular Matlab RGB vector, that can be underlined
%    and bolded: -[0,1,1] means underlined cyan, '*[1,0,0]' is bold red.
%
%    STYLE is case-insensitive and accepts unique partial strings just
%    like handle property names.
%
%    CPRINTF by itself, without any input parameters, displays a demo
%
% Example:
%    cprintf;   % displays the demo
%    cprintf('text',   'regular black text');
%    cprintf('hyper',  'followed %s','by');
%    cprintf('key',    '%d colored', 4);
%    cprintf('-comment','& underlined');
%    cprintf('err',    'elements\n');
%    cprintf('cyan',   'cyan');
%    cprintf('_green', 'underlined green');
%    cprintf(-[1,0,1], 'underlined magenta');
%    cprintf([1,0.5,0],'and multi-\nline orange\n');
%    cprintf('*blue',  'and *bold* (R2011b+ only)\n');
%    cprintf('string');  % same as fprintf('string') and cprintf('text','string')
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab
%    functionality. It works on Matlab 7+, but use at your own risk!
%
%    A technical description of the implementation can be found at:
%    <a href="http://undocumentedmatlab.com/blog/cprintf/">http://UndocumentedMatlab.com/blog/cprintf/</a>
%
% Limitations:
%    1. In R2011a and earlier, a single space char is inserted at the
%       beginning of each CPRINTF text segment (this is ok in R2011b+).
%
%    2. In R2011a and earlier, consecutive differently-colored multi-line
%       CPRINTFs sometimes display incorrectly on the bottom line.
%       As far as I could tell this is due to a Matlab bug. Examples:
%         >> cprintf('-str','under\nline'); cprintf('err','red\n'); % hidden 'red', unhidden '_'
%         >> cprintf('str','regu\nlar'); cprintf('err','red\n'); % underline red (not purple) 'lar'
%
%    3. Sometimes, non newline ('\n')-terminated segments display unstyled
%       (black) when the command prompt chevron ('>>') regains focus on the
%       continuation of that line (I can't pinpoint when this happens). 
%       To fix this, simply newline-terminate all command-prompt messages.
%
%    4. In R2011b and later, the above errors appear to be fixed. However,
%       the last character of an underlined segment is not underlined for
%       some unknown reason (add an extra space character to make it look better)
%
%    5. In old Matlab versions (e.g., Matlab 7.1 R14), multi-line styles
%       only affect the first line. Single-line styles work as expected.
%       R14 also appends a single space after underlined segments.
%
%    6. Bold style is only supported on R2011b+, and cannot also be underlined.
%
% Change log:
%    2015-06-24: Fixed a few discoloration issues (some other issues still remain)
%    2015-03-20: Fix: if command window isn't defined yet (startup) use standard fprintf as suggested by John Marozas
%    2012-08-09: Graceful degradation support for deployed (compiled) and non-desktop applications; minor bug fixes
%    2012-08-06: Fixes for R2012b; added bold style; accept RGB string (non-numeric) style
%    2011-11-27: Fixes for R2011b
%    2011-08-29: Fix by Danilo (FEX comment) for non-default text colors
%    2011-03-04: Performance improvement
%    2010-06-27: Fix for R2010a/b; fixed edge case reported by Sharron; CPRINTF with no args runs the demo
%    2009-09-28: Fixed edge-case problem reported by Swagat K
%    2009-05-28: corrected nargout behavior suggested by Andreas Gäb
%    2009-05-13: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>
%
% See also:
%    sprintf, fprintf

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.10 $  $Date: 2015/06/24 01:29:18 $

  persistent majorVersion minorVersion
  if isempty(majorVersion)
      %v = version; if str2double(v(1:3)) <= 7.1
      %majorVersion = str2double(regexprep(version,'^(\d+).*','$1'));
      %minorVersion = str2double(regexprep(version,'^\d+\.(\d+).*','$1'));
      %[a,b,c,d,versionIdStrs]=regexp(version,'^(\d+)\.(\d+).*');  %#ok unused
      v = sscanf(version, '%d.', 2);
      majorVersion = v(1); %str2double(versionIdStrs{1}{1});
      minorVersion = v(2); %str2double(versionIdStrs{1}{2});
  end

  % The following is for debug use only:
  %global docElement txt el
  if ~exist('el','var') || isempty(el),  el=handle([]);  end 
  if nargin<1, showDemo(majorVersion,minorVersion); return;  end
  if isempty(style),  return;  end
  if all(ishandle(style)) && length(style)~=3
      dumpElement(style);
      return;
  end

  % Process the text string
  if nargin<2, format = style; style='text';  end
  %error(nargchk(2, inf, nargin, 'struct'));
  %str = sprintf(format,varargin{:});

  % In compiled mode
  try useDesktop = usejava('desktop'); catch, useDesktop = false; end
  if isdeployed | ~useDesktop %#ok<OR2> - for Matlab 6 compatibility
      % do not display any formatting - use simple fprintf()
      % See: http://undocumentedmatlab.com/blog/bold-color-text-in-the-command-window/#comment-103035
      % Also see: https://mail.google.com/mail/u/0/?ui=2&shva=1#all/1390a26e7ef4aa4d
      % Also see: https://mail.google.com/mail/u/0/?ui=2&shva=1#all/13a6ed3223333b21
      count1 = fprintf(format,varargin{:});
  else
      % Else (Matlab desktop mode)
      % Get the normalized style name and underlining flag
      [underlineFlag, boldFlag, style, debugFlag] = processStyleInfo(style);

      % Set hyperlinking, if so requested
      if underlineFlag
          format = ['<a href="">' format '</a>'];

          % Matlab 7.1 R14 (possibly a few newer versions as well?)
          % have a bug in rendering consecutive hyperlinks
          % This is fixed by appending a single non-linked space
          if majorVersion < 7 || (majorVersion==7 && minorVersion <= 1)
              format(end+1) = ' ';
          end
      end

      % Set bold, if requested and supported (R2011b+)
      if boldFlag
          if (majorVersion > 7 || minorVersion >= 13)
              format = ['<strong>' format '</strong>'];
          else
              boldFlag = 0;
          end
      end

      % Get the current CW position
      cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
      lastPos = cmdWinDoc.getLength;

      % If not beginning of line
      bolFlag = 0;  %#ok
      %if docElement.getEndOffset - docElement.getStartOffset > 1
          % Display a hyperlink element in order to force element separation
          % (otherwise adjacent elements on the same line will be merged)
          if majorVersion<7 || (majorVersion==7 && minorVersion<13)
              if ~underlineFlag
                  fprintf('<a href=""> </a>');  %fprintf('<a href=""> </a>\b');
              elseif format(end)~=10  % if no newline at end
                  fprintf(' ');  %fprintf(' \b');
              end
          end
          %drawnow;
          bolFlag = 1;
      %end

      % Get a handle to the Command Window component
      mde = com.mathworks.mde.desk.MLDesktop.getInstance;
      cw = mde.getClient('Command Window');

      % Fix: if command window isn't defined yet (startup), use standard fprintf()
      if (isempty(cw))
         count1 = fprintf(format,varargin{:});
         if nargout
             count = count1;
         end
         return;
      end
      
      xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);

      % Store the CW background color as a special color pref
      % This way, if the CW bg color changes (via File/Preferences), 
      % it will also affect existing rendered strs
      com.mathworks.services.Prefs.setColorPref('CW_BG_Color',xCmdWndView.getBackground);

      % Display the text in the Command Window
      % Note: fprintf(2,...) is required in order to add formatting tokens, which
      % ^^^^  can then be updated below (no such tokens when outputting to stdout)
      count1 = fprintf(2,format,varargin{:});

      % Repaint the command window
      %awtinvoke(cmdWinDoc,'remove',lastPos,1);   % TODO: find out how to remove the extra '_'
      drawnow;  % this is necessary for the following to work properly (refer to Evgeny Pr in FEX comment 16/1/2011)
      xCmdWndView.repaint;
      %hListeners = cmdWinDoc.getDocumentListeners; for idx=1:numel(hListeners), try hListeners(idx).repaint; catch, end, end

      docElement = cmdWinDoc.getParagraphElement(lastPos+1);
      if majorVersion<7 || (majorVersion==7 && minorVersion<13)
          if bolFlag && ~underlineFlag
              % Set the leading hyperlink space character ('_') to the bg color, effectively hiding it
              % Note: old Matlab versions have a bug in hyperlinks that need to be accounted for...
              %disp(' '); dumpElement(docElement)
              setElementStyle(docElement,'CW_BG_Color',1+underlineFlag,majorVersion,minorVersion); %+getUrlsFix(docElement));
              %disp(' '); dumpElement(docElement)
              el(end+1) = handle(docElement);  % #ok used in debug only
          end

          % Fix a problem with some hidden hyperlinks becoming unhidden...
          fixHyperlink(docElement);
          %dumpElement(docElement);
      end

      % Get the Document Element(s) corresponding to the latest fprintf operation
      while docElement.getStartOffset < cmdWinDoc.getLength
          % Set the element style according to the current style
          if debugFlag, dumpElement(docElement); end
          specialFlag = underlineFlag | boldFlag;
          setElementStyle(docElement,style,specialFlag,majorVersion,minorVersion);
          if debugFlag, dumpElement(docElement); end
          docElement2 = cmdWinDoc.getParagraphElement(docElement.getEndOffset+1);
          if isequal(docElement,docElement2),  break;  end
          docElement = docElement2;
      end
      if debugFlag, dumpElement(docElement); end

      % Force a Command-Window repaint
      % Note: this is important in case the rendered str was not '\n'-terminated
      xCmdWndView.repaint;

      % The following is for debug use only:
      el(end+1) = handle(docElement);  %#ok used in debug only
      %elementStart  = docElement.getStartOffset;
      %elementLength = docElement.getEndOffset - elementStart;
      %txt = cmdWinDoc.getText(elementStart,elementLength);
  end

  if nargout
      count = count1;
  end
  return;  % debug breakpoint
end

% Process the requested style information
function [underlineFlag,boldFlag,style,debugFlag] = processStyleInfo(style)
  underlineFlag = 0;
  boldFlag = 0;
  debugFlag = 0;

  % First, strip out the underline/bold markers
  if ischar(style)
      % Styles containing '-' or '_' should be underlined (using a no-target hyperlink hack)
      %if style(1)=='-'
      underlineIdx = (style=='-') | (style=='_');
      if any(underlineIdx)
          underlineFlag = 1;
          %style = style(2:end);
          style = style(~underlineIdx);
      end

      % Check for bold style (only if not underlined)
      boldIdx = (style=='*');
      if any(boldIdx)
          boldFlag = 1;
          style = style(~boldIdx);
      end
      if underlineFlag && boldFlag
          warning('YMA:cprintf:BoldUnderline','Matlab does not support both bold & underline')
      end

      % Check for debug mode (style contains '!')
      debugIdx = (style=='!');
      if any(debugIdx)
          debugFlag = 1;
          style = style(~debugIdx);
      end

      % Check if the remaining style sting is a numeric vector
      %styleNum = str2num(style); %#ok<ST2NM>  % not good because style='text' is evaled!
      %if ~isempty(styleNum)
      if any(style==' ' | style==',' | style==';')
          style = str2num(style); %#ok<ST2NM>
      end
  end

  % Style = valid matlab RGB vector
  if isnumeric(style) && length(style)==3 && all(style<=1) && all(abs(style)>=0)
      if any(style<0)
          underlineFlag = 1;
          style = abs(style);
      end
      style = getColorStyle(style);

  elseif ~ischar(style)
      error('YMA:cprintf:InvalidStyle','Invalid style - see help section for a list of valid style values')

  % Style name
  else
      % Try case-insensitive partial/full match with the accepted style names
      matlabStyles = {'Text','Keywords','Comments','Strings','UnterminatedStrings','SystemCommands','Errors'};
      validStyles  = [matlabStyles, ...
                      'Black','Cyan','Magenta','Blue','Green','Red','Yellow','White', ...
                      'Hyperlinks'];
      matches = find(strncmpi(style,validStyles,length(style)));

      % No match - error
      if isempty(matches)
          error('YMA:cprintf:InvalidStyle','Invalid style - see help section for a list of valid style values')

      % Too many matches (ambiguous) - error
      elseif length(matches) > 1
          error('YMA:cprintf:AmbigStyle','Ambiguous style name - supply extra characters for uniqueness')

      % Regular text
      elseif matches == 1
          style = 'ColorsText';  % fixed by Danilo, 29/8/2011

      % Highlight preference style name
      elseif matches <= length(matlabStyles)
          style = ['Colors_M_' validStyles{matches}];

      % Color name
      elseif matches < length(validStyles)
          colors = [0,0,0; 0,1,1; 1,0,1; 0,0,1; 0,1,0; 1,0,0; 1,1,0; 1,1,1];
          requestedColor = colors(matches-length(matlabStyles),:);
          style = getColorStyle(requestedColor);

      % Hyperlink
      else
          style = 'Colors_HTML_HTMLLinks';  % CWLink
          underlineFlag = 1;
      end
  end
end
% Convert a Matlab RGB vector into a known style name (e.g., '[255,37,0]')
function styleName = getColorStyle(rgb)
  intColor = int32(rgb*255);
  javaColor = java.awt.Color(intColor(1), intColor(2), intColor(3));
  styleName = sprintf('[%d,%d,%d]',intColor);
  com.mathworks.services.Prefs.setColorPref(styleName,javaColor);
end
% Fix a bug in some Matlab versions, where the number of URL segments
% is larger than the number of style segments in a doc element
function delta = getUrlsFix(docElement)  %#ok currently unused
  tokens = docElement.getAttribute('SyntaxTokens');
  links  = docElement.getAttribute('LinkStartTokens');
  if length(links) > length(tokens(1))
      delta = length(links) > length(tokens(1));
  else
      delta = 0;
  end
end
% fprintf(2,str) causes all previous '_'s in the line to become red - fix this
function fixHyperlink(docElement)
  try
      tokens = docElement.getAttribute('SyntaxTokens');
      urls   = docElement.getAttribute('HtmlLink');
      urls   = urls(2);
      links  = docElement.getAttribute('LinkStartTokens');
      offsets = tokens(1);
      styles  = tokens(2);
      doc = docElement.getDocument;

      % Loop over all segments in this docElement
      for idx = 1 : length(offsets)-1
          % If this is a hyperlink with no URL target and starts with ' ' and is collored as an error (red)...
          if strcmp(styles(idx).char,'Colors_M_Errors')
              character = char(doc.getText(offsets(idx)+docElement.getStartOffset,1));
              if strcmp(character,' ')
                  if isempty(urls(idx)) && links(idx)==0
                      % Revert the style color to the CW background color (i.e., hide it!)
                      styles(idx) = java.lang.String('CW_BG_Color');
                  end
              end
          end
      end
  catch
      % never mind...
  end
end
% Set an element to a particular style (color)
function setElementStyle(docElement,style,specialFlag, majorVersion,minorVersion)
  %global tokens links urls urlTargets  % for debug only
  global oldStyles
  if nargin<3,  specialFlag=0;  end
  % Set the last Element token to the requested style:
  % Colors:
  tokens = docElement.getAttribute('SyntaxTokens');
  try
      styles = tokens(2);
      oldStyles{end+1} = cell(styles);

      % Correct edge case problem
      extraInd = double(majorVersion>7 || (majorVersion==7 && minorVersion>=13));  % =0 for R2011a-, =1 for R2011b+
      %{
      if ~strcmp('CWLink',char(styles(end-hyperlinkFlag))) && ...
          strcmp('CWLink',char(styles(end-hyperlinkFlag-1)))
         extraInd = 0;%1;
      end
      hyperlinkFlag = ~isempty(strmatch('CWLink',tokens(2)));
      hyperlinkFlag = 0 + any(cellfun(@(c)(~isempty(c)&&strcmp(c,'CWLink')),cell(tokens(2))));
      %}

      jStyle = java.lang.String(style);
      if numel(styles)==4 && isempty(char(styles(2)))
          % Attempt to fix discoloration issues - NOT SURE THAT THIS IS OK! - 24/6/2015
          styles(1) = jStyle;
      end
      styles(end-extraInd) = java.lang.String('');
      styles(end-extraInd-specialFlag) = jStyle;  % #ok apparently unused but in reality used by Java
      if extraInd
          styles(end-specialFlag) = jStyle;
      end

      oldStyles{end} = [oldStyles{end} cell(styles)];
  catch
      % never mind for now
  end
  
  % Underlines (hyperlinks):
  %{
  links = docElement.getAttribute('LinkStartTokens');
  if isempty(links)
      %docElement.addAttribute('LinkStartTokens',repmat(int32(-1),length(tokens(2)),1));
  else
      %TODO: remove hyperlink by setting the value to -1
  end
  %}

  % Correct empty URLs to be un-hyperlinkable (only underlined)
  urls = docElement.getAttribute('HtmlLink');
  if ~isempty(urls)
      urlTargets = urls(2);
      for urlIdx = 1 : length(urlTargets)
          try
              if urlTargets(urlIdx).length < 1
                  urlTargets(urlIdx) = [];  % '' => []
              end
          catch
              % never mind...
              a=1;  %#ok used for debug breakpoint...
          end
      end
  end
  
  % Bold: (currently unused because we cannot modify this immutable int32 numeric array)
  %{
  try
      %hasBold = docElement.isDefined('BoldStartTokens');
      bolds = docElement.getAttribute('BoldStartTokens');
      if ~isempty(bolds)
          %docElement.addAttribute('BoldStartTokens',repmat(int32(1),length(bolds),1));
      end
  catch
      % never mind - ignore...
      a=1;  %#ok used for debug breakpoint...
  end
  %}
  
  return;  % debug breakpoint
end
% Display information about element(s)
function dumpElement(docElements)
  %return;
  disp(' ');
  numElements = length(docElements);
  cmdWinDoc = docElements(1).getDocument;
  for elementIdx = 1 : numElements
      if numElements > 1,  fprintf('Element #%d:\n',elementIdx);  end
      docElement = docElements(elementIdx);
      if ~isjava(docElement),  docElement = docElement.java;  end
      %docElement.dump(java.lang.System.out,1)
      disp(docElement)
      tokens = docElement.getAttribute('SyntaxTokens');
      if isempty(tokens),  continue;  end
      links = docElement.getAttribute('LinkStartTokens');
      urls  = docElement.getAttribute('HtmlLink');
      try bolds = docElement.getAttribute('BoldStartTokens'); catch, bolds = []; end
      txt = {};
      tokenLengths = tokens(1);
      for tokenIdx = 1 : length(tokenLengths)-1
          tokenLength = diff(tokenLengths(tokenIdx+[0,1]));
          if (tokenLength < 0)
              tokenLength = docElement.getEndOffset - docElement.getStartOffset - tokenLengths(tokenIdx);
          end
          txt{tokenIdx} = cmdWinDoc.getText(docElement.getStartOffset+tokenLengths(tokenIdx),tokenLength).char;  %#ok
      end
      lastTokenStartOffset = docElement.getStartOffset + tokenLengths(end);
      try
          txt{end+1} = cmdWinDoc.getText(lastTokenStartOffset, docElement.getEndOffset-lastTokenStartOffset).char; %#ok
      catch
          txt{end+1} = ''; %#ok<AGROW>
      end
      %cmdWinDoc.uiinspect
      %docElement.uiinspect
      txt = strrep(txt',newline,'\n');
      try
          data = [cell(tokens(2)) m2c(tokens(1)) m2c(links) m2c(urls(1)) cell(urls(2)) m2c(bolds) txt];
          if elementIdx==1
              disp('    SyntaxTokens(2,1) - LinkStartTokens - HtmlLink(1,2) - BoldStartTokens - txt');
              disp('    ==============================================================================');
          end
      catch
          try
              data = [cell(tokens(2)) m2c(tokens(1)) m2c(links) txt];
          catch
              disp([cell(tokens(2)) m2c(tokens(1)) txt]);
              try
                  data = [m2c(links) m2c(urls(1)) cell(urls(2))];
              catch
                  % Mtlab 7.1 only has urls(1)...
                  data = [m2c(links) cell(urls)];
              end
          end
      end
      disp(data)
  end
end
% Utility function to convert matrix => cell
function cells = m2c(data)
  %datasize = size(data);  cells = mat2cell(data,ones(1,datasize(1)),ones(1,datasize(2)));
  cells = num2cell(data);
end
% Display the help and demo
function showDemo(majorVersion,minorVersion)
  fprintf('cprintf displays formatted text in the Command Window.\n\n');
  fprintf('Syntax: count = cprintf(style,format,...);  click <a href="matlab:help cprintf">here</a> for details.\n\n');
  url = 'http://UndocumentedMatlab.com/blog/cprintf/';
  fprintf(['Technical description: <a href="' url '">' url '</a>\n\n']);
  fprintf('Demo:\n\n');
  boldFlag = majorVersion>7 || (majorVersion==7 && minorVersion>=13);
  s = ['cprintf(''text'',    ''regular black text'');' 10 ...
       'cprintf(''hyper'',   ''followed %s'',''by'');' 10 ...
       'cprintf(''key'',     ''%d colored'',' num2str(4+boldFlag) ');' 10 ...
       'cprintf(''-comment'',''& underlined'');' 10 ...
       'cprintf(''err'',     ''elements:\n'');' 10 ...
       'cprintf(''cyan'',    ''cyan'');' 10 ...
       'cprintf(''_green'',  ''underlined green'');' 10 ...
       'cprintf(-[1,0,1],  ''underlined magenta'');' 10 ...
       'cprintf([1,0.5,0], ''and multi-\nline orange\n'');' 10];
   if boldFlag
       % In R2011b+ the internal bug that causes the need for an extra space
       % is apparently fixed, so we must insert the sparator spaces manually...
       % On the other hand, 2011b enables *bold* format
       s = [s 'cprintf(''*blue'',   ''and *bold* (R2011b+ only)\n'');' 10];
       s = strrep(s, ''')',' '')');
       s = strrep(s, ''',5)',' '',5)');
       s = strrep(s, '\n ','\n');
   end
   disp(s);
   eval(s);

end
%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%
% - Fix: Remove leading space char (hidden underline '_')
% - Fix: Find workaround for multi-line quirks/limitations
% - Fix: Non-\n-terminated segments are displayed as black
% - Fix: Check whether the hyperlink fix for 7.1 is also needed on 7.2 etc.
% - Enh: Add font support
