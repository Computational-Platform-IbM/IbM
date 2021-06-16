% %%%% bacteria.m -Generates the IbM: Position of the cells and division
% %%%% R.bac is the structure that keeps all the information related with this
% subscript

function [R, StVLiq, Fdiv] = bacteria(L, G, X, R)
St = R.St;
Sxy = R.Sxy;
bac = R.bac;
numStVLiq = St.numStVLiq;
numStVLiq2 = St.numStVLiq2;

U = zeros(St.numSt,1);
for k = 1: St.numSt
    if k <= numStVLiq2
         U(k) = Sxy.Sbc_Dir(k);
    elseif k <= numStVLiq 
        U(k) = G(Sxy.nT*((k-numStVLiq2)-1)+1);
    elseif k > numStVLiq
        U(k) = sum((X.*(bac.atrib(:,5)==(k-numStVLiq)))/Sxy.Vg)/Sxy.nT;
    end
end
% %%%% Update of the variables 
St.StV = U;
St.StVLiq = St.StV(1:numStVLiq);% Liquid and gas variables
St.StVLiq2 = St.StV(1:numStVLiq2);% Liquid variables
St.StVX = St.StV(numStVLiq+1:end);% Biomass
StVLiq = [L; G];
R.bac.atrib(:,3) = X; 
R.bac.atrib(:,6) = ((X*bac.bac_MW/bac.bac_rho)*(3/(4*pi))).^(1/3);

R.Sxy = Sxy; 
R.St = St;
[R.bac, Fdiv] = bac_division(R.bac, Sxy.maxxSys, Sxy.maxySys);
end

function [bac, Fdiv] = bac_division(bac, maxxSys, maxySys)
bac_n = bac.bac_n;

Fdiv = 0;
bac_m = bac.atrib(:,3)*bac.bac_MW; %Using grams from here
j = 0;
iffdiv = (sum(bac_m > bac.bac_mmax) >= 1)+(sum(bac_m < bac.bac_mmin) >= 1);
if ne(iffdiv,0)
    nold = bac_n;
    Fdiv = 1; 
    bac_ns = bac.bac_ns;
    bac_x = bac.atrib(:,1);
    bac_y = bac.atrib(:,2);
    bac_a = bac.atrib(:,4);
    bac_s = bac.atrib(:,5);
    bac_r = bac.atrib(:,6);
    act = bac.atrib(:,8);
%     bac_mu_max = bac.atrib(:,7);
    bac_yield = bac.atrib(:,7);
    bac_Ks = bac.bac_Ks;
%     bac_maint = bac.atrib(:,9);
    while sum(bac_m > bac.bac_mmax) >= 1 %|| sum(bac_m < bac.bac_mmin) >= 1
        for i=1:nold
            if bac_m(i-j) > bac.bac_mmax
                % choose an angle for the division
                fi = (1 + (-1-1).*rand(1))*2*pi;
                % increase number of cells with one
                bac_n = bac_n + 1;
                bac_ns(bac_s(i-j)) = bac_ns(bac_s(i-j)) + 1;

                % coordinates, species and radius of the new cell (bac_n)
                bac_x(bac_n,1) = bac_x(i-j) + bac_r(i-j)*cos(fi); % 0x coordinates for the daughter
                bac_y(bac_n,1) = bac_y(i-j) + bac_r(i-j)*sin(fi); % 0y coordinates for the daughter  
                while bac_x(bac_n,1) < bac.bac_rmax 
                    fi = (1 + (-1-1).*rand(1))*2*pi;
                    bac_x(bac_n,1)  = bac_x(i-j) + 2*bac_r(i-j)*cos(fi);
                end  
                while bac_y(bac_n,1) < bac.bac_rmax 
                    fi = (1 + (-1-1).*rand(1))*2*pi;
                    bac_y(bac_n,1)  = bac_y(i-j) + 2*bac_r(i-j)*cos(fi);
                end  
                bac_s(bac_n,1)  = bac_s(i-j);                    % same cell type  
                bac_a(bac_n,1) = bac_a(i-j);
                act(bac_n,1) = 1;
%                 bac_mu_max(bac_n,1)  = bac_mu_max(i-j);                    % same mu max
                bac_yield(bac_n,1)  = bac_yield(i-j);
                bac_Ks(bac_n,:)  = bac_Ks(i-j, :);                    % same Ks
%                 bac_maint(bac_n,1)  = bac_maint(i-j);                    % same maintenance
                bac_m(bac_n,1)  = (0.45 + 0.1*rand(1,1))*bac_m(i-j);
                bac_r(bac_n,1)  = ((bac_m(bac_n)/bac.bac_rho)*(3/(4*pi)))^(1/3);
                % move old parent cell (ng)
                bac_m(i-j) = bac_m(i-j) - bac_m(bac_n) ;        % remaining mass with the parent
                bac_r(i-j) = ((bac_m(i-j)/bac.bac_rho)*(3/(4*pi)))^(1/3);% new radius of the parent
            end
            if bac_m(i-j) < bac.bac_mmin
                act(i-j) = 0;
            end
            if act(i-j) == 0 && bac_m(i-j) > bac.bac_mmin
                act(i-j) = 1;
            end
        end
     end
end

if Fdiv == 1
    fprintf('..');
    %QuadTree Algorithm
    quadtree = shoving.BiomassQuadtree(0, maxxSys, 0, maxySys);
    r = quadtree.pushing2D(length(bac_x), bac_x, bac_y, bac_r, 0.1, bac.s_dist);
    bac_x = r.bac_x;
    bac_y = r.bac_y;
    
    dxx = bac_x - bac.bac_c;
    dyy = bac_y -  bac.bac_c;
    normd = sqrt(dxx.*dxx + dyy.*dyy);
    while sum(normd > (bac.bac_ymax/2) ) > 0
      [~, I] = max(normd);
      bac_n = bac_n - 1;
      bac_ns(bac_s(I)) = bac_ns(bac_s(I)) - 1;
      % delete of coordinates, species and radius 
      bac_m(I)  = [];
      bac_r(I)  = [];
      bac_x(I)  = [];
      bac_y(I)  = [];
      bac_s(I)  = [];
      bac_a(I)  = [];
      act(I)    = [];
%       bac_mu_max(I) = [];
      bac_yield(I) = [];
      bac_Ks(I, :) = [];
%       bac_maint(I) = [];
      dxx = bac_x -  bac.bac_c;
      dyy = bac_y -  bac.bac_c;
      normd = sqrt(dxx.*dxx+dyy.*dyy);
    end

    bac.bac_n = bac_n;
    bac.bac_ns = bac_ns;
    bac.bac_Ks = bac_Ks;
    bac.atrib = [bac_x, bac_y, bac_m/bac.bac_MW, bac_a, bac_s, bac_r, bac_yield, act];
    bac.bac_rho_bio = sum(bac_m)/((max(bac_y)- min(bac_y))*(max(bac_x)- min(bac_x))*(1e-6));
end
end