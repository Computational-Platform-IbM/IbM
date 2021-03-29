% function [R, NRVliq, Fb_layer, BB1, BB2, LL, UU, CC, DD] = DiffMatrices(R, CC, DD, NRVliq, Fdiv)
% function [R, NRVliq, Fb_layer, a, BBc, LL, CC, DD] = DiffMatrices(R, CC, DD, NRVliq, Df, Fdiv)
% function [R, NRVliq, Fb_layer, BBc, LL, CC, DD] = DiffMatrices(R, CC, DD, NRVliq, Df, Fdiv)
function [R, StVLiq, pH, Fb_layer, Bc, Lxy, CC, DD] = DiffMatrices(Bc, Lxy, R, StVLiq, pH, CC, DD, Fdiv)
    % persistent Lxy Bc

    Fb_layer = 0;

    if Fdiv == 1
        center = R.bac.bac_c;
        bac_x = R.bac.atrib(:,1);
        bac_y = R.bac.atrib(:,2);
        bac_x_c = bac_x - center ;
        bac_y_c = bac_y - center;
        norm_xy = sqrt(bac_x_c.*bac_x_c+bac_y_c.*bac_y_c);

        theta_c = atan2(bac_y_c, bac_x_c);
    %     theta_c(theta_c<0) = +theta_c(theta_c<0) + 2*pi;
        d_theta = (R.Sxy.dx/max(norm_xy));
        v_theta = -pi:d_theta:pi; %v_theta = min(theta_c):d_theta:2*pi
        Mnorm = zeros(length(v_theta)-1,1); Mv_theta = zeros(length(v_theta)-1,1);
        for i = 1:length(v_theta)-1
            ind_theta = (theta_c > v_theta(i)).*(theta_c <= v_theta(i+1)) == 1;
            [M, maxI] = max(norm_xy.*ind_theta);
            Mnorm(i) = M;
            Mv_theta(i) = theta_c(maxI);
        end

        [StVLiq, pH, R.Sxy, Fb_layer, Bc, Lxy, x_pn, y_pn] = b_layer(Bc, Lxy, StVLiq, pH, R.pOp.pH, R.Sxy, R.St, bac_x, bac_y);
        aux_x = x_pn-R.Sxy.dx/2; aux_x1 = x_pn + R.Sxy.dx/2;
        aux_y = y_pn-R.Sxy.dy/2; aux_y1 = y_pn + R.Sxy.dy/2;
        DD = sparse(R.Sxy.nT, 1);
        CC = sparse(R.bac.bac_n, R.Sxy.nT);
        for i = 1:R.Sxy.nT 
            c1 = (bac_x >= aux_x(i)).*(bac_x <= aux_x1(i));
            c2 = (bac_y >= aux_y(i)).*(bac_y <= aux_y1(i));
            theta_nT = atan2(y_pn(i)- center, x_pn(i)- center);
            [~, minI]= min(abs(theta_nT-Mv_theta));
            Mnorm_nT = Mnorm(minI) + R.Sxy.T_blayer;

            vv = (c1.*c2); vv(vv==1) = find(vv);
            ci = sparse(1: R.bac.bac_n, i, vv, R.bac.bac_n, R.Sxy.nT);
            bx = aux_x1(i)-center; by = aux_y1(i)-center;
            di = sparse(i, 1, Mnorm_nT > sqrt(bx.*bx+by.*by), R.Sxy.nT, 1);
            DD = DD + di;
            CC = CC + ci;
        end
        ind = [R.Sxy.nT*((1:R.St.numStVLiq2)'-1)+1,R.Sxy.nT*(1:R.St.numStVLiq2)'];
        for k = 1:R.St.numStVLiq2
            StVLiq(ind(k,1):ind(k,2)) = DD.*StVLiq(ind(k,1):ind(k,2)) + (1-DD)*R.Sxy.Sbc_Dir(k);
        end
    end
end

function [StVLiq, pH, Sxy, Fb_layer, Bc, Lxy, x_pn, y_pn] = b_layer(Bc, Lxy, StVLiq, pH, pHini, Sxy, St, bac_x, bac_y)
% function [ Sxy, kTr, Fb_layer, Bc, Lxy, x_pn, y_pn] = b_layer(Bc, Lxy, Sxy, kTr, St, bac_x, bac_y)

x_biof = [min(bac_x), max(bac_x)];
y_biof = [min(bac_y), max(bac_y)];

xSys = 0:Sxy.dx:Sxy.maxxSys;                       
ySys = 0:Sxy.dy:Sxy.maxySys;

dx = Sxy.dx; dy = Sxy.dy;
%  pos_maxy = length(Sxy.ySys);
maxy =  y_biof(2) + dx + Sxy.T_blayer;
pos_maxy = find((ySys >= maxy),1);
if isempty(pos_maxy)
    pos_maxy = length(ySys);
end

maxx =  x_biof(2) + dx + Sxy.T_blayer;

%  pos_maxx = length(Sxy.xSys);
pos_maxx = find((xSys >= maxx),1);
if isempty(pos_maxx)
    pos_maxx = length(xSys);
end

% pos_minx = 1;
% pos_miny = 1;

minx =  x_biof(1) - dx - Sxy.T_blayer;
pos_minx = find((xSys <= minx),1, 'last');
if isempty(pos_minx)
    pos_minx = 1;
end

miny =  y_biof(1) - dx - Sxy.T_blayer;
pos_miny = find((ySys <= miny),1, 'last');
if isempty(pos_miny)
    pos_miny = 1;
end

x = xSys(pos_minx:pos_maxx);                
y = ySys(pos_miny:pos_maxy);
nx = length(x)-1+2;
ny = length(y)-1+2;

pos_xySys_old = Sxy.pos_xySys;
ax = zeros(Sxy.nxSys-2,1);
ax(pos_minx:pos_minx+nx-3) = 1;
ay = zeros(Sxy.nxSys-2,1);
ay(pos_miny:pos_miny+ny-3) = 1;
Sxy.pos_xySys = sparse(kron(ax, ay));
xp = x+ dx/2;
x_pn = kron(xp(1:end-1)', ones(ny-2,1));
yp = y+ dx/2;
y_pn = kron(ones(nx-2,1),yp(1:end-1)');
Fb_layer = sum(ne(Sxy.pos_xySys, pos_xySys_old)) > 0;

if Fb_layer == 1
    nT=(nx-2)*(ny-2);
    Sxy.nx = nx; Sxy.ny = ny; Sxy.nT = nT; %Sxy.x = x; Sxy.y = y;

%     Sxy.Sbc_DirFull = kron(Sxy.Sbc_Dir, ones(nT,1));

%     S = kron(St.StVIni, ones(Sxy.nTSys,1));

    S = kron(Sxy.Sbc_Dir, ones(Sxy.nTSys,1));
    if ne(sum(pos_xySys_old), 0)
        S(kron(ones(St.numStVLiq,1),pos_xySys_old)==1) = StVLiq;  
    else
         S(kron(ones(St.numStVLiq,1),pos_xySys_old)==1) = StVLiq;
%         S(kron(ones(St.numStVLiq,1),Sxy.pos_xySys)==1) = kron(St.StVIni, ones(nT,1));  
    end
      StVLiq = S(kron(ones(St.numStVLiq,1),Sxy.pos_xySys)==1);

%     Sxy.StVLiq = S(kron(ones(St.numStVLiq,1),Sxy.pos_xySys)==1);
%     Sxy.StVGas = Sxy.StVLiq(St.numStVLiq2*Sxy.nT+1:end); 
%     Sxy.StVLiq2 = Sxy.StVLiq(1:St.numStVLiq2*Sxy.nT); 
    
    S = pHini*ones(Sxy.nTSys,1);
    if ne(sum(pos_xySys_old), 0)
        S(pos_xySys_old==1) = pH;
    end
    pH = S(Sxy.pos_xySys==1);
    

%     kTr.Diffwater = kron(kTr.Diffn, ones(Sxy.nT,1));

%     S = kron(kTr.Diffn, ones(Sxy.nTSys,1));
%     if ne(sum(pos_xySys_old), 0)
%         S(kron(ones(St.numStVLiq2,1),pos_xySys_old)==1) = kTr.Diff;
%     end
%     kTr.Diff = S(kron(ones(St.numStVLiq2,1),Sxy.pos_xySys)==1);
    
    Ix = speye(single(nx-2));
    Iy = speye(single(ny-2));
    %Laplacian Matrix
    Ex1 = sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
    Ax = Ex1+Ex1'-2*Ix;
    Ey1 = sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
    Ay = Ey1+Ey1'-2*Iy;
    %Dirichlet boundary conditions
    %%Wall
%      Ay(1,1)=-1;
%      Ax(1,1)=-1;
%      Ax(nx-2,nx-2)=-1;
    %%Cyclic boundary conditions in the sides
%       Ax(1,nx-2)=1;  Ax(nx-2,1)=1;
    Lxy = kron( Ix, Ay/dy^2 ) + kron( Ax/dx^2, Iy);

    %Boundary conditions Matrix
    bcy = zeros(ny-2, ny-2);
    bcx = zeros((nx-2),(nx-2));
    %%Dirichlet boundary condition
     bcx(nx-2,nx-2)=1;
     bcy(ny-2,ny-2) = 1;
     bcx(1,1)=1;
     bcy(1,1) = 1;
    Bc = kron(Ix, bcy/dy^2 ) + kron( bcx/dx^2, Iy);
end
end