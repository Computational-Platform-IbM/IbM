function [R, StVLiq, pH, Fb_layer, Bc, Lxy, Bce, Ve, Ve2, CC, DD] = DiffMatrices(Bc, Lxy, Bce, Ve, Ve2, R, StVLiq, pH, CC, DD, Fdiv)

    % persistent Lxy Bc

    Fb_layer = 0;

    if Fdiv == 1
        bac_x = R.bac.atrib(:,1);
        bac_y = R.bac.atrib(:,2);
%         bac_r = R.bac.atrib(:,6);
        dx = R.Sxy.dx;
        dy = R.Sxy.dy;
        T_blayer = R.Sxy.T_blayer;
        sM = 0; %Smoothness level

        [StVLiq, pH, R.Sxy, Fb_layer, Bc, Lxy, Bce, Ve, Ve2, x_pn, y_pn] = b_layer(Bc, Lxy, Bce, Ve, Ve2, StVLiq, pH, R.pOp.pH, R.Sxy, R.St, bac_x, bac_y);
        nT = R.Sxy.nT;
        nx = R.Sxy.nxSys - 2;
        ny = R.Sxy.nySys - 2;
        bac_n = R.bac.bac_n;
        
        %Check position of each bacteria in granule and assigment of [S]
        %Detection where is inside or outside of granule (DD matrix)
        aux_x = x_pn - dx/2; aux_x1 = x_pn + dx/2;
        aux_y = y_pn - dy/2; aux_y1 = y_pn + dy/2;
        DD = sparse(R.Sxy.nT, 1);
        CC = sparse(R.bac.bac_n, R.Sxy.nT);
        
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
%             Mnorm_nT = R.Sxy.maxy;
        end
        My = kron(My_aux,ones(ny,1));
        
        %CC & DD construction
        parfor i = 1:nT
            c1 = ((bac_x) > aux_x(i)).*((bac_x) <= aux_x1(i));
            c2 = ((bac_y) > aux_y(i)).*((bac_y) <= aux_y1(i));
            vv = (c1.*c2); vv(vv==1) = find(vv);
            ci = sparse(1: bac_n, i, vv, bac_n, nT); %sparce(position i, position j, value, rows sparse, columns sparse);
            di = sparse(i, 1, My(i,1) > aux_y(i), nT, 1); % %value = 1 or 0 -> 1: inside biofilm + BDL | 0: outside biofilm + BDL
            DD = DD + di;
            CC = CC + ci;
        end
        
        ind = [R.Sxy.nT*((1:R.St.numStVLiq2)'-1)+1,R.Sxy.nT*(1:R.St.numStVLiq2)'];
        for k = 1:R.St.numStVLiq2
            %Outside granule: StVLiq = Sbc_Dir // %1: inside biofilm + BDL | 0: outside biofilm + BDL
            StVLiq(ind(k,1):ind(k,2))= DD.*StVLiq(ind(k,1):ind(k,2)) + (1-DD)*R.Sxy.Sbc_Dir(k);
        end
    end

end

function [StVLiq, pH, Sxy, Fb_layer, Bc, Lxy, Bce, Ve, Ve2, x_pn, y_pn] = b_layer(Bc, Lxy, Bce, Ve, Ve2, StVLiq, pH, pHini, Sxy, St, ~, ~) %All space is included
% function [StVLiq, pH, Sxy, Fb_layer, Bc, Lxy, x_pn, y_pn] = b_layer(Bc, Lxy, StVLiq, pH, pHini, Sxy, St, ~, bac_y) %Only biofilm space is included

    xSys = 0:Sxy.dx:Sxy.maxxSys;                       
    ySys = 0:Sxy.dy:Sxy.maxySys;

    dx = Sxy.dx; dy = Sxy.dy;
    
    pos_minx = 1; pos_miny = 1;
    maxy = Sxy.maxy + Sxy.T_blayer;
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
    pos_xySys_old = Sxy.pos_xySys;
    ax = zeros(Sxy.nxSys-2,1);
    ax(pos_minx:pos_minx+nx-3) = 1;
    ay = zeros(Sxy.nySys-2,1);
    ay(pos_miny:pos_miny+ny-3) = 1;
    Sxy.pos_xySys = sparse(kron(ax, ay));
    xp = x + dx/2;
    x_pn = kron(xp(1:end-1)', ones(ny-2,1));
    yp = y + dy/2;
    y_pn = kron(ones(nx-2,1),yp(1:end-1)');
    Fb_layer = sum(ne(Sxy.pos_xySys, pos_xySys_old)) > 0;

    if Fb_layer == 1
        nT=(nx-2)*(ny-2);
        Sxy.nx = nx; Sxy.ny = ny; Sxy.nT = nT;

        %Substrate concentration in biofilm/granule
        S = kron(Sxy.Sbc_Dir, ones(Sxy.nTSys,1));
        if ne(sum(pos_xySys_old), 0)
            S(kron(ones(St.numStVLiq,1),pos_xySys_old)==1) = StVLiq;  
        else

            S(kron(ones(St.numStVLiq,1),Sxy.pos_xySys)==1) = kron(St.StVIni, ones(nT,1));  
        end
        StVLiq = S(kron(ones(St.numStVLiq,1),Sxy.pos_xySys)==1);

        %pH
        S = pHini*ones(Sxy.nTSys,1);
        if ne(sum(pos_xySys_old), 0)
            S(pos_xySys_old==1) = pH;
        end
        pH = S(Sxy.pos_xySys == 1);

        %Biofilm/granule modifies diffusivity (porosity/tortuosity)
%         kTr.Diffwater = kron(kTr.Diffn, ones(Sxy.nT,1));
%         S = kron(kTr.Diffn, ones(Sxy.nTSys,1));
%         if ne(sum(pos_xySys_old), 0)
%             S(kron(ones(St.numStVLiq2,1),pos_xySys_old)==1) = kTr.Diff;
%         end
%         kTr.Diff = S(kron(ones(St.numStVLiq2,1),Sxy.pos_xySys)==1);

        %Laplacian Matrix
        Ix = speye(single(nx-2));
        Iy = speye(single(ny-2));
        %Laplacian Matrix
        Ex1 = sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
        Ax = Ex1+Ex1'-2*Ix;
        Ey1 = sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
        Ay = Ey1+Ey1'-2*Iy;
        %Dirichlet boundary conditions
        %%Wall
%          Ax(1,1) = -1;          %Left
%          Ax(nx-2,nx-2) = -1;    %Right
%          Ay(1,1) = -1;          %Bottom
%          Ay(ny-2,ny-2) = -1;    %Top
        %%Cyclic boundary conditions in the sides
        Ax(1,nx-2) = 1;  Ax(nx-2,1) = 1;
        
        Lxy = kron(Ix,Ay/dy^2) + kron(Ax/dx^2,Iy);

        %Boundary conditions Matrix
        bcx = zeros(nx-2,nx-2);
        bcy = zeros(ny-2,ny-2);
        %%Dirichlet boundary condition
%         bcx(1,1) = 1;         %Left
%         bcx(nx-2,nx-2) = 1;   %Right
%         bcy(1,1) = 1;         %Bottom
        bcy(ny-2,ny-2) = 1;   %Top
        
        Bc = kron(Ix,bcy/dy^2) + kron(bcx/dx^2,Iy);
        %%Exit boundary condition
        bcx = zeros(nx-2,nx-2);
        bcy = zeros(ny-2,ny-2);
%         bcx(1,1) = 1; ys = ones(ny-2,1); xs = zeros(nx-2,1); = xs(2) = 1;         %Left
%         bcx(nx-2,nx-2) = 1; ys = ones(ny-2,1); xs = zeros(nx-2,1); = xs(end-1) = 1;	%Right
        bcy(1,1) = 1; ys = zeros(ny-2,1); ys(1) = 1; ys2 = zeros(ny-2,1); ys2(2) = 1; xs = ones(nx-2,1);             %Bottom
%         bcy(ny-2,ny-2) = 1; ys = zeros(ny-2,1); ys(end-1) = 1; xs = ones(nx-2,1);   %Top
        Bce = kron(Ix,bcy/dy^2) + kron(bcx/dx^2,Iy);
        Ve = kron(xs,ys);
        Ve2 = kron(xs,ys2);
    end
end