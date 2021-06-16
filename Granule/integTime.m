function R = integTime(R)
Tol_a = R.kTr.Tolabs;
Tol_r = (R.kTr.Tolrel(1:4))'; %percentage (only NH3, NO2, NO3, O2)
dT = R.Sxy.dT;
dT_bac = R.Sxy.dT_bac; 
dT_Print = R.Sxy.dT_Print;
dT_Div = R.Sxy.dT_Div;
nTi =ceil (R.Sxy.maxT/dT + 1);
numStVLiq2 = R.St.numStVLiq2;  

R = boundary(R, dT_bac, 0, 0);
% R.St.StVIni(strcmp(R.St.StNames, 'Na')) = R.Sxy.Sbc_Dir(strcmp(R.St.StNames, 'Na'));
% R.St.StVIni(strcmp(R.St.StNames, 'CO2')) = R.Sxy.Sbc_Dir(strcmp(R.St.StNames, 'CO2'));

[R, StVLiq, pH, ~, BBc, Lxy, CC, DD] = DiffMatrices(0, 0, R, 0, 0, 0, 0, 1);
[R, ~, NRVliq, ~, pH, Df] = massBal(StVLiq(1:numStVLiq2*R.Sxy.nT), StVLiq(numStVLiq2*R.Sxy.nT+1:end), R.bac.atrib(:,3), R, pH, CC);
LL = sparse(R.Sxy.nT, R.Sxy.nT);
for p = 1: R.Sxy.nT
    LLp = sparse(1: R.Sxy.nT, p, Lxy(:,p).*Df, R.Sxy.nT, R.Sxy.nT);
    LL = LL + LLp;
end

Sbc_Dir = 1000*kron(R.Sxy.Sbc_Dir(1:numStVLiq2), ones(R.Sxy.nT,1));
U = StVLiq(1:numStVLiq2*R.Sxy.nT)*1000;

ind = [R.Sxy.nT*((1:numStVLiq2)'-1)+1,R.Sxy.nT*(1:numStVLiq2)'];
for ks = 1:numStVLiq2
    U(ind(ks,1):ind(ks,2))= DD.*U(ind(ks,1):ind(ks,2)) + (1-DD).*Sbc_Dir(ind(ks,1):ind(ks,2));
end
V = StVLiq(numStVLiq2*R.Sxy.nT+1:end);
W = R.bac.atrib(:,3); 
U0 = U;
RES = zeros(R.Sxy.nT,numStVLiq2);
mm = 1; nn = 2; 
j=2; jj=2; i = 2;  l = i + mm; ll = 2;
TimeBac = dT_bac*(j-1); TimePrint = dT_Print*(ll-1);TimeDiv = dT_Div*(jj-1); Time = dT*(i-1); TimeSS = dT*(l-1); %TimePrint2 = dT_Print*(ll-1); %TimePrev = Time;
numSt = numStVLiq2*R.Sxy.nT;

D = speye(R.Sxy.nT);
Sxy = zeros((R.Sxy.nx-2), numStVLiq2);
aux = floor((R.Sxy.ny-2)/2)*R.Sxy.ny;
for ii = 1:numStVLiq2
    S = StVLiq(ind(ii,1):ind(ii,2)); 
    Sxy(:,ii)  = S((aux+1):((R.Sxy.nx-2)+aux));
end
pHxy = pH((aux+1):((R.Sxy.nx-2)+aux));
                                       
out_integTime(0, R, Sxy, 0);

Fn = 2;F = ones(1, Fn); NNN = 0;
SS0 = 1e100*(2:-1:1);
Diff = R.kTr.Diffn*R.Sxy.dT/2;
T0 = R.Inf.St; T0(1:3,1)= R.pOp.NH3sp; T0(6,1)= R.pOp.NH3sp/2;
U_aux = zeros(R.Sxy.nT,numStVLiq2); NRVliq_aux= zeros(R.Sxy.nT,numStVLiq2); Sbc_Dir_aux= zeros(R.Sxy.nT,numStVLiq2);
dx = R.Sxy.dx;Diffn = R.kTr.Diffn;

while i < nTi
    for kn = 1:numStVLiq2
        U_aux(:,kn) = U(ind(kn,1):ind(kn,2));
        NRVliq_aux(:,kn) = NRVliq(ind(kn,1):ind(kn,2));
        Sbc_Dir_aux(:,kn) = Sbc_Dir(ind(kn,1):ind(kn,2));
    end
       
    parfor k = 1:numStVLiq2
        B1 = (D + Diff(k,1)*LL)*U_aux(:,k);
        B2 = (2*Diff(k,1)*BBc)*Sbc_Dir_aux(:,k); 
        B3 = dT*(1000)*NRVliq_aux(:,k);
        B = B1 + B2 + B3;
        U_aux(:,k) = (D - Diff(k,1)*LL)\B;
        U_aux(:,k)= DD.*U_aux(:,k) + (1-DD).*Sbc_Dir_aux(:,k);
        RES(:, k) = ((abs((dx^2)*LL*U_aux(:,k)+ (dx^2)*BBc*Sbc_Dir_aux(:,k)+((dx^2)/Diffn(k,1))*(1000)*NRVliq_aux(:,k))/1000)./(1e-4 + U_aux(:,k)/1000));
    end
    for kn = 1:numStVLiq2
        U(ind(kn,1):ind(kn,2)) = U_aux(:,kn);
    end
    U = (U > 0).*U + (U <= 0).*Tol_a;
    StVLiq = [U/1000; V];

    [R.bac.atrib(:,4), NRVbac, NRV, pH, ~] = my_kinetics(R, StVLiq, pH, CC);
    NRVliq = NRV(1:numSt);
    NRVgas = NRV(numSt+1:end);
    V = V + dT*NRVgas;
    W = W + dT*NRVbac; R.bac.atrib(:,3) = W; 

    if Time >= TimeSS   
        VV = abs(U - U0);
        VVr = max(reshape(VV, R.Sxy.nT,numStVLiq2));
        SS = max(reshape((VV.*(U > Tol_a))./U0 ,R.Sxy.nT, numStVLiq2))*100;
        SS = SS(1:4); %percentage (only NH3, NO2, NO3, O2)
        TolRES = mean(reshape(RES, R.Sxy.nT, numStVLiq2)); %mean() or max()
        TolRES = TolRES(1:4)*100; %perct (only NH3, NO2, NO3, O2)
        
        if max(TolRES) <= 1 %max(SS - Tol_r) < 0
            F = ones(1, Fn); SS0 = 1e100*(2:-1:1); NNN = 0;
            dT = R.Sxy.dT;
            dT_bac = TimeBac - Time;
            i = i + ceil((TimeBac-Time)/dT); Time = TimeBac;
            if i > nTi
                i = nTi; Time = dT*(nTi-1);
            end
        else
            SS0(1:end-1) = SS0(2:end);
            SS0(end) = max(SS);
            F(1:end-1) = F(2:end);
            F(end) = sum(SS0 == sort(SS0,'descend'))/2;
        end
        U0 = U;
        if Time >= TimeBac 
            [R, NRVbac, NRVliq, NRVgas, pH, ~] = massBal( U/1000, V, W, R, pH, CC);
            V = V + dT_bac*NRVgas; 
            W = W + dT_bac*NRVbac;
            
            fprintf('>>>> Current simulation time is %.1f hours.',Time)
            if Time >= TimeDiv
                [R, StVLiq, Fdiv] = bacteria( U/1000, V, W, R);
                fprintf('Number of bacteria: %.0f', R.bac.bac_n)
                save('R.mat', 'R', '-v7.3');
                if Fdiv == 1
                    if R.bac.bac_n > R.bac.bac_nmax
                        error('\n\n Maximum number of cells allowed exceeded !!\n\n Number of cells calculated: %d\n Maximum allowed %d\n\n\n >>> END of SIMULATION\n\n\n', R.bac.bac_n,R.bac.bac_nmax)
                    end
                    W = R.bac.atrib(:,3);
                end
                [R, StVLiq, pH, Fb_layer, BBc, Lxy, CC, DD] = DiffMatrices(BBc, Lxy, R, StVLiq, pH, CC, DD, Fdiv);

                if Fb_layer == 1
                    [R, ~, NRVliq, ~, pH, Df] = massBal(StVLiq(1:numStVLiq2*R.Sxy.nT), StVLiq(numStVLiq2*R.Sxy.nT+1:end), R.bac.atrib(:,3), R, pH, CC);
                    U = StVLiq(1:numStVLiq2*R.Sxy.nT)*1000;
                    U0 = U;
                    V = StVLiq(numStVLiq2*R.Sxy.nT+1:end);
                    numSt = numStVLiq2*R.Sxy.nT;
                    ind = [R.Sxy.nT*((1:numStVLiq2)'-1)+1,R.Sxy.nT*(1:numStVLiq2)'];
                    D = speye(R.Sxy.nT);
                    LL = sparse(R.Sxy.nT, R.Sxy.nT);
                    for p = 1: R.Sxy.nT
                        LLp = sparse(1: R.Sxy.nT, p, Lxy(:,p).*Df, R.Sxy.nT, R.Sxy.nT);
                        LL = LL + LLp;
                    end
                end
                
                if Time >= TimePrint
                    fprintf('.Imprimiendo resultados...')
                    Sxy = zeros((R.Sxy.nx-2), numStVLiq2);
                    aux = floor((R.Sxy.ny-2)/2)*R.Sxy.ny;
                    for ii = 1:numStVLiq2
                        S = StVLiq(ind(ii,1):ind(ii,2)); 
                        Sxy(:,ii)  = S((aux+1):((R.Sxy.nx-2)+aux));
                    end
                    pHxy = pH((aux+1):((R.Sxy.nx-2)+aux));
                                       
                    out_integTime(Time, R, Sxy, pHxy, i == nTi);
                    ll = ll+1; TimePrint = dT_Print*(ll-1);
                    fprintf('.Done!\n')
                end
                jj = jj+1; TimeDiv = dT_Div*(jj-1);
            end
            fprintf('\n')
            dT_bac = R.Sxy.dT_bac;
            R = boundary(R, dT_bac, NRVliq, 1);  
            j = j+1; TimeBac = dT_bac*(j-1);
            l = i + mm; TimeSS = dT*(l-1);
 
            Sbc_Dir = 1000*kron(R.Sxy.Sbc_Dir(1:numStVLiq2), ones(R.Sxy.nT,1));
            U_aux = zeros(R.Sxy.nT,numStVLiq2); NRVliq_aux= zeros(R.Sxy.nT,numStVLiq2); Sbc_Dir_aux= zeros(R.Sxy.nT,numStVLiq2);
            RES = zeros(R.Sxy.nT,numStVLiq2);
        end
    end
    i = i + (dT/R.Sxy.dT)*1; Time = R.Sxy.dT*(i-1);
end
end
function  out_integTime(t, R, Sxy, pHxy, fin)
nTSys = R.Sxy.nTSys;
nxSys = R.Sxy.nxSys;

if t == 0
    tic
    All_StatesVar = [];
    B = [];    
    SSxy = [];
    pHs = [];
    if R.flagGas == 2
        All_StatesVar(end+1, :) = [0; R.St.StV; R.Qgas; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3)); R.bac.bac_rho_bio]'; %#ok<*NASGU>
    else
        All_StatesVar(end+1, :) = [0; R.St.StV; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3)); R.bac.bac_rho_bio]'; 
    end
    
    S = zeros(nxSys-2,R.St.numStVLiq2);
    for k = 1:R.St.numStVLiq2
        A = R.Sxy.Sbc_Dir(k)*ones(nxSys-2,1);
        
        aux = round(((nxSys-2)-(R.Sxy.nx-2))/2);
        A(aux+1:aux+(R.Sxy.nx-2)) = Sxy(:,k);
        S(:, k) = A;
    end
    SSxy(:, :, 1) = S;
    ApH = R.pOp.pH*ones(nxSys-2,1);
    ApH(aux+1:aux+(R.Sxy.nx-2)) = pHxy(:,1);
    pHs(:, 1) = ApH;
    
    bac = R.bac;   
    B1 = t(end)*ones(bac.bac_nmax,1); 
    aux = zeros(bac.bac_nmax-bac.bac_n,1);
    B2 = [bac.atrib(:,1);aux];
    B3 = [bac.atrib(:,2);aux];
    B4 = [bac.atrib(:,5);aux];
    B5 = [bac.atrib(:,6);aux];
    B6 = [bac.atrib(:,4);aux];
    B7 = [bac.atrib(:,8);aux];

    B(:, :, end+1) = [B1 B2 B3 B4 B5 B6 B7];
else 
    load ResultsSim.mat; %#ok<LOAD>
    bac = R.bac;    
    B1 = t(end)*ones(bac.bac_nmax,1); 
    aux = zeros(bac.bac_nmax-bac.bac_n,1);
    B2 = [bac.atrib(:,1);aux];
    B3 = [bac.atrib(:,2);aux];
    B4 = [bac.atrib(:,5);aux];
    B5 = [bac.atrib(:,6);aux];
    B6 = [bac.atrib(:,4);aux];
    B7 = [bac.atrib(:,8);aux];

    B(:, :, end+1) = [B1 B2 B3 B4 B5 B6 B7];
    if R.flagGas == 2
        All_StatesVar(end+1, :) = [t(end); R.St.StV; R.Qgas; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3)); R.bac.bac_rho_bio]';
    else
         All_StatesVar(end+1, :) = [t(end); R.St.StV; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3)); R.bac.bac_rho_bio]';
    end
    
    S = zeros(nxSys-2,R.St.numStVLiq2);
    aux = round(((nxSys-2)-(R.Sxy.nx-2))/2);
    for k = 1:R.St.numStVLiq2
        A = R.Sxy.Sbc_Dir(k)*ones(nxSys-2,1);
        A(aux+1:aux+(R.Sxy.nx-2)) = Sxy(:,k);
        S(:, k) = A;
    end
    SSxy(:, :, end+1) = S;
    ApH = R.pOp.pH*ones(nxSys-2,1);
    ApH(aux+1:aux+(R.Sxy.nx-2)) = pHxy(:,1);
    pHs(:, end+1) = ApH;
    
    if fin == 1
        toc
    end
end
save('ResultsSim.mat', 'All_StatesVar', 'B', 'SSxy', 'pHs', '-v7.3');
end