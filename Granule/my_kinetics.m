function [Mu, NRVbac, NRV, pH, Df] = my_kinetics(R,StVLiq, pH, CC)
Keq = R.pTh.Keq;
chrM = R.pTh.chrM;
MatrixMet = R.rm.MatrixMet;
MatrixDecay = R.rm.MatrixDecay_mod;
Vg = R.Sxy.Vg;
nT = R.Sxy.nT;
T = R.pOp.T;
Sh_ini_all = 10.^(-pH);
SLiqxy_all = (StVLiq > 0).*StVLiq;

numStVLiq = R.St.numStVLiq;
numStVLiq2 = R.St.numStVLiq2;

NH3_pos = strcmp(R.St.StNames(1:numStVLiq), 'NH3');
NO2_pos = strcmp(R.St.StNames(1:numStVLiq), 'NO2');
O2_pos = strcmp(R.St.StNames(1:numStVLiq), 'O2');
rv = R.pTh.react_v;

bac_s = R.bac.atrib(:,5);
bac_m = R.bac.atrib(:,3);
% bac_mu_max = R.bac.atrib(:,7);
% bac_yield = R.bac.atrib(:,8);
bac_Ks = R.bac.bac_Ks;
% bac_maint = R.bac.atrib(:,9);
bac_n = R.bac.bac_n;
act = R.bac.atrib(:,8);

Mu = zeros(bac_n, 1); NRVbac = zeros(bac_n, 1); NRV = zeros(numStVLiq*nT,1); Ngas = NRV; X = zeros(nT,1); pH = X;
NRV_aux = zeros(numStVLiq, nT);
SLiqxy_aux=zeros(numStVLiq,nT);
Sh_ini_all_aux=zeros(nT);
for ii = 1:nT
    indLiq = ii+((1:numStVLiq)-1)*nT;
    SLiqxy_aux(:,ii) = SLiqxy_all(indLiq);
    Sh_ini_all_aux(ii) = Sh_ini_all(indLiq(1));
end
NRVbac_aux = cell(nT);
for i=1:nT
    SLiqxy = SLiqxy_aux(:,i);
    [spcM, Sh] = solve_pH(Sh_ini_all_aux(i), [SLiqxy;1;0], Keq, chrM, 0);
    pH(i) = -log10(Sh);
    c = CC(:,i);
    num_b = nnz(c);
    if  num_b > 0
        c =(c(ne(c,0)));
        NRVbac_aux{i} = zeros(num_b,3);
        for e = 1:1:num_b
            k = bac_s(c(e));
            %- mu_max & maint -%
            if k == 1 %AOB
                mu_max = ((1.28*10^(12)*exp(-8183/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
                maint = (1.651*10^(11)*exp(-8183/T))/24;
            elseif k == 2 %Nitrob
                mu_max = ((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
                maint = (8.626*10^(6)*exp(-5295/T))/24;
            elseif k == 3 %Nitros
                mu_max = 0.63*((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
                maint = 0.63*(8.626*10^(6)*exp(-5295/T))/24;
            elseif k == 4 %AMX
                mu_max = 1.89*10^(8)*exp(-7330/T); %From [Puyol2014], Brocadia spp.
                maint = 0.05*mu_max;
            end
            %- Monod blocks -%
            M = 1;
            if ne(bac_Ks(c(e),1),0)
%                 M = M*(SLiqxy(NH3_pos)/(SLiqxy(NH3_pos) + bac_Ks(c(e),1)));
                M = M*(spcM(NH3_pos,rv(NH3_pos))/(spcM(NH3_pos,rv(NH3_pos)) + bac_Ks(c(e),1)));
            end
            if ne(bac_Ks(c(e),2),0)
%                 M = M*(SLiqxy(NO2_pos)/(SLiqxy(NO2_pos) + bac_Ks(c(e),2)));
                M = M*(spcM(NO2_pos,rv(NO2_pos))/(spcM(NO2_pos,rv(NO2_pos)) + bac_Ks(c(e),2)));
            end
            if ne(bac_Ks(c(e),3),0)
%                 M = M*(SLiqxy(O2_pos)/(SLiqxy(O2_pos) + bac_Ks(c(e),3)));
                M = M*(spcM(O2_pos,rv(O2_pos))/(spcM(O2_pos,rv(O2_pos)) + bac_Ks(c(e),3)));
            end
            if ne(bac_Ks(c(e),4),0)
%                 M = M*(bac_Ks(c(e),4)/(SLiqxy(NH3_pos) + bac_Ks(c(e),4)));
                M = M*(bac_Ks(c(e),4)/(spcM(NH3_pos,rv(NH3_pos)) + bac_Ks(c(e),4)));
            end
            if ne(bac_Ks(c(e),5),0)
%                 M = M*(bac_Ks(c(e),5)/(SLiqxy(NO2_pos) + bac_Ks(c(e),5)));
%                 M = M*(bac_Ks(c(e),5)/(spcM(NO2_pos,rv(NO2_pos)) + bac_Ks(c(e),5)));
                M = M*( (4.43e-3*(1 + SLiqxy(NH3_pos)/5.45e-3)) / (SLiqxy(NO2_pos) + 4.43e-3) );
            end
            if ne(bac_Ks(c(e),6),0)
%                M = M*(bac_Ks(c(e),6)/(SLiqxy(NO2_pos) + bac_Ks(c(e),6)));
                M = M*(bac_Ks(c(e),6)/(spcM(O2_pos,rv(O2_pos)) + bac_Ks(c(e),6)));
            end
            
            mu_aux = mu_max*M;
            mu = mu_aux - maint;
            rg = mu*bac_m(c(e));
            sol_c = MatrixMet(:,k)*mu_aux*bac_m(c(e));
            if mu >= 0
                NRVxy = sol_c;
            else
                if act(c(e)) == 1
                    NRVxy = MatrixDecay(:,k)*(-rg) + sol_c;
                else
                    NRVxy = zeros(size(NRVxy));      %If bac is inactive: no substrate consumption, no change of mass
                    rg = zeros(size(rg));
                end
            end
            NRVbac_aux{i}(e,:) = [mu, rg, c(e)];
            
%             NRVbac(c(e)) = rg;
%             Mu(c(e)) = mu;
            NRVxy(1:numStVLiq2) = NRVxy(1:numStVLiq2)/Vg;
            NRV_aux(:,i) = NRV_aux(:,i) + NRVxy;
%             NRV(indLiq) = NRV(indLiq)+NRVxy;
        end
    end
end
for ii=1:nT
    indLiq = ii+((1:numStVLiq)-1)*nT;
    NRV(indLiq) = NRV(indLiq)+NRV_aux(:,ii);
    ii_aux = NRVbac_aux{ii};
    if ne(size(ii_aux),0)
        for j =1:size(ii_aux,1)
            NRVbac(ii_aux(j,3)) = ii_aux(j,2);
            Mu(ii_aux(j,3)) = ii_aux(j,1);
        end
    end
end

NRV = NRV + Ngas;
if R.flagGas == 2
    nrv_gas = sum((Ngas(1:numStVLiq2*nT)< 0).*Ngas(1:numStVLiq2*nT));
    R.Qgas = nrv_gas*R.pOp.Rg*R.pOp.T/R.pOp.P;
end
%  Df = 1 - (0.43.*X.^0.92)./(11.19 + 0.27.*X.^0.99);
 Df = 1 ;
end

function [spcM, Sh] = solve_pH(Sh_ini, StV, Keq, chrM, flagpH)
w = 1;

if flagpH == 1
    spcM = zeros(size(chrM));
    Sh = Sh_ini;
    Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);

    spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
    spcM(:,2) = (StV * Sh^3)                                    ./Denm;
    spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
    spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
    spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;

else
    % Checking the existence of a zero pool in the function between pH 1 and 14
    a=1e-14;
    b=1;
    Sh_v = [a; b]; F = zeros(2,1);
    for nn = 1:2
        Sh=Sh_v(nn);
        spcM = zeros(size(chrM));
        Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);

        spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
        spcM(:,2) = (StV * Sh^3)                                    ./Denm;
        spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
        spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
        spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;

        % Evaluation of the charge balance for the current Sh value, F(Sh)
        F(nn) = Sh + sum(sum(spcM.*chrM));
    end
    FF = prod(F);
    if FF > 0 || isnan(FF)
        error('ERROR.- The sum of charges returns a wrong value')
    end

    fa = F(1);
    fb = F(2);

    % Newton-Raphson method.-
    Sh = Sh_ini;
    % Counter of convergences
    ipH=1;    Tol = 5.e-15;    maxIter = 20;
    % Inicialization of matrix of species
    spcM = zeros(size(chrM)); dspcM = spcM;
    while ipH <= maxIter
        Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2); 
        spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
        spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
        spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
        spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
        spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;

        % Evaluation of the charge balance for the current Sh value, F(Sh)
        F = Sh + sum(sum(spcM.*chrM));

        % Calculation of all derivated functions
        dDenm = Denm.^2;
        aux = 3*Sh^2*(Keq(:,1)/w + 1) + 2*Sh*Keq(:,2) + Keq(:,2).*Keq(:,3);

        dspcM(:,1) =  (3*Sh^2*Keq(:,1).*StV)./(w*Denm) - ((Keq(:,1).*StV*Sh^3).*aux) ./(w*dDenm);
        dspcM(:,2) =  (3*Sh^2*StV)./Denm - (StV*Sh^3.*aux) ./dDenm;
        dspcM(:,3) = (2*Sh*Keq(:,2).*StV)./Denm - ((Keq(:,2).*StV*Sh^2).*aux) ./dDenm;
        dspcM(:,4) = (Keq(:,2).*Keq(:,3).*StV)./Denm - ((Keq(:,2).*Keq(:,3).*StV*Sh).*aux)./dDenm;
        dspcM(:,5) = -(Keq(:,2).*Keq(:,3).*Keq(:,4).*StV.*aux) ./dDenm;

        % Evaluation of the charge balance for the current Sh value, dF(Sh)
        dF = 1 + sum(sum(dspcM.*chrM));
        %Error
        err = F/dF;
        % Newton-Raphson algorithm
        Sh = Sh - err;

        if (abs(err) < 1e-14) && (abs(F) < Tol)
          % Checking if a valid pH was obtained
            if (Sh > 1e-14) && (Sh < 1)
                ipH = maxIter;
            else
                % Counter of convergence
                ipH = 1; maxIter = 50;
                n1 = 0; n2 = 0;
                while (ipH < maxIter)
                    Sh = (fb*a-fa*b)/(fb-fa);
                    spcM = zeros(size(chrM));
                    Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);

                    spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
                    spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
                    spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
                    spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
                    spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;

                    fc = Sh + sum(sum(spcM.*chrM));
                    if fa*fc > 0
                        n1 = n1+1;
                        if n1 == 2
                            fb = (fc/(fc+fa))*fb;
                            n1 = 0;
                        end
                        a = Sh; fa = fc;
                    elseif fb*fc > 0 % To avoid problems when fc == 0
                        n2 = n2+1;
                        if n2 == 2
                            fa = (fc/(fc+fb))*fa;
                            n2 = 0;
                        end
                        b = Sh; fb = fc;
                    end
                    err1 = abs(fc);
                    err2 = abs(Sh-(fb*a-fa*b)/(fb-fa));
                    if (err1 < Tol) && (err2 < 1e-14)
                        ipH = maxIter;
                    end    
                    ipH = ipH+1;  
                end
            end
        end
        ipH = ipH+1;
    end        
end
end