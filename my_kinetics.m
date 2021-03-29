% function [Mu, NRVbac, NRV, pH, Df] = my_kinetics(R,StVLiq, pH, CC)
function [Mu, NRVbac, NRV, pH, Df] = my_kinetics(R,StVLiq, ~, CC)
% Keq = R.pTh.Keq;
% chrM = R.pTh.chrM;
% flagpH = R.flagpH;
% flagDG = R.flagDG;
MatrixMet = R.rm.MatrixMet;
% MatrixCat = R.rm.MatrixCat_mod;
% MatrixAn = R.rm.MatrixAn_mod;
MatrixDecay = R.rm.MatrixDecay_mod;
% rMatrixFull = R.rm.rMatrixFull';
Vg = R.Sxy.Vg;
nT = R.Sxy.nT;
% T = R.pOp.T;
% Sh_ini_all = 10.^(-pH);
SLiqxy_all = (StVLiq > 0).*StVLiq;
%  SLiqxy_all = (R.Sxy.StVLiq > 0).*R.Sxy.StVLiq;
% Sh_ini_all = 10.^(R.pOp.pH)*ones(nT,1);
% RthT = R.pOp.Rth*R.pOp.T;
% VRgT = R.pOp.Vgas/(R.pOp.Rg*R.pOp.T);
% KhV = R.pTh.Kh.KhV;
% kLa = R.kTr.kLa;
% DG0 = R.pTh.DG0;
% DGr0 = R.pTh.DGr0;
% DGdis = R.pTh.DGdis;
% spcR = R.pTh.spcR;
% lCat_eD = R.pTh.lCat_eD;

numStVLiq = R.St.numStVLiq;
numStVLiq2 = R.St.numStVLiq2;
% numStVGas2Liq = R.St.numStVGas - 1;
% numX = R.St.numX;
% numStFull = R.St.numStFull;
% Liq2Gas = R.St.Liq2Gas;

A_pos = strcmp(R.St.StNames(1:numStVLiq), 'A');
%  rv = R.pTh.react_v;
bac_s = R.bac.atrib(:,5);
bac_m = R.bac.atrib(:,3);
bac_mu_max = R.bac.atrib(:,7);
% bac_yield = R.bac.atrib(:,8);
bac_Ks = R.bac.bac_Ks;
bac_maint = R.bac.atrib(:,9);

% bac_s = R.bac.bac_s;
% bac_m = R.bac.bac_m;
% bac_mu_max = R.bac.bac_mu_max;
% bac_yield = R.bac.bac_yield;
% bac_Ks = R.bac.bac_Ks;
% bac_maint = R.bac.bac_maint;
% CC = R.Sxy.c;
Mu = zeros(R.bac.bac_n, 1); NRVbac = zeros(R.bac.bac_n, 1); NRV = zeros(numStVLiq*nT,1); Ngas = NRV; X = zeros(nT,1); pH = X; 
% DGCatM = zeros(nT*numX,1); DGAnM = DGCatM; DG_mat = ones(numStFull,1);

for i=1:nT
    indLiq = i+((1:numStVLiq)-1)*nT;
    SLiqxy = SLiqxy_all(indLiq); 
%     [spcM, Sh] =solve_pH(Sh_ini_all(indLiq(1)), [SLiqxy;1;0], Keq, chrM, 0);
%     pH(i) = -log10(Sh);
    pH(i) = R.pOp.pH;
    spcM = SLiqxy;
    c = CC(:,i);
    num_b = nnz(c);
    if  num_b > 0
        c =(c(ne(c,0)));
        X(i) = sum(bac_m(c))/Vg; 
        for e = 1:num_b
            k = bac_s(c(e));
            M = 1;             
            if ne(bac_Ks(c(e),1),0)
                   M = M*(spcM(A_pos)/(spcM(A_pos) + bac_Ks(c(e),1)));
            end
            mu_aux = bac_mu_max(c(e))*M;
            mu = mu_aux - bac_maint(c(e));
            rg = mu*bac_m(c(e));
            sol_c = MatrixMet(:,k)*mu_aux*bac_m(c(e));
            if mu >= 0
                NRVxy = sol_c;
            else
                NRVxy = MatrixDecay(:,k)*(-rg) + sol_c;
            end
            NRVbac(c(e)) = rg;
            Mu(c(e)) = mu;
            NRVxy(1:numStVLiq2) = NRVxy(1:numStVLiq2)/Vg;
            NRV(indLiq) = NRV(indLiq)+NRVxy;
        end
    end
end

NRV = NRV + Ngas;
%  X = X*R.bac.bac_MW; 
if R.flagGas == 2
    nrv_gas = sum((Ngas(1:numStVLiq2*nT)< 0).*Ngas(1:numStVLiq2*nT));
    R.Qgas = nrv_gas*R.pOp.Rg*R.pOp.T/R.pOp.P;
end
%  Df = 1 - (0.43.*X.^0.92)./(11.19 + 0.27.*X.^0.99);
 Df = 1 ;
end



% function [spcM, Sh] = solve_pH(Sh_ini, StV, Keq, chrM, flagpH)
% w = 1;
% 
% if flagpH == 1
%     spcM = zeros(size(chrM));
%     Sh = Sh_ini;
%     Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
% 
%     spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
%     spcM(:,2) = (StV * Sh^3)                                    ./Denm;
%     spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
%     spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
%     spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
% 
% else
%     % Checking the existence of a zero pool in the function between pH 1 and 14
%     a=1e-14;
%     b=1;
%     Sh_v = [a; b]; F = zeros(2,1);
%     for nn = 1:2
%         Sh=Sh_v(nn);
%         spcM = zeros(size(chrM));
%         Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
% 
%         spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
%         spcM(:,2) = (StV * Sh^3)                                    ./Denm;
%         spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
%         spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
%         spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
% 
%         % Evaluation of the charge balance for the current Sh value, F(Sh)
%         F(nn) = Sh + sum(sum(spcM.*chrM));
%     end
%     FF = prod(F);
%     if FF > 0 || isnan(FF)
%         error('ERROR.- The sum of charges returns a wrong value')
%     end
% 
%     fa = F(1);
%     fb = F(2);
% 
%     % Newton-Raphson method.-
%     Sh = Sh_ini;
%     % Counter of convergences
%     ipH=1;    Tol = 5.e-15;    maxIter = 20;
%     % Inicialization of matrix of species
%     spcM = zeros(size(chrM)); dspcM = spcM;
%     while ipH <= maxIter
%         Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2); 
%         spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
%         spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
%         spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
%         spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
%         spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
% 
%         % Evaluation of the charge balance for the current Sh value, F(Sh)
%         F = Sh + sum(sum(spcM.*chrM));
% 
%         % Calculation of all derivated functions
%         dDenm = Denm.^2;
%         aux = 3*Sh^2*(Keq(:,1)/w + 1) + 2*Sh*Keq(:,2) + Keq(:,2).*Keq(:,3);
% 
%         dspcM(:,1) =  (3*Sh^2*Keq(:,1).*StV)./(w*Denm) - ((Keq(:,1).*StV*Sh^3).*aux) ./(w*dDenm);
%         dspcM(:,2) =  (3*Sh^2*StV)./Denm - (StV*Sh^3.*aux) ./dDenm;
%         dspcM(:,3) = (2*Sh*Keq(:,2).*StV)./Denm - ((Keq(:,2).*StV*Sh^2).*aux) ./dDenm;
%         dspcM(:,4) = (Keq(:,2).*Keq(:,3).*StV)./Denm - ((Keq(:,2).*Keq(:,3).*StV*Sh).*aux)./dDenm;
%         dspcM(:,5) = -(Keq(:,2).*Keq(:,3).*Keq(:,4).*StV.*aux) ./dDenm;
% 
%         % Evaluation of the charge balance for the current Sh value, dF(Sh)
%         dF = 1 + sum(sum(dspcM.*chrM));
%         %Error
%         err = F/dF;
%         % Newton-Raphson algorithm
%         Sh = Sh - err;
% 
%         if (abs(err) < 1e-14) && (abs(F) < Tol)
%           % Checking if a valid pH was obtained
%             if (Sh > 1e-14) && (Sh < 1)
%                 ipH = maxIter;
%             else
%                 % Counter of convergence
%                 ipH = 1; maxIter = 50;
%                 n1 = 0; n2 = 0;
%                 while (ipH < maxIter)
%                     Sh = (fb*a-fa*b)/(fb-fa);
%                     spcM = zeros(size(chrM));
%                     Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
% 
%                     spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
%                     spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
%                     spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
%                     spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
%                     spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
% 
%                     fc = Sh + sum(sum(spcM.*chrM));
%                     if fa*fc > 0
%                         n1 = n1+1;
%                         if n1 == 2
%                             fb = (fc/(fc+fa))*fb;
%                             n1 = 0;
%                         end
%                         a = Sh; fa = fc;
%                     elseif fb*fc > 0 % To avoid problems when fc == 0
%                         n2 = n2+1;
%                         if n2 == 2
%                             fa = (fc/(fc+fb))*fa;
%                             n2 = 0;
%                         end
%                         b = Sh; fb = fc;
%                     end
%                     err1 = abs(fc);
%                     err2 = abs(Sh-(fb*a-fa*b)/(fb-fa));
%                     if (err1 < Tol) && (err2 < 1e-14)
%                         ipH = maxIter;
%                     end    
%                     ipH = ipH+1;  
%                 end
%             end
%         end
%         ipH = ipH+1;
%     end        
% end
% end