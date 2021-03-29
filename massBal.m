% function R = massBal( L, G, X, R)
% function [NRVbac, NRVliq, NRVgas, Df] = massBal( L, G, X, R, CC)
function [R, NRVbac, NRVliq, NRVgas, pH, Df] = massBal( L, G, X, R, pH, CC)
St = R.St;
Sxy = R.Sxy;
numStVLiq2 = St.numStVLiq2;
numStVLiq = St.numStVLiq;
StVLiq = [L; G];
% Sxy.StVLiq = [L; G];
% Sxy.StVLiq2 = L; 
% Sxy.StVGas = G; 
% R.Sxy = Sxy;
R.bac.atrib(:,3) = X;

[R.bac.atrib(:,4), NRVbac, NRV, pH, Df] = my_kinetics(R,StVLiq, pH, CC);

% [R.bac.atrib(:,4), NRVbac, NRV, Df] = my_kinetics(R, CC);
% [R.bac.bac_a, NRVbac, NRV, R.Sxy.pH, R.Sxy.DGcat, R.Sxy.DGan, Df] = my_kinetics(R);
% rm = R.rm;

% Diff = kron(R.kTr.Diffn, ones(R.Sxy.nT,1)).*kron(Df, ones(St.numStVLiq2,1));%%%OJO
% R.kTr.a = repmat(Sxy.dT*Diff/2, 1, Sxy.nT*numStVLiq2);
% R.kTr.Diff = Diff;

if R.flagGas == 1
    %%%%% Gas correction
    ind = [Sxy.nT*((numStVLiq2+1:numStVLiq)'-1)+1,Sxy.nT*(numStVLiq2+1:numStVLiq)'];
    StVGas = StVLiq(numStVLiq2*Sxy.nT+1:end);

%     StVGas = Sxy.StVLiq(numStVLiq2*Sxy.nT+1:end);
    for k = 1:(numStVLiq - numStVLiq2)
        if R.flagGas == 2
            nrv_gas = sum(NRV(ind(k,1):ind(k,2)))/Sxy.nT - (R.Qgas/R.Sxy.Vgas).*StVGas(1+(k-1)*Sxy.nT:k*Sxy.nT);
            NRV(ind(k,1):ind(k,2)) = nrv_gas;
        else
            NRV(ind(k,1):ind(k,2)) = 0;
        end
    end
end
aux = numStVLiq2*Sxy.nT;
NRVgas = NRV(aux+1:end);
NRVliq = NRV(1:aux);
% rm.NRVgas = NRV(aux+1:end);
% rm.NRVliq = NRV(1:aux);
% rm.NRVbac = NRVbac;
% R.rm = rm;
end