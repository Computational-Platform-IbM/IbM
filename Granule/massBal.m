function [R, NRVbac, NRVliq, NRVgas, pH, Df] = massBal( L, G, X, R, pH, CC)
St = R.St;
Sxy = R.Sxy;
numStVLiq2 = St.numStVLiq2;
numStVLiq = St.numStVLiq;
StVLiq = [L; G];
R.bac.atrib(:,3) = X;
[R.bac.atrib(:,4), NRVbac, NRV, pH, Df] = my_kinetics(R,StVLiq, pH, CC);
if R.flagGas == 1
    %%%%% Gas correction
    ind = [Sxy.nT*((numStVLiq2+1:numStVLiq)'-1)+1,Sxy.nT*(numStVLiq2+1:numStVLiq)'];
    StVGas = StVLiq(numStVLiq2*Sxy.nT+1:end);
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
end