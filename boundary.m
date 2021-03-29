% function R = boundary(R, dT, Flag)
function R = boundary(R, dT, NRVliq, Flag)
Sxy = R.Sxy;
Dir_k = R.Inf.Dir_k;

% numStVLiq = R.St.numStVLiq;
numStVLiq2 = R.St.numStVLiq2;
% numStVLiq3 = R.St.numStVLiq3;

invHRT = R.pOp.invHRT;

if ne(Flag,0)
ind = [Sxy.nT*((1:numStVLiq2)'-1)+1,Sxy.nT*(1:numStVLiq2)'];    
G = zeros(numStVLiq2, 1);
Vr = R.pOp.Vr;
for k = 1:numStVLiq2
    G(k) = sum(NRVliq(ind(k,1):ind(k,2))*Sxy.Vg)/Vr;
end
options = odeset('RelTol',1e-8,'AbsTol',1e-20,'NonNegative',ones(numStVLiq2,1));

[~,Y] = ode45(@(t,y) massbal(t,y,G, R.Inf.St, R.pOp.Ssp, R.flagN),[0 dT],Sxy.Sbc_Dir(1:numStVLiq2), options);
Reactor = Y(end, :)';

if ne(sum(Reactor < 0), 0)
    ap = find(Reactor <= 0);
    for k = 1:length(ap)
        if ap(k) == 1
            Reactor(2) = Reactor(2) + Reactor(ap(k));
            Reactor(1) = 0;
        elseif ap(k) == 2
            Reactor(3) = Reactor(3) + Reactor(ap(k));
            Reactor(2) = 0;
        end
    end
end

end
aux = zeros(numStVLiq2, 1);
for k = 1:numStVLiq2
    if Dir_k(k) %Dirichlet  
        aux(k) = Sxy.Sbc_Dir(k);        
    else %Neumann  
        if ne(Flag,0)
            aux(k) = Reactor(k);
        else
            aux(k) = Sxy.Sbc_Dir(k);   
        end
    end
end

% aux(strcmp(R.St.StNames,'SO4')) = aux(strcmp(R.St.StNames,'NH3'))/2;
% controlpH( R.pTh.Keq, R.pTh.chrM, R.St.StNames, R.pOp.pH);

Sbc_Dir_new = aux;

if sum(isnan(Sbc_Dir_new(k))) > 0
    fprintf('Some Sbc_Dir NaN!!\n')
end
for k = 1:numStVLiq2
    if isnan(Sbc_Dir_new(k))
        R.Sxy.Sbc_Dir(k) = Sxy.Sbc_Dir(k);
    else
        R.Sxy.Sbc_Dir(k) = Sbc_Dir_new(k);
    end
end
R.pOp.invHRT = invHRT;

function dy = massbal(~,y,G, Xi, Ssp, flagN)
    dy = zeros(length(y),1);
    if flagN == 1 && y(1) >= Ssp
        if G(1) > 0
            invHRTaux = 0;
        else
            invHRTaux = -G(1)/(Xi(1) - y(1));
        end
        if ne(invHRTaux, 0)
            invHRT = invHRTaux;
        end
    else
        dy(1) = invHRT*(Xi(1) - y(1)) + 0;
    end
end

% function controlpH(Keq, chrM, StNames, pH) %#ok<DEFNU>
% % u = [aux(1:numStVLiq3); zeros((numStVLiq-numStVLiq3),1);1;0];
% Tol = 1e-14;
% Tp = 1;
% u = [aux;0;1;0];
% 
% NaHCO3 = aux(strcmp(StNames,'CO2'));
% w = 1;
% 
% spcM = zeros(size(chrM));
% Sh = 10^(-pH);
% 
% while abs(Tp) > Tol
% u(strcmp(StNames,'Na')) = NaHCO3;
% u(strcmp(StNames,'CO2')) = NaHCO3;
% 
% Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
% 
% spcM(:,1) = ((Keq(:,1)/w).*u*Sh^3)                        ./Denm;        
% spcM(:,2) = (u * Sh^3)                                    ./Denm;
% spcM(:,3) = (u * Sh^2 .* Keq(:,2))                        ./Denm;
% spcM(:,4) = (u * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
% spcM(:,5) = (u      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
% Tp = Sh + sum(sum(spcM.*chrM));
% 
% NaHCO3 = NaHCO3 - Tp;
% 
% end
% % % posCspcM = CspcM > 0;
% % % Cl = (sum(sum(CspcM.*posCspcM))+Sh)/(-chrM(find(strcmp(StNames,'Cl'),1),2));
% % % Na = sum(sum(CspcM.*(1-posCspcM)))/(-chrM(find(strcmp(StNames,'Na'),1),2));
% % Na = sum(sum(CspcM))/(-chrM(find(strcmp(StNames,'Na'),1),2))
% % if Na < 0
% %     error('Correct IC initial concentration -- Na < 0 ')
% % end
% aux(strcmp(StNames,'Na')) = NaHCO3;
% aux(strcmp(StNames,'CO2')) = NaHCO3;
% 
% end
end