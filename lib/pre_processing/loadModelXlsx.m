% %%%% loadModelXls.m -Load all the parameters from Excel sheet and
% generates the initial structure

function R = loadModelXlsx(filename)
R = []; % R initialization
rng(082021);


% GENERAL MODEL PARAMETERS from the Excel file
fprintf('> LOADING AND CREATING MODEL STRUCTURE AND PARAMETERS...\n')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters of the simulation (Temperature, pH, ...)
[OperatParam, pOpNames]   = xlsread(filename, strcat('Parameters'));
aux = '';
for i=1:length(OperatParam)
    % Building the string command for the structure, last without ', '
    if i<length(OperatParam)
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i)),', ');
    else
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i)));
    end
end
eval(strcat('pOp = struct(', char(aux), ');'));
pOp.Vgas = pOp.Vgas*1000; %#ok<NODEF> %L
pOp.Vr = pOp.Vr*1000; %L
pOp.invHRT = 1/pOp.HRT;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE VARIABLES Structure and initial values
[States, StNames] = xlsread(filename, strcat('States'));
St.StV = States(:,1);
St.Phase = StNames(:,end);
St.StNames = StNames(:,1);
St.StVLiq = St.StV(1:(find(strcmp(St.Phase,'S'),1)-1));
St.StVLiq2 = St.StV(1:(find(strcmp(St.Phase,'G'),1)-1));
St.StVLiq3 = St.StV(1:(find(strcmp(St.Phase,'P'),1)-1));
St.StVIni = St.StVLiq;
St.StVX = St.StV(find(strcmp(St.Phase,'S'),1):end);
St.numX = length(St.StVX);
St.numStVLiq = length(St.StVLiq);
St.numStVLiq2 = length(St.StVLiq2);
St.numStVLiq3 = length(St.StVLiq3);
St.numStVGas = St.numStVLiq - St.numStVLiq2;
St.numSt = length(St.StNames);

name = char(St.StNames(St.numStVLiq2 + 1:St.numStVLiq));
name = name(:,2:end);
Liq2Gas = zeros(St.numSt, 1);
for k = 1 : St.numStVGas
    Liq2Gas(strcmp(strcat(name(k, :)),St.StNames)) = k;
end
St.Liq2Gas = Liq2Gas;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the microbial species (qmax, Ks, Yield)
[SpecParam, SpecNames] = xlsread(filename, strcat('SpecParam'));

% aux = strcmp('mu_max',SpecNames(1,2:end));
% pTh.mu_max = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_NH3',SpecNames(1,2:end));
Ks_NH3 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_NO2',SpecNames(1,2:end));
Ks_NO2 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_O2',SpecNames(1,2:end));
Ks_O2 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ki_NH3',SpecNames(1,2:end));
Ki_NH3 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ki_HNO2',SpecNames(1,2:end));
Ki_HNO2 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ki_O2',SpecNames(1,2:end));
Ki_O2 = SpecParam((isfinite(SpecParam(:,aux))),aux);

pTh.Ks = [Ks_NH3, Ks_NO2, Ks_O2, Ki_NH3, Ki_HNO2, Ki_O2];

aux = strcmp('Yield',SpecNames(1,2:end));
pTh.Y = SpecParam((isfinite(SpecParam(:,aux))),aux);

% aux = strcmp('DGdis',SpecNames(1,2:end));
% pTh.DGdis = SpecParam((isfinite(SpecParam(:,aux))),aux);

% aux = strcmp('Maintenance',SpecNames(1,2:end));
% pTh.maint = SpecParam((isfinite(SpecParam(:,aux))),aux);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOICHIOMETRY/REACTION MATRIX
[rMatrixFull, rNames] = xlsread(filename, strcat('ReactionMatrix'));
rm.rMatrixFull = rMatrixFull(:,[1:3:end;2:3:end]);
rm.rNamesX = rNames(1,2:3:end)';
rm.rNamesS = rNames(3:end,1);
MatrixCat = rMatrixFull(:,1:3:end);
MatrixAn = rMatrixFull(:,2:3:end);
MatrixDecay = rMatrixFull(:,3:3:end);
St.numStFull = length(MatrixCat);

MatrixCat_mod = [MatrixCat(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); MatrixCat((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)];
namesAux = [rm.rNamesS(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); rm.rNamesS((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)];
MatrixCat_mod = [MatrixCat_mod(1:(find(strcmp(namesAux, 'H'))-1),:);MatrixCat_mod((find(strcmp(namesAux, 'H'))+1):end,:)];
MatrixCat_mod = MatrixCat_mod(1:St.numStVLiq,:);
MatrixAn_mod = [MatrixAn(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); MatrixAn((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)];
MatrixAn_mod = [MatrixAn_mod(1:(find(strcmp(namesAux, 'H'))-1),:);MatrixAn_mod((find(strcmp(namesAux, 'H'))+1):end,:)];
MatrixAn_mod = MatrixAn_mod(1:St.numStVLiq,:);
MatrixDecay_mod = [MatrixDecay(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); MatrixDecay((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)];
MatrixDecay_mod = [MatrixDecay_mod(1:(find(strcmp(namesAux, 'H'))-1),:);MatrixDecay_mod((find(strcmp(namesAux, 'H'))+1):end,:)];
MatrixDecay_mod = MatrixDecay_mod(1:St.numStVLiq,:);

for k = 1:St.numX
    rm.MatrixMet(:,k) = MatrixCat_mod(:,k)*(1/pTh.Y(k)) + MatrixAn_mod(:,k);
end

rm.MatrixDecay_mod = MatrixDecay_mod;

eD = SpecNames(2:end,end);
lCat_eD = zeros(St.numX,1);
for i=1:St.numX
    lCat_eD(i) = -MatrixAn_mod(strcmp(rm.rNamesS, eD(i,:)), i);
end
pTh.lCat_eD = lCat_eD;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure with the THERMODYNAMIC PARAMETERS
[pThValues, pThNames]   = xlsread(filename, strcat('ThermoParam'));
pThValues = pThValues(2:end,:);
pThNames = pThNames(2:end,:);
DG0 = pThValues(1:(size(pThValues,1)/2), 1:(end-1));
nan_locations = isnan(DG0);
DG0(nan_locations) = 1e10;
chrM = pThValues((size(pThValues,1)/2 +1):size(pThValues,1), 1:(end-1));
nan_locations = isnan(chrM);
chrM(nan_locations) = 0;

sDG0col = size(DG0,2);

Keq = zeros(length(DG0),sDG0col-2);
DG0H2O = -237.18; % Water kJ/mol
Keq(:,1) = exp((DG0H2O + DG0(:,1) - DG0(:,2))/(-pOp.Rth*pOp.T));
for i=2:(sDG0col - 1)
    Keq(:,i) = exp((DG0(:,i+1)+(i*1e10*(DG0(:,i+1) == 1e10))-(DG0(:,i)))/(-pOp.Rth*pOp.T));
end
spcR = accumarray([(1:(St.numStFull))', pThValues(1:length(pThNames)/2,end)],1,[St.numStFull, sDG0col]);

pTh.Phase = pThNames(1:length(pThNames)/2,end);
pTh.chrM = chrM(1:(find(strcmp(pTh.Phase,'S'),1)-1), :);
pTh.Keq = Keq(1:(find(strcmp(pTh.Phase,'S'),1)-1), :);
pTh.DGr0 = rm.rMatrixFull'*sum(spcR.*DG0,2);
pTh.DG0 = DG0;
spcR = spcR(1:(find(strcmp(pTh.Phase,'S'),1)-1), :);
pTh.spcR = spcR == 1;
pTh.react_v = pThValues(1:length(pThNames)/2,end);
vaux = find(strcmp(St.Phase,'G'));
vaux = vaux(ne(vaux, find(strcmp(St.StNames,'gN2'))));
KhV = zeros(length(St.StV),1);
for i=1:length(vaux)
    name = char(St.StNames(vaux(i)));
    name = name(2:end);
    w = strcmp(strcat(name),St.StNames);
    if DG0(w,1)== 1e10
        Kh.(name) = exp((DG0(w,2) - DG0(vaux(i),2))  / (-pOp.Rth*pOp.T));% M/atm
    else
        Kh.(name) = exp((DG0(w,1) - DG0(vaux(i),2))  / (-pOp.Rth*pOp.T));% M/atm
    end
    KhV(i) = Kh.(name);
end

pTh.Kh = 0;
% pTh.Kh.KhV = KhV;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the IbModel (maximum cell's weight, parameters for shoving...)
[BacteriaParam, BacteriaNames]   = xlsread(filename, strcat('Bacteria'));
aux = '';
for i=1:length(BacteriaParam)
    % Building the string command for the structure, last without ', '
    if i<length(BacteriaParam)
        aux = strcat(aux, char(39), BacteriaNames(i), char(39), ', ',num2str(BacteriaParam(i)),', ');
    else
        aux = strcat(aux, char(39), BacteriaNames(i), char(39), ', ',num2str(BacteriaParam(i)));
    end
end
eval(strcat('bac = struct(', char(aux), ');'));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Space discretization
[DiscretizationParam, DiscretizationNames] = xlsread(filename, strcat('Discretization'));
aux = '';
for i=1:length(DiscretizationParam)
    % Building the string command for the structure, last without ', '
    if i<length(DiscretizationParam)
        aux = strcat(aux, char(39), DiscretizationNames(i), char(39), ', ',num2str(DiscretizationParam(i)),', ');
    else
        aux = strcat(aux, char(39), DiscretizationNames(i), char(39), ', ',num2str(DiscretizationParam(i)));
    end
end
eval(strcat('Sxy = struct(', char(aux), ');'));

Sxy.nxSys = Sxy.nx +2; %#ok<NODEF>
Sxy.nySys = Sxy.ny +2;
Sxy.nTSys=(Sxy.nxSys-2)*(Sxy.nySys-2);
Sxy.maxxSys = Sxy.maxx;
Sxy.maxySys = Sxy.maxy;
% Sxy.xSys = 0:Sxy.dx:Sxy.maxxSys;
% Sxy.ySys = 0:Sxy.dy:Sxy.maxySys;
% xp = Sxy.xSys+ Sxy.dx/2;
% Sxy.x_pnSys = xp(1:end-1);
% yp = Sxy.ySys+ Sxy.dy/2;
% Sxy.y_pnSys = yp(1:end-1);
Sxy.Vg = (Sxy.dx*Sxy.dy*Sxy.dz)*1000; %L


[Influent, BLayer] = xlsread(filename, strcat('Influent'));
R.Inf.St = Influent(:,1);
BLayer = BLayer(:,end);
Dir_k = zeros(St.numStVLiq2,1);
for k = 1: St.numStVLiq2
    if strcmp(BLayer(k),'D')
        Dir_k(k) = 1;
    else
        Dir_k(k) = 0;
    end
end
R.Inf.Dir_k = Dir_k;
Sxy.Sbc_Dir = St.StVIni;
for k = 1: St.numStVLiq2
    if Dir_k(k) %Dirichlet
        Sxy.Sbc_Dir(k) = R.Inf.St(k);
    end
end

% Sxy.pHnT = pOp.pH*ones(Sxy.nTSys,1);%%%%OJO

bac.bac_rho_bio = 0;
bac.bac_c =(Sxy.maxySys)/2;

xaux = ((Sxy.dx+Sxy.T_blayer):(2*bac.bac_rmax):(Sxy.maxxSys-(Sxy.dx+Sxy.T_blayer)));
yaux = ((Sxy.dy+Sxy.T_blayer):(2*bac.bac_rmax):(Sxy.maxxSys-(Sxy.dy+Sxy.T_blayer)));
x_pn = kron(xaux, ones(size(yaux)));
y_pn = kron(ones(size(xaux)),yaux);

bac_x = []; bac_y = []; bac_s = [];
Vlot = 0:2*bac.bac_rmax:(bac.bac_ystart/2);
%- AOB-NOB seeds -%
% for i = 1: length(x_pn)
%     dxx = x_pn(i)-bac.bac_c;
%     dyy = y_pn(i)-bac.bac_c;
%     normd = sqrt(dxx*dxx+dyy*dyy);
%      if normd <= (bac.bac_ystart/2)
%         bac_x(end+1,1) = x_pn(i); %#ok<AGROW>
%         bac_y(end+1,1) = y_pn(i); %#ok<AGROW>
%         norm = sqrt((x_pn(i)-bac.bac_c)*(x_pn(i)-bac.bac_c)+(y_pn(i)-bac.bac_c)*(y_pn(i)-bac.bac_c));
%         bs_aux = (round(round(norm*10^(6))/2) == (round(norm*10^(6))/2))+1;
%         if bs_aux == 2
%             bs_aux = bs_aux + floor(rand()*2); % Nitrobacter/NOB = Nitrospira/NOB = 50%
%         end
%         bac_s(end+1,1) = bs_aux; %#ok<AGROW>
%      end
% end

%- 100% random distribution -%
for i = 1: length(x_pn)
    dxx = x_pn(i)-bac.bac_c;
    dyy = y_pn(i)-bac.bac_c;
    normd = sqrt(dxx*dxx+dyy*dyy);
    if normd < Vlot(end)
        for j = 1: length(Vlot)-1
            SSaux = floor(rand*(St.numX-1))+1;
        end
        if SSaux == 3
            SSaux = 4; %AMX
        end
        if SSaux == 2
            SSaux = SSaux + 1; % Nitrospira/NOB = 100% ... double(rand>0.5)
        end
        bac_s(end+1,1) = SSaux; %#ok<AGROW>
        bac_x(end+1,1) = x_pn(i); %#ok<AGROW>
        bac_y(end+1,1) = y_pn(i); %#ok<AGROW>
        norm_bac = sqrt((bac_x-bac.bac_c).*(bac_x-bac.bac_c)+(bac_y-bac.bac_c).*(bac_y-bac.bac_c)); %#ok<NASGU>
    end
end

auxSize = size(bac_y);

bac_mm = 0.9*bac.bac_mmax*ones(auxSize);
bac_m = bac_mm/bac.bac_MW; %mol
St.StVX = bac_m;
bac_a = ones(size(bac_m));
bac_r = ((bac_mm/bac.bac_rho)*(3/(4*pi))).^(1/3);%radius
bac.bac_n = auxSize(1); % number of cells
% bac_mu_max = pTh.mu_max(bac_s);
bac_yield = pTh.Y(bac_s);
bac.bac_Ks = pTh.Ks(bac_s, :);
% bac_maint = pTh.maint(bac_s);
bac_ns = zeros(St.numX, 1);
for i  = 1:St.numX
    for j = 1:auxSize(1)
        bac_ns(i) =  (i == bac_s(j)) + bac_ns(i); % number of cells of each
    end
end
bac.bac_ns = bac_ns;
bac_act = ones(size(bac_x));
bac.atrib = [bac_x, bac_y, bac_m, bac_a, bac_s, bac_r, bac_yield, bac_act];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DiffusionParam, ~] = xlsread(filename, strcat('Diffusion'));

kTr.Diffn = DiffusionParam(1:St.numStVLiq2,:);
% kTr.DiffNames = DiffusionNames(1:St.numStVLiq2,1);
kTr.kLa = DiffusionParam(St.numStVLiq2+1:end,:);
% kTr.kLaNames = DiffusionNames(St.numStVLiq2+1:end,1);
% kTr.DiffwaternT = kron(kTr.Diffn, ones(Sxy.nTSys,1));%%OJO

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sxy.pos_xySys = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[TolParam] = xlsread(filename, strcat('Tolerances'));

kTr.Tolabs = TolParam(1,:);
kTr.Tolrel = TolParam(2:end,:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Keq = pTh.Keq; chrM = pTh.chrM; spcR = pTh.spcR; w = 1; numStVLiq = St.numStVLiq; %#ok<NASGU>
spcM = zeros(size(chrM));
StV = [St.StVLiq;1;0];
Sh = 10^(-pOp.pH);
Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);

spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
spcM(:,2) = (StV * Sh^3)                                    ./Denm;
spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;

spcMc = (spcM == 0)*1e-20 + spcM; %#ok<NASGU>
% DGr = pTh.DGr0 + pOp.Rth*pOp.T*(rm.rMatrixFull'*[nansum(spcR.*log(spcMc),2);ones(St.numStFull-(numStVLiq+2),1)]);
% R.DGCat_bulk =  DGr(1:2:length(DGr));
% R.DGAn_bulk =  DGr(2:2:length(DGr));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Update of the structure R
R.pOp = pOp;
R.pTh = pTh;
R.rm = rm;
R.St = St;
R.bac = bac;
R.Sxy = Sxy;
R.kTr = kTr;


R.flagGas = 1;

R.flagDG = 1;

fprintf('\n>>>>>>>>>>>>>>>>>>>>>>> VERSION WITH pH VARIABLE >>>>>>>>>>>>>>>>>>>>>>>>>>\n\n')

prompt = '>> Do you want to keep the NH3 fixed (variable HRT)? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
x = input(prompt);
if isempty(x)
    x = 1;
end
if x == 1
    fprintf('--> YES\n')
    %         if R.Inf.St(1) ==  R.St.StVIni(1)
    %             error('NH3 influent = NH3 initial - Change the NH3 influent!')
    %         end
else
    fprintf('--> NO\n')
    x = 2;
    if ne(R.Inf.St(1),  R.St.StVIni(1))
        prompt = '>> NH3 influent is different than NH3 initial. Do you want to continue? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
        y = input(prompt);
        if isempty(y)
            y = 1;
        end
        if y == 1
            fprintf('--> YES\n')
        else
            fprintf('--> NO\n')
            error('NH3 influent is different than NH3 initial')
        end
    end
    
end
R.flagN = x;
assignin('base', 'R', R) % Writing R in the Workspace
clear x prompt

fprintf('> ... MODEL STRUCTURE AND PARAMETERS LOADED.')
fprintf('\n------------ooo\n')
VN = Sxy.dT*kTr.Diffn/(Sxy.dx*Sxy.dy);
for i =1:St.numStVLiq2
    fprintf('Von Neumann coefficient of stability for %s: %.2e\n', char(St.StNames(i)),VN(i))
end
fprintf('>>> If Von Neumann coefficient larger than 0.5, consider Backward Euler instead of Crank-Nicolson')
fprintf('\n------------ooo\n')
for i = 1:(R.St.numX)
    fprintf('Number of cells of type %d Name %s : %d\n', i,char(R.rm.rNamesX(i)),R.bac.bac_ns(i))
end
fprintf('>>> INITIAL TOTAL Number of cells: %d', R.bac.bac_n)
fprintf('\n------------ooo\n')
fprintf('\n>>> Push any key to start the simulation...\n')
pause()