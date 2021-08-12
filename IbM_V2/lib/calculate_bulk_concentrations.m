function [bulk_concentrations, invHRT] = calculate_bulk_concentrations(constants, prev_conc, invHRT, reactionMatrix, dT)
    % Function to calculate the bulk layer concentrations. Assumes that the
    % simulated bio-aggregate is representative of the entire reactor. 
    %
    % Includes pH correction using NaHCO3
    % 
    % Has a correction for when accumulated reaction in the aggregate
    % exceeds the available supply/concentration of compounds.
    %
    % constants: struct containing all constants
    % prev_conc: vector with the previous bulk concentration per compound
    % invHRT: 1 / (hydrolic retention time)
    % reactionMatrix: matrix with per grid cell and per compound the
    %   change [h-1] due to bacterial activity 
    % dT: time of integration with diffusion
    %
    % -> bulk_concentrations: vector with the bulk concentration per
    %   compound
    % -> invHRT: new 1 / (hydrolic retention time) [h-1]

    %% unpack constants for easy use
    Keq = constants.Keq;
    chrM = constants.chrM;
    StNames = constants.StNames;
    pH = constants.pHsetpoint;
    Vr = constants.Vr;
    Vg = constants.Vg;
    Dir_k = constants.Dir_k;
    isLiquid = constants.isLiquid;
    influent = constants.influent_concentrations;
    NH3sp = constants.pOp.NH3sp;
    keepNH3fixed = constants.constantN;
    
    %% set bulk_concentrations 
    if isscalar(reactionMatrix) % if no reaction matrix formed, then init with previous concentrations
        bulk_concentrations = prev_conc;
    else
        cumulative_reacted = squeeze(sum(reactionMatrix(:,:,isLiquid), [1,2])) * Vg / Vr;
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-20, 'NonNegative', ones(size(cumulative_reacted)));
        [~, Y] = ode45(@(t, y) massbal(t, y, cumulative_reacted(isLiquid), influent, NH3sp, keepNH3fixed), [0 dT], prev_conc(isLiquid), options);
        bulk_conc_temp = Y(end, :)';
        bulk_concentrations = correct_negative_concentrations(bulk_conc_temp);
        bulk_concentrations(Dir_k) = prev_conc;
    end
    
    %% apply pH correction to bulk_concentrations
    bulk_concentrations(strcmp(StNames, 'SO4')) = bulk_concentrations(strcmp(StNames, 'NH3')) / 2;        
    bulk_concentrations(1:sum(isLiquid)) = controlpH(Keq, chrM, StNames, pH, bulk_concentrations(isLiquid));
    
    if any(bulk_concentrations < 0)
        warning('DEBUG:actionRequired', 'debug: negative bulk concentration encountered... correction required?')
    end
    
    %% helper function
    function dy = massbal(~, bulk_conc, cumulative_reacted, reactor_influx, NH3_reactor_influx, keepNH3fixed)
        % differential equation for the mass balance over the entire reactor
        dy = zeros(length(bulk_conc), 1);

        if keepNH3fixed == 1 && bulk_conc(1) >= NH3_reactor_influx
            if cumulative_reacted(1) < 0
                invHRT = -cumulative_reacted(1) / (reactor_influx(1) - bulk_conc(1));
            end
        else
            dy(1) = invHRT * (reactor_influx(1) - bulk_conc(1)) + 0;
        end

        dy(2:end) = invHRT * (reactor_influx(2:end) - bulk_conc(2:end)) + cumulative_reacted(2:end);
        dy(4:end) = 0;
    end
end


%%
function conc = controlpH(Keq, chrM, StNames, pH, conc)
    % Add NaHCO3 to control the pH at a certain point
    %
    % Keq:
    % chrM:
    % StNames: names of the compounds
    % pH: setpoint of the pH
    % conc: bulk concentration per compound before NaHCO3 addition
    %
    % -> conc: bulk concentration per compound after pH control
    
    Tol = 1e-14;
    Tp = 1;
    u = [conc; 0; 1; 0]; % C: why [..., 0, 1, 0] ???

    NaHCO3 = conc(strcmp(StNames, 'CO2'));
    w = 1;

    spcM = zeros(size(chrM));
    Sh = 10^(-pH);

    while abs(Tp) > Tol
        u(strcmp(StNames, 'Na')) = NaHCO3;
        u(strcmp(StNames, 'CO2')) = NaHCO3;

        Denm = (1 + Keq(:, 1) / w) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);

        spcM(:, 1) = ((Keq(:, 1) / w) .* u * Sh^3) ./ Denm;
        spcM(:, 2) = (u * Sh^3) ./ Denm;
        spcM(:, 3) = (u * Sh^2 .* Keq(:, 2)) ./ Denm;
        spcM(:, 4) = (u * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
        spcM(:, 5) = (u .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;
        Tp = Sh + sum(sum(spcM .* chrM));

        NaHCO3 = NaHCO3 - Tp;
    end

    conc(strcmp(StNames, 'Na')) = NaHCO3;
    conc(strcmp(StNames, 'CO2')) = NaHCO3;
end


%%
function conc = correct_negative_concentrations(conc)
    % Perform a correction to get rid of negative concentrations.
    %
    % conc: bulk concentrations before correction
    %
    % -> conc: bulk concentrations after correction
    
    negative_indices = find(conc < 0);

    if negative_indices
        warning('DEBUG:noActionRequired', 'debug: negative concentration encountered... correction applied')
        for k = 1:length(negative_indices)
            % if [NH3] < 0, then remove excess consumption from [NO2]
            if negative_indices(k) == 1 
                conc(2) = conc(2) + conc(1);
                conc(1) = 0;
                
            % if [NO2] < 0, then remove excess consumption from [NO3]
            elseif negative_indices(k) == 2
                conc(3) = conc(3) + conc(2);
                conc(2) = 0;
            end
        end
    end
end

