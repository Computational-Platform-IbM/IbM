function [bulk_concentrations, invHRT] = calculate_bulk_concentrations(constants, prev_conc, invHRT, reactionMatrix, dT, settings)
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
    compoundNames = constants.compoundNames;
    pH = constants.pHsetpoint;
    Vr = constants.Vr;
    Vg = constants.Vg;
    Dir_k = constants.Dir_k;
    influent = constants.influent_concentrations;
    % NH3sp = constants.pOp.NH3sp;
    % keepNH3fixed = constants.constantN;
    variableHRT = settings.variableHRT;
    if variableHRT
        bulk_setpoint = constants.bulk_setpoint;
        setpoint_index = constants.setpoint_index;
    else
        bulk_setpoint = 0;
        setpoint_index = 0;
    end
    
    %% set bulk_concentrations 
    if isscalar(reactionMatrix) % if no reaction matrix formed, then init with previous concentrations
        bulk_concentrations = prev_conc;
    else
        cumulative_reacted = squeeze(sum(reactionMatrix, [1,2])) * Vg / Vr; % [mol/h] / [L]
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-20, 'NonNegative', ones(size(cumulative_reacted)));
        try
            [~, Y] = ode45(@(t, y) massbal(t, y, cumulative_reacted, influent, variableHRT, bulk_setpoint, setpoint_index, Dir_k, settings), [0 dT], prev_conc, options);
            bulk_conc_temp = Y(end, :)';
            bulk_concentrations = correct_negative_concentrations(bulk_conc_temp); %<E: Negative concentration from mass balance of reactor. />
        catch e % really bad practice....
            fprintf(2, '\n\nSomething went wrong...\nreason: %s\n', e.identifier)
%             bulk_concentrations = prev_conc;
%             warning(e)
            rethrow(e);
        end
        
        % set dirichlet boundary condition
        bulk_concentrations(Dir_k) = prev_conc(Dir_k);

    end
    
    %% apply pH correction to bulk_concentrations
    if settings.pHbulkCorrection
        bulk_concentrations(strcmp(compoundNames, 'SO4')) = bulk_concentrations(strcmp(compoundNames, 'NH3')) / 2;        
        bulk_concentrations = controlpH(Keq, chrM, compoundNames, pH, bulk_concentrations);
    
        if any(bulk_concentrations < 0) %<E: Negative concentration from control pH of reactor. />
            warning('DEBUG:actionRequired', 'debug: negative bulk concentration encountered after pH control... correction required?')
            bulk_concentrations = bulk_concentrations .* (bulk_concentrations > 0);
        end
    end
    
    %% helper function
    function dy = massbal(~, bulk_conc, cumulative_reacted, reactor_influx, variableHRT, bulk_setpoint, setpoint_index, Dir_k, settings)
        % Differential equation for the mass balance over the entire reactor
        % Will modify the HRT to match the setpoint of NH3 iff the outflow
        % concentration is larger than the setpoint.
        %
        % bulk_conc: concentration in the bulk liquid [mol/L]
        % cumulative_reacted: cumulative reaction rate [mol/L/h] per
        %   compound
        % reactor_influx: influx concentrations [mol/L] per compound
        % bulk_setpoint: setpoint for outflux of compound <setpoint_index>
        % variableHRT: boolean representing whether we want to control the
        %   the concentration in the outflow of the reactor
        % (invHRT): 1/(Hydrolic retention rate) [1/h]
        %
        % -> dy: derivative of bulk concentration
        
        dy = zeros(length(bulk_conc), 1);

        if settings.structure_model
            switch settings.type
                case 'Neut'
                    if variableHRT == 1 && (bulk_conc(1) > NH3sp || bulk_conc(1) < NH3sp)
                        if cumulative_reacted(1) < 0
                            invHRT = -cumulative_reacted(1) / (reactor_influx(1) - NH3sp);
                        end
                    else
                        dy(1) = invHRT * (reactor_influx(1) - bulk_conc(1)) + cumulative_reacted(1);
                    end
                    
                    if variableHRT == 1 && (bulk_conc(2) > NH3sp || bulk_conc(2) < NH3sp)
                        if cumulative_reacted(2) < 0
                            invHRT = -cumulative_reacted(2) / (reactor_influx(2) - NH3sp);
                        end
                    else
                        dy(2) = invHRT * (reactor_influx(2) - bulk_conc(2)) + cumulative_reacted(2);
                    end
                    
                    if variableHRT == 1 && (bulk_conc(3) > NH3sp || bulk_conc(3) < NH3sp)
                        if cumulative_reacted(3) < 0
                            invHRT = -cumulative_reacted(3) / (reactor_influx(3) - NH3sp);
                        end
                    else
                        dy(3) = invHRT * (reactor_influx(3) - bulk_conc(3)) + cumulative_reacted(3);
                    end
                    dy(4:end) = invHRT * (reactor_influx(4:end) - bulk_conc(4:end)) + cumulative_reacted(4:end);
                    dy(5:end) = 0;
                case {'Comp','Comm','Copr'}
                    if variableHRT == 1 && (bulk_conc(1) > NH3sp || bulk_conc(1) < NH3sp)
                        if cumulative_reacted(1) < 0
                            invHRT = -cumulative_reacted(1) / (reactor_influx(1) - NH3sp);
                        end
                    else
                        dy(1) = invHRT * (reactor_influx(1) - bulk_conc(1)) + cumulative_reacted(1);
                    end

                    dy(2:end) = invHRT * (reactor_influx(2:end) - bulk_conc(2:end)) + cumulative_reacted(2:end);
                    dy(5:end) = 0;
                otherwise
                    error(['Type <', settings.type, '> is not a registered set of simulations.'])
            end
            
        else
            if variableHRT && bulk_conc(setpoint_index) > bulk_setpoint
                % if bulk conc of set compound is larger than setpoint, decrease HRT
                if cumulative_reacted(setpoint_index) < 0
                    invHRT = -cumulative_reacted(setpoint_index) / (reactor_influx(setpoint_index) - bulk_setpoint); 
                end
            end

            % always calculate all non-dirichlet bulk concentrations
            dy(~Dir_k) = invHRT * (reactor_influx(~Dir_k) - bulk_conc(~Dir_k)) + cumulative_reacted(~Dir_k);
        end
    end
end


%%
function conc = controlpH(Keq, chrM, compoundNames, pH, conc)
    % Add NaHCO3 to control the pH at a certain point
    %
    % Keq:
    % chrM:
    % compoundNames: names of the compounds
    % pH: setpoint of the pH
    % conc: bulk concentration per compound before NaHCO3 addition
    %
    % -> conc: bulk concentration per compound after pH control
    
    Tol = 1e-14;
    Tp = 1;
    u = [conc; 1; 0]; % for H2O and H concentrations

    NaHCO3 = conc(strcmp(compoundNames, 'CO2'));
    w = 1;

    spcM = zeros(size(chrM));
    Sh = 10^(-pH);

    while abs(Tp) > Tol
        u(strcmp(compoundNames, 'Na')) = NaHCO3;
        u(strcmp(compoundNames, 'CO2')) = NaHCO3;

        Denm = (1 + Keq(:, 1) / w) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);

        spcM(:, 1) = ((Keq(:, 1) / w) .* u * Sh^3) ./ Denm;
        spcM(:, 2) = (u * Sh^3) ./ Denm;
        spcM(:, 3) = (u * Sh^2 .* Keq(:, 2)) ./ Denm;
        spcM(:, 4) = (u * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
        spcM(:, 5) = (u .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;
        Tp = Sh + sum(sum(spcM .* chrM));

        NaHCO3 = NaHCO3 - Tp;
    end

    conc(strcmp(compoundNames, 'Na')) = NaHCO3;
    conc(strcmp(compoundNames, 'CO2')) = NaHCO3;
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
        warning('DEBUG:noActionRequired', 'debug: negative concentration encountered and corrected')
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

