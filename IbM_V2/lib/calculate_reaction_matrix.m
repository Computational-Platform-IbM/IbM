function [reaction_matrix, mu, pH] = calculate_reaction_matrix(grid2bac, grid2nBacs, bac, grid, conc, constants, pH)
    % Calculate how much compounds are consumed/produced per grid cell due
    % to bacterial activity. Also updates the growth rate of the respective
    % bacteria.
    % 
    % grid2bac: matrix with per grid cell which bacteria reside there
    % grid2nBacs: matrix with how many bacteria per grid cell
    % bac: struct containing all information regarding the bacteria
    % grid: struct containing all information regarding the grid
    % conc: matrix containing all concentrations per grid cell as of
    %   (ix, iy, compound)
    % constants: struct containing all simulation constants
    % pH_old: matrix with per grid cell the pH
    %
    % -> reaction_matrix: matrix with per grid cell and per compound the
    %   change [h-1] due to bacterial activity
    % -> mu: vector with updated growth rates per bacterium
    % -> pH: matrix with per grid cell the pH

    % convert init pH to pH matrix
    if isscalar(pH)
        pH = ones(size(grid2nBacs))*pH;
    end
    
    % init
    T = constants.T;
    Keq = constants.Keq;
    chrM = constants.chrM;
    Vg = constants.Vg;
    compounds = constants.StNames;
    reactive_form = constants.react_v;
    Ks = constants.Ks;
    Ki = constants.Ki;
    mMetabolism = constants.MatrixMet;
    mDecay = constants.MatrixDecay;
    
    iNH3 = strcmp(compounds, 'NH3');
    iNO2 = strcmp(compounds, 'NO2');
    iO2 = strcmp(compounds, 'O2');
    
    reaction_matrix = zeros(size(conc));
    mu = zeros(size(bac.x));
    
    % for each gridcell
    for ix = 1:grid.nX % parfor?
        for iy = 1:grid.nY
            % calculate pH & speciation
            Sh_old = 10^-pH(ix, iy);
            [spcM, Sh] = solve_pH(Sh_old, [squeeze(conc(ix, iy, :)); 1; 0], Keq, chrM, 0); % <C: remove flagpH???/> <C: why [...; 1; 0]? />
            pH(ix, iy) = -log10(Sh);
                
            if grid2nBacs(ix, iy)
                % get bacteria in this grid cell
                iBacs = squeeze(grid2bac(ix, iy, 1:grid2nBacs(ix, iy))); % n-by-1 vector of bacterial indices
                
                speciesGrid = bac.species(iBacs);
                unique_species = unique(speciesGrid); % which species
                for i = 1:length(unique_species)
                    curr_species = unique_species(i);
                    
                    % update mu_max per species based on pH
                    [mu_max, maint] = determine_max_growth_rate_and_maint(curr_species, T, Sh);

                    % get reactive concentrations for Ks & Ki
                    reactive_conc = [spcM(iNH3, reactive_form(iNH3)), ...
                        spcM(iNO2, reactive_form(iNO2)), ...
                        spcM(iO2, reactive_form(iO2))];
                    
                    % set mu for all bacteria of same species in that gridcell
                    M = calculate_monod(Ks(curr_species,:), Ki(curr_species, :), reactive_conc);
                    mu_noMaintenance = M * mu_max;
                    mu_withMaintenance = mu_noMaintenance - maint;
                    mu(iBacs(speciesGrid == curr_species)) = mu_withMaintenance;

                    % update reaction_matrix element for this grid cell
                    concentrationChange = mMetabolism(:, curr_species) * mu_noMaintenance;
                    if mu_withMaintenance < 0
                        concentrationChange = concentrationChange - mDecay(:, curr_species) * mu_withMaintenance;
                    end
                    concentrationChange = concentrationChange * sum(bac.molarMass(iBacs(speciesGrid == curr_species)));
                    
                    reaction_matrix(ix, iy, :) = reaction_matrix(ix, iy, :) + reshape(concentrationChange, 1,1,[]);
                end

                % divide by Vg?
                reaction_matrix = reaction_matrix / Vg;
            end
        end
    end
end

function [spcM, Sh] = solve_pH(Sh_ini, StV, Keq, chrM, flagpH)  % <C: remove flagpH???/>
    % Solve the pH and speciation per grid cell
    %
    % <>
    % original as copied from latest version, to be improved?
    % </>
    %
    % -> spcM: species matrix for the respective pH
    % -> Sh: proton concentration
    
    w = 1;

    if flagpH == 1
        spcM = zeros(size(chrM));
        Sh = Sh_ini;
        Denm = (1 + Keq(:, 1) / w) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);

        spcM(:, 1) = ((Keq(:, 1) / w) .* StV * Sh^3) ./ Denm;
        spcM(:, 2) = (StV * Sh^3) ./ Denm;
        spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
        spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
        spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;

    else
        % Checking the existence of a zero pool in the function between pH 1 and 14
        a = 1e-14;
        b = 1;
        Sh_v = [a; b]; F = zeros(2, 1);

        for nn = 1:2
            Sh = Sh_v(nn);
            spcM = zeros(size(chrM));
            Denm = (1 + Keq(:, 1) / w) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);

            spcM(:, 1) = ((Keq(:, 1) / w) .* StV * Sh^3) ./ Denm;
            spcM(:, 2) = (StV * Sh^3) ./ Denm;
            spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
            spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
            spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;

            % Evaluation of the charge balance for the current Sh value, F(Sh)
            F(nn) = Sh + sum(sum(spcM .* chrM));
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
        ipH = 1; Tol = 5.e-15; maxIter = 20;
        % Inicialization of matrix of species
        spcM = zeros(size(chrM)); dspcM = spcM;

        while ipH <= maxIter
            Denm = (1 + Keq(:, 1) / w) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);
            spcM(:, 1) = ((Keq(:, 1) / w) .* StV * Sh^3) ./ Denm;
            spcM(:, 2) = (StV .* Sh^3) ./ Denm;
            spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
            spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
            spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;

            % Evaluation of the charge balance for the current Sh value, F(Sh)
            F = Sh + sum(sum(spcM .* chrM));

            % Calculation of all derivated functions
            dDenm = Denm.^2;
            aux = 3 * Sh^2 * (Keq(:, 1) / w + 1) + 2 * Sh * Keq(:, 2) + Keq(:, 2) .* Keq(:, 3);

            dspcM(:, 1) = (3 * Sh^2 * Keq(:, 1) .* StV) ./ (w * Denm) - ((Keq(:, 1) .* StV * Sh^3) .* aux) ./ (w * dDenm);
            dspcM(:, 2) = (3 * Sh^2 * StV) ./ Denm - (StV * Sh^3 .* aux) ./ dDenm;
            dspcM(:, 3) = (2 * Sh * Keq(:, 2) .* StV) ./ Denm - ((Keq(:, 2) .* StV * Sh^2) .* aux) ./ dDenm;
            dspcM(:, 4) = (Keq(:, 2) .* Keq(:, 3) .* StV) ./ Denm - ((Keq(:, 2) .* Keq(:, 3) .* StV * Sh) .* aux) ./ dDenm;
            dspcM(:, 5) = -(Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4) .* StV .* aux) ./ dDenm;

            % Evaluation of the charge balance for the current Sh value, dF(Sh)
            dF = 1 + sum(sum(dspcM .* chrM));
            %Error
            err = F / dF;
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
                        Sh = (fb * a - fa * b) / (fb - fa);
                        spcM = zeros(size(chrM));
                        Denm = (1 + Keq(:, 1) / w) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);

                        spcM(:, 1) = ((Keq(:, 1) / w) .* StV * Sh^3) ./ Denm;
                        spcM(:, 2) = (StV .* Sh^3) ./ Denm;
                        spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
                        spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
                        spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;

                        fc = Sh + sum(sum(spcM .* chrM));

                        if fa * fc > 0
                            n1 = n1 + 1;

                            if n1 == 2
                                fb = (fc / (fc + fa)) * fb;
                                n1 = 0;
                            end

                            a = Sh; fa = fc;
                        elseif fb * fc > 0 % To avoid problems when fc == 0
                            n2 = n2 + 1;

                            if n2 == 2
                                fa = (fc / (fc + fb)) * fa;
                                n2 = 0;
                            end

                            b = Sh; fb = fc;
                        end

                        err1 = abs(fc);
                        err2 = abs(Sh - (fb * a - fa * b) / (fb - fa));

                        if (err1 < Tol) && (err2 < 1e-14)
                            ipH = maxIter;
                        end

                        ipH = ipH + 1;
                    end

                end

            end

            ipH = ipH + 1;
        end

    end

end

function [mu_max, maint] = determine_max_growth_rate_and_maint(species, T, Sh)
    % Determine the maximum growth rate and maintenance [h-1] for a
    % specific species under certain conditions
    %
    % species: species of bacterium
    % T: temperature (in Kelvin)
    % Sh: concentration of protons
    %
    % -> mu_max: max growth rate
    % -> maint: maintenance requirement

    switch species
        case 1 % AOB
            mu_max = ((1.28*10^(12)*exp(-8183/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
            maint = (1.651*10^(11)*exp(-8183/T))/24;

        case 2 % NOB; Nitrobacter
            mu_max = ((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
            maint = (8.626*10^(6)*exp(-5295/T))/24;

        case 3 % NOB; Nitrospira
            mu_max = 0.63 * ((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
            maint = 0.63 * (8.626*10^(6)*exp(-5295/T))/24;

        case 4 % AMX; Brocadia spp [Puyol 2014]
            mu_max = 1.89*10^(8)*exp(-7330/T);
            maint = 0.05 * mu_max;
            
        otherwise
            error('Bacterial species not implemented: %d', species);
    end
end

function M = calculate_monod(Ks, Ki, conc)
    % Calculate the Monod coefficient given the Ks's, Ki's and the
    % corresponding concentrations
    %
    % Ks: vector with all Ks values (1-by-n)
    % Ki: vector with all Ki values (1-by-n)
    % conc: vector with all concentrations corresponding to the Ks and Ki
    %   values (1-by-n)
    %
    % -> M: Monod coefficient
    
    % apply Ks
    M = prod(conc ./ (conc + Ks));
    
    % apply Ki
    for i = 1:length(Ki)
        if Ki ~= 0
            M = M * Ki(i) / (Ki(i) + conc(i));
        end
    end
end
