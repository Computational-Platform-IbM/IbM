function [reaction_matrix, mu, pH_new] = rMatrix_chunk(pH, conc, grid2bac, grid2nBacs, diffRegion, ...
    grouped_bac, nBacs, bacOffset, ...
    reactive_form, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constants, structure_model, pHincluded)
    % Calculate reaction matrix, mu and pH in one part of the simulation
    % domain

    pH_bulk = constants(1);
    constantpH = constants(2);
    pHtolerance = constants(3);
    T = constants(4);
    if structure_model
        iA = constants(5);
        iB = constants(6);
        iC = constants(7);
        iO2 = constants(8);   
    else
        iNH3 = constants(5);
        iNO2 = constants(6);
        iO2 = constants(7);
    end

    bac_species = grouped_bac(:,1);
    bac_molarMass = grouped_bac(:,2);

    mu = zeros(nBacs, 1);
    reaction_matrix = zeros(size(conc));
    pH_new = zeros(size(pH));

    % for each gridcell
    for ix = 1:size(conc, 1)
        for iy = 1:size(conc, 2)
            if ~diffRegion(ix, iy) % bulk layer
                pH_new(ix, iy) = pH_bulk;
            
            else % in diffusion layer, thus pH calculation needs to be performed
                % calculate pH & speciation
                if pHincluded
                    Sh_old = 10^-pH(ix, iy);
                    [spcM, Sh] = solve_pH(Sh_old, [reshape(conc(ix,iy,:), [], 1, 1); 1; 0], Keq, chrM, constantpH, pHtolerance); % <C: why [...; 1; 0]? />
                    pH_new(ix, iy) = -log10(Sh);
                else
                    pH(ix, iy) = pH_bulk;
                    Sh = 10^-pH(ix, iy);
                    spcM = reshape(conc(ix,iy,:), [], 1, 1);
                end

                if grid2nBacs(ix, iy) % if also cells are found here, then update reaction matrix
                    % get bacteria in this grid cell
                    iBacs = reshape(grid2bac(ix, iy, 1:grid2nBacs(ix, iy)), [], 1); % n-by-1 vector of bacterial indices

                    % correct for chunk indexing...
                    iBacs = iBacs - bacOffset;
                    
                    speciesGrid = bac_species(iBacs);
                    unique_species = unique(speciesGrid); % which species
                    for i = 1:length(unique_species)
                        curr_species = unique_species(i);

                        % update mu_max per species based on pH
                        [mu_max, maint] = determine_max_growth_rate_and_maint(curr_species, T, Sh, structure_model);

                        % get reactive concentrations for soluble components
                        if structure_model
                            if pHincluded
                                reactive_conc = [spcM(iA, reactive_form(iA)), ...
                                    spcM(iB, reactive_form(iB)), ...
                                    spcM(iC, reactive_form(iC)), ...
                                    spcM(iO2, reactive_form(iO2))];
                            else
                                reactive_conc = [spcM(iA), ...
                                                spcM(iB), ...
                                                spcM(iC), ...
                                                spcM(iO2)];  
                            end
                        else
                            reactive_conc = [spcM(iNH3, reactive_form(iNH3)), ...
                                spcM(iNO2, reactive_form(iNO2)), ...
                                spcM(iO2, reactive_form(iO2))];
                        end

                        % set mu for all bacteria of same species in that gridcell
                        M = calculate_monod(Ks(curr_species,:), Ki(curr_species, :), reactive_conc);                        % [DN]
                        mu_noMaintenance = mu_max * M;                                                                      % [1/h]
                        mu_withMaintenance = mu_noMaintenance - maint;                                                      % [1/h]
                        mu(iBacs(speciesGrid == curr_species)) = mu_withMaintenance;                                        % [1/h]

                        % update reaction_matrix element for this grid cell
                        concentrationChange = mMetabolism(:, curr_species) * mu_noMaintenance;                              % [molS/molX/h]
                        if mu_withMaintenance < 0
                            concentrationChange = concentrationChange - mDecay(:, curr_species) * mu_withMaintenance;       % [molS/molX/h]
                        end
                        concentrationChange = concentrationChange * sum(bac_molarMass(iBacs(speciesGrid == curr_species))); % [molS/h]

                        reaction_matrix(ix, iy, :) = reaction_matrix(ix, iy, :) + reshape(concentrationChange, 1,1,[]);     % [molS/h]
                    end
                end
            end
        end
    end
end