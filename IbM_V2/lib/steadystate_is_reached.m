function [isReached, max_RES_value] = steadystate_is_reached(conc, reaction_matrix, dx, bulk_concentrations, diffRegion, constants)
    % Check whether a steady state is reached during the diffusion
    % 
    % conc0: last concentration values
    % conc: latest concentration values
    %
    % -> isReached: boolean signifying whether steady state is reached

    % unpack/declare variables
    correction_concentration_steadystate = constants.correction_concentration_steady_state; % [mol/L]
    steadystate_tolerance = constants.steadystate_tolerance;
    method = constants.RESmethod; % {'mean', 'max', 'norm'}
    
    L = [0 1 0; 1 -4 1; 0 1 0];  % 2D laplacian stencil base
    nCompounds = sum(constants.isLiquid);
    characteristic_time = dx^2 ./ constants.diffusion_rates;
    compound_steadystate = zeros(nCompounds, 1, 'logical');
    
    
    % DEBUG
    max_RES_value = zeros(nCompounds, 1);
    % END DEBUG
    
    
    % per compound, calculate delta-concentration
    for iCompound = 1:nCompounds % parfor?
        padded_conc = create_dirichlet_boundary(conc(:,:,iCompound), bulk_concentrations(iCompound));
        delta_conc = convn(padded_conc, L, 'valid');
        delta_conc = delta_conc + characteristic_time(iCompound) * reaction_matrix(:,:,iCompound);
        RES = delta_conc ./ (correction_concentration_steadystate + conc(:,:,iCompound)); % ==> RES [mol/L] ./ ([mol/L] + [mol/L])
        RES = RES(diffRegion); % only look at RES values in diffusion region, edge in bulk-layer is by definition not zero in SS
        compound_steadystate(iCompound) = isReached_compound(RES, method, steadystate_tolerance);
        
        % DEBUG
%         max_RES_value(iCompound) = max(abs(RES), [], 'all');
%         if ~compound_steadystate(iCompound)
%             % plot RES over total domain
%             figure(iCompound + 3)
%             heatmap(RES');
%             max_val = max(abs(caxis()));
%             caxis([-max_val, max_val])
%             colormap(redblue());
%             title(constants.StNames(iCompound))
%             fprintf('%s not in steady state yet\n', constants.StNames{iCompound});
%         end
        % END DEBUG
        
    end
    
    isReached = all(compound_steadystate);
end

function SSreached = isReached_compound(RES, method, steadystate_tolerance)
    % Determine for one compound whether the steady state is reached
    %
    % RES: matrix with per grid cell the delta-concentration, corrected
    %   as per <conc / (correction + conc)>. 
    % method: steady state method applied on the entire grid.
    % steadystate_tolerance: tolerance when it comes to steadystate
    %
    % -> SSreached: boolean whether steady state is reached for this compound
    
    switch method
        case 'mean'
            %{
            C: When starting with small granules in large domain, this
            method doesn't make any sense any more, because by definition
            bulk layer has 0 RES. And in this model, the domain isn't
            shrunk to only the diffusion region, thus artificially has a
            low RES when using the mean. Instead now taking the mean over
            all non-zero entries.
            %}
            SSdif = mean(abs(nonzeros(RES)), 'all');
            if isnan(SSdif) && all(RES == 0, 'all') % if delta-concentration == 0 everywhere, mean returns NaN, but still in steady state...
                SSdif = 0;
            end
        case 'max'
            SSdif = max(abs(RES), [], 'all');
        case 'norm'
            SSdif = sqrt(sum(RES.^2, 'all'));
            %{
            C: in testing this was too harsh of a constraint, because all
            little errors accumulate and make it so that it takes
            (figuratively) a million iterations before steady state...
            %}
            
        otherwise
            error(['RES method <', method, '> is not a valid method.'])
    end
    
    SSreached = SSdif <= steadystate_tolerance;
end



