function conc = diffusionMG(conc, reaction_matrix, bulk_concentrations, diffRegion, grid, constants, Time)
    % Solve diffusion for all molecules in the liquid phase using the 
    % multigrid method. IMPORTANT: only runs for all dirichlet conditions 
    % as of now. Future versions should include variable conditions per 
    % boundary.
    %
    % conc: concentration of each molecule in the grid [mol/L]
    %   (ix, iy, compound)
    % reaction_matrix: matrix with per grid cell and per compound the
    %   change [mol/L/h] due to bacterial activity
    % bulk_concentrations: vector with the bulk concentration per compound
    % diffRegion: matrix with per gridcell whether the cell is in the
    %   diffusion region
    % grid: struct containing all information regarding the grid
    % constants: struct containing all simulation constants
    % dT: time over which to solve the diffusion equations
    %
    % -> conc: concentrations after solving the diffusion equations [mol/L]
    
    % variable declarations/unpacking
    diffusion_coef = constants.diffusion_rates; % [m2/h]
    accuracy = constants.diffusion_accuracy;
    nCompounds = length(diffusion_coef);
    dx = grid.dx;
    dT = Time.dT;

    % set parameters for V-cycle 
    iter_pre = 6;
    iter_post = 7;
    iter_final = 5;
    
    % convert concentration to mol/m3
    conc = conc * 1000; 
    
    % stencil initialisation
    L = [0 1 0; 1 -4 1; 0 1 0];                 % laplacian stencil
    I = zeros(3,3);
    I(2,2) = 1;
    
    base = [1 2 1; 2 4 2; 1 2 1];
    L_restriction = base/16;
    L_prolongation = base/4;
    
    for iCompound = 1:nCompounds                % parfor?
        % stencil updates/declarations
        alpha = dT * diffusion_coef(iCompound) / (2*dx^2);
        L_lhs = I - alpha*L;                    % lefthand-side stencil 
        L_0 = alpha*L;                          % basis stencil laplacian
        L_rhs = I + alpha*L;                    % righthand-side stencil

        % create right hand side
        % - boundary conditions
        rhs_bc = calculate_rhs_dirichlet(conc(:,:,iCompound), L_rhs, bulk_concentrations(iCompound)*1000, diffRegion); % (nX, nY, nCompounds)

        % - reaction matrix
        rhs_react = dT*1000*reaction_matrix(:,:,iCompound);
        
        rhs = rhs_bc + rhs_react;

        % solve using multigrid
        isSolution = false;             % without running, no solution yet
        while ~isSolution
            conc(:,:,iCompound) = V_Cycle(conc(:,:,iCompound), diffRegion, bulk_concentrations(iCompound)*1000, rhs, L_0, L_restriction, L_prolongation, 9, 0, iter_pre, iter_post, iter_final);
            residual_diffRegion = residual(conc(:,:,iCompound), rhs, L_lhs);
            residual_diffRegion = diffRegion.*residual_diffRegion;
            isSolution = sum(residual_diffRegion.^2, 'all') < accuracy^2;
        end
        
        % check for negative values, raise error
        negative_concentration = conc(:,:,iCompound) < 0;
        if any(negative_concentration, 'all')
            if Time.dT ~= Time.minDT && abs(min(conc(negative_concentration))) > accuracy^2 / 100 % dT can be reduced and significant negative value
                error('Diffusion:NegativeConcentration', 'Negative concentration encountered in diffusion solution of compound %s', constants.compoundNames{iCompound})
            else % dT cannot be reduced, thus return corrected concentration or insignificant negative value
                temp = conc(:,:,iCompound);
                warning('Diffusion:NegativeConcentration', 'Negative concentration encountered in diffusion solution of compound %s, but cannot correct dT value thus corrected %d value(s) (smallest number %g) to 0', constants.compoundNames{iCompound}, sum(negative_concentration, 'all'), min(temp(negative_concentration)))
                conc(:,:,iCompound) = (conc(:,:,iCompound) > 0) .* conc(:,:,iCompound);
            end
        end
    end
    
    % convert concentrations back to mol/L
    conc = conc / 1000;
end


function rhs = calculate_rhs_dirichlet(phi, L_rhs, value, diffRegion)
    % Calculate (part of) the RHS of the diffusion equation using
    % convolution for 1 compound only, and only in the diffusion region.
    % Assumes that there is at least one layer of bulk liquid around the
    % diffusion region.
    %
    % phi: matrix with starting values (:, :)
    % L_rhs: kernel used for the rhs convolution
    % value: boundary value for this compound
    % diffRegion: matrix with per gridcell whether it is in the diffusion
    %   region
    % 
    % Only required if else the concentration would be 0 (outside of
    % domain)
    %
    % -> rhs: right-hand-side of the diffusion equation due to diffusion,
    %   valid only in the diffusion region. Outside an artificial value of
    %   <bulk_concentration> is set.
    
    rhs_diffRegion = convn(diffRegion.*phi + ~diffRegion*value, L_rhs, 'same');
    rhs = diffRegion .* rhs_diffRegion + ~diffRegion .* ones(size(phi)) * value;
end


