function conc = diffusion(conc, reaction_matrix, bulk_concentrations, diffRegion, grid, constants)
    % Solve diffusion for all molecules in the liquid phase using the 
    % multigrid method. IMPORTANT: only runs for all dirichlet conditions 
    % as of now. Future versions should include variable conditions per 
    % boundary.
    %
    % conc: concentration of each molecule in the grid [mol/L]
    %   (ix, iy, compound)
    % reaction_matrix: matrix with per grid cell and per compound the
    %   change [h-1] due to bacterial activity
    % bulk_concentrations: vector with the bulk concentration per compound
    % diffRegion: matrix with per gridcell whether the cell is in the
    %   diffusion region
    % grid: struct containing all information regarding the grid
    % constants: struct containing all simulation constants
    % dT: time over which to solve the diffusion equations
    %
    % -> conc: concentrations after solving the diffusion equations [mol/L]
    
    % variable declarations/unpacking
    diffusion_coef = constants.diffusion_rates; % [m2/h] <E: diffusion_rates -> diffusion_coef. /> 
    accuracy = constants.diffusion_accuracy;
    absolute_tolerance = constants.Tol_a;
    nCompounds = length(diffusion_coef);
    dT = constants.dT;
    dx = grid.dx;

    % set parameters for V-cycle 
    % <TODO: optimize under realistic conditions/>
    iter_pre = 3;
    iter_post = 3;
    iter_final = 4;
    
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
            residual_diffRegion = residual_diffRegion(diffRegion);
            isSolution = sum(residual_diffRegion.^2, 'all') < accuracy^2;
        end
        
        % apply correction for negative values
        negative_concentration = conc(:,:,iCompound) < 0;
        if any(negative_concentration, 'all')
            warning('DEBUG:noActionRequired', 'debug: negative concentration encountered in diffusion solution of compound %s... correction applied', constants.StNames{iCompound})
            conc(:,:,iCompound) = ~negative_concentration.*conc(:,:,iCompound) + negative_concentration*absolute_tolerance; % set negative concentrations to very small number, not to 0 because of divide-by-0 in other parts of the code
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
    % '2*value' in create_dirichlet_boundary arguments is from
    % Crank-Nicholson method (Dirichlet Boundary of both t=t and t=t+dt).
    %
    % -> rhs: right-hand-side of the diffusion equation due to diffusion,
    %   valid only in the diffusion region. Outside an artificial value of
    %   <bulk_concentration> is set.
    
    rhs_diffRegion = convn(diffRegion.*phi + ~diffRegion*value, L_rhs, 'same'); % The factor 2 is not required, because we already have the correct conc at t=t, and t=t+dt?
    rhs = diffRegion .* rhs_diffRegion + ~diffRegion .* ones(size(phi)) * value;
    
%     phi = create_dirichlet_boundary(phi, 2*value);
%     rhs = convn(phi, L_rhs, 'valid');
end


