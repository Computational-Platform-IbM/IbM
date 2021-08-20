function conc = diffusion(conc, reaction_matrix, bulk_concentrations, grid, constants)
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
    % grid: struct containing all information regarding the grid
    % constants: struct containing all simulation constants
    % dT: time over which to solve the diffusion equations
    %
    % -> conc: concentrations after solving the diffusion equations [mol/L]
    
    % variable declarations/unpacking
    diffusion_rates = constants.diffusion_rates;
    accuracy = constants.diffusion_accuracy;
    absolute_tolerance = constants.Tol_a;
    nCompounds = length(diffusion_rates);
    dT = constants.dT;
    
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
        alpha = constants.dT * diffusion_rates(iCompound) / (2*grid.dx^2);
        L_lhs = I - alpha*L;                    % lefthand-side stencil 
        L_0 = alpha*L;                          % basis stencil laplacian
        L_rhs = I + alpha*L;                    % righthand-side stencil

        % create right hand side
        % - boundary conditions
        rhs = calculate_rhs_dirichlet(conc(:,:,iCompound), L_rhs, bulk_concentrations(iCompound)); % (nX, nY, nCompounds)

        % - reaction matrix
        rhs = rhs + dT*1000*reaction_matrix(:,:,iCompound);

        % solve using multigrid
        while sum(residual(conc(:,:,iCompound), rhs, L_lhs).^2, 'all') > accuracy^2 % absolute norm of residual > accuracy
            conc(:,:,iCompound) = V_Cycle(conc(:,:,iCompound), rhs, L_0, L_restriction, L_prolongation, 9, 0, iter_pre, iter_post, iter_final);
        end
    end
    
    % apply correction for negative values
    negative_concentration = conc < 0;
    if any(negative_concentration)
        warning('DEBUG:noActionRequired', 'debug: negative concentration encountered in diffusion solution... correction applied')
        conc = ~negative_concentration.*conc + negative_concentration*absolute_tolerance; % set negative concentrations to very small number, not to 0 because of divide-by-0 in other parts of the code
    end    
    % convert concentrations back to mol/L
    conc = conc / 1000;
end


function rhs = calculate_rhs_dirichlet(phi, L_rhs, value)
    % Calculate (part of) the RHS of the diffusion equation using
    % convolution for 1 compound only.
    %
    % phi: matrix with starting values (:, :)
    % L_rhs: kernel used for the rhs convolution
    % value: boundary value for this compound
    %
    % rhs: right-hand-side of the diffusion equation due to diffusion
    
    phi = create_dirichlet_boundary(phi, value);
    rhs = convn(phi, L_rhs, 'valid');
end


