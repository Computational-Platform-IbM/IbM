function conc = diffusion(conc, reaction_matrix, bulk_concentrations, grid, constants, dT)
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
    nCompounds = length(diffusion_rates);
    accuracy = constants.diffusion_accuracy;
    isLiquid = constants.isLiquid;
    
    % set parameters for V-cycle
    iter_pre = 3;
    iter_post = 3;
    iter_final = 4;
    
    % convert concentration to mol/m3
    conc = conc * 1000; 
    
    % stencil initialisation
    alpha = dT * diffusion_rates ./ (2*grid.dx^2);
    alpha = reshape(alpha, 1, 1, nCompounds);
    
    L = [0 1 0; 1 -4 1; 0 1 0]; % laplacian stencil
    I = zeros(3,3);
    I(2,2) = 1;
    L_lhs = I - alpha.*L; % lefthand-side stencil (stacked)
    L_0 = alpha.*L; % basis stencil laplacian (stacked)
    L_rhs = I + alpha.*L; % righthand-side stencil (stacked)
    
    base = [1 2 1; 2 4 2; 1 2 1];
    L_restriction = base/16;
    L_prolongation = base/4;
    
    
    % create right hand side
    % - boundary conditions
    rhs = calculate_rhs_dirichlet(conc, L_rhs, bulk_concentrations); % (nX, nY, nCompounds)
    
    % - reaction matrix
    rhs = rhs + dT*1000*reaction_matrix;
    
    % solve using multigrid
    while sum(residual(conc(:,:,isLiquid), rhs, L_lhs).^2, 'all') > accuracy^2 % absolute norm of residual > accuracy
        conc(:,:,isLiquid) = V_Cycle(conc(:,:,isLiquid), rhs, L_0, L_restriction, L_prolongation, 9, 0, iter_pre, iter_post, iter_final);
    end

    % convert concentrations back to mol/L
    conc = conc / 1000;
end


function rhs = calculate_rhs_dirichlet(phi, L_rhs, values)
    % Calculate (part of) the RHS of the diffusion equation using
    % convolution.
    %
    % phi: matrix with starting values (:, :, nCompounds)
    % L_rhs: kernel used for the rhs convolution
    % values: vector with boundary value for each compound
    %
    % rhs: right-hand-side of the diffusion equation due to diffusion
    
    phi = create_dirichlet_boundary(phi, values);
    rhs = convn(phi, L_rhs, 'valid');
end

function x_dirichlet = create_dirichlet_boundary(x, values)
    % Construct dirichlet boundary around the given matrix x
    %
    % x: matrix of (n, n, nCompounds) around which to construct boundary
    % values: vector with boundary value for each compound
    %
    % x_dirichlet: padded matrix x with the boundary values
    
    x_dirichlet = zeros(size(x, [1,2]) + 2, size(x, 3));
    x_dirichlet(2:end-1, 2:end-1, :) = x;

    mask = zeros(size(x, [1,2])+2, 'logical');
    mask([1, end], :) = 1;
    mask(:, [1, end]) = 1;
    
    x_dirichlet = set_concentrations(x_dirichlet, values, mask);
end


%{

conc_og = reshape(grid.cells.conc(1:nCells), grid.nX, grid.nX);
rhs = reshape(rhs, grid.nX, grid.nX);

%% time MG cpu V-cycle & F-cycle

% V-cycle presets
if randomnessFactor
    iter_pre = 7;
    iter_post = 7;
    iter_final = 10;
else
    iter_pre = 1;
    iter_post = 2;
    iter_final = 5;
end

for i = 1:n_timing
    conc = conc_og;
    tic
%     while norm(residual(conc, rhs, L_lhs)) > accuracy % absolute residual
    while sum(residual(conc, rhs, L_lhs).^2, 'all') > accuracy^2 % absolute residual
        conc = V_Cycle(conc, rhs, L_0, L_restriction, L_prolongation, 5, 0, iter_pre, iter_post, iter_final);
    end
    timings(2, i) = toc;
end
fprintf('norm V-cycle: %e\n', norm(residual(conc, rhs, L_lhs)));
fprintf('norm2 V-cycle: %e\n\n', sqrt(sum(residual(conc, rhs, L_lhs).^2, 'all')));

% F-cycle presets
if randomnessFactor
    iter_pre = 4;
    iter_resmoothing = 5;
    iter_post = 5;
    iter_final = 10;
else
    iter_pre = 1;
    iter_resmoothing = 1;
    iter_post = 1;
    iter_final = 5;
end

for i = 1:n_timing
    conc = conc_og;
    tic
%     while norm(residual(conc, rhs, L_lhs)) > accuracy % absolute residual
    while sum(residual(conc, rhs, L_lhs).^2, 'all') > accuracy^2 % absolute residual
        conc = F_Cycle(conc, rhs, L_0, L_restriction, L_prolongation, 5, 0, iter_pre, iter_resmoothing, iter_post, iter_final);
    end
    timings(3, i) = toc;
end
fprintf('norm F-cycle: %e\n', norm(residual(conc, rhs, L_lhs)));
fprintf('norm2 F-cycle: %e\n\n', sqrt(sum(residual(conc, rhs, L_lhs).^2, 'all')));




%% time smoothing

for i = 1:n_timing
    conc = conc_og;
    tic
%     while norm(residual(conc, rhs, L_lhs)) > accuracy % absolute residual
    while sum(residual(conc, rhs, L_lhs).^2, 'all') > accuracy^2 % absolute residual
        conc = smoothing(conc, rhs, L_lhs);
    end
    conc_final = gather(conc);
    timings(4, i) = toc;
end
fprintf('norm smoothing cpu: %e\n', norm(residual(conc, rhs, L_lhs)));
fprintf('norm2 smoothing cpu: %e\n\n', sqrt(sum(residual(conc, rhs, L_lhs).^2, 'all')));


%% time gpu_smoothing

rhs_gpu = gpuArray(rhs);
for i = 1:n_timing
    conc = gpuArray(conc_og);
    tic
%     while norm(residual(conc, rhs_gpu, L_lhs)) > accuracy
    while sum(residual(conc, rhs, L_lhs).^2, 'all') > accuracy^2 % absolute residual
        conc = smoothing(conc, rhs_gpu, L_lhs);
    end
    conc_final = gather(conc);
    timings(5, i) = toc;
end
fprintf('norm smoothing gpu: %e\n', norm(residual(conc, rhs, L_lhs)));
fprintf('norm2 smoothing gpu: %e\n\n', sqrt(sum(residual(conc, rhs, L_lhs).^2, 'all')));


%% time gpu_MG V-cycle & F-cycle

% V-cycle presets
if randomnessFactor
    iter_pre = 7;
    iter_post = 7;
    iter_final = 10;
else
    iter_pre = 1;
    iter_post = 2;
    iter_final = 5;
end

for i = 1:n_timing
    conc = gpuArray(conc_og);
    tic
%     while norm(residual(conc, rhs, L_lhs)) > accuracy % absolute residual
    while sum(residual(conc, rhs, L_lhs).^2, 'all') > accuracy^2 % absolute residual
        conc = V_Cycle_gpu(conc, rhs, L_0, L_restriction, L_prolongation, 5, 0, iter_pre, iter_post, iter_final);
    end
    conc_final = gather(conc);
    timings(6, i) = toc;
end
fprintf('norm V-cycle gpu: %e\n', norm(residual(conc, rhs, L_lhs)));
fprintf('norm2 V-cycle gpu: %e\n\n', sqrt(sum(residual(conc, rhs, L_lhs).^2, 'all')));


% F-cycle presets
if randomnessFactor
    iter_pre = 4;
    iter_resmoothing = 5;
    iter_post = 5;
    iter_final = 10;
else
    iter_pre = 1;
    iter_resmoothing = 1;
    iter_post = 1;
    iter_final = 5;
end

for i = 1:n_timing
    conc = gpuArray(conc_og);
    tic
%     while norm(residual(conc, rhs, L_lhs)) > accuracy % absolute residual
    while sum(residual(conc, rhs, L_lhs).^2, 'all') > accuracy^2 % absolute residual
        conc = F_Cycle_gpu(conc, rhs, L_0, L_restriction, L_prolongation, 5, 0, iter_pre, iter_resmoothing, iter_post, iter_final);
    end
    conc_final = gather(conc);
    timings(7, i) = toc;
end
fprintf('norm F-cycle gpu: %e\n', norm(residual(conc, rhs, L_lhs)));
fprintf('norm2 F-cycle gpu: %e\n\n', sqrt(sum(residual(conc, rhs, L_lhs).^2, 'all')));


% %%
% 
% conc = gpuArray(conc_og);
% fprintf('\ngpu euclidean norm\n');
% tic
% for i = 1:1000
%     sqrt(sum(residual(conc, rhs, L_lhs).^2, 'all')) > accuracy;
% end
% toc
% 
% fprintf('gpu sum(x^2) > acc^2\n');
% tic
% for i = 1:1000
%     sum(residual(conc, rhs, L_lhs).^2, 'all') > accuracy^2; % winner!!!
% end
% toc
% 
% fprintf('\n');




%% plot timings

names = {'mldivide', 'MG: V-cycle', 'MG: F-cycle', 'Jacobi smoothing', 'Jacobi smoothing (gpu)', 'MG: V-cycle (gpu)', 'MG: F-cycle (gpu)'};
means = mean(timings, 2);
err = std(timings, 0, 2);

[means_sorted, sortIndex] = sort(means, 'descend');
names_cat = categorical(names(sortIndex));
names_cat = reordercats(names_cat, names(sortIndex));

figure(10); clf;
bar(names_cat, means_sorted, 0.7); hold on;
er = errorbar(names_cat, means_sorted, err(sortIndex));
er.Color = [0, 0, 0];
er.LineStyle = 'none'; hold off;

title([{'--- Dirichlet boundary (start w/ central concentration) ---'}, {'|R| \leq 5e-3, accuracy \leq 1e-7,'}, {'nx=ny=1025 (n=10), Neumann=0.2'}])
ylabel('Time per iteration \pm 1 std [s]')


%% F_cycle function

function phi = F_Cycle(phi, f, L_0, L_restriction, L_prolongation, maxDepth, depth, iter_pre, iter_resmoothing, iter_post, iter_final)
	% Recursive F-cycle multigrid for solving the Poisson equation (\nabla^2 phi = f) on a uniform grid of spacing h
    L_lhs = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*depth))*L_0;
    
    % Pre-smoothing
    for i = 1:iter_pre
        phi = smoothing(phi,f,L_lhs);
    end
	
	% Compute Residual Errors
    r = residual(phi,f,L_lhs);
    
	% Restriction
	rhs = restriction(r, L_restriction);

	eps = zeros(size(rhs));

	% stop recursion at smallest grid size, otherwise continue recursion
	if ceil(sqrt(numel(phi))) <= maxDepth
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*(depth+1)))*L_0;
        for i = 1:iter_final
            eps = smoothing(eps,rhs,L_lhs_deeper);
        end
    else
        eps = F_Cycle(eps, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, iter_pre, iter_resmoothing, iter_post, iter_final);        
    end
	
	% Prolongation and Correction
	phi = phi + prolongation(eps, L_prolongation, size(phi));
    
	% Re-smoothing
    for i=1:iter_resmoothing
        phi = smoothing(phi,f,L_lhs);
    end

	% Compute residual errors
	r = residual(phi,f,L_lhs);
	
	% Restriction
	rhs = restriction(r, L_restriction);

	% stop recursion at smallest grid size, otherwise continue recursion
	if ceil(sqrt(numel(phi))) <= maxDepth
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*(depth+1)))*L_0;
        for i = 1:iter_final
            eps = smoothing(eps,rhs,L_lhs_deeper);
        end
    else
        eps = V_Cycle(eps, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, 5, 5, 5); % iter_pre, iter_post, iter_final);        
    end
    
	% Prolongation and Correction
	phi = phi + prolongation(eps, L_prolongation, size(phi));
    
	% Post-smoothing
    for i = 1:iter_post
        phi = smoothing(phi,f,L_lhs);
    end
end

function phi = F_Cycle_gpu(phi, f, L_0, L_restriction, L_prolongation, maxDepth, depth, iter_pre, iter_resmoothing, iter_post, iter_final)
	% Recursive F-cycle multigrid for solving the Poisson equation (\nabla^2 phi = f) on a uniform grid of spacing h
    L_lhs = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*depth))*L_0;
    
    % Pre-smoothing
    for i = 1:iter_pre
        phi = smoothing(phi,f,L_lhs);
    end
	
	% Compute Residual Errors
    r = residual(phi,f,L_lhs);
    
	% Restriction
	rhs = restriction(r, L_restriction);

	eps = zeros(size(rhs), 'gpuArray');

	% stop recursion at smallest grid size, otherwise continue recursion
	if ceil(sqrt(numel(phi))) <= maxDepth
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*(depth+1)))*L_0;
        for i = 1:iter_final
            eps = smoothing(eps,rhs,L_lhs_deeper);
        end
    else
        eps = F_Cycle_gpu(eps, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, iter_pre, iter_resmoothing, iter_post, iter_final);        
    end
	
	% Prolongation and Correction
	phi = phi + prolongation_gpu(eps, L_prolongation, size(phi));
    
	% Re-smoothing
    for i=1:iter_resmoothing
        phi = smoothing(phi,f,L_lhs);
    end

	% Compute residual errors
	r = residual(phi,f,L_lhs);
	
	% Restriction
	rhs = restriction(r, L_restriction);

	% stop recursion at smallest grid size, otherwise continue recursion
	if ceil(sqrt(numel(phi))) <= maxDepth
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*(depth+1)))*L_0;
        for i = 1:iter_final
            eps = smoothing(eps,rhs,L_lhs_deeper);
        end
    else
        eps = V_Cycle_gpu(eps, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, 5, 5, 5); % iter_pre, iter_post, iter_final);        
    end
    
	% Prolongation and Correction
	phi = phi + prolongation_gpu(eps, L_prolongation, size(phi));
    
	% Post-smoothing
    for i = 1:iter_post
        phi = smoothing(phi,f,L_lhs);
    end
end

%% V_cycle function

function phi = V_Cycle(phi, f, L_0, L_restriction, L_prolongation, maxDepth, depth, iter_pre, iter_post, iter_final)
	% Recursive V-Cycle Multigrid for solving the Diffusion equation on a uniform grid

    % Create correct left-hand side stencil
    L_lhs = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*depth))*L_0;
    
	% Pre-Smoothing
    for i = 1:iter_pre
        phi = smoothing(phi,f,L_lhs);
    end
	
	% Compute Residual Errors
    r = residual(phi,f,L_lhs);
    
	% Restriction
	rhs = restriction(r, L_restriction);

	eps = zeros(size(rhs));

	% stop recursion at smallest grid size, otherwise continue recursion
	if ceil(sqrt(numel(phi))) <= maxDepth
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*(depth+1)))*L_0;
        for i = 1:iter_final
            eps = smoothing(eps,rhs,L_lhs_deeper);
        end
    else
        eps = V_Cycle(eps, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, iter_pre, iter_post, iter_final);        
    end
    
	% Prolongation and Correction
	phi = phi + prolongation(eps, L_prolongation, size(phi));
    
	% Post-Smoothing
    for i = 1:iter_post
        phi = smoothing(phi,f,L_lhs);
    end
end

function phi = V_Cycle_gpu(phi, f, L_0, L_restriction, L_prolongation, maxDepth, depth, iter_pre, iter_post, iter_final)
	% Recursive V-Cycle Multigrid for solving the Diffusion equation on a uniform grid

    % Create correct left-hand side stencil
    L_lhs = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*depth))*L_0;
    
	% Pre-Smoothing
    for i = 1:iter_pre
        phi = smoothing(phi,f,L_lhs);
    end
	
	% Compute Residual Errors
    r = residual(phi,f,L_lhs);
    
	% Restriction
	rhs = restriction(r, L_restriction);

	eps = zeros(size(rhs), 'gpuArray');

	% stop recursion at smallest grid size, otherwise continue recursion
	if ceil(sqrt(numel(phi))) <= maxDepth
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] - (1/2^(2*(depth+1)))*L_0;
        for i = 1:iter_final
            eps = smoothing(eps,rhs,L_lhs_deeper);
        end
    else
        eps = V_Cycle_gpu(eps, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, iter_pre, iter_post, iter_final);        
    end
    
	% Prolongation and Correction
	phi = phi + prolongation_gpu(eps, L_prolongation, size(phi));
    
	% Post-Smoothing
    for i = 1:iter_post
        phi = smoothing(phi,f,L_lhs);
    end
end
%}

%{
%% convenience functions

function phi = smoothing(phi, rhs, L_lhs)
    % apply a Jacobi iteration on the system using convolution
    L_sm = L_lhs;
    L_sm(2,2) = 0;
    phi = (rhs - convn(phi, L_sm, 'same'))/L_lhs(2,2);
end

function r = residual(phi, rhs, L_lhs)
    % calculate the residual, given a concentration matrix and a lhs
    % stencil.
    r = rhs - convn(phi, L_lhs, 'same');
end

function rhs = restriction(r, L_restriction)
    % creates new rhs matrix for a coarser mesh, based on the residuals.
    rhs = convn(r, L_restriction, 'same');
    rhs = rhs(1:2:end, 1:2:end); % *4?
end

function phi_fine = prolongation(phi_coarse, L_prolongation, sz)
    % linearly interpolates the coarse values to the fine values
    phi_fine = zeros(sz);
    phi_fine(1:2:end, 1:2:end) = phi_coarse;
    phi_fine = convn(phi_fine, L_prolongation, 'same');
end

function phi_fine = prolongation_gpu(phi_coarse, L_prolongation, sz)
    % linearly interpolates the coarse values to the fine values
    phi_fine = zeros(sz, 'gpuArray');
    phi_fine(1:2:end, 1:2:end) = phi_coarse;
    phi_fine = convn(phi_fine, L_prolongation, 'same');
end



%}
