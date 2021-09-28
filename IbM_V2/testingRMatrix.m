% compare parallel reaction matrix calculation with sequential version

load('testingRMatrix_large.mat')

[diffusion_region, focus_region] = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid);
xRange = focus_region.x0:focus_region.x1;
yRange = focus_region.y0:focus_region.y1;

%% parallel compute (4 chunks)
[par_rMatrix, par_mu, par_pH] = par_calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
                grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), ...
                conc(xRange, yRange, :), constants, pH(xRange, yRange), chunks, nChunks_dir);


%% parallel compute (1 chunk) 
chunks = create_chunks(1, focus_region);

[par1_rMatrix, par1_mu, par1_pH] = par_calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
                grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), ...
                conc(xRange, yRange, :), constants, pH(xRange, yRange), chunks, 1);            
            

%% sequential compute            
[res_rMatrix, res_mu, res_pH] = calculate_reaction_matrix(grid2bac(xRange, yRange, :), ...
            grid2nBacs(xRange, yRange), bac, diffusion_region(xRange, yRange, :), conc(xRange, yRange, :), constants, constants.pHsetpoint);            

        
%% compute differences
diff_rMatrix = sqrt(sum((res_rMatrix - par_rMatrix).^2, 'all'));
diff_mu = sqrt(sum((res_mu - par_mu).^2, 'all'));
diff_pH = sqrt(sum((res_pH - par_pH).^2, 'all'));

fprintf('Norm difference 4 chunks vs sequential\n')
fprintf('\trMatrix: %.4g \n\tmu: %.4g \n\tpH: %.4g\n\n', diff_rMatrix, diff_mu, diff_pH)

diff_rMatrix = sqrt(sum((par1_rMatrix - par_rMatrix).^2, 'all'));
diff_mu = sqrt(sum((par1_mu - par_mu).^2, 'all'));
diff_pH = sqrt(sum((par1_pH - par_pH).^2, 'all'));

fprintf('Norm difference 1 chunk vs 4 chunks\n')
fprintf('\trMatrix: %.4g \n\tmu: %.4g \n\tpH: %.4g\n\n', diff_rMatrix, diff_mu, diff_pH)