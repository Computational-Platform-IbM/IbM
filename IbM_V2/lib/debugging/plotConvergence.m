function plotConvergence(RESvalues, iRES, constants, Time)
    % debug visualisation of convergence to steady state
    %
    % RESvalues: array (nCompounds-by-n) with maximum RES value per
    %   compound
    % iRES: index for which to plot the convergence
    % constants: struct with constants of the system
    % Time: struct with information regarding the simulation times
    %
    % -> plot with maximum RES value per compound over the different
    % diffusion cycles until it is declared steady state
    
    figure('Name', 'Convergence of RES values'); 
    plot((1:iRES-1)*constants.nDiffusion_per_SScheck, RESvalues(:,2:iRES)'*100, 'LineWidth', 2); 
%     ylim([0, 3]); 
    legend(constants.StNames{1:7}); 
    ylabel('Maximum error [%]'); 
    xlabel('Diffusion iteration')
    title(sprintf('Convergence plot at t=%.1f', Time))
end