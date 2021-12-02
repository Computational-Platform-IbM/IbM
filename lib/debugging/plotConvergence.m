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
    
    figure(4); 
    plot((1:iRES)*constants.nDiffusion_per_SScheck, RESvalues(:,1:iRES)'*100, 'LineWidth', 2); 
%     ylim([0, 3]); 
    legend(constants.compoundNames{1:size(RESvalues, 1)}); 
    ylabel('Maximum error [%]'); 
    xlabel('Diffusion iteration')
    title(sprintf('Convergence plot at t=%.1f', Time))
end