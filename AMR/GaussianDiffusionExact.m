function conc = GaussianDiffusionExact(X, Y, t)
    % Gaussian diffusion concentration at time t (exact solution)
    % set parameters
    t0 = 1; eta = 1;
    scaling = 1;
    
    % calculate concentration
    r2 = X.^2 + Y.^2;
    conc = scaling / (4 * pi * eta * (t + t0)) * exp(-r2 ./ (4 * eta * (t + t0)));
end