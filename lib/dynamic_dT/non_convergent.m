function no_convergence = non_convergent(iRES, RESvalues, tol)
    % Detect whether the diffusion is no longer converging after a certain
    % number of diffusion iteration in the steady-state cycle

    if iRES > 5
        no_conv_positive = abs(max(RESvalues(:, iRES)) - max(RESvalues(:, iRES-1))) < tol && ...
        abs(max(RESvalues(:, iRES)) - max(RESvalues(:, iRES-2))) < 2*tol;
        no_conv_negative = max(RESvalues(:, iRES)) - max(RESvalues(:, iRES-1)) > 0 && ...
        abs(max(RESvalues(:, iRES)) - max(RESvalues(:, iRES-1))) < 1e-5;
        no_convergence = no_conv_positive || no_conv_negative;
    else
        no_convergence = false;
    end
    
end
