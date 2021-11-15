function no_convergence = non_convergent(iRES, RESvalues, constants)
    % Detect whether the diffusion is no longer converging after a certain
    % number of diffusion iteration in the steady-state cycle

    if iRES > 5
        no_convergence = abs(max(RESvalues(:, iRES)) - max(RESvalues(:, iRES-1))) < constants.dynamicDT.tolerance_no_convergence;
    else
        no_convergence = false;
    end
end
