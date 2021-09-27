function [spcM, Sh] = solve_pH(Sh_ini, StV, Keq, chrM, keepConstantpH, Tol)
    % Solve the pH and speciation per grid cell using a modified
    % Newton-Raphson algorithm
    %
    % Sh_ini: initial guess for H+ concentration
    % StV: concentration vector
    % Keq: equilibrium constants per compound in the StV vector
    % chrM: charge per subcompound of the StV vector
    % keepConstantpH: (bool) if true, then assume pH is set value,
    %   otherwise, calculate pH from speciation
    % Tol: Newton-Raphson tolerance
    %
    % -> spcM: species matrix for the respective pH
    % -> Sh: proton concentration
    
    if keepConstantpH
        % assume steady pH value, only calculate speciation
        spcM = zeros(size(chrM));
        Sh = Sh_ini;
        Denm = (1 + Keq(:, 1)) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);

        spcM(:, 1) = ((Keq(:, 1)) .* StV * Sh^3) ./ Denm;
        spcM(:, 2) = (StV * Sh^3) ./ Denm;
        spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
        spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
        spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;

    else
        % Newton-Raphson method
        Sh = Sh_ini;
        % Counter of convergences
        ipH = 1; 
        maxIter = 20;
        err = 1; % initial error
        F = 1; % initial electroneutrality error

        % Inicialization of matrix of species
        spcM = zeros(size(chrM)); 
        dspcM = zeros(size(chrM));

        while (abs(err) > Tol) && (abs(F) > Tol) && ipH <= maxIter
            Denm = (1 + Keq(:, 1)) * Sh^3 + Keq(:, 2) * Sh^2 + Keq(:, 3) .* Keq(:, 2) * Sh + Keq(:, 4) .* Keq(:, 3) .* Keq(:, 2);
            spcM(:, 1) = ((Keq(:, 1)) .* StV * Sh^3) ./ Denm;
            spcM(:, 2) = (StV .* Sh^3) ./ Denm;
            spcM(:, 3) = (StV * Sh^2 .* Keq(:, 2)) ./ Denm;
            spcM(:, 4) = (StV * Sh .* Keq(:, 2) .* Keq(:, 3)) ./ Denm;
            spcM(:, 5) = (StV .* Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4)) ./ Denm;

            % Evaluation of the charge balance for the current Sh value, F(Sh)
            F = Sh + sum(spcM .* chrM, 'all');

            % Calculation of all derivated functions
            dDenm = Denm.^2;
            aux = 3 * Sh^2 * (Keq(:, 1) + 1) + 2 * Sh * Keq(:, 2) + Keq(:, 2) .* Keq(:, 3);

            dspcM(:, 1) = (3 * Sh^2 * Keq(:, 1) .* StV) ./ (Denm) - ((Keq(:, 1) .* StV * Sh^3) .* aux) ./ (dDenm);
            dspcM(:, 2) = (3 * Sh^2 * StV) ./ Denm - (StV * Sh^3 .* aux) ./ dDenm;
            dspcM(:, 3) = (2 * Sh * Keq(:, 2) .* StV) ./ Denm - ((Keq(:, 2) .* StV * Sh^2) .* aux) ./ dDenm;
            dspcM(:, 4) = (Keq(:, 2) .* Keq(:, 3) .* StV) ./ Denm - ((Keq(:, 2) .* Keq(:, 3) .* StV * Sh) .* aux) ./ dDenm;
            dspcM(:, 5) = -(Keq(:, 2) .* Keq(:, 3) .* Keq(:, 4) .* StV .* aux) ./ dDenm;

            % Evaluation of the charge balance for the current Sh value, dF(Sh)
            dF = 1 + sum(dspcM .* chrM, 'all');
            %Error
            err = F / dF;
            % Newton-Raphson algorithm
            Sh = max(Sh - err, 1e-10);

            ipH = ipH + 1;
        end
        if any(spcM < 0)
            warning('DEBUG:actionRequired', 'debug: negative concentration encountered after pH calculation...');
        end
        if (Sh < 1e-14) || (Sh > 1)
            warning('DEBUG:actionRequired', 'debug: pH found outside of 1-14 range...');
        end
        if ipH == maxIter
            error('pH solver did not converge to a solution within given number of iterations');
        end
    end
end

