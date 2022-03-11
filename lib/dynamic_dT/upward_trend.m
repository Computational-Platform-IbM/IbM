function bool = upward_trend(iDiffusion, iRES, RESvalues, constants)
    % Detect whether there is an upward trend in the RES values
    
%     bool = mod(iDiffusion, constants.dynamicDT.nItersCycle) > 50 && ...
%         max(RESvalues(:, iRES)) - max(RESvalues(:, iRES - 2)) > 0.01;
%     bool = iRES >= 3 && ...
%         any(RESvalues(:, iRES) - RESvalues(:, iRES - 2) > 1e-4);
    bool = iRES >= 3 && ...
        max(RESvalues(:, iRES)) - max(RESvalues(:, iRES - 2)) > 1e-4;

end

