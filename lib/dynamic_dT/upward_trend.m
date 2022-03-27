function bool = upward_trend(iRES, RESvalues)
    % Detect whether there is an upward trend in the RES values
    
%     bool = mod(iDiffusion, constants.dynamicDT.nItersCycle) > 50 && ...
%         max(RESvalues(:, iRES)) - max(RESvalues(:, iRES - 2)) > 0.01;
%     bool = iRES >= 3 && ...
%         any(RESvalues(:, iRES) - RESvalues(:, iRES - 2) > 1e-4);
%     changed_this_ss = false; %Time.current - Time.changed_dT < 2*Time.dT;
%     if iRES >=3
%         trend_up = max(RESvalues(:, iRES)) - max(RESvalues(:, iRES - 2)) > 1e-4 || any(RESvalues(:, iRES) - RESvalues(:, iRES - ) > max(RESvalues(:, iRES)) / 10);
%     else
%         trend_up = false;
%     end
%     bool = ~changed_this_ss && trend_up;
    bool = iRES > 2 && any(RESvalues(:,iRES) - RESvalues(:,iRES-1) > 0.1);

end

