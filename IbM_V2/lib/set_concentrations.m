function conc = set_concentrations(conc_old, set_concs, mask)
    % Set the initial concentration in the bioaggregate
    %
    % conc: concentration of each compound per grid cell (nX, nY, compound)
    % init_concs: vector of initial concentration per compound in the 
    %   bioaggregate (nCompounds-by-1)
    % mask: logical matrix (nX, nY) with cells marked for which
    %   concentration needs to be set
    %
    % -> conc: concentration of each compound per grid cell (nX, nY, comp)
    
    set_concs = reshape(set_concs, 1,1,[]);
    conc = ~mask .* conc_old + mask .* set_concs;
end