function bool = bulk_conc_diff_within_limit(new_bulk_concs, bulk_concs, constants)
    % Determine whether all bulk concentrations only differ at most the
    % threshold with the previous bulk concentrations

    bool = all(abs(new_bulk_concs - bulk_concs) ./ (bulk_concs + 1e-20) < constants.dynamicDT.maxRelDiffBulkConc);
end