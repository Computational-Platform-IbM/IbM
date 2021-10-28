function [mu_max, maint] = determine_max_growth_rate_and_maint(species, T, Sh, structure_model)
    % Determine the maximum growth rate and maintenance [h-1] for a
    % specific species under certain conditions
    %
    % species: species of bacterium
    % T: temperature (in Kelvin)
    % Sh: concentration of protons
    %
    % -> mu_max: max growth rate
    % -> maint: maintenance requirement
    
    if structure_model
        mu_max = 1/24;
        maint = 0.25*mu_max;
    else
        switch species
            case 1 % AOB
                mu_max = ((1.28*10^(12)*exp(-8183/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
                maint = (1.651*10^(11)*exp(-8183/T))/24;

            case 2 % NOB; Nitrobacter
                mu_max = ((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
                maint = (8.626*10^(6)*exp(-5295/T))/24;

            case 3 % NOB; Nitrospira
                mu_max = 0.63 * ((6.69*10^(7)*exp(-5295/T))/(1+((2.05*10^(-9))/Sh)+ (Sh/(1.66*10^(-7)))))/24;
                maint = 0.63 * (8.626*10^(6)*exp(-5295/T))/24;

            case 4 % AMX; Brocadia spp [Puyol 2014]
                mu_max = 1.89*10^(8)*exp(-7330/T);
                maint = 0.05 * mu_max;

            otherwise
                error('Bacterial species not implemented: %d', species);
        end
    end
end