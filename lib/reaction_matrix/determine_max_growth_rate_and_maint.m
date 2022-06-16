function [mu_max, maint] = determine_max_growth_rate_and_maint(species, T, Sh)
    % Determine the maximum growth rate and maintenance [h-1] for a
    % specific species under certain conditions
    %
    % species: species of bacterium
    % T: temperature (in Kelvin)
    % Sh: concentration of protons [10^(-pH)]
    %
    % -> mu_max: max growth rate
    % -> maint: maintenance requirement

    switch species
        case 1 % AOB [1/h]
            mu_max = (5.333*10^(10)*exp(-8183/T))/(1+((10^(-8.688))/Sh)+ (Sh/(10^(-6.780))));
            maint = 6.879*10^(9)*exp(-8183/T);

        case 2 % NOB; Nitrobacter [1/h]
            mu_max = (2.788*10^(6)*exp(-5295/T))/(1+((10^(-8.688))/Sh)+ (Sh/(10^(-6.780))));
            maint = 3.594*10^(5)*exp(-5295/T);

        case 3 % NOB; Nitrospira [1/h]
            mu_max = (1.756*10^(6)*exp(-5295/T))/(1+((10^(-8.688))/Sh)+ (Sh/(10^(-6.780))));
            maint = 2.264*10^(5)*exp(-5295/T);
            
        case 4 % NOB; Nitrotoga [1/h]
            mu_max = (0.025*((305.6-T)/4.67)*(T/300.93)^(64.439))/(1+((10^(-8.688))/Sh)+ (Sh/(10^(-6.780))));
            maint = 3.225*10^(-3)*((305.6-T)/4.67)*(T/300.93)^(64.439);
            
        case 5 % AMX; Brocadia spp. [1/h]
            mu_max = 1.888*10^(8)*exp(-7330/T);
            maint = 0.05 * mu_max;

        otherwise
            error('Bacterial species not implemented: %d', species);
    end
end