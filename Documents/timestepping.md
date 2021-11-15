# Timestepping documentation
The size of timesteps, both in diffusion and bacterial activity are taken dynamically in order to speed up the simulation and circumvent infinite oscillations.
The size of these timesteps is dynamically changed based on certain conditions.

## Diffusion timesteps
The timestep for diffusion (`dT`) is dynamically set between a minimum and a maximum which are based on the Neumann stability condition.
By default the maximum is set using a Neumann stability coefficient of 0.5, while the minimum is set with a value of 0.01.
The minimum is only set to prevent infinite decreasing cycles reducing the dT to a point where it will not increase again naturally in a reasonable time.

### Decrease dT for diffusion
There are three situations in which the dT is reduced dynamically:
1) If there is an upward trend in the RES values detected, then the dT is reduced. An increase in the RES values means that the timestep is too large, resulting in diverging oscillations (the RES value is the absolute value).
2) If after a certain number of diffusion iterations (`constants.dynamicDT.iterThresholdDecrease`) the RES values are no longer converging. In this case, the non-convergence is fixed by reducing the dT.
3) If during diffusion a negative concentration is corrected back to 0. When negative concentrations (resulting from the diffusion equation), are corrected back to a value of zero, the RES values will remain the same/high, as from the reaction_matrix different values are expected in the steady state.

### Increase dT for diffusion
A reduced dT is sometimes necessary to get over some stiff situations. However, during most of the simulation the conditions hardly change, thus a larger dT can be tolerated.
In order to increase the dT after stiff situations, the dT is increased under the following condition:
1) If there have been multiple (`constants.dynamicDT.nIterThresholdIncrease`) steady states reached with a high number (`> constants.dynamicDT.iterThresholdIncrease`) of diffusion iterations.

For all increase and decrease actions that are taken the increase/decrease is only called if there has not been a recent dT change. The only exception to this is if there is an upward trend in the RES values, that will always decreases the dT.
This reduces the possibility for infinite loops and ensures that a decrease cannot be followed by an immediate increase (reducing its effectiveness).

## Bacterial activity timesteps
The timestep for bacterial activity (`dT_bac') is dynamically set between a minimum and maximum. These are determined based on the dT_bac that has been set by the user.
The maximum dT_bac is set as the dT_bac that has been set by the user. The minimum dT_bac is 20 times less than the maximum (artificial number). 

### Decrease dT for bacterial activity
The decrease of dT_bac is used as a tool to reduce the change in bulk concentrations between consecutive steady state calculations. 
After calculating the new bulk concentrations, the relative difference between this and the previous bulk concentrations is checked.
If the relative difference is larger than a certain threshold (`constants.dynamicDT.maxRelDiffBulkConc`), the dT_bac is reduced.
This is repeated untill either the minimum dT_bac is reached or the relative concentration change is within the tolerance.
Also, the dT_bac is reduced if in the bacterial division multiple division rounds are required before all bacteria are below the division threshold.

### Increase dT for bacterial activity
Similar to the dT, the dT_bac can be increased in smooth regions. This is determined in a similar fashion as with the dT increase.
If there is multiple steady states (`constants.dynamicDT.nIterThresholdIncrease`) reached with an initial RES value below a certain threshold (`constants.dynamicDT.initRESThresholdIncrease`), the dT_bac in increased.
This is only done if the dT_bac has not recently been changed.
