# Individual-based Model framework.

From: Gonzalez-Cabaleiro group.

*Contributors: Chiel van Amstel, Eloi Martinez-Rabert, Rebeca Gonzalez-Cabaliero

Individual-based Models (IbM) are characterized for simulating the microbes in an aggregate (granule or biofilm) as a discrete entity with unique traits. Each of the microbe is a small point of reaction that shapes its kinetics function of the local environment and therefore, its growth is affected by the activity of the microbes in the surroundings. Thus, IbM consider the intrinsic interaction between diverse microbes, either a positive (mutualism, syntrophism…) or negative (competition, ammensalism…) interaction. Diffusion of components across the simulation domain is modelled together with the biotic reactions. The substrates and products of this activity are diffusing, and these dynamics are resolved function the boundary conditions imposed over the limits of the domain. This type of models is used to describe microbial growth assuming that cellular division occurs a certain microbe’s age, cycle stage or size. Therefore, division or death occurs independently for each of the cells and only function of local conditions, capturing the heterogeneity of the aggregate and supporting non-linear growth models

## Before having fun



## Granule version

Lorem ipsum

## Suspension version

**Comming soon...**

## References

[1] 

Individual based model on bacterial colonies

Mathemathical model simulating the diffusion of compounds and bacterial growth.

To execute the model (Granule versions):

1. Create/modify Excel with all parameters
2. Create preset-file using the 'create_mat' function
3. Execute the model with a call to 'IbM(sim_number)'
