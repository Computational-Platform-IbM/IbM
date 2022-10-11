# Individual-based Model framework.

From: Gonzalez-Cabaleiro group.

*Contributors: Chiel van Amstel, Eloi Martinez-Rabert, Rebeca Gonzalez-Cabaliero

Individual-based Models (IbM) are characterized for simulating the microbes in an 
aggregate (granule or biofilm) as a discrete entity with unique traits. Each of the microbe is 
a small point of reaction that shapes its kinetics function of the local environment and 
therefore, its growth is affected by the activity of the microbes in the surroundings. Thus, 
IbM consider the intrinsic interaction between diverse microbes, either a positive 
(mutualism, syntrophism…) or negative (competition, ammensalism…) interaction. Diffusion 
of components across the simulation domain is modelled together with the biotic reactions. 
The substrates and products of this activity are diffusing, and these dynamics are solved 
function the boundary conditions imposed over the limits of the domain. This type of models 
is used to describe microbial growth assuming that cellular division occurs a certain 
microbe’s age, cycle stage or size. Therefore, division or death occurs independently for 
each of the cells and only function of local conditions, capturing the heterogeneity of the aggregate 
and supporting non-linear growth models [1].


## Before having fun

This IbM framework is build up in MATLAB. Therefore, MATLAB must be installed in your computer.

You can use the [free 30-day trial](https://mathworks.com/campaigns/products/trials.html?ef_id=CjwKCAjwqJSaBhBUEiwAg5W9p96Y1NtC8BCa4Pw_wm3sswXR27ZkvuHZtWMOMUntOrmDSc1Ib3MGCRoCILQQAvD_BwE:G:s&s_kwcid=AL!8664!3!463011314378!p!!g!!matlab%20downlaod&s_eid=ppc_6588247642&q=matlab%20downlaod&gclid=CjwKCAjwqJSaBhBUEiwAg5W9p96Y1NtC8BCa4Pw_wm3sswXR27ZkvuHZtWMOMUntOrmDSc1Ib3MGCRoCILQQAvD_BwE){:target="_blank"}

## Granule version

Lorem ipsum

## Suspension version

**Comming soon...**

## References

[1] Hellweger, F.L., et al., (2016). doi: 10.1038/nrmicro.2016.62

Individual based model on bacterial colonies

Mathemathical model simulating the diffusion of compounds and bacterial growth.

To execute the model (Granule versions):

1. Create/modify Excel with all parameters
2. Create preset-file using the 'create_mat' function
3. Execute the model with a call to 'IbM(sim_number)'
