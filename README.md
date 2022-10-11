# Individual-based Model framework.

From: Gonzalez-Cabaleiro group.

*Contributors: Chiel van Amstel, Eloi Martinez-Rabert, Rebeca Gonzalez-Cabaliero

**:warning: To open the links in a new tab: right click on the link + "Open link in new tab". :warning:**

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


## Before having fun...

### :gear: MATLAB installation

This IbM framework is build up in MATLAB. Therefore, MATLAB must be installed in your 
computer.You can use the [free 30-day trial](https://www.mathworks.com/campaigns/products/trials.html?ef_id=CjwKCAjwqJSaBhBUEiwAg5W9p96Y1NtC8BCa4Pw_wm3sswXR27ZkvuHZtWMOMUntOrmDSc1Ib3MGCRoCILQQAvD_BwE:G:s&s_kwcid=AL!8664!3!463011314378!p!!g!!matlab%20downlaod&s_eid=ppc_6588247642&q=matlab%20downlaod&gclid=CjwKCAjwqJSaBhBUEiwAg5W9p96Y1NtC8BCa4Pw_wm3sswXR27ZkvuHZtWMOMUntOrmDSc1Ib3MGCRoCILQQAvD_BwE) or "Sign in"/"Create account" if your university already has a MATLAB campus license (see ["Campus-Wide License" Search](https://www.mathworks.com/academia/tah-support-program/eligibility.html)).

Click [here](https://www.mathworks.com/help/install/) for more information about MATLAB Installation and Licensing.

### :arrow_down: Code download and seed-file<sup>1</sup> creation
<sup><sup>1</sup>Seed-file: `.mat` file with all information of simulation. This file is used to execute the code.</sup><br>
**  <sup>All seed-files must be saved as 'sim_xxxx.mat'. More info in step 6 (Create the seed-file)</sup>**

Instructions to Download and run IbM:
1. Download .zip code. Last version: `v2.4.1`. [Download code](https://github.com/Computational-Platform-IbM/IbM/archive/refs/tags/v2.4.0.zip).
2. Extract files to a destination (:bulb: recommendation - Desktop).
3. Open MATLAB.
    - For more information about MATLAB Layout and how to change it, click [here](https://www.mathworks.com/help/matlab/matlab_env/change-the-desktop-layout.html).
4. Go to the **Code folder<sup>2</sup>**...
    &#09;<br><sup><sup>2</sup>Code folder: folder with `IbM.m` file </sup>.
    &#09;<br>→ writing `cd newFolder` to *Command Window* - more info about [cd](https://www.mathworks.com/help/matlab/ref/cd.html).
    &#09;<br>→ using *Folder Toolbar* - paste the folder name where the code was extracted.
5. Generate the path to the main code:
    &#09;<br>→ writing `addpath(genpath('lib'))` to *Command Window*.
6. Create the seed-file:
    1. Modify main Excel (lib\planning\Excels\main.xlsx) with all parameters. 
    <br><sup>Instructions on how to use main.xlsx in *Information* sheet, and [Granule version](https://github.com/Computational-Platform-IbM/IbM#granule-version) | [Suspension version](https://github.com/Computational-Platform-IbM/IbM#suspension-version).</sup>
    2. Create seed-file writing `create_mat` to *Command Window*.
    3. Save seed-file writing <br> `save('planning/sim_xxxx.mat','grid','bac','constants','init_params','settings','-v7.3')`</br> to *Command Window* (:bulb: `sim_xxxx.mat`, where xxxx is the simulation number [from 0001 to 9999]).
7. Execute IbM code:
    1. Copy the desired seed-file to Code folder (folder with `IbM.m` file)
    2. Call to `IbM(sim_xxxx)` (:bulb: `sim_xxxx`, where xxxx is the chosen simulation number).

## Granule version

Lorem ipsum

## Suspension version

:building_construction: **Comming soon...**

## References

[1] Hellweger, F.L., et al., (2016). *Nature Reviews Microbiology*. doi: 10.1038/nrmicro.2016.62

________________
Individual based model on bacterial colonies

Mathemathical model simulating the diffusion of compounds and bacterial growth.

To execute the model (Granule versions):

1. Create/modify Excel with all parameters
2. Create preset-file using the 'create_mat' function
3. Execute the model with a call to 'IbM(sim_number)'
