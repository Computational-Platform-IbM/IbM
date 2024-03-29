# Individual-based Model framework

*Contributors: Chiel van Amstel, Eloi Martinez-Rabert, Rebeca Gonzalez-Cabaleiro

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
and supporting non-linear growth models[^1].
____________________________
**Methods of IbM can be downloaded [here](https://github.com/Computational-Platform-IbM/IbM/raw/main/Documents/Methods.pdf).**<br>
Detailed description of [detachment](https://github.com/Computational-Platform-IbM/IbM/blob/main/Documents/detachment.md) and [dynamic time step](https://github.com/Computational-Platform-IbM/IbM/blob/main/Documents/timestepping.md) is also available.

## Before having fun...

**:warning: To open the links in a new tab: right click on the link + "Open link in new tab". :warning:**

<!-- ### :information_source: Relevant information for reviewers

There is no peer-review process at this time. 

**Thanks for your help in improving the quality of our manuscript.** -->

### :gear: MATLAB installation

This IbM framework is build up in MATLAB. Therefore, MATLAB must be installed in your computer.
<br>You can use the [free 30-day trial](https://www.mathworks.com/campaigns/products/trials.html?ef_id=CjwKCAjwqJSaBhBUEiwAg5W9p96Y1NtC8BCa4Pw_wm3sswXR27ZkvuHZtWMOMUntOrmDSc1Ib3MGCRoCILQQAvD_BwE:G:s&s_kwcid=AL!8664!3!463011314378!p!!g!!matlab%20downlaod&s_eid=ppc_6588247642&q=matlab%20downlaod&gclid=CjwKCAjwqJSaBhBUEiwAg5W9p96Y1NtC8BCa4Pw_wm3sswXR27ZkvuHZtWMOMUntOrmDSc1Ib3MGCRoCILQQAvD_BwE) or "Sign in"/"Create account" if your university already has a MATLAB campus license (see [Campus-Wide License](https://www.mathworks.com/academia/tah-support-program/eligibility.html) search).

Click [here](https://www.mathworks.com/help/install/) for more information about MATLAB Installation and Licensing.

### :clipboard: Instructions to Download and Run IbM

1. Download .zip code. Last version: `v2.4.3`. [Download code](https://github.com/Computational-Platform-IbM/IbM/archive/refs/tags/v2.4.3.zip).
2. Extract files to a destination (:bulb: recommendation - Desktop).
3. Open MATLAB.
    - For more information about MATLAB Layout and how to change it, click [here](https://www.mathworks.com/help/matlab/matlab_env/change-the-desktop-layout.html).
4. Go to the **Code folder<sup>2</sup>**...
    &#09;<br><sup><sup>2</sup>Code folder: folder with `IbM.m` file. </sup>
    &#09;<br>→ writing `cd newFolder` to *Command Window* (more info about [cd](https://www.mathworks.com/help/matlab/ref/cd.html)).
    &#09;<br>→ or using *Folder Toolbar* - paste the folder name where the code was extracted.
5. Generate the path to the main code writing `addpath(genpath('lib'))` to *Command Window*.
6. Create the seed-file<sup>1</sup>:
<br><sup><sup>1</sup>Seed-file: `.mat` file with all information of simulation. This file is used to execute the code.</sup><br>
    1. Modify main Excel (lib\planning\Excels\main.xlsx) with all parameters. 
    <br><sup>Instructions on how to use main.xlsx in *Information* sheet, and [Granule version](https://github.com/Computational-Platform-IbM/IbM#granule-version) | [Suspension version](https://github.com/Computational-Platform-IbM/IbM#suspension-version). (**Excel setup** section)</sup>
    2. Create seed-file writing `create_mat` to *Command Window*.
    3. Save seed-file writing <br> `save('planning/sim_xxxx.mat','grid','bac','constants','init_params','settings','-v7.3')`</br> to *Command Window* (:bulb: `sim_xxxx.mat`, where xxxx is the simulation number [from 0001 to 9999]).
    4. Remove items from workspace, freeing up system memory - writing `clear all` to *Command Window* (more info about [clear](https://www.mathworks.com/help/matlab/ref/clear.html)).
7. Execute IbM code:
    1. Copy the desired seed-file to Code folder (folder with `IbM.m` file)
    2. Call to `IbM(xxxx)` (:bulb: where xxxx is the chosen simulation number).<br>
    Once simulation has been run, `sim_xxxx.mat` (i.e., seed-file) is moved to the corresponding **Results** folder.
8. Get Data (see [Get Data](https://github.com/Computational-Platform-IbM/IbM#get-data) section) or Visualization of Results (see [Data visualization](https://github.com/Computational-Platform-IbM/IbM#data-visualization) section).

### :card_file_box: Supplementary Information

· The seed-file can be reviewed writing `load('sim_xxxx.mat')` (or `load('planning\sim_xxxx.mat')` if this isn't in **Code folder**) to *Command Window*, or double click to `sim_xxxx.mat` in *Current folder*. All data is on *Workspace*.
>**Seed-file structure** (`sim_xxxx.mat`)<p>
`bac`               - Information of initial position (*x,y*), molar mass, radius, species and which cell is active/inactive.
<br>`constants`     - Summary of all parameters of model (Discretization, Kinetic parameters, Metabolic matrix...).
<br>`grid`          - Information of grid.
<br>`init_params`	- Initial conditions of the system (HRT, concentration in aggregate and bulk liquid).
<br>`settings`      - Settings of the model (discretization, HRT, detachment, pH, model version, parallization).

· During the execution of the `IbM()` , *Warning messages* may be displaied to *Command Window*. In general, *Warning messages* are merely informative about the simulation progress. See list of [Warning messages](https://github.com/Computational-Platform-IbM/IbM#warning-list).
<br>A MATLAB error looks like:
> Error using `function` (line x)<br>*Details of error*

## Granule version

### :bar_chart: Excel setup

#### · Microbial spatial distribution study[^2]
:bulb: You can find an Excel template (`granule_template.xlsx`) in `planning\Excels\Templates`. Remember to copy this Excel to `planning\Excels\`, and rename it to `main.xlsx`.<p>
In order to setup the IbM for this study, the following key settings must be selected on the main Excel (see *Settings* sheet in `main.xlsx`):
- ***Dynamic dT***                      - *TRUE*
- ***Variable HRT***                    - *TRUE*
- ***Initialisation method***           - *granule*
- ***kDist***                           - *1.0*
- ***Inactivation enabled***            - *TRUE*
- ***pH bulk concentration corrected*** - *FALSE*
- ***pH solving included***             - *FALSE*
- ***Speciation included***             - *FALSE*
- ***Structure model***                 - *TRUE*
- ***Structure model type***            - *Neut*, *Comp*, *Comm* or *Copr* → changes in the feeding regime.

:warning: *Settings* sheet is a summary of all seetings in the model. **DO NOT CHANGE THEM IN THIS SHEET**, but in the respective tab.

### :newspaper: Publications

```
E. Martinez-Rabert, C. van Amstel, C. Smith, W. T. Sloang, R. Gonzalez-Cabaleiro (2022), 
"Environmental and ecological controls of the spatial distribution of microbial populations in aggregates". 
PLOS Computational Biology. doi: 10.1371/journal.pcbi.1010807
```

## Suspension version

### :bar_chart: Excel setup

#### · Metabolic heterogeneity of Comammox *Nitrospira* [^3]
:bulb: You can find an Excel template (`CmxNitrospira_template.xlsx`) in `planning\Excels\Templates`. Remember to copy this Excel to `planning\Excels\`, and rename it to `main.xlsx`.<p>
In order to setup the IbM for this study, the following key settings must be selected on the main Excel (see *Settings* sheet in `main.xlsx`):
- ***Dynamic dT***                      - *TRUE*
- ***Variable HRT***                    - *FALSE*
- ***Initialisation method***           - *suspension*
- ***kDist***                           - *1.50*
- ***Inactivation enabled***            - *TRUE*
- ***pH bulk concentration corrected*** - *FALSE*
- ***pH solving included***             - *FALSE*
- ***Speciation included***             - *FALSE*
- ***Structure model***                 - *FALSE*
- ***Structure model type***            - *"-"*

:warning: *Settings* sheet is a summary of all seetings in the model. **DO NOT CHANGE THEM IN THIS SHEET**, but in the respective tab.

### :newspaper: Publications

```
Martinez-Rabert, E., Smith, C.J., Sloan, W.T. et al. (2023)
Competitive and substrate limited environments drive metabolic heterogeneity for comammox Nitrospira.
ISME COMMUN. 3, 91 (2023). doi: 10.1038/s43705-023-00288-8
```

## Get Data
Relevant results are saved every **dt saved** time (specified in Excel, *Discretization* sheet). Results are saved in the corresponding **Result folder** (`\Results\xxxx`, where `xxxx` is the simulation number) in `results1D.mat` file.
>**Results-file structure** (`\xxxx\results1D.mat`)<p>
`bac_saved`               - Number of cells (*nBacs*), position (*x,y*), radius (*radius*), species/metabolism (*species*), state (*active*), actual growth rate (*mu*) [*columns* - cell iD; *rows* - time].
<br>`conc_saved`          - 3D Matrix of concentration profiles on the transverse plane of aggregate for all substrates [*columns* - x position; *rows* - time; *z* - substrate/product]
<br>`pH_saved`            - Matrix of local pH on the transverse plane of aggregate for all substrates [*columns* - x position; *rows* - time].
<br>`reactor_saved`       - Bulk liquid concentration (*bulk_conc*), Hydraulic Retention Time (*HRT*), density of granule (*granule_density*) [*columns* - substrate/products; *rows* - time].

The data that is saved to the **Results-file** can be modified in the `lib\post_processing\save_slice.m` script.

### :clipboard: Instructions to Get Data
1. Go to the **Result folder** of simulation xxxx (`\Results\xxxx\`) through *Current Folder* browser.
2. Load **Result-file** (`results1D.mat`) - double-click on `.mat` file or writing `load('results1D.mat')` on *Command Window*.
3. Saved data to a file (`.txt`,`.dat`,`.csv` or `.xls*`) using `writematrix(input_data,filename,__)` function (more info about [writematrix](https://www.mathworks.com/help/matlab/ref/writematrix.html)). The file is automatically generated in the **Result folder**.
    &#09;<br><sup>:warning: `writematrix` cannot save *struct*-typed data. Then, `input_data` must start with *struct* name (`bac_saved`, `conc_saved`, `pH_saved` or `reactor_saved`) followed by a dot `.` and the name of variable (e.g., `writematrix(bac_saved.nBacs,'test.xlsx')`).</sup>
    &#09;<br><sup>:bulb: Data can be saved even if the simulation is still running.

## Data visualization

### :clipboard: Instructions for Data visualization
1. Go to the **Code folder** (folder with `IbM.m` file).
2. Generate the path to the main code writing `addpath(genpath('lib'))` to *Command Window* (if not already generated).
3. Run Data visualization → `function()` listed below. 

### · Draw aggregate/s (Granule & Suspension version)
Function to plot cells for a given simulation number and time → `plotBacsMatlab(simulation_number, Time)`.<br>
<sub>Inputs: `simulation_number` is the chosen simulation, and `Time` in days.</sub><br>
<br>**To modify colours, open `lib\debugging\plotBacs.m` script. (Lines 15-17)**<br>
:warning: Make sure that *species_per_color* vector corresponds to number of species/metabolisms.
> If the simulation has not finishied yet, copy the seed-file (`simulation_xxxx.mat`) to the corresponding **Result folder** (`\Results\xxxx`).

### · Draw substrate/product profiles (Granule version)
Function to draw 2D substrate/product profiles → `plot2d_S(simulation_number, substrates_number, Time)`.<br>
<sub>Inputs: `simulation_number` is the chosen simulation, `substrate_number` following substrate order of Excel, and `Time` in days.</sub><br>
<br>**To modify colours, open `lib\lib\post_processing.m` script. (Line 6)**<br>

### · Analyse results over time (Granule & Suspension version)
Function to plot over time (_i_) number of bacteria (active and inactive) in the simulation domain, (_ii_) active mass per species, (_iii_) density of the reactor, (_iv_) oxygen uptake rate, and (_v_) bulk liquid concentrations over time for a given simulation number and its replicates → `analyseGranule_replicates(simulation_number, nReplicates, finished)`.<br>
<sub>Inputs: `simulation_number` is the chosen simulation, `nReplicates` is the number of replicates (subsequent numbers), and `finished` indicates if simulation is completed (1) or not (0) .</sub><br>

### · Draw spatial distribution of microbial colonies (Granule version)
Function to plot radial distribution (*layered stratification*) and angular distribution (*columned stratification*) of microbial colonies → `plotStratification(simulation_number, finished, Time)`.<br>
<sub>Inputs: `simulation_number` is the chosen simulation, `finished` Simulation finished: 0 or 1 [N/Y], `Time` in days.</sub><br>
<br>Images `.tif` of plots are automaticatly created to **Result folder** (`\Results\xxxx`). 
> Resolution of Spatial distribution plots can be modified → `lib\post_processing\plotStratification.m` script.
 <br>· **Line 59** for radial distribution (*layered stratification*).
 <br>· **Line 107** for angular distribution (*columned stratification*).

### · Ecological analysis (Suspension version)
Function to plot pairwise interactions → `flocbyfloc(simulation_number, nReplicates, finished)`.<br>
<sub>Inputs: `simulation_number` is the chosen simulation, `finished` Simulation finished: 0 or 1 [N/Y], `Time` in days.</sub><br>
1. Run `flocbyfloc(simulation_number, substrates_number, Time)`. An Excel file (`RawData_flocbyfloc.xlsx`) is automatically generated with the mass of all species for each floc at time 0 (*t_0* sheet) and at the end of simulation (*t_max* sheet).<br>
<sub>Inputs: `simulation_number` is the chosen simulation, `substrate_number` following substrate order of Excel, and `Time` in days.</sub><br>
2. Drop .xlsx file (name `data_2OEN.xlsx`) with the list of flocs that you want analyze to `lib\post_processing\Data` folder, following the same structure as `RawData_flocbyfloc.xlsx` (*t_0* and *t_max* sheets).
3. Run `f2OEN` function. This will give the corresponding ecological anlaysis with the sample size of each pirewise interaction and the representation of the data in which Kendall's τ has been calculated.

## Warning List

All possible *Warning messages* are listed here:
<br><sub>Format: *WarningTag* (`script.m` source) - Action/Message</sub>
- *DEBUG:noactionRequired* (`lib\calculate_bulk_concentration.m`) - Negative concentration encountered and corrected.
- *DEBUG:noactionRequired* (`lib\detachment\recalculateT.m`) - Detachment speed equals 0, thus infinite time of crossing.
- *dt_diffusion:Limits* (`lib\integTime.m`) - Maximum/Minimum **dt_diffusion** value is reached.
- *dt_bac:Limits* (`lib\integTime.m`) - Maximum/Minimum **dt_bac** value is reached.
- *Diffusion:SlowConvergence* (`lib\integTime.m`) - **dt_diffusion** is increased. 
- *Diffusion:SlowConvergence* (`lib\integTime.m`) - **dt_bac** is increased. 
- *Diffusion:ConvergenceStuck* (`lib\integTime.m`) - **dt_diffusion** is decreased. 
- *Diffusion:UpwardRES* (`lib\integTime.m`) - **dt_diffusion** is decreased. 
- *Diffusion:NegativeConcentration* (`lib\diffusion\diffusionMG.m`) - **dt_diffusion** is reduced.
- *CellDivision:TooFastDivision* (`lib\integTime.m`) - **dt_bac** is decreased. 
- *BulkConcentration:ValueJump* (`lib\integTime.m`) - **dt_bac** is decreased.

:warning: Special *Warning messages*, **POSSIBLE ACTION REQUIRED**: 
- *DEBUG:actionRequired* (`lib\determine_diffusion_region.m`) - Higher simulation domain must be set.
- *DEBUG:actionRequired* (`lib\calculate_bulk_concentration.m`) - Negative bulk concentration encountered after pH control.
- *DEBUG:actionRequired* (`lib\detachment\computeRoot.m`) - ValueError:Should always be above 0...
- *DEBUG:actionRequired* (`lib\detachment\recalculateT.m`) - All neighbours have infinite time of crossing.

>**Legend**
<br>dt_diffusion - step time applied in solving the diffusion-reaction equation 
<br>dt_bac - step time applied in computation ofbateria mass change, division and inactivation

## Contact

**Eloi Martinez-Rabert**. :envelope: 2424069M@student.gla.ac.uk or eloi.mrp@gmail.com

### References

[^1]: Hellweger, F.L. *et al.*, (2016). *Nature Reviews Microbiology*. doi: 10.1038/nrmicro.2016.62<br>
[^2]: Martinez-Rabert, E. *et al.*, (2022). *PLOS Computational Biology*. doi: 10.1371/journal.pcbi.1010807<br>
[^3]: Martinez-Rabert, E. *et al.*, (2023). *ISME Communications Journal*. doi: 10.1038/s43705-023-00288-8<br>
