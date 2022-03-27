# Detachment

Detachment of biomass can be done in multiple ways in the model.

## Naive

When the detachment method is set to naive, the granule is set to have a maximum radius. Every bacterium that is outside of this radius is removed.

## Mechanistic

In the mechanistic detachment, the general idea is that due to friction on the outside of the granule cells are removed continuously. In the discretised model, this means that not only cells on the outside have to be removed. To account for this temporal discretisation, cells that are not directly removed (i.e. on the outside of the granule) also have a mass decreasing correction.

To determine which cells need to be detached (i.e. removed) and which cells need mass reduction due to erosion and how much, the time of detachment (`Tdet`) is calculated, stating the time it will take untill all biomass is removed per gridcell.

The `Tdet` can be calculated by moving the time of detachment front (starting a zero at the edge of the granule) inwards. This can be done by the fast-marching method, which does exactly this on a spacially discretised domain.

### Fast marching method

In order to move a front forward, the direction and speed of the front at each point is required. The local detachment speed is calculated as a function of the distance from the centre of the granule: `Fdet = kDet * r^2`.

The fast marching method sequentially finds the `Tdet` values for the moving front. After initiating the front at the edge of the granule, the algorithm is continuously moving the front forward using three arrays that are continuously updated:

- Visited: all gridcells that the front has already passed, the value for `Tdet` has been determined in these cells
- Narrow band: all gridcells that are close to being passed by the moving front
- Far: all gridcells that yet have to be passed by the moving front that are not in the narrow band

Initialisation is done as follows:

- All gridcells outside of the granule are assigned `Tdet = 0` and are placed in Visited
- The gridcells at the edge of the granule are assigned a `Tdet` value, based on the local detachment speed and are placed in the Narrow band
- All other gridcells are assigned a `Tdet` value of infinity and are placed in Far

The fast marching algorithm then consists of only a few steps, which are repeated until all gridcells are placed in Visited:

1. Pick the gridcell with the lowest value of `Tdet`
2. Place this gridcell in Visited
3. Update all neighbouring gridcells of this gridcell
   1. All neighbours in Far get placed in Narrow band
   2. The `Tdet` values of all neighbours in Narrow band are updated based on their respective Visited neighbours

Updating the `Tdet` values of neighbours is only done using the Visited neighbours. This updating can be seen as using the new position of the front (i.e. the Narrow band) to find out when the front will reach the gridcells that have to be updated.

### Detachment ratio

Using the time of detachment per gridcell, the ratio of mass that is lost due to detachment/erosion is calculated as: `r = dT / Tdet`. `dT` is the timestep of detachment, i.e. `dT_bac`.

If this ratio is greater than 1, this means that all cells in these gridcells are to be removed. For all gridcells in which this ratio is smaller than 1, the mass of bacteria is updated as `mass_new = mass_old * (1 - r)`.

If cells, after removal of mass, are below the minimum size of bacteria they are removed. However, if inactivation is enabled in the model, then cells in the middle of the granule are assumed to be protected from removal by the outside bacteria. In the mechanistic detachment algorithm cells in the middle of the granule also lose mass, thus are at a certain point below this threshold. In order to keep these cells in the system, only cells within a certain radius (determined by the furthest active bacterium in the system) are removed.

## Suspension

When the model is running in suspension modus, bacteria are removed based on the dilution rate. All bacteria that have a growth rate lower than the dilution rate (set via the HRT), are removed.
