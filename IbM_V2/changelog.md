# Changelog v0.1-alpha
List all changes and/or improvements from the original model of R. Gonzalez-Cabaleiro as performed by C. van Amstel during the refactoring of the code.


## General
*Improved readability*
- Changed folder structure
- Separation of functional components into files
- Added function docstrings
- Improved variable names
- Separation of constants and variables
- Vectorisation of bacterial activity using masks

*Profiling*
- Add build-in profiling capability grouped per major functionality of the model

*Debugging/plotting*
- Build-in functions to plot:
    - Profiling data
    - Bacteria (species)
    - Diffusion region (with/without radius of boundary layer)
    - Convergence to steady state during diffusion iterations
    - Maximum error over simulation time
    - Norm between new and previous concentrations over time
    - Bulk concentrations over time
    - 2D Concentration profiles
- Introduction DEBUG warnings (togglable)



## Computational improvements
*General*
- Focus region based on diffusion region

*Calculation of reaction matrix*
- Compute bulk layer pH only once
- Complete parallelisation possible using chunks

*pH and speciation algorithm*
- Removed check to see whether root is found in 1 < pH < 14 (Should always hold)
- Removed bracketing whenever Sh < 1e-14
- Introduced bounded correction, so that always Sh > `setValue`

*Diffusion*
- Change from `mldivide` to Multigrid solver using V-cycles
- Use of convolutions for potential further speedup (gpu or FFT method)
- Due to use of convolutions, only Jacobi smoothing is possible, even though Gauss-Seidell would be faster (in theory)
- Integrated diffusion-region-only diffusion solver using the `diffusion_region`

*Bacterial activity*
- Complete vectorisation of bacterial activity using boolean masks

*Shoving*
- Rewritten in Java
- Quadtree implementation for shoving
- Only check neighbours for overlap instead of each bacterium


## Algorithmic changes
*Calculation of reaction matrix*
- Calculate kinetics only once per species per gridcell. Reaction matrix uses cumulative mass per species to calculate concentration change.

*pH and speciation*
- Add check for pH outside of [1, 14]-range
- Add check for non-convergence
- Add check for negative concentrations in speciation 

*Steady-state check*
- Change number of diffusion iterations per steady-state check
- Check maximum RES value (default, but changeable)
- Add non-converging check to continue even if steady-state is not reached based on a set convergence-tolerance
- Use convolutions for checking of steady-state


## Error corrections
*General*
- Remove redundant function calls
- Cleanup/correct function call order in model flow

*Diffusion*
- Diffusion was previously performed over entire (focus) domain, disregarding the bulk layer that was directly around the granule

*Steady-state check*
- Check only the RES in the `diffusion_region`. By definition, the RES value in the edge of the bulk layer will be nonzero, thus the model never converged.
