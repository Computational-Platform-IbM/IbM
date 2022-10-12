# Changelog
All notable changes to this project will be documented in this file.

## [2.4.1]
### Added
- Proper README file
### Fixed
- Bugs in *Structure model*
### Removed
- Unnecessary scripts

## [2.4.0]
### Added
- Mechanistic detachment implemented
- Function to plot detachment time per grid cell
- Function to plot bacteria with overlay of detachment times
- Excel template for AOB/NOB/CMX/NRMX/AMX simulation
- Excel template for AOB/NOB/AMX simulation
- Setting for suspension initialisation
- Height correction for single aggregate simulations reaction rate in simulation domain
### Changed
- Loading of model settings from Excel completely refactored
- Also reduce dT_bac when multiple round of division are needed per timestep (overcrowding)
- HRT not only changed when bulk concentrations are higher than setpoint, but also when below the setpoint
- Paramters of adaptive step size optimized
- Diffusion region determined based on center of the grid cell
### Fixed
- Only removing outside bacteria due to erosion, no longer erroneous removal of inside bacteria
- Inactivation method now correctly takes into account the growth rate
- Correctly calculating bulk liquid concentrations (gradient was not calculated when HRT was changed)
- Increase maximum possible number of bacteria per gridcell from 255 to 65535 (uint8 to uint16)
### Removed
- Unnecessary diffusion of non-diffusive compounds (e.g. N2 was in concentration & reaction matrix, but does not diffuse)

## [2.3.0] - 2021-11-15
### Added
- Dynamic dT for diffusion
- Dynamic dT_bac for timestepping
- Automatic switch to parallelisation based on the number of bacteria
- Automatic initialisation of the parallel pool with the correct number of workers
- Automatized .mat creator and registration in simulation_log.md
- kDist constant implementation to shoving algorithm
- New setting: structure model and types (Neutralism, competition, commensalism, co-protection)
- New setting: consideration of pH or not
- Documentation for dynamic dT
### Changed
- Calculus of "correction_concentration_steady_state" based on chosen TolAbs by user
- Independent scripts for mu_max/decay and Monod terms calculus
- All reaction matrix calculations are now in one function, which changes behaviour based on settings.parallelized
- Only structure_model: HRT is recalculated when bulk_conc is higher or lower than setpoint


## [2.2.0] - 2021-10-20
### Added
- Plotting function for reaction matrix (2D profile).
- Improved workflow for simulation using a single file (IbM.m).
- Saving of results.
- Functions for analysing the result files.


## [2.1.0] - 2021-10-19
### Added
- Plotting function for concentrations (2D profile).
- Plotting function for convergence to steady state (norm of delta-concentration).
- Parallelisation of reaction matrix calculations.

### Changed
- Updated the Materials & Methods document.

### Fixed
- Minor bug where incorrect diffusion region was taken in the final smoothing step (multigrid method).


## [2.0.0] - 2021-10-18
### Added
- Docstrings for every function and file
- Check pH after `pH_solve` to be in range [1, 14].
- Check for negative concentrations after `pH_solve`.
- Check for non-convergence using consecutive RES values, continue with steady-state if this is detected.
- Variable for the number of diffusion iterations per steady-state check.
- Shoving algorithm in Java.
- Quadtree datastructure for shoving.
- Complete vectorisation of the bacterial activity (division, death, inactivation, etc.).
- Multigrid solver using V-cycles for solving the diffusion.
- Use of convolutions in the multigrid method for potential further speedup (gpu or FFT method).
- Adaptation in standard multigrid method to only calculate diffusion in the diffusion region.
- Clamp the error correction in Newton-Raphson within the valid range, so that it will always converge.
- Reimplemented the focus-region (only solve smaller domain) for significant speed-up.
- Build-in profiling option.
- Debug warnings for negative concentrations, non-convergence, etc.
- Plotting function for bacteria.
- Plotting function for bulk concentration.
- Plotting function for convergence to steady state (RES value).
- Plotting function for diffusion region.
- Plotting function for maximum error over simulation time (max RES value & delta-concentration).
- Plotting function for profiling.

### Changed
- Complete refactorisation of the code for improved readability and ease of development.
- Restructuring the R.mat struct into separate structs/objects with different functionality.
- Folder structure with better organisation.
- Separation of functionality in different files.
- Calculate kinetics only once per species per gridcell, then using the cumulative mass per species to calculate concentration change.
- Check steady-state based on the maximum RES value instead of the average RES value.
- Use different calculation method (convolutions) to determine the RES values.
- Only check the close neighbours in the shoving algorithms instead of all other bacteria.
- Only compute bulk-layer pH once per reaction matrix calculation.
- 
### Removed
- Check in `pH_solve` to see whether root is present in range [1, 14] (will always hold).
- Bracketing after Newton-Raphson jumps outside of valid pH range.

### Fixed
- Removed redundant function calls.
- Correct function call order in model flow.
- Only perform diffusion in diffusion region, not in bulk layer too.
- Only check RES values in the diffusion region. 


## [1.0.0] - 2021-08-28
Direct import of code from RGC with minor bug fixes.

### Fixed
- Bugs regarding reaction matrix, diffusion and bulk concentration.



[2.4.0]: https://github.com/Computational-Platform-IbM/IbM/compare/v2.3.0...v2.4.0
[2.3.0]: https://github.com/Computational-Platform-IbM/IbM/compare/v2.2.0...v2.3.0
[2.2.0]: https://github.com/Computational-Platform-IbM/IbM/compare/v2.1.0...v2.2.0
[2.1.0]: https://github.com/Computational-Platform-IbM/IbM/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/Computational-Platform-IbM/IbM/releases/tag/v2.0.0
[1.0.0]: https://github.com/Computational-Platform-IbM/IbM/releases/tag/v1.0.0
