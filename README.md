# Particle Disaggregation Code
This repository contains Matlab code to calculate the evolution of the 
particle size distribution in a single layer of the water column. The 
model allows for particle aggregation, disaggregation, and sinking, 
and also changes in aggregate size from cell growth (see `SetUpCoag.m`).

## What the code does.
The code numerically solves the aggregation-disaggregation mass 
balance equations using a so-called *sectional approach* developed
by Gelbard and Seinfeld (J. Colloid and Interface Sci., 68:363-382, 1979).

All the code is documented within the code itself. What follows is a brief 
description of the code to get things started. 

The main code driver is the function `coag_driver.m` and this runs the setup, numerical 
solution, and plotting of the simulation. The code `run_driver.m` is a script which
shows how to run the model for a range of different parameter values.  

At the simplest level, the only file that needs to be changed is `SetUpCoag.m`. 
This documented file contains the parameters for any simulation. It is currently set
up to run a simulation for a 65 m thick layer of water. 

### Coagulation kernels
The model can run simulations with multiple coagulation kernels and these are
calculated in `CalcBetas.m`. These include 

 * Brownian, Shear, and Differential Sedimentation collision processes.
 * Rectilinear kernels (i.e. no fluid interactions between colliding particles) and 
   curvilinear kernels (i.e. accounting for fluid interactions) for shear and 
   differential sedimentation
 * A fractal kernel based on the work of Bruce Logan.
 
The kernels are calculated one at a time via calls from `coag_driver.m` and can then be
added (see the example in `coag_driver.m`). 

There are two versions of `CalcBetas.m`

 * `CalcBetas.m` assumes a fractal-modified Stokes' Law for the sinking speed of the aggregates
   when calculating the kernels for differential sedimentation. 
 * `CalcBetas_vs.m` allows for other settling speed relationships using a user-supplied function
   `SettlingVelocity.m`
   
### Disaggregation
Disaggregation of particles by turbulent shear is modeled using results of Jackson 
(Deep-Sea Research II, 42: 159-184, 1995). This is a very simple model that does
prevents the formation of very large particles. It includes both erosion and fragmentation,
though the former process is thought not to be important in the oceans. Disaggregation
can be switched off by setting the constants `c3` and `c4` to zero in `SetupCoag.m`. 

## Things to watch for as you play
There are some things to watch for as you vary parameters in the simulations. Key amongst
which is the possibility that you may create particles larger than the largest one allowed
in the simulation (this is set by the variable `n_sections` in `SetUpCoag.m` which sets
the number of logarithmically defined size bins). If these coagulation losses are greater
than the sinking losses (i.e. the loss due to sinking out of the bottom of the layer), then
you should increase the number sections until this ratio is acceptable. These quantities
are calculated in the functions `SectionalMassBalance.m` (which calculates all processes
adding and removing particles to a size bin) and `TotalMassBalance.m` which calculates
the system-wide gains and losses (these should balance in steady state). 
