# upic-emma-2.0

## PSTAD Branch

(December 13, 2018 -- make a version for PSTAD solver.)

##


This is the repo for UPIC-EMMA 2.0, a UPIC-EMMA based on UPIC 2.0 by Viktor Decyk (which includes MPI/OpenMP parallelization). 

New features include:

* non-square cells
* FDTD operators for benchmarking
* PSATD operator for boosted-frame simulation
* additional boundary conditions for fields and particles to allow the study of laser plasma interactions at different intensities.
* OpenPMD support

Upcoming features include:

* binary collision
* Galilean transformation
* OpenACC/OpenMP with offloading support



It also supports Doxygen documentation, to use it, type:

doxygen Doxyfile



## Input parameters:

* N_threads (integer): number of OpenMP threads
* FDTD (integer): FDTD=1 -> finite difference operator follows FDTD FDTD=0 -> use traditional difference opreator ( i k )
* cfl (real): timestep in terms of the CFL condition, if cfl=0.95 then dt is 0.95 * CFL limit
* laserpulse (logical): if laserpulse = .TRUE. then a laser is initialized
* relativistic (logical): if relativistic = .TRUE., then relativity is turned on for orbit integration
* moving_ions (logical): if moving_ions = .TRUE., then ions are mobile and have mass, if .FALSE., then ions are treated as immobile fluid
* 




Revision History:

1/8/2018:  Initial README (F. S. Tsung)
11/19/2020: Adding documentation for input deck

