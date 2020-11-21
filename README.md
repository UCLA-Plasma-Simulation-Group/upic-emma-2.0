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
* FDTD (integer): FDTD=1 -> finite difference operator of N-th order (N is even), or 0 -> use traditional difference opreator ( i k )
* cfl (real): timestep in terms of the CFL condition, if cfl=0.95 then dt is 0.95 * CFL limit
* laserpulse (logical): if laserpulse = .TRUE. then a laser is initialized
* relativistic (logical): if relativistic = .TRUE., then relativity is turned on for orbit integration
* moving_ions (logical): if moving_ions = .TRUE., then ions are mobile and have mass, if .FALSE., then ions are treated as immobile fluid
* Delta_x, Delta_y   = spatial cells dimension in units of (Delta)
* indx, indy = spatial dimensions of the simulation box   
                      along the x and y axis. Pay attention :   
                      The number of gridpoints Nx and Ny            
                      along x and y axis must be a power of 2 due   
                      to the FFT algorithm. We impose so :          
                      Nx = 2^indx and Ny = 2^indy where             
                      indx = 1 + floor(ln(Lx/Delta_x)/ln(2)) and
                      indy = 1 + floor(ln(Ly/Delta_y)/ln(2))
                      such that Lx = 2^indx and Ly = 2^indy
                      For example, with Delta = 1. :
                      Lx = 2, 4, 16 , 32, 64, 128, 256, 512, 1024,  
                           2048, 4096, 8192, 16384, 32768, 65536,   
                           131072, ... with                         
                      indx = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                             13, 14, 15, 16, 17, ... respectively  

* t_end = time at end of simulation, in units of              
                      (/ omega_p) or (/ omega_0) depending on units



#######################################################################
###             Plasma properties (useful only if plasma = 1) :         #
#######################################################################
* atomic_weight      = plasma atomic weight A (g/mol)                 #
* atomic_number      = plasma atomic number Z ()                      #
* ionization_state   = plasma ionization state Z* ()                  #
* ax/ay          = smoothed particle size in x/y direction            #
* den_me         = Number of electrons per macro electron             #
* Npic           = Number of macro electrons per spatial cell         #
* param          = vTe / c with c the velocity of light in vacuum     #
* ptx, pty, ptz  = root mean square electron momentum (relativistic=1)
                  or velocity (relativistic = 0) components in the x,
                  components in the x, y and z directions resp.     
                  in units of ( me ) Delta omega_p                 
                  Be careful, since Te = (Te,x + Te,y + Te,z) / 3, 
                  vTe = sqrt( (vtx^2 + vty^2 + vtz^2) / 3 )        
                  must be respected                                 

* px0, py0, pz0  = mean plasma momentum (relativistic=1) or velocity  
                  (relativistic = 0) components in the x, y and z    
                  directions resp.                                   
                                                                     


Revision History:

1/8/2018:  Initial README (F. S. Tsung)
11/19/2020: Adding documentation for input deck

