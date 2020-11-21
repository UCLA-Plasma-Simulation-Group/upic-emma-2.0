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



###             Plasma properties (useful only if plasma = 1) :         #
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
                                                                     


###         Laser pulse properties (useful only if laser = 1) :         #
* propdir  = 1 (x) or 2 (y) = main propagation direction              #
            (MPI parallelization is done along the y-direction)      #
                                                                     #
* polardir = polarisation parameter such that -1. <= polardir < 1. :  #
          _ polardir =     0.  corresponds to a linear polarization  #
                                           along : _ z (propdir = 1) #
                                                   _ x (propdir = 2) #
          _ polardir =   - 1.  correspond to a linear polarization   #
                                          along : _ y (propdir = 1)  #
                                                  _ z (propdir = 2)  #
          _ polardir = +/- 0.5 correspond to a circular polarization #
          _ polardir =    else correspond to elliptic polarizations  #
                       polardir > 0 -> clockwise                     #
                       polardir < 0 -> counterclockwise              #
                       
* shape    = 1 (Gaussian wave shape) or 0 (Plane wave)                #
* theta    = angle of incidence relative to the 'propdir'-axis        #
            assuming that shape = 1 (for a plane wave, rotate the    #
            plasma!!!)                                               #
            Be careful : theta must be between -45. and 45. so       #
                                                                     #
* tlaunch  = time at which the laser pulse starts to be launched      #
            in units of omega^-1                                     #
 * FWHMt    = time duration Full Width at Half Maximum of the laser    #
            pulse envelop in units of omega^-1                       #
                                                                     #
 * xfocal/yfocal = focal point coordinates along x and y in units of   #
                 (Delta) for Gaussian waves                          #
 * FWHMs    = laser pulse spot size Full Width at Half Maximum at the  #
            focal point in units of (Delta) for Gaussian waves       #
            (the laser pulse waist at focal point :                  #
                           W0 = FWHMx / sqrt(4ln2) )                 #
                                                                     #
* omega0   = Laser pulse angular frequency in units of (omega)        #
            -> omega0 = 1. is imposed if units = LPI                 #
                                                                     #
* a0       = Maximum amplitude of the laser pulse potential vector    #
            in units of (me Delta omega / e)                         #
            -> E0 = omega0 * a0 is the maximum amplitude of the      #
                    laser pulse electric/magnetic field              #
                                                                     #

###                       Boundary conditions                           #
                                                                    
 * BCx = 1 (periodic) or 0 (absorbing)                                 #
       condition along the x-axis for both particles and fields      #
                                                                     #
 * BCy = 1 (periodic) or 0 (absorbing)                                 #
       condition along the y-axis for both particles and fields      #
                                                                     #
 * For absorbing cond., the Perfectly Matched Layer (PML) technique    #
                      is used for the E-M fields (Berenger, 1996) :  #
                      _ PML_scheme = 0 : Yee scheme                  #
                        PML_scheme = 1 : Implicit scheme             #
                      _ L_PML is the PML layer thickness in units    #
                        of Delta                                     #
                      _ n_cond_PML is the integer n such that the PML#
                        electrical and magnetic conductivity read :  #
                        sigma(x) = sigma_m * (|x-x_I|/L_PML)^n       #
                        where x_I is the position of the interface   #
                        between the PML and the effective sim. box.  #
                        3 < n_cond_PML < 4 is optimal                #
                                                                     #

###                            Diagnostics :                            #
* Delta_t_diag = results are saved in text files in results/          #
                every Delta_t_diag in unit of (omega_p or omega0)^-1 #

###                        density profile :                            #
* density_x/y, write down expressions for density along x/y           #
              directions. The names of coordinate of x is x1         #
              and y is x2. many math functions and branching         #
              function are supported.                                #
                                                                     #
######################################################################

Revision History:

* 1/8/2018:  Initial README (F. S. Tsung)
* 11/19/2020: Adding documentation for input deck

