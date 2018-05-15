!-----------------------------------------------------------------------
! 2-1/2D Electromagnetic MPI/OpenMP PIC code
! written by Viktor K. Decyk, UCLA
! Non Squared Cells implemented by Michael J. Touati, UCLA, 2016, March
! -> New subroutine       : GRID
! -> Modified subroutines : mpdcomp2, mppois2_init, mpdistr2h, mpdblkp2
!                           mpmovin2, mpcheck2, wmpdjpost2, ompmove2,
!                           mppost2, mpcuperp2, mpibpois2, wmpgmjpost2,
!                           mpdcuperp2, mpmaxwel2, mppois2, wmpbpush2
! Simple diagnostic tool for code testing added by Michael J. Touati, UCLA, 2016, March
! -> New module diag.f90 with new subroutines : DIAG_INIT, DIAG_REAL_MOMENTS, 
!                                               DIAG_FOURIER_MOMENTS, DIAG_FOURIER_FIELDS, 
!                                               DIAG_REAL_FIELDS, DIAG_PART, DIAG_ENERGY
! Moving ions implemented by Michael J. Touati, UCLA, 2016, July 
! Input deck implemented by Michael J. Touati, UCLA, 2016, August
! -> New module :: input.f90
! Finite Time and Arbitrary Order Spatial Finite Difference Scheme implemented
! by Michael J. Touati, UCLA, 2016, October
! -> New module FTFDfield2.f90
      program upic_emma
      use upic_m_parameters
      use modmpinit2
      use modmpbpush2
      use modmpcurd2
      use modmpfield2
      use modmpfieldpml2
      use mppmod2
      use omplib
      use ompplib2
      use diag
      use antenna
      use input
! HDF5/OpenPMD modules
      use parallel_class
      use hdf5io_class
! 
      use mpi

      
      implicit none
! npx/npy = number of macro electrons distributed in x/y direction.
! de = number of electrons per macro electron
      real            :: de
      integer         :: npx
      integer         :: npy
! Cells dimensions in x/y directions in spatial unit choosen by the user ( M Touati)
      real, dimension(2) :: delta
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
      real :: tend
      real :: dt
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
	  real :: vtx, vty, vx0, vy0
! vtz/vz0 = thermal/drift velocity of electrons in z direction
!           which is not spatially resolved in this 2D-3P code
      real :: vtz, vz0
! ci = reciprocal of velocity of light.
	  real :: ci
! idimp = dimension of phase space = 5
! ipbc = particle boundary condition: 1 = periodic 
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 5
      integer :: ipbc = 1
      integer :: relativity
! omx/omy/omz = external magnetic field in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! idps = number of partition boundaries
      integer, parameter :: idps = 2
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
! wsrc  = time integrated source of electrical energy in the simulation box
! wpsrc = time integrated source of electrical energy in the simulation box per MPI node
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0, wsrc = 0., wpsrc = 0., wl =0.0
! sorting tiles, should be less than or equal to 32
      integer :: mx = 16, my = 16
! fraction of extra particles needed for particle management
      real :: xtras = 0.5 !xtras = 0.2
! list = (true,false) = list of particles leaving tiles found in push
      logical :: list = .true.
! declare scalars for standard code
      integer :: n
      integer :: nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy
      integer :: mx1, ntime, nloop, isign, ierr
      real :: affp, dth, omt, wp0
      double precision :: np
!
! declare scalars for MPI code
      integer :: nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn
      integer :: nyp, noff, npp, nps, myp1, mxyp1
      integer :: nterf
!
! declare scalars for OpenMP code
      integer :: nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc
! nvpp = number of shared memory nodes (0=default)
      integer :: nvpp
!
! declare arrays for standard code
! part = macro electron array
      real, dimension(:,:), allocatable :: part
! qe = electron charge density with guard cells
      real, dimension(:,:), allocatable :: qe
! cue = electron current density with guard cells
	  real, dimension(:,:,:), allocatable :: cue
! fxyze/bxyze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:,:), allocatable :: fxyze, bxyze
! pixyze = smoothed Poynting vector with guard cells 
!        = (c / 4 pi) * fxyze x bxyze in non-normalized units
      real, dimension(:,:,:), allocatable :: pixyze
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:), allocatable :: exyz, bxyz
! exyz/bxyz = laser electric/magnetic field in fourier space
      complex, dimension(:,:,:), allocatable :: exyzl, bxyzl
! qt = scalar charge density field array in fourier space
      complex, dimension(:,:), allocatable :: qt
! cut = vector current density field array in fourier space
! fxyt/bxyt = vector electric/magnetic field in fourier space
      complex, dimension(:,:,:), allocatable :: cut, fxyt, bxyt
! ffc = form factor array for poisson solver
      complex, dimension(:,:), allocatable :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), allocatable :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), allocatable :: sct
      double precision, dimension(15) :: wtot, work
! arrays required for darwin initial fields
      real, dimension(:,:,:), allocatable :: amu
      complex, dimension(:,:,:), allocatable :: dcut, amut
! x, y, kx, ky = spatial grid points/spatial frequency gridpoint arrays ( M Touati )
! yp, kpx = spatial grid points (y-dir) and spatial frequency (x-dir) arrays 
!           for each MPI node ( M Touati )
! moff = lowermost global spatial freq. array kpx in particle partition ( M Touati )
      real, dimension(:), allocatable :: x, y, kx, ky, yp, kpx
      integer :: moff
!
! declare arrays for MPI code
! edges(1:2) = lower:upper y boundaries of particle partition
      real, dimension(:), allocatable  :: edges
!
! declare arrays for OpenMP code
! ppart = tiled macro electron array
      real, dimension(:,:,:), allocatable :: ppart
! kpic = number of macro electrons in each tile
      integer, dimension(:), allocatable :: kpic
! ncl = number of macro electrons departing tile in each direction
      integer, dimension(:,:), allocatable :: ncl
! ihole = location/destination of each macro electron departing tile
      integer, dimension(:,:,:), allocatable :: ihole
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime, ltime
      real :: tinit = 0.0, tloop = 0.0
      real :: tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0, tmov = 0.0
      real :: tfmov = 0.0
      real, dimension(2) :: tfft = 0.0
      double precision :: dtime
!
! declare and initialize macro ions quantities 
! movion = .true. allows for simulating a fully ionized plasma with one moving ion species
      logical :: movion
! Z_ato
! Z_ion = ionization state (fully ionized ->   ionization state Z_ion
!                                                 = Atomic number Z)
! A_ion = Standard atomic weight (Atomic_mass_unit is in unit of the electron mass!)
      real, parameter    :: Atomic_mass_unit = 1836.2
      real               :: Z_ato
      real               :: A_ion
      integer            :: Z_ion
! me = macro electron mass ( = - qme in such a way to conserve the correct ratio 
!                              charge/mass)
! qmi = macro ions charge (= - qme because we choose the same number of macro ions 
!                            than macro electrons and we want to respect electrical
!                            neutrality)
! me_real  = real electron mass in units of electron mass (me_real = 1.)
! me       = macro electron mass in units of me_real
! qme_real = real electron charge in units of e
! qme      = macro electron charge in units of e
! qbme     = macro (or real) electron charge / mass ratio
! di       = number of real ion per macro ion
!            such that there is the same number of macro electrons and macro ions
!            in the simulation
! mi_real  = real ion mass in units of me_real
! mi       = macro ion mass in units of me_real
! qmi_real = real ion charge in units of e
! qmi      = macro ion charge in units of e
! qbmi     = macro (or real) ion charge / mass ration
      real :: me_real, me, qme_real, qme, qbme
      real :: di, mi_real, mi, qmi_real, qmi, qbmi
      real :: vtxi, vtyi, vtzi
      real :: vx0i, vy0i, vz0i
      real :: wki = 0.
! parti = macro ion array
      real, dimension(:,:), allocatable :: parti
! declare scalars for MPI code
      integer :: nppi, npsi
! declare scalars for OpenMP code
      integer :: nppmxi, nppmx0i, nbmaxpi, ntmaxpi, npbmxi
! pparti = tiled macro ion array
      real, dimension(:,:,:), allocatable :: pparti
! kpici = number of macro ions in each tile
      integer, dimension(:), allocatable :: kpici
! ncli = number of macro ions departing tile in each direction
      integer, dimension(:,:), allocatable :: ncli
! iholei = location/destination of each macro ion departing tile
      integer, dimension(:,:,:), allocatable :: iholei
! qi = ion charge density with guard cells
      real, dimension(:,:), allocatable :: qi
! cui = ion current density with guard cells
	  real, dimension(:,:,:), allocatable :: cui 
! qtot = total charge density with guard cells
      real, dimension(:,:), allocatable :: qtot
! cutot = total current density with guard cells
	  real, dimension(:,:,:), allocatable :: cutot
! declare PML variables : 
! sigma = 2D PML conductivity. sigma(1,i,j) = sigma_x(x(i)) and sigma(2,i,j) = sigma_y(yp(j))
! vpml_exyz, vpml_bxyz = PML EM fields in the Fourier space
! For example, in 2D : vpml_exyz(1,q,p) = IDFT{E_xy^{i,j}}^{q,p},
!                      vpml_exyz(2,q,p) = IDFT{E_yx^{i,j}}^{q,p},
!                      vpml_exyz(3,q,p) = IDFT{E_zx^{i,j}}^{q,p} and 
!                      vpml_exyz(4,q,p) = IDFT{E_zy^{i,j}}^{q,p}.
! wmppml/wfppml     = instantaneous magnetic field/transverse electric field 
!                     energy in the pml layers per MPI node
! wmpml/wfpml/wtpml = instantaneous magnetic field/transverse electric field/total 
!                     energy in the pml layers
! wmsim/wfsim/wtsim = instantaneous magnetic field/transverse electric field/total 
!                     energy in the effective simulation box
! wmpdump/wfpdump      = time integrated dumped magnetic/transverse electric
!                        energy in the pml layers per MPI node
! wmdump/wfdump/wtdump = time integrated dumped magnetic/transverse electric/total 
!                        energy in the pml layers
! wesc  = instantaneous escaping/entering energy from/to the effective simulation box / pml layers
! wpesc = same but per MPI node
      real, dimension(:,:,:), allocatable    :: sigma
      complex, dimension(:,:,:), allocatable :: vpml_exyz, vpml_bxyz
      complex, dimension(:,:), allocatable   :: vpml_exyz_corr, vpml_bxyz_corr
      real :: wfppml = 0.0, wmppml = 0.0
      real :: wfpml = 0.0, wmpml = 0.0, wtpml = 0.0
      real :: wfsim = 0.0, wmsim = 0.0, wtsim = 0.0
      real :: wfpdump = 0., wmpdump = 0.
      real :: wfdump = 0., wmdump = 0., wtdump = 0.
      real :: wescx  = 0., wescy  = 0., wesc  = 0.
      real :: wpescx = 0., wpescy = 0., wpesc = 0.
! declare variables for emulating a Finite Time Finite Difference Yee Scheme
! stencil             = order of acuracy / 2
! exyz_corr/bxyz_corr = array containing the 4 particular modes of transverse 
!                       electric/magnetic field in fourier space. For example :
!                       exyz_corr(1,1) = Ex(kx=0,ky=0),
!                       exyz_corr(2,1) = Ex(kx=pi/dx,ky=0),
!                       exyz_corr(3,2) = Ey(kx=0,ky=pi/dy) and 
!                       exyz_corr(4,3) = Ez(kx=pi/dx,ky=pi/dy).
	  integer :: stencil
      complex, dimension(:,:), allocatable :: exyz_corr, bxyz_corr
! declare and initialize diagnostic tools data ( M Touati)
! phtime  = physical time in / omega_p
! dphtime = time interval in / omega_p at which the sim. results are stored
! Ndiag   = only one macro-particle orbit out of Ndiag is stored 
      real    :: phtime
      real    :: dphtime
      integer :: Ndiag
      logical :: store_cond
      real, dimension(:,:), allocatable :: part_temp
! data needed to diagnose physical quantities in Fourier space
      integer :: modesx, modesy, modesxpd, modesyd
      complex, dimension(:,:), allocatable :: qf
      complex, dimension(:,:,:), allocatable :: cuf, exyzf, bxyzf
! tdiag = timing data of the diagnostic tools
      real :: tdiag = 0.0
! volume = volume of a spatial cell in (Delta^3)
      real :: volume
! HDF5 output
! sfield => 2d array for HDF5 dump
      real*4, allocatable  :: sfield(:,:)
! p, pp => parallel configuration needed for HDF5 dumps
      type(parallel), target :: p
      class(parallel), pointer :: pp => null()
      type(hdf5file) :: file
      integer ix,iy 
! HDF5
      
!
! start timing initialization
      call dtimer(dtime,itime,-1)
! read the input deck ( M Touati)
      call read_input_deck()
      relativity = 0
      if (relativistic) relativity = 1
      movion     = moving_ions
      delta(1)   = delta_x
      delta(2)   = delta_y
      tend       = t_end
      Z_ato      = atomic_number
      A_ion      = atomic_weight
      Z_ion      = ionization_state
      de         = den_me
      npx        = int(sqrt(Npic)*(2**indx))
      npy        = int(sqrt(Npic)*(2**indy)) 
      dphtime    = Delta_t_diag
      store_cond = .false.
      select case (units)
      	case ('ES')
      		ci = real(param)
      		vtx = ptx
      		vty = pty
      		vtz = ptz
      		if (((vtx**2.)+(vty**2.)+(vtz**2.)) .ne. 3. ) then
      			vtx = vtx * sqrt(3. / ((vtx**2.)+(vty**2.)+(vtz**2.)))
      			vty = vty * sqrt(3. / ((vtx**2.)+(vty**2.)+(vtz**2.)))
      			vtz = vtz * sqrt(3. / ((vtx**2.)+(vty**2.)+(vtz**2.)))
      		else
      			vtx = vtx
      			vty = vty
      			vtz = vtz
      		end if
      		vx0 = px0
      		vy0 = py0
      		vz0 = pz0
      	case ('EM')
      		ci  = 1.
      		vtx = ptx
      		vty = pty
      		vtz = ptz
      		if (((vtx**2.)+(vty**2.)+(vtz**2.)) .ne. 3. ) then
      			vtx = vtx * param * sqrt(3. / ((vtx**2.)+(vty**2.)+(vtz**2.)))
      			vty = vty * param * sqrt(3. / ((vtx**2.)+(vty**2.)+(vtz**2.)))
      			vtz = vtz * param * sqrt(3. / ((vtx**2.)+(vty**2.)+(vtz**2.)))
      		else
      			vtx = vtx * param
      			vty = vty * param
      			vtz = vtz * param
      		end if
      		vx0 = px0
      		vy0 = py0
      		vz0 = pz0
      end select
!     
      Ndiag = (4**indx) - 1
! initialize scalars for standard code
! np = total number of particles in simulation
      np =  dble(npx)*dble(npy)
! nx/ny = number of grid points in x/y direction
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = ndim*nxe
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny)
! mx1 = number of tiles in x direction
      mx1 = (nx - 1)/mx + 1
      affp = dble(nx)*dble(ny)/np
! Modification of affp for non square cells ( M Touati )
      volume = delta(1) * delta(2) * 1.
      affp = affp*volume
! Modification of affp for higher than 1 electron per macro electron ( M Touati )
	  affp = affp / de
! Real and macro particles masses and charges ( M. Touati )
	  qme_real = - 1.
	  qme      = de * qme_real
	  me_real  = 1.
	  me       = de * me_real
	  qbme     = qme / me
	  di       = de / Z_ion
	  qmi_real = Z_ion
	  qmi      = di * qmi_real
	  mi_real  = A_ion * Atomic_mass_unit
	  mi       = di * mi_real
      qbmi     = qmi / mi
      if (movion) then
        	vtxi = vtx / sqrt(mi_real) ! vtx,y,z are thermal velocities if relativity == 0
         	vtyi = vty / sqrt(mi_real) ! ( v_Ti = (kB Ti / mi)^1/2 = v_Te * (me/mi)^1/2 )
         	vtzi = vtz / sqrt(mi_real) 
         	vx0i = vx0
         	vy0i = vy0
         	vz0i = vz0
      else
         vtxi = 0.
         vtyi = 0. 
         vtzi = 0.
         vx0i = 0.
         vy0i = 0.
         vz0i = 0.
      endif
!
! nvp = number of distributed memory nodes
! initialize for distributed memory parallel processing
      call PPINIT2(idproc,nvp)
      kstrt = idproc + 1
! check if too many processors
      if (nvp > ny) then
         if (kstrt==1) then
         write (*,*) 'Too many processors requested: ny, nvp=', ny, nvp
         endif
         call PPEXIT()
         stop
      endif
! initialize for shared memory parallel processing
	  nvpp = N_threads
      call INIT_OMP(kstrt, nvp, nvpp)
!
! initialize data for MPI code
      allocate(edges(idps))
! calculate partition variables: edges, nyp, noff, nypmx
! edges(1:2) = lower:upper boundary of particle partition
! nyp = number of primary (complete) gridpoints in particle partition
! noff = lowermost global gridpoint in particle partition
! nypmx = maximum size of particle partition, including guard cells
! nypmn = minimum value of nyp
      call mpdcomp2(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,delta)
!
! check for unimplemented features
      if ((list).and.(ipbc.ne.1)) then
         if (kstrt==1) then
            write (*,*) 'ipbc /= 1 and list = .true. not yet supported'
            write (*,*) 'list reset to .false.'
            list = .false.
         endif
      endif
!
! initialize additional scalars for MPI code
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvp + 1
! npmax = maximum number of electrons in each partition
      npmax = (np/nvp) * part_factor
! myp1 = number of tiles in y direction
      myp1 = (nyp - 1)/my + 1; mxyp1 = mx1*myp1
! nterf = number of shifts required by field manager (0=search)
      nterf = 0
!
	  if (kstrt==1) then
      	print*,'=================='
      	print*,'| Initialization |'
      	print*,'=================='
      end if
!
! allocate and initialize data for standard code
      allocate(part(idimp,npmax))
      allocate(qe(nxe,nypmx),fxyze(ndim,nxe,nypmx))
      allocate(cue(ndim,nxe,nypmx),bxyze(ndim,nxe,nypmx))
      allocate(pixyze(ndim,nxe,nypmx))
      allocate(exyz(ndim,nye,kxp),bxyz(ndim,nye,kxp))
      allocate(exyz_corr(4,ndim), bxyz_corr(4,ndim) )
      allocate(qt(nye,kxp),fxyt(ndim,nye,kxp))
      allocate(cut(ndim,nye,kxp),bxyt(ndim,nye,kxp))
      allocate(ffc(nyh,kxp),mixup(nxhy),sct(nxyh))
      allocate(kpic(mxyp1),kpici(mxyp1))
      allocate(sfield(nx,nyp))
      if (laserpulse) then
	      allocate(exyzl(ndim,nye,kxp),bxyzl(ndim,nye,kxp))
      	  exyzl = 0.
          bxyzl = 0.
      end if
! allocate and initialize data for ions ( M. Touati )
      if (movion) then 
      	allocate(parti(idimp,npmax))
      	allocate(qi(nxe,nypmx),cui(ndim,nxe,nypmx))
      end if
      allocate(qtot(nxe,nypmx),cutot(ndim,nxe,nypmx))
! Prepare spatial gridpoint and spatial frequency arrays ( M. Touati )
      allocate(x(nxe),y(nye),kx(nxe),ky(nye),yp(nypmx),kpx(kxp))
      moff = kxp * (kstrt - 1)
      call GRID(x,y,kx,ky,yp,kpx,nxe,nye,nypmx,kxp,noff,moff,nx,ny,delta)
!
! open diagnostic files to write results during the simulation ( M. Touati )
	  call DIAG_INIT(movion, kstrt, N_rho, N_rho_ele, N_rho_ion,&
           			 N_jx, N_jy, N_jz,&  
           			 N_jex, N_jey, N_jez,&
           			 N_jix, N_jiy, N_jiz,&
           			 N_rho_Fourier, N_rho_ele_Fourier, N_rho_ion_Fourier,&
           			 N_jx_Fourier, N_jy_Fourier,N_jz_Fourier,&
           			 N_histogram_ele, N_histogram_ion,N_jez_Fourier,&
           		 	 N_jix_Fourier, N_jiy_Fourier,N_jiz_Fourier,&
           			 N_rho_Fourier_arg, N_rho_ele_Fourier_arg, N_rho_ion_Fourier_arg,&
				     N_jx_Fourier_arg, N_jy_Fourier_arg, N_jz_Fourier_arg,&
				     N_jex_Fourier_arg, N_jey_Fourier_arg, N_jez_Fourier_arg,&
           			 N_jix_Fourier_arg, N_jiy_Fourier_arg, N_jiz_Fourier_arg,&
           			 N_Ex_Fourier, N_Ey_Fourier,N_Ez_Fourier,&
           			 N_Bx_Fourier, N_By_Fourier, N_Bz_Fourier,&
           			 N_Ex_Fourier_arg, N_Ey_Fourier_arg, N_Ez_Fourier_arg,&
           			 N_Bx_Fourier_arg, N_By_Fourier_arg, N_Bz_Fourier_arg,&
           			 N_Ex, N_Ey, N_Ez, N_Bx, N_By, N_Bz,N_Pix, N_Piy, N_Piz,&
           			 N_ele_px_vs_xy, N_ele_py_vs_xy, N_ele_pz_vs_xy,&
           			 N_ion_px_vs_xy, N_ion_py_vs_xy, N_ion_pz_vs_xy,&
           			 N_es_fields, N_el_fields, N_ma_fields,&
           			 N_kinetic_energy, N_ele_kinetic_energy, N_ion_kinetic_energy,&
           			 N_el_fields_pml, N_ma_fields_pml, N_fields_src,&
           			 N_fields_esc, N_el_dumped, N_ma_dumped)
! Initialize Diagnostics for HDF5/OpenPMD
      call p%new()
      pp => p
! OpenPMD
!	
! Variables needed to unpack array containing Fourier quantities ( M. Touati )
      modesx = nx/4; modesy = ny/4
      modesxpd = min(modesx,kxp); modesyd = min(2*modesy-1,ny)
      if (modesx==(nx/2+1)) modesxpd = modesxpd + 1
      allocate(qf(modesyd,modesxpd))
      allocate(cuf(ndim,modesyd,modesxpd))
      allocate(exyzf(ndim,modesyd,modesxpd))
      allocate(bxyzf(ndim,modesyd,modesxpd))
! Time step : dt < 2 / (k_max c) with k_max = pi * sqrt( (1/delta)^2 + (1/delta)^2 )
      dt = cfl*0.45*ci ! (sqrt(2)/pi = 0.45...)      
! Modification of the CFL condition for non square cells ( M Touati )
      dt = dt / sqrt(0.5*((1./delta(1)**2.)+(1./delta(2)**2.)))
! Allocate pml arrays and set up coefficients needed to computed the PML equations and
! modification of the CFL condition for solving the PML equations
	  if ((BCx == 0) .or. (BCy == 0)) then
	  	  call setup_pml(sigma,vpml_exyz,vpml_bxyz,FTFD,vpml_exyz_corr,vpml_bxyz_corr,&
	  	                 ax,ay,dt,x,yp,delta,L_PML,n_cond_PML,ci,&
	                     nx,ny,nxe,nye,BCx,BCy,nyp,nypmx,kxp,nvp,kstrt)
	  end if  
!
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      dth = 0.0
!
! prepare fft tables
      call mpfft2_init(mixup,sct,indx,indy)
! Emulation of a Finite Time Finite Difference Yee Scheme
      stencil = 0
      if (FTFD .ge. 2) then
      	stencil = floor(FTFD/2.)
      	call setup_bracket(stencil,nx,ny,delta,kx,ky)
      end if
! calculate form factors
      call mppois2_init(ffc,ax,ay,affp,nx,ny,kstrt,kx,ky)
! initialize electrons
      nps = 1
      npp = 0
      call mpdistr2h(relativity,me_real,ci,part,edges,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,&
     		        &nx,ny,ipbc,x,y,density_x,density_y,ierr)
! check for macro electrons initialization error
      if (ierr /= 0) then
        if (kstrt==1) then
            write (*,*) 'particle initialization error: ierr=', ierr
         endif
         call PPEXIT()
         stop
      endif
! initialize ions if moving ions
      if (movion) then 
      	npsi = 1
        nppi = 0
      	call mpdistr2h(relativity,mi_real,ci,parti,edges,nppi,npsi,vtxi,vtyi,vtzi,vx0i,vy0i,vz0i,npx,npy,&
                      &nx,ny,ipbc,x,y,density_x,density_y,ierr)
! check for macro ions initialization error
     	if (ierr /= 0) then
        	if (kstrt==1) then
           		write (*,*) 'particle initialization error: ierr=', ierr
        	end if
         	call PPEXIT()
         	stop
      	endif
      end if
!
! find number of macro electrons in each of mx, my tiles: updates kpic, nppmx
      irc = 0
      call mpdblkp2(part,kpic,npp,noff,nppmx,mx,my,mx1,delta,irc)
! allocate vector macro electrons data
      nppmx0 = (1.0 + xtras)*nppmx
      ntmaxp = xtras*nppmx
      npbmx = xtras*nppmx
      nbmaxp = 0.25*mx1*npbmx
      allocate(ppart(idimp,nppmx0,mxyp1))
      allocate(ncl(8,mxyp1),ihole(2,ntmaxp+1,mxyp1))
!
! copy ordered macro electrons data for OpenMP
      call mpmovin2(part,ppart,kpic,npp,noff,mx,my,mx1,delta,irc)
! sanity check for macro electrons
      call mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,nypmx,x,yp,irc)
!
! find number of macro ions in each of mx, my tiles: updates kpici, nppmxi
	  if (movion) then
      	call mpdblkp2(parti,kpici,nppi,noff,nppmxi,mx,my,mx1,delta,irc)
! allocate vector macro ions data
      	nppmx0i = (1.0 + xtras)*nppmxi
      	ntmaxpi = xtras*nppmxi
      	npbmxi  = xtras*nppmxi
      	nbmaxpi = 0.25*mx1*npbmxi
      	allocate(pparti(idimp,nppmx0i,mxyp1))
      	allocate(ncli(8,mxyp1),iholei(2,ntmaxpi+1,mxyp1))
!
! copy ordered macro ions data for OpenMP
      	call mpmovin2(parti,pparti,kpici,nppi,noff,mx,my,mx1,delta,irc)
! sanity check for macro ions
      	call mpcheck2(pparti,kpici,noff,nyp,nx,mx,my,mx1,nypmx,x,yp,irc)
	  end if
!
	  if (plasma) then
! Particle loading diagnostics for checking particles momentum/velocity initialization
! in some particular spatial cells
 	  	allocate(part_temp(idimp,nppmx0))
 	  	part_temp(:,:) =  ppart(:,:,1)
      	call DIAG_HISTOGRAM(ci,vtx,vty,vtz,N_histogram_ele,nppmx0,idimp,nppmx0,nxe,nypmx,&
                            part_temp,x,yp,delta,tdiag)
      	deallocate(part_temp)
      	if (movion) then
      		allocate(part_temp(idimp,nppmx0))
			part_temp(:,:) = pparti(:,:,1)      
			call DIAG_HISTOGRAM(ci,vtxi,vtyi,vtzi,N_histogram_ion,nppmx0i,idimp,nppmx0i,&
		                        nxe,nypmx,part_temp,x,yp,delta,tdiag)
      		deallocate(part_temp)
      	end if
      end if
!
! initialize transverse electromagnetic fields
      exyz = cmplx(0.0,0.0)
      bxyz = cmplx(0.0,0.0)
! set magnitude of external magnetic field
      omt = sqrt(omx*omx + omy*omy + omz*omz)
!
      if (dt > 0.45*ci) then
         if (kstrt==1) then
            write (*,*) 'Warning: Courant condition may be exceeded!'
         endif
      endif
!
! initialization time
      	call dtimer(dtime,itime,1)
      	tinit = tinit + real(dtime)
!
! start timing loop
      call dtimer(dtime,ltime,-1)   
!
! * * * start main iteration loop * * *
!
! Initialization of physical time
      phtime = 0.
      do n = 1, nloop 
      ntime = n - 1
! Simulation results are saved every dphtime :
	  store_cond = (mod(phtime,dphtime) < dt).or.(ntime == nloop - 1)
!
! Just for information during the simulation run
	  if ( (kstrt==1) .and. (store_cond))then
		write (*,*) '============================='
		write (*,*) '| time = ', phtime,        ' |'
		write (*,*) '============================='
	  end if
!
      if (plasma) then
! deposit current with OpenMP:
! updates ppart and cue, and possibly ncl, ihole, irc
      	call dtimer(dtime,itime,-1)
      	call dtimer(dtime,itime,1)
      	tdjpost = tdjpost + real(dtime)
		cue = 0.
		call wmpdjpost2(ppart,cue,kpic,ncl,ihole,noff,nyp,qme,dth,ci,     &
		                &tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,list,x,y,delta,irc)
!
! updates pparti and cui, and possibly ncli, iholei, irc ( M. Touati )
	  	if (movion) then
	  	    cui = 0.
			call wmpdjpost2(pparti,cui,kpici,ncli,iholei,noff,nyp,qmi,dth,ci,     &
		                  &tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,list,x,y,delta,irc)
	  	end if
!
! reorder particles by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
		 call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,tsort,tmov,kstrt,nvp, &
		              &nx,ny,mx,my,npbmx,nbmaxp,mx1,list,x,y,irc)
! updates: pparti, kpici, and irc and possibly ncli and iholei ( M. Touati )
	  	if (movion) then
			call ompmove2(pparti,kpici,ncli,iholei,noff,nyp,tsort,tmov,kstrt,nvp, &
			             &nx,ny,mx,my,npbmx,nbmaxp,mx1,list,x,y,irc)
	  	end if
!
! deposit electron charge with OpenMP: updates qe
      	call dtimer(dtime,itime,-1)
      	call dtimer(dtime,itime,1)
      	tdpost = tdpost + real(dtime)
      	qe = 0.
     	call mppost2(ppart,qe,kpic,noff,qme,tdpost,mx,my,mx1,x,yp,delta)
!
! add guard cells with OpenMP: updates cue, qe
     	call wmpnacguard2(cue,nyp,tguard,nx,kstrt,nvp)
     	call wmpaguard2(qe,nyp,tguard,nx,kstrt,nvp)
! deposit ion charge with OpenMP: updates qi ( M. Touati )
	  	if (movion) then
	  	    qi = 0.
	  		call mppost2(pparti,qi,kpici,noff,qmi,tdpost,mx,my,mx1,x,yp,delta)
! add guard cells with OpenMP: updates cui, qi ( M. Touati )
	    	call wmpnacguard2(cui,nyp,tguard,nx,kstrt,nvp)
      		call wmpaguard2(qi,nyp,tguard,nx,kstrt,nvp)
	  	end if
! Diagnostic of the two first particle distribution moments in real space ( M Touati )
	  cutot = cue
      qtot  = qe
      if (movion) then
        cutot = cutot + cui
        qtot  = qtot  + qi
      end if
      	if (store_cond) then
        	call DIAG_REAL_MOMENTS(N_rho, N_jx, N_jy , N_jz, nyp, nx, nxe, nypmx,&
                               de, phtime, yp, x, qtot, cutot, tdiag)
        	if (movion) then
        		call DIAG_REAL_MOMENTS(N_rho_ele, N_jex, N_jey , N_jez, nyp, nx, nxe, nypmx,&
                               de, phtime, yp, x, qe, cue, tdiag)
            	call DIAG_REAL_MOMENTS(N_rho_ion, N_jix, N_jiy , N_jiz, nyp, nx, nxe, nypmx,&
                               de, phtime, yp, x, qi, cui, tdiag)
        	end if
      	end if
	  end if
!
! transform charge to fourier space with OpenMP:
! updates qt, nterf, and ierr, modifies qtot
      isign = -1
      call wmpfft2r(qtot,qt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,indy,&
     &kstrt,nvp,kyp,ny,nterf,ierr)
!
! transform current to fourier space with OpenMP:
! updates cut, nterf, and ierr, modifies cut
      isign = -1
      call wmpfft2rn(cutot,cut,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,  &
     &indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! Diagnostic of the two first particle distribution moments in Fourier space ( M Touati )
      if (store_cond) then
        call mprdvmodes2(cut,cuf,tdiag,nx,ny,modesx,modesy,kstrt)
        call mprdmodes2(qt,qf,tdiag,nx,ny,modesx,modesy,kstrt)
        call DIAG_FOURIER_MOMENTS(N_rho_Fourier, N_rho_Fourier_arg,&
                                  N_jx_Fourier, N_jx_Fourier_arg,&
                                  N_jy_Fourier, N_jy_Fourier_arg,&
                                  N_jz_Fourier, N_jz_Fourier_arg,&
                                  modesyd, modesxpd, kxp, nye,&
                                  de, phtime, kpx, ky, qf, cuf, tdiag)
      end if
!
! take transverse part of current with OpenMP: updates cut
      call mpcuperp2(cut,FTFD,tfield,nx,ny,kstrt,kx,ky)
!
! calculate electromagnetic fields in fourier space with OpenMP:
! updates exyz, bxyz, wf, wm
      if (ntime==0) then
! initialize electromagnetic fields from darwin fields
! calculate initial darwin magnetic field
         call mpibpois2(cut,bxyz,ffc,ci,wm,tfield,nx,ny,kstrt,kx,ky)
         wf = 0.0
! calculate initial darwin electric field
         allocate(amu(4,nxe,nypmx),amut(4,nye,kxp),dcut(ndim,nye,kxp))
         amu = 0.0
         call wmpgmjpost2(ppart,amu,kpic,noff,qme,ci,tdjpost,mx,my,mx1, &
     &relativity,x,yp,delta)
         call wmpnacguard2(amu,nyp,tguard,nx,kstrt,nvp)
         isign = -1
         call wmpfft2rn(amu,amut,noff,nyp,isign,mixup,sct,tfft,tfmov,   &
     &indx,indy,kstrt,nvp,kyp,ny,nterf,ierr)
         deallocate(amu)
         call mpdcuperp2(dcut,amut,tfield,nx,ny,kstrt,kx,ky)
         deallocate(amut)
         call mpetfield2(dcut,exyz,ffc,affp,ci,wf,tfield,nx,ny,kstrt)
         deallocate(dcut)
         dth = 0.5*dt
! initialize pml fields in case of absorbing boundary conditions
         if (((BCx == 0) .or. (BCy == 0))) then
			 call mpmaxwel2pml_init(vpml_exyz,vpml_bxyz,vpml_exyz_corr,vpml_bxyz_corr,&
	  	                            exyz,bxyz,nye,kxp,nyh)
	  	 else if (FTFD > 0) then
	  	 	call mpmaxwel2_init(exyz,bxyz,exyz_corr,bxyz_corr,nye,kxp,nyh)
         end if
! update electromagnetic fields
      else
! absorbing boundary conditions case
      	 if (((BCx == 0) .or. (BCy == 0))) then
      	 	call mpmaxwel2pml(PML_scheme,ndim,indx,indy,nx,nxe,ny,nye,nyp,nypmx,nxhy,nxyh,nyh,kxp,kyp,&
	                          nvp,kstrt,noff,BCx,BCy,L_PML,ax,ay,affp,ci,dt,kx,ky,x,yp,delta,sigma,&
	                          ffc,cut,nterf,ierr,mixup,wm,wf,tfield,tfmov,tguard,&
	                          tfft,vpml_exyz,vpml_bxyz,FTFD,vpml_exyz_corr,vpml_bxyz_corr,&
	  	                      sct,exyz,bxyz,wescx,wescy,wfppml,wmppml,wfpdump,wmpdump)
       	 else
! periodic boundary conditions case
         	call mpmaxwel2(exyz,bxyz,FTFD,exyz_corr,bxyz_corr,cut,ffc,&
         	               affp,ci,dt,wf,wm,tfield,nx,ny,kstrt,kx,ky,ax,ay)
     	 endif
      endif
!
! Modification of cells volume in case of Non Square Cells
      we = we * volume
      wf = wf * volume
      wm = wm * volume
!
! calculate longitudinal force/charge in fourier space with OpenMP:
! updates fxyt, we
      call mppois2(qt,fxyt,ffc,we,tfield,nx,ny,kstrt,kx,ky)
!
! Add laser pulse fields in real space
	  if ((phtime >= tlaunch) .and. laserpulse) then
	   !  isign = 1
!         call wmpfft2rn(fxyze,exyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
!                     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!         call wmpfft2rn(bxyze,bxyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
!                      &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
	  	call LaunchLaser(fxyze,bxyze,wl,x,yp,affp,delta,ci,dt,phtime,L_PML,theta,polardir,tlaunch,FWHMt,&
	                     FWHMs,xfocal,yfocal,omega0,a0,BCx,BCy,nx,nyp,nxe,nypmx,ntime,propdir,shape,laserpulse)
	  	! copy guard cells with OpenMP: updates fxyze, bxyze
      	call wmpncguard2(fxyze,nyp,tguard,nx,kstrt,nvp)
      	call wmpncguard2(bxyze,nyp,tguard,nx,kstrt,nvp)
	  	! transform E-M fields to Fourier space with OpenMP:
	  	isign = -1
      	call wmpfft2rn(fxyze,exyzl,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,  &
     	    &indy,kstrt,nvp,kyp,ny,nterf,ierr)
     	call wmpfft2rn(bxyze,bxyzl,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,  &
     	    &indy,kstrt,nvp,kyp,ny,nterf,ierr)
     	! Filter laser fields and add it to the electromagnetic fields with OpenMP:
	    isign = 1
      	call mpemfield2(exyz,exyzl,FTFD,ffc,isign,kx,ky,ax,ay,tfield,nx,ny,kstrt)
      	call mpemfield2(bxyz,bxyzl,FTFD,ffc,isign,kx,ky,ax,ay,tfield,nx,ny,kstrt)
! 	  	exyz = exyz + exyzl
! 	  	bxyz = bxyz + bxyzl
	  end if
!
! add longitudinal and transverse electric fields with OpenMP:
! updates fxyt
      isign = 1
      call mpemfield2(fxyt,exyz,FTFD,ffc,isign,kx,ky,ax,ay,tfield,nx,ny,kstrt)
! copy magnetic field with OpenMP: updates bxyt
      isign = -1
      call mpemfield2(bxyt,bxyz,FTFD,ffc,isign,kx,ky,ax,ay,tfield,nx,ny,kstrt)
! Diagnostic of electromagnetic fields in Fourier space ( M Touati )
      if (store_cond) then
        call mprdvmodes2(exyz,exyzf,tdiag,nx,ny,modesx,modesy,kstrt)
        call mprdvmodes2(bxyz,bxyzf,tdiag,nx,ny,modesx,modesy,kstrt)
        call DIAG_FOURIER_FIELDS(N_Ex_Fourier, N_Ey_Fourier, N_Ez_Fourier,&
                                 N_Bx_Fourier, N_By_Fourier, N_Bz_Fourier,&
                                 N_Ex_Fourier_arg, N_Ey_Fourier_arg, N_Ez_Fourier_arg,&
                                 N_Bx_Fourier_arg, N_By_Fourier_arg, N_Bz_Fourier_arg,&
                                 modesxpd, modesyd, kxp, nye,&
                                 de, phtime, kpx, ky, exyzf, bxyzf, tdiag)
      end if
!
! transform force to real space with OpenMP:
! updates fxyze, nterf, and ierr, modifies fxyt
      isign = 1
      call wmpfft2rn(fxyze,fxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! transform magnetic field to real space with OpenMP:
! updates bxyze, nterf, and ierr, modifies bxyt
      isign = 1
      call wmpfft2rn(bxyze,bxyt,noff,nyp,isign,mixup,sct,tfft,tfmov,indx&
     &,indy,kstrt,nvp,kyp,ny,nterf,ierr)
!
! add constant to magnetic field with OpenMP: updates bxyze
      if (omt > 0.0) call mpbaddext2(bxyze,nyp,tfield,omx,omy,omz,nx)
!
! copy guard cells with OpenMP: updates fxyze, bxyze
      call wmpncguard2(fxyze,nyp,tguard,nx,kstrt,nvp)
      call wmpncguard2(bxyze,nyp,tguard,nx,kstrt,nvp)
! Update the Poynting vector
	  call UPDATE_POYNTING(ci,affp,dt,L_PML,delta,x,yp,fxyze,bxyze,&
	                       BCx,BCy,nx,ny,nyp,nxe,nypmx,pixyze,wpesc,wpescx,wpescy,tdiag)
	  call wmpncguard2(pixyze,nyp,tguard,nx,kstrt,nvp)
! Diagnostic of electromagnetic fields in real space ( M Touati )
      if (store_cond) then
!        call DIAG_REAL_FIELD(N_Ex, N_Ey, N_Ez, N_Bx, N_By, N_Bz,&
!                             nxe, nypmx, nyp, nx,&
!                             de, phtime, yp, x, fxyze, bxyze, tdiag)
         do ix=1,nyp
             do iy=1,nx
                 sfield(ix,iy) = fxyze(1,ix,iy)
              end do
         end do
         call file%new(iter=ntime, basePath='MS', axisLabels=(/'x','y'/), &
     &   gridSpacing = real(delta,4), records='E1', filenamebase = 'EandB',  &
     &   filepath='EMF/' )
!    &   gridGlobalOffset=(/0.0d0, 0.0d0/),&
!    &   position=(/0.0,0.0/))
!        call file%new(iter=ntime, axisLabels = (/'x','y'/), &
!    & gridSpacing=delta, gridGlobalOffset=(/ 0.0d0, 0.0d0 /), &
!    & basePath='MS',  records='E1')
         call pwfield(pp,file,sfield,(/nx,ny/),(/nx,nyp/),(/0,noff/),ierr)
         do ix=1,nyp
             do iy=1,nx
                 sfield(ix,iy) = fxyze(2,ix,iy)
              end do
         end do
         call file%new(iter=ntime, basePath='MS', axisLabels=(/'x','y'/), &
     &   gridSpacing = real(delta,4), records='E2',filenamebase ='EandB',filepath='EMF/')
!        call file%new(iter=ntime,axisLabels = (/'x','y'/), &
!    & gridSpacing=delta, gridGlobalOffset=(/ 0.0d0, 0.0d0 /),&
!    & basePath='MS', meshesPath='FLD', records='E2')
         call pwfield(pp,file,sfield,(/nx,ny/),(/nx,nyp/),(/0,noff/),ierr)
         do ix=1,nyp
             do iy=1,nx
                 sfield(ix,iy) = fxyze(3,ix,iy)
              end do
         end do
         call file%new(iter=ntime, basePath='MS', axisLabels=(/'x','y'/), &
     &   gridSpacing = real(delta,4), records='E3',filenamebase ='EandB',filepath='EMF/')
!        call file%new(iter=ntime,axisLabels = (/'x','y'/), &
!    & gridSpacing=delta, gridGlobalOffset=(/ 0.0d0, 0.0d0 /),&
!    & basePath='MS', meshesPath='FLD', records='E3')
         call pwfield(pp,file,sfield,(/nx,ny/),(/nx,nyp/),(/0,noff/),ierr)
        call DIAG_POYNTING(N_Pix,N_Piy,N_Piz,&
                           nxe, nypmx, nyp, nx,&
                           de, phtime, yp, x, pixyze, tdiag)
      end if
!
	 if (plasma) then
! push particles with OpenMP:
! updates ppart and wke, and possibly ncl, ihole, irc
    	wke = 0.0
     	call wmpbpush2(ppart,fxyze,bxyze,kpic,ncl,ihole,noff,nyp,qbme,dt, &
    				  &dth,ci,wke,tpush,nx,ny,mx,my,mx1,ipbc,relativity,list,x,y,delta,  &
    				  &irc)   
    	wke = me_real * wke
! updates pparti and wki, and possibly ncli, iholei, irc ( M. Touati )
	  	wki = 0.
    	if (movion) then
   			call wmpbpush2(pparti,fxyze,bxyze,kpici,ncli,iholei,noff,nyp,qbmi,dt, &
   						  &dth,ci,wki,tpush,nx,ny,mx,my,mx1,ipbc,relativity,list,x,y,delta,  &
   						  &irc)
    		wki = mi_real * wki
    	end if
! reorder particles by tile with OpenMP and MPI
! updates: ppart, kpic, and irc and possibly ncl and ihole
     	call ompmove2(ppart,kpic,ncl,ihole,noff,nyp,tsort,tmov,kstrt,nvp, &
    				 &nx,ny,mx,my,npbmx,nbmaxp,mx1,list,x,y,irc)
! updates: pparti, kpici, and irc and possibly ncli and iholei ( M. Touati )
  		if (movion) then
    		call ompmove2(pparti,kpici,ncli,iholei,noff,nyp,tsort,tmov,kstrt,nvp, &
   						 &nx,ny,mx,my,npbmx,nbmaxp,mx1,list,x,y,irc)     
     	end if
!
! Diagnostic of macro-particle phase space ( M Touati )
      	if (store_cond) then
        	call DIAG_PART(N_ele_px_vs_xy, N_ele_py_vs_xy, N_ele_pz_vs_xy, &
                           npp, idimp, nppmx0, mxyp1, Ndiag, de, phtime, ppart, tdiag)
        	if (movion) then
        		call DIAG_PART(N_ion_px_vs_xy, N_ion_py_vs_xy, N_ion_pz_vs_xy, &
                               nppi, idimp, nppmx0i, mxyp1, Ndiag, de, phtime, pparti, tdiag)
        	end if
      	end if
     end if
!
! energy diagnostic                   
      wt = we + wf + wm
      wtot(1) = wt
      wtot(2) = wke
      wtot(3) = wki
      wtot(4) = wfpdump
      wtot(5) = wmpdump
      wtot(6) = we
      wtot(7) = wf
      wtot(8) = wm
      wtot(9) = wpsrc
      wtot(10) = wfppml
      wtot(11) = wmppml
      wtot(12) = wpesc
      wtot(13) = wpescx
      wtot(14) = wpescy
      wtot(15) = wl
      call PPDSUM(wtot,work,15)
      wt   = wtot(1)
      wke  = wtot(2)
      wki  = wtot(3)
      wfdump = wtot(4)
      wmdump = wtot(5)
      wtdump = wfdump + wmdump
      we   = wtot(6)
      wf   = wtot(7)
      wm   = wtot(8)
      wsrc = wtot(9)
      wfpml = wtot(10)
      wfsim = wf - wfpml
      wmpml = wtot(11) 
      wmsim = wm - wmpml
      wtpml = wfpml + wmpml
      wtsim = wt - wtpml
      wesc  = wtot(12)
      wescx = wtot(13)
      wescy = wtot(14)
      wl    = wtot(15)
      if ( (kstrt==1) .and. (store_cond)) then
      	write (*,*)
        write (*,*) 'e- total particle instantaneous kinetic energy :' 
        write (*,*) '------------------------------------------------'
        write (*,'(3e14.7)') wke
        write (*,*)
        write (*,*) 'ions total particle instantaneous kinetic energy :' 
        write (*,*) '--------------------------------------------------'
        write (*,'(3e14.7)') wki
        write (*,*)
        write (*,*) 'instantaneous laser energy :' 
        write (*,*) '--------------------------------------------------'
        write (*,'(3e14.7)') wl
        write (*,*)
        if ((BCx == 1.) .and. (BCy == 1.)) then
			write (*,*) 'Instantaneous Electrostatic, Transverse Electric, Magnetic Field and Total Field Energies :'
			write (*,*) '-------------------------------------------------------------------------------------------'
			write (*,'(4e14.7)') we, wf, wm, wt
			write (*,*)
			write (*,*) 'Total Instantaneous Energy in the simulation box:'
			write (*,*) '-------------------------------------------------'
			write (*,'(3e14.7)') wki + wke + wt
			write (*,*)
        end if
        if ((BCx == 0.) .or. (BCy == 0.)) then
			write (*,*) 'Time integrated source of Electric energy in the effective simulation box :'
			write (*,*) '---------------------------------------------------------------------------'
			write (*,'(3e14.7)') wsrc
			write (*,*)
			write (*,*) 'Instantaneous Magnetic, Transverse Electric and Total energy in the effective simulation box :'
			write (*,*) '----------------------------------------------------------------------------------------------'
			write (*,'(3e14.7)') wmsim, wfsim, wtsim
			write (*,*)
			write (*,*) 'Instantaneous escaping energy from the effective simulation box to the pml layers:'
			write (*,*) '---------------------------------------------------------------------------------'
			write (*,'(3e14.7)') wesc
			write (*,*)
			write (*,*) 'Instantaneous Magnetic, Transverse Electric and Total energy in the PML layers :'
			write (*,*) '--------------------------------------------------------------------------------'
			write (*,'(3e14.7)') wmpml, wfpml, wtpml
			write (*,*)
			write (*,*) 'Time integrated Magnetic, Transverse Electric and Total Energy dumped in the PML layers :' 
			write (*,*) '-----------------------------------------------------------------------------------------'
			write (*,'(3e14.7)') wfdump, wmdump, wtdump
			write (*,*)
			write (*,*) 'Total Instantaneous Energy in the effective simulation box:'
        	write (*,*) '-----------------------------------------------------------'
        	write (*,'(3e14.7)') wki + wke + wtsim
        	write (*,*)
        end if 
      endif
      if (kstrt.eq.1) then
        call DIAG_ENERGY(BCx,BCy,movion,phtime,tdiag, &
                         N_es_fields, N_el_fields, N_ma_fields, &
                         N_kinetic_energy,N_ele_kinetic_energy, N_ion_kinetic_energy, &
                         N_el_fields_pml, N_ma_fields_pml, N_fields_src, &
           				 N_fields_esc, N_el_dumped, N_ma_dumped, &
                         we, wf, wm, wki+wke, wke, wki, wsrc, wfsim, wmsim, &
                         wesc,wfpml,wmpml,wfdump,wmdump)
      end if
! Physical time increment
      phtime = phtime + dt
      end do
!
! * * * end main iteration loop * * *
!
      ntime = ntime + 1
!
! loop time
      call dtimer(dtime,ltime,1)
      tloop = tloop + real(dtime)
!
      if (kstrt.eq.1) then
         write (*,*)
         write (*,*) 'ntime, relativity = ', ntime, relativity
         write (*,*) 'MPI nodes nvp = ', nvp
         wt = we + wf + wm
         write (*,*) 'Final Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
         write (*,*) 'Final Electrostatic, Transverse Electric and Magne&
     &tic Field Energies:'
         write (*,'(3e14.7)') we, wf, wm
!
         write (*,*)
         write (*,*) 'initialization time = ', tinit
         write (*,*) 'deposit time = ', tdpost
         write (*,*) 'current deposit time = ', tdjpost
         tdpost = tdpost + tdjpost
         write (*,*) 'total deposit time = ', tdpost
         write (*,*) 'guard time = ', tguard
         write (*,*) 'solver time = ', tfield
         write (*,*) 'field move time = ', tfmov
         write (*,*) 'fft and transpose time = ', tfft(1), tfft(2)
         write (*,*) 'push time = ', tpush
         write (*,*) 'particle move time = ', tmov
         write (*,*) 'sort time = ', tsort
         tfield = tfield + tguard + tfft(1) + tfmov
         write (*,*) 'total solver time = ', tfield
         tsort = tsort + tmov
         time = tdpost + tpush + tsort
         write (*,*) 'total particle time = ', time
         write (*,*) 'diagnostic time = ', tdiag
         wt = time + tfield
         tloop = tloop - wt
         write (*,*) 'total and additional time = ', wt, tloop
         write (*,*)
!
         wt = 1.0e+09/(real(nloop)*real(np))
         write (*,*) 'Push Time (nsec) = ', tpush*wt
         write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
         write (*,*) 'Sort Time (nsec) = ', tsort*wt
         write (*,*) 'Total Particle Time (nsec) = ', time*wt
      endif
!
      call PPEXIT()
      end program
