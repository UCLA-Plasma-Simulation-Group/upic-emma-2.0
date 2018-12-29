module modmpfieldpml2
! 
	use libmpfieldpml2_h
    use ompplib2
      
	implicit none
	public  :: setup_pml
	public  :: mpmaxwel2pml_init, mpmaxwel2pml
	private :: mpmaxwel2pml_yee, mpmaxwel2pml_implicit
	private :: mpmaxwelb2pml, mpmaxwele2pml
	contains
!
	subroutine setup_pml(sigma,vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
	  	                 ax,ay,dt,x,yp,delta,L_PML,n_cond_PML,ci,&
	                     nx,ny,nxe,nye,BCx,BCy,nyp,nypmx,kxp,nvp,kstrt)
! This subroutine allocates the Perfectly Matched Layer (PML) elctromagnetic  
! field arrays vpmlx_exyz / vpml_bxyz and conductivity sigma
		implicit none
! input :  _ BCx (BCy) = 0 if PML along x (y)
!          _ L_PML = PML layer thickness
!          _ n_cond_PML = real n such that the PMLelectrical and magnetic conductivity reads :
!            sigma(x) = sigma_m * (|x-x_I|/L_PML)^n where x_I is the position of the 
!            interface between the PML and the effective sim. box
!          _ nvp   = number of MPI node used in the simulation
!          _ kstrt = index of the MPI node from 1 to nvp
!          _ nx = physical x size of the simulation box partition
!          _ ny = physical y size of the simulation box
!          _ nxe = nx + 2 = x size of the simulation box counting guard cells
!          _ nyp = physical y size of the simulation box partition
!          _ nypmx = nyp + 2 = y size of the simulation box partition counting guard cells
!          _ ci = inverse speed of light
!          _ delta(1), delta(2) = cell size along x , cell size along y
!          _ x(i) and yp(j) = y(joff+j) are the values of space grid points positions
		integer, intent(in)                                 :: nx, ny, nxe, nye
		integer, intent(in)                                 :: Bcx, Bcy, FDTD
		integer, intent(in)                                 :: kxp, nyp, nypmx
		integer, intent(in)                                 :: nvp, kstrt
		real, intent(in)                                    :: L_PML
		real, intent(in)                                    :: n_cond_PML
		real, intent(in)                                    :: ci
		real, dimension(2), intent(in)                      :: delta
		real, dimension(nxe), intent(in)                    :: x
        real, dimension(nypmx), intent(in)                  :: yp
! output : _ sigma PML conductivity array : sigma(1,i,j) = sigma_x(x(i))
!                                           sigma(2,i,j) = sigma_y(y(j))
!          _ allocate the PML E / B field arrays vpml_exyz and vpml_bxyz in the Fourier space
!            exemple : vpml_exyz(1,:,:) = Exy(:,:),
!                      vpml_exyz(2,:,:) = Eyx(:,:),
!                      vpml_exyz(2,:,:) = Ezx(:,:) and
!                      vpml_exyz(2,:,:) = Ezy(:,:).
		real, intent(inout)                                 :: ax,ay,dt
		real, dimension(:,:,:), allocatable, intent(out)    :: sigma
		complex, dimension(:,:,:), allocatable, intent(out) :: vpml_exyz, vpml_bxyz
		complex, dimension(:,:), allocatable, intent(out)   :: vpml_exyz_corr, vpml_bxyz_corr
! local variables
		real               :: sigma_m
!         
! Find the optimal maximum value of the PML magnetic and electric conductivity
! according to Computational electrodynamics, the FDTD method, A. Taflove and 
! S C Hagness, Eq. (7.62) and (7.66)
		sigma_m = 0.8 * real(n_cond_PML + 1) / ci
! More general formula for getting R=10^-15 in theory at normal incidence
		!sigma_m = 15.*log(10.)*real(n_cond_PML + 1) / (2.*ci*L_PML)
! Modify the time step
		dt = min(dt,dt/sigma_m)
! Setup PML conductivity
        allocate(sigma(2,nxe,nypmx))
        call setup_sigma(sigma,ax,ay,dt,x,yp,delta,sigma_m,L_PML,n_cond_PML,&
	                     nx,ny,nxe,nyp,nypmx,BCx,BCy,FDTD,nvp,kstrt)
! Allocate PML EM fields
        allocate(vpml_exyz(4,nye,kxp),vpml_bxyz(4,nye,kxp))
        allocate(vpml_exyz_corr(4,4),vpml_bxyz_corr(4,4))
!
    end subroutine setup_pml
	
	elemental function filter(n_cond_PML,L_PML,side,x0,x1,x2,x3,x)
		implicit none
		integer, intent(in) :: side
		real, intent(in)    :: n_cond_PML, L_PML, x0, x1, x2, x3, x
		real                :: filter
! locals 
		real :: pi, dx, xm
!
		pi = 3.14159265359
		dx = L_PML/2.
		if (side == 0) then
! 			filter = sin(0.5*pi*(x-x0)/x1)
!			filter = 1. - exp(-(4.**n_cond_PML)*log(2.)*(     x/dx)**(2.*n_cond_PML))
			xm = 0.5*L_PML
			filter = 1. - exp(-log(2.)*(2.*(x-xm)/dx)**(2.*n_cond_PML))
		else if (side == 1) then
! 			filter = cos(0.5*pi*(x-x2)/x1)
! 			filter = 1. - exp(-(4.**n_cond_PML)*log(2.)*((x-x3)/dx)**(2.*n_cond_PML))
			xm = x3 - (0.5*L_PML)
			filter = 1. - exp(-log(2.)*(2.*(x-xm)/dx)**(2.*n_cond_PML))
		end if
	end function filter
	
	subroutine setup_sigma(sigma,ax,ay,dt,x,yp,delta,sigma_m,L_PML,n_cond_PML,&
	                       nx,ny,nxe,nyp,nypmx,BCx,BCy,FDTD,nvp,kstrt)
! this subroutine setup the PNL conductivity sigma in the real space
! sigma(1,i,j) = sigma_x(x(i)) and sigma(2,i,j) = sigma_y(y(j))
		implicit none
! inputs
		integer, intent(in)                :: nx, ny, nxe, nyp, nypmx
		integer, intent(in)                :: BCx, BCy, FDTD, nvp, kstrt
		real, intent(in)                   :: ax, ay, dt, sigma_m, L_PML, n_cond_PML
		real, dimension(nxe), intent(in)   :: x
		real, dimension(nypmx), intent(in) :: yp
		real, dimension(2), intent(in)     :: delta
!outputs
        real, dimension(2,nxe,nypmx), intent(out) :: sigma
!locals
		integer                      :: i, j, dim, filtre, filter_new, polynomial, ii, jj
		real                         :: x0, y0, x1, x2, x3, y1, y2, y3
		real                         :: tguard, S0, norm, pi, bx, by, ft
		real, dimension(2,nxe,nypmx) :: sigma_s
!
		sigma  = 0.
		x0 = 0.
		y0 = 0.
        x1 = L_PML
        y1 = L_PML
        x3 = real(nx)*delta(1)
        y3 = real(ny)*delta(2)
        x2 = x3 - L_PML
        y2 = y3 - L_PML
        filter_new = 0
        polynomial = 1
        filtre     = polynomial
        select case (filtre)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case (0)
! Case PML along x only
		if ((BCx == 0) .and. (BCy ==1)) then
			!$OMP PARALLEL DO PRIVATE(i,j,ft)
			do j=1,nyp,1
				do i=1,nx,1
					if      ( (x(i) .le. x1) .and. (x(i) .ne. x0) ) then
						ft = filter(n_cond_PML,L_PML,0,x0,x1,x2,x3,x(i))
						sigma(1,i,j) = ( 1. - ft ) / ( dt * ft )
					else if ( (x(i) .ge. x2) .and. (x(i) .ne. x3) ) then
						ft = filter(n_cond_PML,L_PML,1,x0,x1,x2,x3,x(i))
						sigma(1,i,j) = ( 1. - ft ) / ( dt * ft )
					else
						sigma(1,i,j) = 0.
					end if
				end do
			end do
			!$END OMP PARALLEL DO
		end if
! Case PML along y only
		if ((BCx == 1) .and. (BCy ==0)) then
			!$OMP PARALLEL DO PRIVATE(i,j,ft)
			do j=1,nyp,1
				do i=1,nx,1
					if      ( (yp(j) .le. y1) .and. (yp(j) .ne. y0) ) then
						ft = filter(n_cond_PML,L_PML,0,y0,y1,y2,y3,yp(j))
						sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
					else if ( (yp(j) .ge. y2) .and. (yp(j) .ne. y3) ) then
						ft = filter(n_cond_PML,L_PML,1,y0,y1,y2,y3,yp(j))
						sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
					else
						sigma(2,i,j) = 0.
					end if
				end do
			end do
			!$END OMP PARALLEL DO
		end if
! Case PML along x and y
		if ((BCx == 0) .and. (BCy ==0)) then
			!$OMP PARALLEL DO PRIVATE(i,j,ft)
			do j=1,nyp,1
				do i=1,nx+1,1
! 					if      ( (x(i) .le. x1) .and. (x(i) .ne. x0) ) then
					if      ( (x(i) .le. x3/2.) ) then
						ft = filter(n_cond_PML,L_PML,0,x0,x1,x2,x3,x(i))
						sigma(1,i,j)     = ( 1. - ft ) / ( dt * ft )
! 						if      ( (yp(j) .le. y1) .and. (yp(j) .ne. y0) )then
						if      ( (yp(j) .le. y3/2.) )then
							ft = filter(n_cond_PML,L_PML,0,y0,y1,y2,y3,yp(j))
							sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
! 						else if ( (yp(j) .ge. y2) .and. (yp(j) .ne. y3) )then
						else if ( (yp(j) .gt. y3/2.) )then
							ft = filter(n_cond_PML,L_PML,1,y0,y1,y2,y3,yp(j))
							sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
						end if
! 					else if ( (x(i) .ge. x2) .and. (x(i) .ne. x3) ) then
					else if ( ( x(i) .gt. x3/2.) ) then
						ft = filter(n_cond_PML,L_PML,1,x0,x1,x2,x3,x(i))
						sigma(1,i,j)     = ( 1. - ft ) / ( dt * ft )
! 						if      ( (yp(j) .le. y1) .and. (yp(j) .ne. y0) )then
						if      ( (yp(j) .le. y3/2.) )then
							ft = filter(n_cond_PML,L_PML,0,y0,y1,y2,y3,yp(j))
							sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
! 						else if ( (yp(j) .ge. y2) .and. (yp(j) .ne. y3) )then
						else if ( (yp(j) .gt. y3/2.) )then
							ft = filter(n_cond_PML,L_PML,1,y0,y1,y2,y3,yp(j))
							sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
						end if
					else
! 						if      ( (yp(j) .le. y1) .and. (yp(j) .ne. y0) )then
						if      ( (yp(j) .le. y3/2.) )then
							ft = filter(n_cond_PML,L_PML,0,y0,y1,y2,y3,yp(j))
							sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
! 						else if ( (yp(j) .ge. y2) .and. (yp(j) .ne. y3) )then
						else if ( (yp(j) .gt. y3/2.) )then
							ft = filter(n_cond_PML,L_PML,1,y0,y1,y2,y3,yp(j))
							sigma(2,i,j) = ( 1. - ft ) / ( dt * ft )
						end if
					end if
				end do
			end do
			!$END OMP PARALLEL DO
		end if
		!$OMP PARALLEL DO PRIVATE(i)
		do i=1,nx,1
			if (kstrt == 1  ) sigma(2,i,1  ) = sigma(2,i,2    )
			if (kstrt == nvp) sigma(2,i,nyp) = sigma(2,i,nyp-1)
		end do
		!$OMP END PARALLEL DO
		!$OMP PARALLEL DO PRIVATE(j)
		do j=1,nyp,1
			sigma(1,1 ,j) = sigma(1,2   ,j)
			sigma(1,nx,j) = sigma(1,nx-1,j)
		end do
		!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case (1)
! Case PML along x only
		if ((BCx == 0) .and. (BCy ==1)) then
			!$OMP PARALLEL DO PRIVATE(i,j)
			do j=1,nyp,1
				do i=1,nx,1
					if (x(i) .le. x1) then
						sigma(1,i,j) = sigma_m * (((x1-x(i))/L_PML)**n_cond_PML)
					else if (x(i) .ge. x2) then
						sigma(1,i,j) = sigma_m * (((x(i)-x2)/L_PML)**n_cond_PML)
					else
						sigma(1,i,j) = 0.
					end if
				end do
			end do
			!$END OMP PARALLEL DO
		end if
! Case PML along y only
		if ((BCx == 1) .and. (BCy ==0)) then
			!$OMP PARALLEL DO PRIVATE(i,j)
			do j=1,nyp,1
				do i=1,nx,1
					if (yp(j) .le. y1) then
						sigma(2,i,j) = sigma_m * (((y1-yp(j))/L_PML)**n_cond_PML)
					else if (yp(j) .ge. y2) then
						sigma(2,i,j) = sigma_m * (((yp(j)-y2)/L_PML)**n_cond_PML)
					else
						sigma(2,i,j) = 0.
					end if
				end do
			end do
			!$END OMP PARALLEL DO
		end if
! Case PML along x and y
		if ((BCx == 0) .and. (BCy ==0)) then
			!$OMP PARALLEL DO PRIVATE(i,j)
			do j=1,nyp,1
				do i=1,nx+1,1
					if (x(i) .le. x1) then
						sigma(1,i,j) = sigma_m * (((x1-x(i))/L_PML)**n_cond_PML)
						if (yp(j) .le. y1) then
							sigma(2,i,j) = sigma_m * (((y1-yp(j))/L_PML)**n_cond_PML)
						else if (yp(j) .ge. y2) then
							sigma(2,i,j) = sigma_m * (((yp(j)-y2)/L_PML)**n_cond_PML)
						end if
					else if (x(i) .ge. x2) then
						sigma(1,i,j) = sigma_m * (((x(i)-x2)/L_PML)**n_cond_PML)
						if (yp(j) .le. y1) then
							sigma(2,i,j) = sigma_m * (((y1-yp(j))/L_PML)**n_cond_PML)
						else if (yp(j) .ge. y2) then
							sigma(2,i,j) = sigma_m * (((yp(j)-y2)/L_PML)**n_cond_PML)
						end if
					else
						if (yp(j) .le. y1) then
							sigma(2,i,j) = sigma_m * (((y1-yp(j))/L_PML)**n_cond_PML)
						else if (yp(j) .ge. y2) then
							sigma(2,i,j) = sigma_m * (((yp(j)-y2)/L_PML)**n_cond_PML)
						end if
					end if
				end do
			end do
			!$END OMP PARALLEL DO
		end if
        end select
        ! Filtering
!         pi = 3.14159265359
!         bx = .912871
!         by = .912871
!         norm =        ((2.*pi)**(-0.5)) / bx
!         norm = norm * ((2.*pi)**(-0.5)) / by
!         do dim=1,2,1
! 			do j=1,nyp,1
! 				do i=1,nx,1
! 						sigma_s(dim,i,j) = 0.
! 						do jj=1,nyp,1
! 							do ii=1,nx,1
! 								S0 = norm * exp(-0.5 * ( (( x(ii)- x(i))/bx)**2. ) )
! 								S0 = S0   * exp(-0.5 * ( ((yp(jj)-yp(j))/by)**2. ) )
! 		! 						print*,'S0=',S0
! 								sigma_s(dim,i,j) = sigma_s(dim,i,j) +  ( sigma(dim,ii,j) * S0 * delta(1) * delta(2) )
! 		! 						print*,'sigma_s',sigma_s
! 							end do
! 						end do
! 						sigma(dim,i,j) = sigma_s(dim,i,j)
! 				end do
! 			end do
!         end do
        do j=1,nyp,1
			do dim = 1,2,1
				sigma(dim,nx+1,j) = sigma(dim,1,j)
				sigma(dim,nx+2,j) = sigma(dim,1,j)
			end do
		end do
        call wmpncguard2(sigma,nyp,tguard,nx,kstrt,nvp)
	end subroutine setup_sigma
    
    subroutine mpmaxwel2pml_init(vpml_exyz,vpml_bxyz,vpml_exyz_corr,vpml_bxyz_corr,&
                                 exyz,bxyz,nye,kxp,nyh)
! This subroutine copy the initial Darwin electromagnetic fields at t = 0 in the 
! PML field arrays
		implicit none
! input
		integer, intent(in)                        :: nye, kxp, nyh
		complex, dimension(3,nye,kxp), intent(in)  :: exyz, bxyz
! output 
		complex, dimension(4,nye,kxp), intent(out) :: vpml_exyz, vpml_bxyz
		complex, dimension(4,4), intent(out)       :: vpml_exyz_corr, vpml_bxyz_corr
! locals
		integer :: q, p, q1, dim
!
		!$OMP PARALLEL DO PRIVATE(q,p)
		do p =1,kxp,1
			do q=1,nye,1
				! Exy / Bxy component
				vpml_exyz(1,q,p) =     exyz(1,q,p)
				vpml_bxyz(1,q,p) =     bxyz(1,q,p)
				! Eyx / Byx component
				vpml_exyz(2,q,p) =     exyz(2,q,p)
				vpml_bxyz(2,q,p) =     bxyz(2,q,p)
				! Ezx / Bzx component
				vpml_exyz(3,q,p) = 0.5*exyz(3,q,p)
				vpml_bxyz(3,q,p) = 0.5*bxyz(3,q,p)
				! Ezy / Bzy component
				vpml_exyz(4,q,p) = 0.5*exyz(3,q,p)
				vpml_bxyz(4,q,p) = 0.5*bxyz(3,q,p)
			end do
		end do
		!$OMP END PARALLEL DO
		q1 = nyh + 1
		do dim = 1,4,1
			! mode (kx,ky) = (0,0) :
			vpml_exyz_corr(1,dim) = cmplx(real(vpml_exyz(dim,1,1)),0.0)
			vpml_bxyz_corr(1,dim) = cmplx(real(vpml_bxyz(dim,1,1)),0.0)
			! mode (kx,ky) = (pi/delta_x,0) :
			vpml_exyz_corr(2,dim) = cmplx(aimag(vpml_exyz(dim,1,1)),0.0)
			vpml_bxyz_corr(2,dim) = cmplx(aimag(vpml_bxyz(dim,1,1)),0.0)
			! mode (kx,ky) = (0,pi/delta_y) :
			vpml_exyz_corr(3,dim) = cmplx(real(vpml_exyz(dim,q1,1)),0.0)
			vpml_bxyz_corr(3,dim) = cmplx(real(vpml_bxyz(dim,q1,1)),0.0)
			! mode (kx,ky) = (pi/delta_x,pi/delta_y) :
			vpml_exyz_corr(4,dim) = cmplx(aimag(vpml_exyz(dim,q1,1)),0.0)
			vpml_bxyz_corr(4,dim) = cmplx(aimag(vpml_bxyz(dim,q1,1)),0.0)
		end do
    end subroutine mpmaxwel2pml_init
    
    subroutine mpmaxwel2pml(pml_scheme,ndim,indx,indy,nx,nxe,ny,nye,nyp,nypmx,nxhy,nxyh,nyh,kxp,kyp,&
	                        nvp,kstrt,noff,BCx,BCy,L_PML,ax,ay,affp,ci,dt,kx,ky,x,yp,delta,sigma,&
	                        ffc,cut,nterf,ierr,mixup,wm,wf,tfield,tfmov,tguard,&
	                        tfft,vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
	  	                    sct,exyz,bxyz,wescx,wescy,wfppml,wmppml,wfpdump,wmpdump)
	    implicit none
! input :
	    integer, intent(in) 				     		:: pml_scheme, ndim, indx, indy
	    integer, intent(in)			    	     		:: nx, nxe
	    integer, intent(in) 				     		:: ny, nye, nyp, nypmx 
	    integer, intent(in) 				     		:: nxhy, nxyh, nyh
	    integer, intent(in) 				     		:: kxp, kyp
	    integer, intent(in) 				     		:: nvp, kstrt, noff
	    integer, intent(in)                      		:: BCx, BCy, FDTD
	    real, intent(in)                         		:: L_PML, ax, ay, affp, ci, dt   
	    real, dimension(nxe), intent(in)           		:: kx
	    real, dimension(nye), intent(in)           		:: ky
	    real, dimension(nxe), intent(in)                :: x
        real, dimension(nypmx), intent(in)              :: yp
	    real, dimension(2), intent(in)           		:: delta
	    real, dimension(2,nxe,nypmx), intent(in) 		:: sigma
        complex, dimension(nyh,kxp), intent(in)         :: ffc
        complex, dimension(ndim,nye,kxp), intent(in)    :: cut
! in/output
		integer, intent(inout)                          :: nterf, ierr 
		integer, dimension(nxhy), intent(inout)         :: mixup
		complex, dimension(nxyh), intent(inout)         :: sct
		real, intent(inout)                             :: wm, wf 
        real, intent(inout)                             :: tfield, tfmov, tguard
	    real, dimension(2), intent(inout)               :: tfft
	    complex, dimension(4,nye,kxp), intent(inout)    :: vpml_exyz, vpml_bxyz 
	    complex, dimension(4,4), intent(inout)          :: vpml_exyz_corr, vpml_bxyz_corr 
	    complex, dimension(ndim,nye,kxp), intent(inout) :: exyz, bxyz
	    real, intent(in)                                :: wescx,wescy
	    real, intent(inout)                             :: wfppml, wmppml
	    real, intent(inout)                             :: wfpdump, wmpdump
! locals
		integer, parameter :: fully_implicit = 1
!
    	if (pml_scheme == fully_implicit) then
    		call mpmaxwel2pml_implicit(ndim,indx,indy,nx,nxe,ny,nye,nyp,nypmx,nxhy,nxyh,nyh,kxp,kyp,&
	                                   nvp,kstrt,noff,BCx,BCy,L_PML,ax,ay,affp,ci,dt,kx,ky,x,yp,delta,sigma,&
	                                   ffc,cut,nterf,ierr,mixup,wm,wf,tfield,tfmov,tguard,&
	                                   tfft,vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
	                                   sct,exyz,bxyz,wescx,wescy,wfppml,wmppml,wfpdump,wmpdump)
    	else
    		call mpmaxwel2pml_yee(ndim,indx,indy,nx,nxe,ny,nye,nyp,nypmx,nxhy,nxyh,nyh,kxp,kyp,&
	                              nvp,kstrt,noff,BCx,BCy,L_PML,ax,ay,affp,ci,dt,kx,ky,x,yp,delta,sigma,&
	                              ffc,cut,nterf,ierr,mixup,wm,wf,tfield,tfmov,tguard,&
	                              tfft,vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
	                              sct,exyz,bxyz,wescx,wescy,wfppml,wmppml,wfpdump,wmpdump)
    	end if
    end subroutine mpmaxwel2pml
    
    subroutine mpmaxwel2pml_implicit(ndim,indx,indy,nx,nxe,ny,nye,nyp,nypmx,nxhy,nxyh,nyh,kxp,kyp,&
	                                 nvp,kstrt,noff,BCx,BCy,L_PML,ax,ay,affp,ci,dt,kx,ky,x,yp,delta,sigma,&
	                                 ffc,cut,nterf,ierr,mixup,wm,wf,tfield,tfmov,tguard,&
	                                 tfft,vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
	                                 sct,exyz,bxyz,wescx,wescy,wfppml,wmppml,wfpdump,wmpdump)
	    implicit none
! input :
	    integer, intent(in) 				     		:: ndim, indx, indy
	    integer, intent(in)			    	     		:: nx, nxe
	    integer, intent(in) 				     		:: ny, nye, nyp, nypmx 
	    integer, intent(in) 				     		:: nxhy, nxyh, nyh
	    integer, intent(in) 				     		:: kxp, kyp
	    integer, intent(in) 				     		:: nvp, kstrt, noff
	    integer, intent(in)                      		:: BCx, BCy, FDTD
	    real, intent(in)                         		:: L_PML, ax, ay, affp, ci, dt   
	    real, dimension(nxe), intent(in)           		:: kx
	    real, dimension(nye), intent(in)           		:: ky
	    real, dimension(nxe), intent(in)                :: x
        real, dimension(nypmx), intent(in)              :: yp
	    real, dimension(2), intent(in)           		:: delta
	    real, dimension(2,nxe,nypmx), intent(in) 		:: sigma
        complex, dimension(nyh,kxp), intent(in)         :: ffc
        complex, dimension(ndim,nye,kxp), intent(in)    :: cut
! in/output
		integer, intent(inout)                          :: nterf, ierr 
		integer, dimension(nxhy), intent(inout)         :: mixup
		complex, dimension(nxyh), intent(inout)         :: sct
		real, intent(inout)                             :: wm, wf 
        real, intent(inout)                             :: tfield, tfmov, tguard
	    real, dimension(2), intent(inout)               :: tfft
	    complex, dimension(4,nye,kxp), intent(inout)    :: vpml_exyz, vpml_bxyz
	    complex, dimension(4,4), intent(inout)          :: vpml_exyz_corr, vpml_bxyz_corr 
	    complex, dimension(ndim,nye,kxp), intent(inout) :: exyz, bxyz
	    real, intent(in)                         		:: wescx,wescy
	    real, intent(inout)                             :: wfppml, wmppml
	    real, intent(inout)                             :: wfpdump, wmpdump 
! locals
		integer                                         :: isign
        real, dimension(4,nxe,nypmx)                    :: vpml_exyze, vpml_bxyze
        complex, dimension(4,nye,kxp)                   :: vpml_bxyz_np1
        real, dimension(4,nxe,nypmx)                    :: vpml_bxyze_np1
!
		! Update the PML B field from n to n+1/2'
		call mpmaxwelb2pml(vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
		                   cut,ffc,affp,ci,dt,wm,tfield,nx,ny,&
						   kstrt,kx,ky,ax,ay)
		! Transform the PML B field at n+1/2 to real space'
		isign = 1
		call wmpfft2rn(vpml_bxyze,vpml_bxyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_bxyze,nyp,tguard,nx,kstrt,nvp)
		! Dump the B-field at n+1/2
		if ((wescx > 0.) .or. (wescy > 0.)) then
			call mpdumppml_implicit(vpml_bxyze,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		end if
		! Transform the PML B-field at n+1/2 fields back to Fourier space'
		isign = -1
		call wmpfft2rn(vpml_bxyze,vpml_bxyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		! Update the PML E field from n to n+1'
		call mpmaxwele2pml(vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
		                   cut,ffc,affp,ci,dt,wf,tfield,nx,ny,&
						   kstrt,kx,ky,ax,ay)
		! Transform the PML E field at n+1 to real space'
		isign = 1
		call wmpfft2rn(vpml_exyze,vpml_exyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_exyze,nyp,tguard,nx,kstrt,nvp)
		! Dump the E field at n+1
		if ((wescx > 0.) .or. (wescy > 0.)) then
			call mpdumppml_implicit(vpml_exyze,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		end if
		! Update the dumped electric energy wfpdump between n and n+1 and
		!        the instantaneous electric energy wfppml at n+1 in the PML layers 
		call update_energy(sigma,x,yp,vpml_exyze,BCx,BCy,L_PML,dt,delta,affp,ci,&
		                   wescx,wescy,wfppml,wfpdump,nx,ny,nyp,nxe,nypmx)
		! Transform the PML E field at n+1 back to Fourier space'
		isign = -1
		call wmpfft2rn(vpml_exyze,vpml_exyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		! Update the PML B field from n+1/2 to n+1'
		call mpmaxwelb2pml(vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
		                   cut,ffc,affp,ci,dt,wm,tfield,nx,ny,&
						   kstrt,kx,ky,ax,ay)
		! Copy the B field at n+1
		call mpcopypml(vpml_bxyz,vpml_bxyz_np1,4,nye,kxp)
		! Transform the PML B field at n+1 to real space'
		isign = 1
		call wmpfft2rn(vpml_bxyze_np1,vpml_bxyz_np1,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_bxyze_np1,nyp,tguard,nx,kstrt,nvp)
		! Dump the B field at n+1
		if ((wescx > 0.) .or. (wescy > 0.)) then
			call mpdumppml_implicit(vpml_bxyze_np1,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		end if
		! Update the dumped magnetic energy wmpdump between n and n+1 and
		!        the instantaneous magnetic energy wmppml at n+1 in the PML layers 
		call update_energy(sigma,x,yp,vpml_bxyze_np1,BCx,BCy,L_PML,dt,delta,affp,ci,&
		                   wescx,wescy,wmppml,wmpdump,nx,ny,nyp,nxe,nypmx)
		! Transform the PML B field at n+1 back to Fourier space'
		isign = -1
		call wmpfft2rn(vpml_bxyze_np1,vpml_bxyz_np1,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		! Update the box E and B fields from the PML E and B fields
		call mpupdatefields(exyz,bxyz,vpml_exyz,vpml_bxyz_np1,ny,nye,kxp,kstrt)
    end subroutine mpmaxwel2pml_implicit
    
    subroutine mpmaxwel2pml_yee(ndim,indx,indy,nx,nxe,ny,nye,nyp,nypmx,nxhy,nxyh,nyh,kxp,kyp,&
	                            nvp,kstrt,noff,BCx,BCy,L_PML,ax,ay,affp,ci,dt,kx,ky,x,yp,delta,sigma,&
	                            ffc,cut,nterf,ierr,mixup,wm,wf,tfield,tfmov,tguard,&
	                            tfft,vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
	                            sct,exyz,bxyz,wescx,wescy,wfppml,wmppml,wfpdump,wmpdump)
	    implicit none
! input :
	    integer, intent(in) 				     		:: ndim, indx, indy
	    integer, intent(in)			    	     		:: nx, nxe
	    integer, intent(in) 				     		:: ny, nye, nyp, nypmx 
	    integer, intent(in) 				     		:: nxhy, nxyh, nyh
	    integer, intent(in) 				     		:: kxp, kyp
	    integer, intent(in) 				     		:: nvp, kstrt, noff
	    integer, intent(in)                      		:: BCx, BCy, FDTD
	    real, intent(in)                         		:: L_PML, ax, ay, affp, ci, dt   
	    real, dimension(nxe), intent(in)           		:: kx
	    real, dimension(nye), intent(in)           		:: ky
	    real, dimension(nxe), intent(in)                :: x
        real, dimension(nypmx), intent(in)              :: yp
	    real, dimension(2), intent(in)           		:: delta
	    real, dimension(2,nxe,nypmx), intent(in) 		:: sigma
        complex, dimension(nyh,kxp), intent(in)         :: ffc
        complex, dimension(ndim,nye,kxp), intent(in)    :: cut
! in/output
		integer, intent(inout)                          :: nterf, ierr 
		integer, dimension(nxhy), intent(inout)         :: mixup
		complex, dimension(nxyh), intent(inout)         :: sct
		real, intent(inout)                             :: wm, wf 
        real, intent(inout)                             :: tfield, tfmov, tguard
	    real, dimension(2), intent(inout)               :: tfft
	    complex, dimension(4,nye,kxp), intent(inout)    :: vpml_exyz, vpml_bxyz
	    complex, dimension(4,4), intent(inout)          :: vpml_exyz_corr, vpml_bxyz_corr
	    complex, dimension(ndim,nye,kxp), intent(inout) :: exyz, bxyz
	    real, intent(in)                         		:: wescx,wescy
	    real, intent(inout)                             :: wfppml, wmppml
	    real, intent(inout)                             :: wfpdump, wmpdump 
! locals
		integer                                         :: isign
        real, dimension(4,nxe,nypmx)                    :: vpml_exyze, vpml_bxyze
        complex, dimension(4,nye,kxp)                   :: vpml_bxyz_np1
        real, dimension(4,nxe,nypmx)                    :: vpml_bxyze_np1
!
		! Update the PML B field from n to n+1/2'
		call mpmaxwelb2pml(vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
		                   cut,ffc,affp,ci,dt,wm,tfield,nx,ny,&
						   kstrt,kx,ky,ax,ay)
		! Transform the PML B field at n+1/2 and the E field at n to real space'
		isign = 1
		call wmpfft2rn(vpml_bxyze,vpml_bxyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_bxyze,nyp,tguard,nx,kstrt,nvp)
		call wmpfft2rn(vpml_exyze,vpml_exyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_exyze,nyp,tguard,nx,kstrt,nvp)
		! Dump the B-field at n+1/2
		if ((wescx > 0.) .or. (wescy > 0.)) then
			call mpdumppml_yee(vpml_bxyze,vpml_exyze,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		end if
		! Transform the PML fields to Fourier space'
		isign = -1
		call wmpfft2rn(vpml_bxyze,vpml_bxyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpfft2rn(vpml_exyze,vpml_exyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		! Update the PML E field from n to n+1'
		call mpmaxwele2pml(vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
		                   cut,ffc,affp,ci,dt,wf,tfield,nx,ny,&
						   kstrt,kx,ky,ax,ay)
		! Transform the PML E field at n+1 and the B field at n+1/2 to real space'
		isign = 1
		call wmpfft2rn(vpml_bxyze,vpml_bxyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_bxyze,nyp,tguard,nx,kstrt,nvp)
		call wmpfft2rn(vpml_exyze,vpml_exyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_exyze,nyp,tguard,nx,kstrt,nvp)
		! Dump the E field at n+1
		if ((wescx > 0.) .or. (wescy > 0.)) then
			call mpdumppml_yee(vpml_exyze,vpml_bxyze,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		end if
		! Update the dumped electric energy wfpdump between n and n+1 and
		!        the instantaneous electric energy wfppml at n+1 in the PML layers 
		call update_energy(sigma,x,yp,vpml_exyze,BCx,BCy,L_PML,dt,delta,affp,ci,&
		                   wescx,wescy,wfppml,wfpdump,nx,ny,nyp,nxe,nypmx)
		! Transform the PML fields to Fourier space'
		isign = -1
		call wmpfft2rn(vpml_bxyze,vpml_bxyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpfft2rn(vpml_exyze,vpml_exyz,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		! Update the PML B field from n+1/2 to n+1'
		call mpmaxwelb2pml(vpml_exyz,vpml_bxyz,FDTD,vpml_exyz_corr,vpml_bxyz_corr,&
		                   cut,ffc,affp,ci,dt,wm,tfield,nx,ny,&
						   kstrt,kx,ky,ax,ay)
		! Copy the B field at n+1
		call mpcopypml(vpml_bxyz,vpml_bxyz_np1,4,nye,kxp)
		! Transform the PML B field at n+1 to real space'
		isign = 1
		call wmpfft2rn(vpml_bxyze_np1,vpml_bxyz_np1,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		call wmpncguard2(vpml_bxyze_np1,nyp,tguard,nx,kstrt,nvp)
		! Dump the B field at n+1
		if ((wescx > 0.) .or. (wescy > 0.)) then
			call mpdumppml_yee2(vpml_bxyze_np1,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		end if
		! Update the dumped magnetic energy wmpdump between n and n+1 and
		!        the instantaneous magnetic energy wmppml at n+1 in the PML layers
		call update_energy(sigma,x,yp,vpml_bxyze_np1,BCx,BCy,L_PML,dt,delta,affp,ci,&
		                   wescx,wescy,wmppml,wmpdump,nx,ny,nyp,nxe,nypmx)
		! Transform the PML B field at n+1 to Fourier space'
		isign = -1
		call wmpfft2rn(vpml_bxyze_np1,vpml_bxyz_np1,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
					   indy,kstrt,nvp,kyp,ny,nterf,ierr)
		! Update the box E and B fields
		call mpupdatefields(exyz,bxyz,vpml_exyz,vpml_bxyz_np1,ny,nye,kxp,kstrt)
    end subroutine mpmaxwel2pml_yee

	subroutine mpdumppml_implicit(f,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		implicit none
! inputs
		integer, intent(in)                      :: nx,nxe,nyp,nypmx
		integer, intent(in)                      :: BCx,BCy,kstrt,nvp
		real, intent(in)                         :: wescx,wescy
		real, dimension(2,nxe,nypmx), intent(in) :: sigma
		real, intent(in)                         :: dt
		real, intent(inout)                      :: tguard
! outputs
		real, dimension(4,nxe,nypmx), intent(inout) :: f
! locals
		integer :: i, j, dim
		!$OMP PARALLEL DO PRIVATE(i,j)
		do j=1,nyp,1
			do i = 1,nx,1
				if (wescx > 0.) then
					! yx components
					f(2,i,j) = f(2,i,j) / (1. + (sigma(1,i,j)*dt))
					! zx components
					f(3,i,j) = f(3,i,j) / (1. + (sigma(1,i,j)*dt))
				end if
				if (wescy > 0.) then
					! xy components
					f(1,i,j) = f(1,i,j) / (1. + (sigma(2,i,j)*dt))
					! zy components
					f(4,i,j) = f(4,i,j) / (1. + (sigma(2,i,j)*dt))
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!$OMP PARALLEL DO PRIVATE(dim,j)
		do j=1,nyp,1
			do dim = 1,4,1
				f(dim,nx+1,j) = f(dim,1,j)
				f(dim,nx+2,j) = f(dim,1,j)
			end do
		end do
		!$OMP END PARALLEL DO
        call wmpncguard2(f,nyp,tguard,nx,kstrt,nvp)
	end subroutine mpdumppml_implicit

	subroutine mpdumppml_yee(f,g,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		implicit none
! inputs
		integer, intent(in)                      :: nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp
		real, intent(in)                         :: wescx,wescy
		real, dimension(2,nxe,nypmx), intent(in) :: sigma
		real, intent(in)                         :: dt
		real, intent(inout)                      :: tguard
! outputs
		real, dimension(4,nxe,nypmx), intent(inout) :: f, g
! locals
		integer :: i, j, dim
    	!$OMP PARALLEL DO PRIVATE(i,j)
		do j=1,nyp,1
			do i = 1,nx,1
				if (wescx > 0.) then
					! yx components
					f(2,i,j) = f(2,i,j) / (1. + 0.5 * (sigma(1,i,j)*dt))
					g(2,i,j) = g(2,i,j) * (1. - 0.5 * (sigma(1,i,j)*dt))
					! zx components
					f(3,i,j) = f(3,i,j) / (1. + 0.5 * (sigma(1,i,j)*dt))
					g(3,i,j) = g(3,i,j) * (1. - 0.5 * (sigma(1,i,j)*dt))
				end if
				if (wescy > 0.) then
					! xy components
					f(1,i,j) = f(1,i,j) / (1. + 0.5 * (sigma(2,i,j)*dt))
					g(1,i,j) = g(1,i,j) * (1. - 0.5 * (sigma(2,i,j)*dt))
					! zy components
					f(4,i,j) = f(4,i,j) / (1. + 0.5 * (sigma(2,i,j)*dt))
					g(4,i,j) = g(4,i,j) * (1. - 0.5 * (sigma(2,i,j)*dt))
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!$OMP PARALLEL DO PRIVATE(dim,j)
		do j=1,nyp,1
			do dim = 1,4,1
				f(dim,nx+1,j) = f(dim,1,j)
				f(dim,nx+2,j) = f(dim,1,j)
				g(dim,nx+1,j) = g(dim,1,j)
				g(dim,nx+2,j) = g(dim,1,j)
			end do
		end do
		!$OMP END PARALLEL DO
        call wmpncguard2(f,nyp,tguard,nx,kstrt,nvp)
        call wmpncguard2(g,nyp,tguard,nx,kstrt,nvp)
	end subroutine mpdumppml_yee

	subroutine mpdumppml_yee2(f,wescx,wescy,sigma,dt,nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp,tguard)
		implicit none
! inputs
		integer, intent(in)                      :: nx,nxe,nyp,nypmx,BCx,BCy,kstrt,nvp
		real, intent(in)                         :: wescx,wescy
		real, dimension(2,nxe,nypmx), intent(in) :: sigma
		real, intent(in)                         :: dt
		real, intent(inout)                      :: tguard
! outputs
		real, dimension(4,nxe,nypmx), intent(inout) :: f
! locals
		integer :: i, j, dim
		!$OMP PARALLEL DO PRIVATE(i,j)
		do j=1,nyp,1
			do i = 1,nx,1
				if (wescx > 0.) then
					! yx components
					f(2,i,j) = f(2,i,j) / (1. + 0.5 * (sigma(1,i,j)*dt))
					! zx components
					f(3,i,j) = f(3,i,j) / (1. + 0.5 * (sigma(1,i,j)*dt))
				end if
				if (wescy > 0.) then
					! xy components
					f(1,i,j) = f(1,i,j) / (1. + 0.5 * (sigma(2,i,j)*dt))
					! zy components
					f(4,i,j) = f(4,i,j) / (1. + 0.5 * (sigma(2,i,j)*dt))
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!$OMP PARALLEL DO PRIVATE(dim,j)
		do j=1,nyp,1
			do dim = 1,4,1
				f(dim,nx+1,j) = f(dim,1,j)
				f(dim,nx+2,j) = f(dim,1,j)
			end do
		end do
		!$OMP END PARALLEL DO
        call wmpncguard2(f,nyp,tguard,nx,kstrt,nvp)
	end subroutine mpdumppml_yee2

	subroutine update_energy(sigma,x,yp,vpml_fxyz,BCx,BCy,L_PML,dt,delta,affp,ci,&
		                     wescx,wescy,wpml,wdump,nx,ny,nyp,nxe,nypmx)
    	implicit none
! input
		integer, intent(in)                             :: BCx, BCy, nx, ny, nyp, nxe, nypmx
		real, intent(in)                                :: wescx,wescy
		real, dimension(2,nxe,nypmx), intent(in) 		:: sigma
    	real, dimension(4,nxe,nypmx), intent(in)        :: vpml_fxyz
    	real, intent(in)                                :: L_PML, dt, affp, ci
    	real, dimension(nxe), intent(in)                :: x
        real, dimension(nypmx), intent(in)              :: yp
    	real, dimension(2), intent(in)                  :: delta
! in/output
		real, intent(inout)                             :: wpml,wdump
! locals
		integer :: i,j
		real    :: anorm, anorm2, x1, x2, y1, y2
		logical :: inside_pml_layers
!
		x1 = L_PML
		x2 = real(nx)*delta(1)-L_PML
		y1 = L_PML
		y2 = real(ny)*delta(2)-L_PML
		anorm  = dt * delta(1) * delta(2)
		anorm2 = 0.5 *delta(1) * delta(2) / affp
		wpml = 0.
		!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:wdump,wpml)
		do j=1,nyp,1
			do i=1,nx,1
				if (wescx > 0.) then
					! y component : Fy = Fyx + 0.
					wdump = wdump + ( anorm * ( (sigma(1,i,j)*vpml_fxyz(2,i,j)) ) &
									 * ( vpml_fxyz(2,i,j) + 0.              ) )
					! z component : Fz = Fzx + Fzy
					wdump = wdump + ( anorm * ( (sigma(1,i,j)*vpml_fxyz(3,i,j)) ) &
									 * ( vpml_fxyz(3,i,j) + vpml_fxyz(4,i,j)) )
				end if
				if (wescy > 0.) then
					! x component : Fx = Fxy + 0.
					wdump = wdump + ( anorm * ( (sigma(2,i,j)*vpml_fxyz(1,i,j)) ) &
									 * ( vpml_fxyz(1,i,j) + 0.              ) )
					! z component : Fz = Fzx + Fzy
					wdump = wdump + ( anorm * ( (sigma(2,i,j)*vpml_fxyz(4,i,j)) ) &
									 * ( vpml_fxyz(3,i,j) + vpml_fxyz(4,i,j)) )
				end if
				if      ((BCx == 0) .and. (BCy == 1)) then
					inside_pml_layers = ( (x(i)  .le. x1) .or. (x(i)  .ge. x2) )
				else if ((BCx == 1) .and. (BCy == 0)) then
					inside_pml_layers = ( (yp(j) .le. y1) .or. (yp(j) .ge. y2) )
				else if ((BCx == 0) .and. (BCy == 0)) then
					inside_pml_layers = ( (x(i)  .le. x1) .or. (x(i)  .ge. x2) )
					if (inside_pml_layers .eqv. .false.) then
				    	inside_pml_layers = (yp(j) .le. y1) .or. (yp(j) .ge. y2)
				    end if
				end if
				if (inside_pml_layers) then
					wpml = wpml + ( anorm2 * (  ((vpml_fxyz(1,i,j)+0.              )**2.)&
					                          + ((vpml_fxyz(2,i,j)+0.              )**2.)& 
					                          + ((vpml_fxyz(3,i,j)+vpml_fxyz(4,i,j))**2.)&
					                         ) )	
				end if
			end do
		end do
		!$OMP END PARALLEL DO
    end subroutine

    subroutine mpcopypml(vpml_bxyz,vpml_bxyz_np1,ndim,nye,kxp)
		implicit none
! input
		integer, intent(in) 			              :: ndim, nye, kxp
		complex, dimension(ndim,nye,kxp), intent(in)  :: vpml_bxyz
! output
		complex, dimension(ndim,nye,kxp), intent(out) :: vpml_bxyz_np1
! local
        integer                       :: p, q, dim
!        
        !$OMP PARALLEL DO PRIVATE(p,q,dim)
        do p=1,kxp,1
	        do q=1,nye,1
        		do dim=1,ndim,1
        			vpml_bxyz_np1(dim,q,p) = vpml_bxyz(dim,q,p)
        		end do
        	end do
        end do
	    !$OMP END PARALLEL DO
	end subroutine mpcopypml

	subroutine mpupdatefields(exyz,bxyz,vpml_exyz,vpml_bxyz,ny,nye,kxp,kstrt)
		implicit none
! input
		integer, intent(in)                          :: ny, nye, kxp, kstrt
		complex, dimension(4,nye,kxp), intent(inout) :: vpml_exyz, vpml_bxyz
! output
		complex, dimension(3,nye,kxp), intent(out)   :: exyz, bxyz
! local
		integer :: j, k
!
        !$OMP PARALLEL DO PRIVATE(j,k)
        do k=1,nye,1
        	do j=1,kxp,1
        	    ! Ex = Exy
        		exyz(1,k,j) = vpml_exyz(1,k,j)
        		! Ey = Eyx
				exyz(2,k,j) = vpml_exyz(2,k,j)
				! Ez = Ezx + Ezy
				exyz(3,k,j) = vpml_exyz(3,k,j) + vpml_exyz(4,k,j)
				! Bx = Bxy
				bxyz(1,k,j) = vpml_bxyz(1,k,j)
				! By = Byx
				bxyz(2,k,j) = vpml_bxyz(2,k,j)
				! Bz = Bzx + Bzy
				bxyz(3,k,j) = vpml_bxyz(3,k,j) + vpml_bxyz(4,k,j)
        	end do
        end do
        !$OMP END PARALLEL DO
	end subroutine mpupdatefields
	
	subroutine mpmaxwelb2pml(exy,bxy,FDTD,exy_corr,bxy_corr,&
	                         cu,ffc,affp,ci,dt,wm,tfield,nx,ny,&
     					     kstrt,kx,ky,ax,ay)
		! solves 2d maxwell's equation in fourier space for
		! transverse magnetic field for a half time step
      	implicit none
      	integer, intent(in)                      :: nx, ny, kstrt, FDTD
      	real, intent(in)                         :: affp, ci, dt, ax, ay
     	real, intent(inout)                      :: wm, tfield
      	complex, dimension(:,:,:), intent(in)    :: cu
      	complex, dimension(:,:), intent(in)      :: ffc
      	complex, dimension(:,:,:), intent(inout) :: exy, bxy
      	complex, dimension(:,:), intent(inout)   :: exy_corr, bxy_corr
      	real, dimension(:), intent(in) :: kx, ky
		! local data
      	integer :: nyv, kxp, nyhd
      	integer, dimension(4) :: itime
      	double precision :: dtime
		! extract dimensions
      	nyv = size(cu,2); kxp = size(cu,3)
      	nyhd = size(ffc,1)
		! initialize timer
      	call dtimer(dtime,itime,-1)
		! call low level procedure
		if (FDTD > 0) then
      		call MPPMAXWELB2YEEPML(exy,bxy,exy_corr,bxy_corr,&
	          	  	               cu,ffc,affp,ci,dt,wm,nx,ny,kstrt,nyv,&
     	       	                   kxp,nyhd,kx,ky,ax,ay)
     	else if (FDTD == 0) then
     		call MPPMAXWELB2SPECTRALPML(exy,bxy,&
	                                    cu,ffc,affp,ci,dt,wm,nx,ny,kstrt,nyv,&
     	                                kxp,nyhd,kx,ky,ax,ay)
     	else
     		print*,'FDTD < 0 NOT SUPPORTED'
     		stop
     	end if
		! record time
      	call dtimer(dtime,itime,1)
		tfield = tfield + real(dtime)
	end subroutine
	
	subroutine mpmaxwele2pml(exy,bxy,FDTD,exy_corr,bxy_corr,&
	                         cu,ffc,affp,ci,dt,wf,tfield,nx,ny,&
     					     kstrt,kx,ky,ax,ay)
		! solves 2d maxwell's equation in fourier space for
		! transverse magnetic field for a half time step
      	implicit none
      	integer, intent(in)                      :: nx, ny, kstrt, FDTD
      	real, intent(in)                         :: affp, ci, dt, ax, ay
     	real, intent(inout)                      :: wf, tfield
      	complex, dimension(:,:,:), intent(in)    :: cu
      	complex, dimension(:,:), intent(in)      :: ffc
      	complex, dimension(:,:,:), intent(inout) :: exy, bxy
      	complex, dimension(:,:), intent(inout)   :: exy_corr, bxy_corr
      	real, dimension(:), intent(in)           :: kx, ky
		! local data
      	integer :: nyv, kxp, nyhd
      	integer, dimension(4) :: itime
      	double precision :: dtime
		! extract dimensions
      	nyv = size(cu,2); kxp = size(cu,3)
      	nyhd = size(ffc,1)
		! initialize timer
      	call dtimer(dtime,itime,-1)
		! call low level procedure
		if (FDTD > 0) then
      		call MPPMAXWELE2YEEPML(exy,bxy,exy_corr,bxy_corr,&
	                      	    cu,ffc,affp,ci,dt,wf,nx,ny,kstrt,nyv,&
     	                    	kxp,nyhd,kx,ky,ax,ay)
     	else if (FDTD == 0) then
     		call MPPMAXWELE2SPECTRALPML(exy,bxy,&
	                                    cu,ffc,affp,ci,dt,wf,nx,ny,kstrt,nyv,&
     	                                kxp,nyhd,kx,ky,ax,ay)
     	else
     		print*,'FDTD < 0 NOT SUPPORTED'
     		stop
     	end if
		! record time
      	call dtimer(dtime,itime,1)
		tfield = tfield + real(dtime)
	end subroutine

end module
