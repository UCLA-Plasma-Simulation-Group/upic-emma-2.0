!-----------------------------------------------------------------------
!
      module modmpinit2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpinit2.f
! mpdcomp2 determines spatial decomposition for uniform distribution
!          calls PDNICOMP2L
! mpdistr2 calculates initial particle co-ordinates and velocities
!          with uniform density and maxwellian velocity with drift
!          for 2d code
!          calls PDISTR2
! mpdistr2h calculates initial particle co-ordinates and velocities
!           with uniform density and maxwellian velocity with drift
!           for 2-1/2d code
!           calls PDISTR2H
! mpdblkp2 finds the maximum number of particles in each tile
!          calls PPDBLKP2L
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 9, 2016
!
      use libmpinit2_h
      use libprofile_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
!      subroutine mpdcomp2(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp)            ! M. Touati
      subroutine mpdcomp2(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,delta)       ! M. Touati
! determines spatial decomposition for uniform density distribution
      implicit none
      integer, intent(in) :: ny, kstrt, nvp
      integer, intent(inout) :: nyp, noff, nypmx, nypmn
      real, dimension(:), intent(in) :: delta                                  ! M. Touati
      real, dimension(:), intent(inout) :: edges
! local data
      integer :: idps
! extract dimensions
      idps = size(edges,1)
! call low level procedure
!	  call PDNICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps)            ! M. Touati
      call PDNICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,idps,delta)      ! M. Touati
! check for errors
      if (nypmn < 1) then
         if (kstrt==1) then
            write (*,*) 'combination not supported nvp, ny =',nvp,ny
         endif
         call PPEXIT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpdistr2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,&
!     &ny,ipbc,ierr)                                                           ! M. Touati
      subroutine mpdistr2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,&
     &ny,ipbc,x,y,ierr)                                                        ! M. Touati
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
      implicit none
      integer, intent(in) :: nps, npx, npy, nx, ny, ipbc
      integer, intent(inout) :: npp, ierr
      real, intent(in) :: vtx, vty, vdx, vdy
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
      real, dimension(:), intent(in) :: x, y                                   ! M. Touati
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
!      call PDISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,    &
!     &idimp,npmax,idps,ipbc,ierr)                                             ! M. Touati
      call PDISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,    &
     &idimp,npmax,idps,ipbc,x,y,ierr)                                          ! M. Touati
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpdistr2h(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,  &
!     &npx,npy,nx,ny,ipbc,ierr)                                                ! M. Touati
	  subroutine mpdistr2h(relativity, mass_real,ci, part, edges, npp, nps,&
	 &vtx, vty, vtz, vdx, vdy, vdz, npx, npy, nx, ny, ipbc, x, y, densx, densy, ierr)     ! M. Touati
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
      implicit none
      integer, intent(in) :: relativity
      real, intent(in)    :: mass_real, ci
      integer, intent(in) :: nps, npx, npy, nx, ny, ipbc
      integer, intent(inout) :: npp, ierr
      real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:), intent(in) :: edges
      real, dimension(nx+2), intent(in) :: x
      real, dimension(ny+2), intent(in) :: y
      character (len=*), intent(in)  :: densx, densy
! local data
      integer :: idimp, npmax, idps
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      idps = size(edges,1)
! call low level procedure
      if (relativity == 1) then
      	call PDISTR2HR(ci,mass_real,part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
       &npx,npy,nx,ny,idimp,npmax,idps,ipbc,x,y,densx,densy,ierr)                          ! M. Touati
      else
!      call PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy, &
!     &nx,ny,idimp,npmax,idps,ipbc,ierr)                                       ! M. Touati
      	call PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
       &npx,npy,nx,ny,idimp,npmax,idps,ipbc,x,y,ierr)                          ! M. Touati
      end if
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpdblkp2(part,kpic,npp,noff,nppmx,mx,my,mx1,irc)             ! M. Touati
      subroutine mpdblkp2(part,kpic,npp,noff,nppmx,mx,my,mx1,delta,irc)        ! M. Touati
! finds the maximum number of particles in each tile
      implicit none
      integer, intent(in) :: mx, my, mx1, npp
      integer, intent(in) :: noff
      integer, intent(inout) :: nppmx, irc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:), intent(in) :: delta
      integer, dimension(:), intent(inout) :: kpic
! local data
      integer :: idimp, npmax, mxyp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      mxyp1 = size(kpic,1)
! call low level procedure
!      call PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mx1,    &
!     &mxyp1,irc)                                                              ! M. Touati
      call PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mx1,    &
     &mxyp1,delta,irc)                                                         ! M. Touati
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpdblkp2 error, irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
	  subroutine PDISTR2HR(ci,mass_real,part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
     &npx,npy,nx,ny,idimp,npmax,idps,ipbc,x,y,densx,densy,ierr)
     	implicit none
     	real, intent(in)                            :: mass_real,ci
	    integer, intent(in)                         :: nps, npx, npy, nx, ny
	    integer, intent(in)                         :: idimp, npmax, idps, ipbc
	    integer, intent(inout)                      :: npp, ierr
      	real, intent(in)                            :: vtx, vty, vtz, vdx, vdy, vdz
        real, dimension(idimp,npmax), intent(inout) :: part
        real, dimension(idps), intent(in)           :: edges
        real, dimension(nx+2), intent(in)           :: x
        real, dimension(ny+2), intent(in)           :: y
        character (len=*), intent(in)               :: densx, densy
! local data
      	integer :: j, k, npt, npxyp
      	real    :: edgelx, edgely, at1, at2, xt, yt, vxt, vyt, vzt, mu
      	real(8) :: dnpx, dnpxy, dt1, pi, theta, phi, p
        integer, dimension(1) :: ierr1, iwork1
        real(8), dimension(4) :: sum4, work4
        integer :: N
        real(8), dimension(:,:), allocatable :: Fp
        real(8)  :: rand
        real, dimension(:), allocatable :: xp, yp
        pi = 3.141592653589793
      	ierr = 0
! particle distribution constant
      	dnpx = real(npx,8)
! set boundary values
      	edgelx = x(1)
      	edgely = y(1)
      	at1 = (x(nx+1)-x(1))/real(npx)
      	at2 = (y(ny+1)-y(1))/real(npy)
      	if (ipbc.eq.2) then
        	edgelx = x(2)
        	edgely = y(2)
        	at1 = (x(nx)-x(2))/real(npx)
        	at2 = (y(ny)-y(2))/real(npy)
      	else if (ipbc.eq.3) then
        	edgelx = x(2)
        	at1 = (x(nx)-x(2))/real(npx)
      	endif
	mu = (ci*vtx)**(-2.)
        allocate(xp(npx), yp(npy))
        call getdensity1d(densx, 'x1', edgelx, x(2)-x(1), nx, npx, 32, xp)
        call getdensity1d(densy, 'x2', edgely, y(2)-y(1), ny, npy, 32, yp)
        write(*,*)'calculate constants -- done',nx,ny,npx,npy
        if (Bessel_Kn(2,mu)==0.) then
! uniform density profile

      		do k = 1, npy
      			yt = yp(k) !edgely + at2*(real(k) - 0.5)
      			do j = 1, npx
     				xt = xp(j) !edgelx + at1*(real(j) - 0.5)
! maxwellian velocity distribution
				rand  = ranorm()
      				vxt   = vtx*rand
      				rand  = ranorm()
      				vyt   = vty*rand
      				rand  = ranorm()
      				vzt   = vtz*rand
! tests :
      		                if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
        				npt = npp + 1
         				if (npt.le.npmax) then
            				part(1,npt) = xt
            				part(2,npt) = yt
            				part(3,npt) = vxt
            				part(4,npt) = vyt
            				part(5,npt) = vzt
            				npp = npt
         				else
            				ierr = ierr + 1
         				endif
      				endif
   				end do
   			end do
		else
			N = max(100000,int(npp/(10*nx*ny)))
			allocate(Fp(2,N))
			call equipartition_function(mu, ci, N, Fp)
! uniform density profile
			do k = 1, npy, 1
      			yt = yp(k) !edgely + at2*(real(k) - 0.5)
      			do j = 1, npx, 1
      				xt = xp(j) !edgelx + at1*(real(j) - 0.5)
! Maxwell-Juttner distribution :
	  				if ( (vtx .ne. vty) .or. (vtx .ne. vtz) .or. (vty .ne. vtz) ) then
	  					stop
	  				else if ( (vtx .eq. 0.) .and. (vty .eq. 0.) .and. (vtz .eq. 0.) ) then 
	  					vxt = vdx
      					vyt = vdy
      					vzt = vdz	
	  				else if ( (vdx .ne. 0.) .or. (vdy .ne. 0.) .or. (vdz .ne. 0.) ) then
	  					stop
	  				else
! - isotropic and without drift
						rand = randum()
	  					theta = Inverse_Ftheta(rand)
	  					rand = randum()
	  					phi   = 2. * pi * rand
	  					rand = randum()
	  					call Inverse_Fp(N, Fp,rand, p)
      					vxt = p * sin(theta) * cos(phi)
      					vyt = p * sin(theta) * sin(phi)
      					vzt = p * cos(theta)
      				end if
! tests :
      				if ((yt.ge.edges(1)).and.(yt.lt.edges(2))) then
         				npt = npp + 1
         				if (npt.le.npmax) then
            				part(1,npt) = xt
            				part(2,npt) = yt
            				part(3,npt) = vxt
            				part(4,npt) = vyt
            				part(5,npt) = vzt
            				npp = npt
         				else
            				ierr = ierr + 1
         				endif
      				endif
      			end do
      		end do
      		deallocate(Fp)
		end if
        deallocate(xp)
        deallocate(yp)
      	npxyp = 0
! add correct drift
      	sum4(1) = 0.0d0
      	sum4(2) = 0.0d0
      	sum4(3) = 0.0d0
     	do j = nps, npp
      		npxyp = npxyp + 1
      		sum4(1) = sum4(1) + part(3,j)
      		sum4(2) = sum4(2) + part(4,j)
      		sum4(3) = sum4(3) + part(5,j)
      	end do
      	sum4(4) = npxyp
      	call PPDSUM(sum4,work4,4)
      	dnpxy = sum4(4)
      	ierr1(1) = ierr
      	call PPIMAX(ierr1,iwork1,1)
      	ierr = ierr1(1)
        write(*,*)'dnpxy=',dnpxy
      	dt1 = 1.0e0/dnpxy
      	sum4(1) = sum4(1)/dnpxy
      	sum4(2) = sum4(2)/dnpxy
      	sum4(3) = sum4(3)/dnpxy
      	do j = nps, npp
      		part(3,j) = part(3,j) - sum4(1)
      		part(4,j) = part(4,j) - sum4(2)
      		part(5,j) = part(5,j) - sum4(3)
   		end do
! process errors
      	dnpxy = dnpxy - dnpx*dble(npy)
      	if (abs(dnpxy).ge.0.1d0) ierr = dnpxy
	  end subroutine
!
!-----------------------------------------------------------------------
	  elemental function Inverse_Ftheta(u)
	  	implicit none
	  	real(8), intent(in) :: u
	  	real(8)             :: Inverse_Ftheta
	  	Inverse_Ftheta = acos(1.-(2.*u))
	  end function Inverse_Ftheta
!
!-----------------------------------------------------------------------
	  elemental function Maxwell_Juttner_p(mu, ci, p)
! This function computes the Maxwell-Juttner distribution function amplitude-component
! f_p(p) = (mu / K2(mu)) * (p**2) * exp(- mu * sqrt( 1 + (p/mc)**2 ) ) / (mc)**3
! << f(px, py, px) dpx dpy dpz = f_p(p) dp f_theta(theta) dtheta f_phi(phi) dphi
! mu = m c^2 / kB T is the inverse temperature normalized to the rest mass energy
! mass_real is the real mass of particles contained in each macro particle
! ci is the inverse velocity normalized to Delta omega_p
! p is in unit of (m_part Delta omega_p) where m_part is the real mass of particles
! p_norm is in unit of (m_part c)
! Maxwell_Juttner_p is in unit of (1 / m_part Delta omega_p)
		implicit none
		real, intent(in)    :: mu, ci
		real(8), intent(in) :: p
		real(8)             :: Maxwell_Juttner_p
! local
		real(8) :: p_norm
		p_norm = p * ci
		Maxwell_Juttner_p = real(mu,8) / Bessel_Kn(2,mu)
		Maxwell_Juttner_p = Maxwell_Juttner_p * exp(-real(mu,8)*dsqrt(1.e0 + (p_norm**2.e0)))
		Maxwell_Juttner_p = Maxwell_Juttner_p * (p_norm**2.e0) * ci
	  end function Maxwell_Juttner_p
!
!-----------------------------------------------------------------------	  
	  elemental function p_max(mu, ci)
! This function computes the normalized value of the momentum that maximizes 
! f_p(p) = (mu / K2(mu)) * (p**2) * exp(- mu * sqrt( 1 + (p/mc)**2 ) ) / (mc)**3
! << f(px, py, px) dpx dpy dpz = f_p(p) dp f_theta(theta) dtheta f_phi(phi) dphi
! -> p_max is the most probable momentum amplitude in units of mass_real Delta omega_p
! mu = m c^2 / kB T is the inverse temperature normalized to the rest mass energy
! mass_real is the real mass of particles contained in each macro particle
! ci is the inverse velocity normalized to Delta omega_p
! p_max is in unit of (m_e \Delta \omega_p)
	  	implicit none
	  	real, intent(in)    :: mu, ci
	  	real(8) :: p_max
		p_max = dsqrt( 2.e0 + ( 2.e0 * dsqrt( 1.e0 + (real(mu,8)**2.e0) ) ) )
		p_max = p_max / ( mu * ci )	  	
	  end function p_max
!
!-----------------------------------------------------------------------
	  pure subroutine equipartition_function(mu, ci, N, Fp)
! This function computes array containing the (p(i),Fp(i)) where 
!                     /p(i)
! Fp(i) = F_p(p(i)) = | f_p(p') dp'
!                     /0
! where f_p is the Maxwell-Juttner distribution function amplitude-component
! f_p(p) = (mu / K2(mu)) * (p**2) * exp(- mu * sqrt( 1 + (p/mc)**2 ) ) / (mc)**3
! << f(px, py, px) dpx dpy dpz = f_p(p) dp f_theta(theta) dtheta f_phi(phi) dphi.
! p(i) goes from 0 to 10. p_max where f_p(p_max) = max{f_p(p), p > 0}
		real, intent(in)                     :: mu, ci
		integer, intent(in)                  :: N
	  	real(8), dimension(2,N), intent(out) :: Fp
! local
		real(8) :: pmax, dp
		integer :: i
		pmax = 5.e0 * p_max(mu, ci) 
		dp = pmax / real(N,8)
		Fp(1,1) = 0.e0
		Fp(2,1) = Maxwell_Juttner_p(mu, ci, Fp(1,1)) * dp
		do i=1,N-1,1
		    Fp(1,i+1) = Fp(1,i) + dp
			Fp(2,i+1) = Fp(2,i) + ( Maxwell_Juttner_p(mu, ci, Fp(1,i+1)) * dp)
		end do
	  end subroutine equipartition_function 
!
!-----------------------------------------------------------------------
	  pure subroutine Inverse_Fp(N,Fp,u,Fpm1_of_u)
! This function computes the inverse function F_p^-1(u) of 
!                   /p
!          F_p(p) = | f_p(p') dp' = u in [0..1]
!                   /0
! where f_p is the Maxwell-Juttner distribution function amplitude-component
! f_p(p) = (mu / K2(mu)) * (p**2) * exp(- mu * sqrt( 1 + (p/mc)**2 ) ) / (mc)**3
! << f(px, py, px) dpx dpy dpz = f_p(p) dp f_theta(theta) dtheta f_phi(phi) dphi.
! mu = m c^2 / kB T is the inverse temperature normalized to the rest mass energy
! mass_real is the real mass of particles contained in each macro particle
! ci is the inverse velocity normalized to Delta omega_p
! p is the momentum normalized to mass_real Delta omega_p
		implicit none
		integer, intent(in)                 :: N
		real(8), dimension(2,N), intent(in) :: Fp
		real(8), intent(in)                 :: u
		real(8), intent(out)                :: Fpm1_of_u
! local
		integer :: i, il, ir, im
		logical :: not_found
		not_found = .true.
		il = 1
		ir = N
		im = (il+ir)/2
		do while (not_found)
			if (Fp(2,im) <= u) then
				il = im
			else
				ir = im
			end if
			im = (il+ir)/2
			not_found = (ir - il) > 1   
		end do
		Fpm1_of_u = Fp(1,im)
	  end subroutine Inverse_Fp
!
!-----------------------------------------------------------------------
	  elemental function Bessel_I0(x)
! This function computes the modified Bessel function of the first kind I_0(x)
! according to an algorithm given in Numerical Recipes.
! written by Michael J Touati, UCLA, October 2016
      	implicit none
      	real(8), intent(in)   :: x
      	real(8)               :: Bessel_I0
! local
		real(8)               :: y, a, b
		real(8), dimension(7) :: P
		real(8), dimension(9) :: Q
      	P(:) = (/1.e0, 3.5156229e0, 3.0899424e0, 1.2067429e0, &
      			 0.2659732e0, 0.360768e-1, 0.45813e-2/)
      	Q(:) = (/0.39894228e0, 0.1328592e-1,  &
      			 0.225319e-2, -0.157565e-2, 0.916281e-2, -0.2057706e-1, &
      			 0.2635537e-1, -0.1647633e-1, 0.392377e-2/)
      	if (abs(x).lt.3.75e0) then
      		y = (x/3.75e0)**2
      		Bessel_I0 = P(1) + y * ( P(2) +&
      		                   y * ( P(3) +&
      		                   y * ( P(4) +&
      		                   y * ( P(5) +&
      		                   y * ( P(6) + y * P(7) ) ) ) ) )
      	else
      		a = abs(x)
      		y = 3.75e0 / a
     		b = exp(a) / dsqrt(a)
      		a = Q(1) + y * ( Q(2) +&
      		           y * ( Q(3) +&
      		           y * ( Q(4) +&
      		           y * ( Q(5) +&
      		           y * ( Q(6) +&
      		           y * ( Q(7) +&
      		           y * ( Q(8) + y * Q(9) ) ) ) ) ) ) )
      		Bessel_I0 = a * b
      	end if
      end function Bessel_I0
!
!-----------------------------------------------------------------------	
	  elemental function Bessel_I1(x)
! This function computes the modified Bessel function of the first kind I_1(x)
! according to an algorithm given in Numerical Recipes.
! written by Michael J Touati, UCLA, October 2016
      	implicit none
      	real(8), intent(in)  :: x
      	real(8)              :: Bessel_I1
! local
		real(8)               :: y, a, b
      	real(8), dimension(7) :: P 
      	real(8), dimension(9) :: Q
       	P(:) = (/0.5e0, 0.87890594e0, 0.51498869e0,  &
      			 0.15084934e0, 0.2658733e-1, 0.301532e-2, 0.32411e-3/)
      	Q(:) = (/0.39894228e0, -0.3988024e-1, &
      			 -0.362018e-2, 0.163801e-2, -0.1031555e-1, 0.2282967e-1, &
      			 -0.2895312e-1, 0.1787654e-1, -0.420059e-2/)
      	if (abs(x) .lt. 3.75e0) then
      		y = (x/3.75e0)**2
      		Bessel_I1 = x * ( P(1) + y * ( P(2) +&
      		                         y * ( P(3) +&
      		                         y * ( P(4) +&
      		                         y * ( P(5) +&
      		                         y * ( P(6) + y * P(7) ) ) ) ) ) )
      	else
      		a = abs(x)
      		y = 3.75e0 / a
      		b = exp(a) / dsqrt(a)
      		a = Q(1) + y * ( Q(2) +&
      		           y * ( Q(3) +&
      		           y * ( Q(4) +&
      		           y * ( Q(5) +&
      		           y * ( Q(6) +&
      		           y * ( Q(7) +&
      		           y * ( Q(8) + y * Q(9) ) ) ) ) ) ) )
      		Bessel_I1 = a * b
      	end if
      end function Bessel_I1
!
!-----------------------------------------------------------------------	
	  elemental function Bessel_In(n,x)
! This function computes the modified Bessel function of the first kind 
! I_n(x) with an error eps < 10^-15 according to the formula :
! I_n(x) = (x/2)^n * sum( (x^2/4)^k / (Gamma(k+1) * Gamma(n+k+1)), 
!                         k=0..infinity )       
! where Gamma is the Gamma function ( Gamma(k+1) = k! for integers k)
! see for ex. :
! http://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
! written by Michael J Touati, UCLA, October 2016
      	implicit none
      	integer, intent(in) :: n
      	real, intent(in)    :: x
      	real(8)             :: Bessel_In
! local
      	real(8)             :: z, nu, epsilon, a, b, c, k
        logical             :: no_converge
      	z  = real(x,8)
      	nu = real(n,8)
      	if (n==0) then
			Bessel_In = Bessel_I0(z)
		else if (n==1) then
			Bessel_In = Bessel_I1(z)
		else
      		Bessel_In = ((z/2._8)**nu) / Gamma(nu+1._8)
      		k = 1._8
      		no_converge = .true.
      		do while (no_converge)
      			a = nu + (2._8*k)
      			b = k + 1._8
      			c = n + k + 1._8
      			epsilon = ((z/2._8)**a) / ( Gamma(b) * Gamma(c) )
      			Bessel_In = Bessel_In + epsilon 
      			k = k + 1._8
      			no_converge = abs(epsilon) .gt. 1.e-15_8
      			if (Bessel_In .gt. 1.e32_8) then
      				no_converge = .false.
      				Bessel_In = 1.e32_8
      			else
      				no_converge = abs(epsilon) .gt. 1.e-15_8
      			end if
      		end do
      	end if
      end function Bessel_In
!
!-----------------------------------------------------------------------      
      elemental function Bessel_K0(x)
! This function computes the modified Bessel function of the second kind K_0(x)
! according to an algorithm given in Numerical Recipes.
! written by Michael J Touati, UCLA, October 2016
      	implicit none
      	real(8), intent(in) :: x
      	real(8)             :: Bessel_K0
! local
      	real(8) :: y, a
      	real(8), dimension(7) :: P, Q 
      	P(:) = (/-0.57721566e0,0.42278420e0,0.23069756e0, &
                 0.3488590e-1,0.262698e-2,0.10750e-3,0.74e-5/)
        Q(:) = (/1.25331414e0,-0.7832358e-1,0.2189568e-1, & 
                -0.1062446e-1,0.587872e-2,-0.251540e-2,0.53208e-3/)
      	if (x.eq.0.e0) then
      		Bessel_K0 =1.e32
      	else if (x.le.2.e0) then
      		y         = x * x / 4.e0
      		a         = -log(x/2.e0) * Bessel_I0(x)
        	Bessel_K0 = a + ( P(1) + y * &
        	                ( P(2) + y * &
        	                ( P(3) + y * &
        	                ( P(4) + y * &
        	                ( P(5) + y * &
        	                ( P(6) + y * P(7) ) ) ) ) ) )
      	else
      		y  = 2.e0/x
      		a  = exp(-x) / dsqrt(x)
      		Bessel_K0 = a * ( Q(1) + y * &
      		                ( Q(2) + y * &
      		                ( Q(3) + y * &
      		                ( Q(4) + y * &
      		                ( Q(5) + y * &
      		                ( Q(6) + y * Q(7) ) ) ) ) ) )
      	end if
	  end function Bessel_K0
!
!-----------------------------------------------------------------------      
      elemental function Bessel_K1(x)
! This function computes the modified Bessel function of the second kind K_1(x)
! according to an algorithm given in Numerical Recipes.
! written by Michael J Touati, UCLA, October 2016
      	implicit none
      	real(8), intent(in) :: x
      	real(8)             :: Bessel_K1
! local
		real(8) ::y, a
		real(8), dimension(7) :: P, Q
     	P(:) = (/1.e0, 0.15443144e0, -0.67278579e0,  &
      			 -0.18156897e0, -0.1919402e-1, -0.110404e-2, -0.4686e-4/)
      	Q(:) = (/1.25331414e0, 0.23498619e0, -0.3655620e-1, &
      			 0.1504268e-1, -0.780353e-2, 0.325614e-2, -0.68245e-3/)
      	if (x == 0.e0) then	
      		Bessel_K1 = 1.e32
      	else if (x .le. 2.e0) then
      		y = x * x / 4.e0
      		a = log(x/2.e0) * Bessel_I1(x)
      		Bessel_K1 = a + (1.e0/x) * ( P(1) + y * ( P(2) +&
      		                                    y * ( P(3) +&
      		                                    y * ( P(4) +&
      		                                    y * ( P(5) +&
      		                                    y * ( P(6) + y * P(7) ) ) ) ) ) )
      	else
      		y = 2.e0 / x
      		a = exp(-x) / dsqrt(x)
      		Bessel_K1 = a * ( Q(1) + y * ( Q(2) +&
      		                         y * ( Q(3) +&
      		                         y * ( Q(4) +&
      		                         y * ( Q(5) +&
      		                         y * ( Q(6) + y * Q(7) ) ) ) ) ) )
      	end if
      end function Bessel_K1
!
!-----------------------------------------------------------------------      
	  elemental function Bessel_Kn(n,x)
! This function computes the modified Bessel function of the second kind 
! K_n(x) = ( pi / 2 ) * ( I_{-n}(x) - I_n(x) ) / sin( n * pi )
! where I_n(x) is the modified Bessel function of the first kind Bessel_In(n,x)
! see for ex. :
! http://mathworld.wolfram.com/ModifiedBesselFunctionoftheSecondKind.html
! Since the above formula can only be computed with sufficient accuracy if z 
! is relatively small when n is close to an integer, we use the recursive relation 
! K_{n+1}(x) = (2 n / x) K_n(x) + K_{n-1}(x), knowing K_0(x) and K_1(x) given by
! algorithms provided by Numerical Recipes. 
! written by Michael J Touati, UCLA, October 2016
      	implicit none
      	integer, intent(in) :: n
      	real, intent(in)    :: x
      	real(8)             :: Bessel_Kn
! local
      	integer             :: k
      	real(8)             :: z, Bessel_Kp, Bessel_K, Bessel_Km
		z = real(x, 8)
		if (n==0) then
			Bessel_Kn = Bessel_K0(z)
		else if (n==1) then
			Bessel_Kn = Bessel_K1(z)
		else
			if (x .eq. 0.e0) then
      			Bessel_Kn = 1.e32
      		else
      			Bessel_K  = Bessel_K1(z)
      			Bessel_Km = Bessel_K0(z)
      			do k=1,n-1,1
      				Bessel_Kp = Bessel_Km + ( real(k) * 2.e0 * Bessel_K / z )
      				Bessel_Km = Bessel_K
      				Bessel_K = Bessel_Kp
      			end do
      			Bessel_Kn = Bessel_K
      		end if
		end if
      end function Bessel_Kn
!-----------------------------------------------------------------------      
      end module
