!-----------------------------------------------------------------------
!
      module ompplib2
!
! ompmove2 reorder particles by tile with OpenMP and MPI
!          calls mporderf2a or mporder2a, and mpmove2, mporder2b
! wmpfft2r performs 2d real/complex FFT for scalar data,
!          moving data between uniform and non-uniform partitions
!          calls mpfmove2 and mpfft2r
! wmpfft2rn performs 2d real/complex FFT for n component vector data,
!           moving data between uniform and non-uniform partitions
!           calls mpfnmove2 and mpfft2rn
! wmpcguard2 copy scalar guard cells to local and remote partitions
! wmpncguard2 copy vector guard cells to local and remote partitions
!             calls mpncguard2, mpcguard2x
! wmpaguard2 add scalar guard cells from local and remote partitions
!            calls mpaguard2x, mpnaguard2
! wmpnacguard2 add vector guard cells from local and remote partitions
!              calls mpacguard2x, mpnacguard2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 13, 2016
!
! GRID calculates spatial and spatial frequency grids :
! _ for the full simulation box : x, y, kx and ky
! _ for each MPI partition : x, yp, kpx and ky
! written by Michael J. Touati, UCLA, 2016, March

      use modmpsort2
      use modmpfft2
      use modmpgard2
      use mppmod2, only: mpmove2, mpfmove2, mpfnmove2, mpcguard2,       &
     &mpncguard2, mpnaguard2, mpnacguard2
      implicit none
!
! ppbuff = buffer array for reordering tiled particle array
      real, dimension(:,:,:), allocatable :: ppbuff
      integer :: szpbuf = 0
! sbufl/sbufr = particle buffers sent to nearby processors
! rbufl/rbufr = particle buffers received from nearby processors
      real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer :: szbufs = 0
! ncll/nclr/mcll/mclr = number offsets send/received from processors
      integer, dimension(:,:), allocatable :: ncll, nclr, mcll, mclr
      integer :: sznbufs = 0
      save
!
      private :: ppbuff, szpbuf
      private :: sbufl, sbufr, rbufl, rbufr, szbufs
      private :: ncll, nclr, mcll, mclr, sznbufs
!
      contains
!
!-----------------------------------------------------------------------
      subroutine GRID(x,y,kx,ky,yp,kpx,nxe,nye,nypmx,kxp,noff,moff,nx,ny,delta)
! This subroutine calculates spatial and spatial frequency grids :
! _ for the full simulation box : x, y, kx and ky
! _ for each MPI partition : x, yp, kpx and ky
! depending on the cell size in the x and y-direction (delta(1) and delta(2) resp.)
! input  : nxe, nye, nypmx, kxp, noff, moff, nx, ny, delta
! output : x, y, kx, ky, yp, kpx
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      integer, intent(in)                 :: nxe, nye, nypmx, kxp, noff, moff, nx, ny
      real, dimension(2), intent(in)      :: delta
      real, dimension(nxe), intent(out)   :: x, kx
      real, dimension(nye), intent(out)   :: y, ky
      real, dimension(nypmx), intent(out) :: yp
      real, dimension(kxp), intent(out)   :: kpx
! local data
      integer :: i, j
      real dnx, dny, dnz
      dnx = 6.28318530717959/real(nx)
      dnx = dnx / delta(1)
      do i = 1, nxe
      	x(i) = real(i-1)*delta(1)
      	kx(i) = dnx*real(i-1)
      end do
      do i=1, kxp
      	kpx(i) = kx(i+moff)
      end do
      dny = 6.28318530717959/real(ny)
      dny = dny / delta(2)
      do j = 1, nye
      	y(j) = real(j-1)*delta(2)
      	ky(j) = dny*real(j-1)
      end do
      do j=1,nypmx
      	yp(j) = y(j+noff)
      end do
      end subroutine
!
!-----------------------------------------------------------------------
	subroutine setup_coef_center_scheme(s,C)
! This subroutine setup the centered scheme coefficients array for discretizing 
! a first order derivative with an order of acuracy 2*s :
! (dB/dx)^i = (1/dx)*sum(C(m)*(B^{i+(m-1/2)}-B^{i-(m-1/2)}),m = 1..s) + O[dx^(2s)]
		implicit none
		integer, intent(in)               :: s
		real, dimension(1:s), intent(out) :: C
! locals
		integer :: m, n
!
		select case (s)
			case (1) ! order 2
				C(1)=1.
			case (2) ! order 4
				C(1)=9./8.       ; C(2)=-1./24.
			case (3) ! order 6
				C(1)=75./64.     ; C(2)=-25./384. 
				C(3)=3./640.
			case (4) ! order 8
				C(1)=1225./1024. ; C(2)=-245./3072.
				C(3)=  49./5120. ; C(4)=-  5./7168.
			case (5) ! order 10
				C(1)=7281626625./6072165376. ; C(2)=-495978735./6072165376.
				C(3)=315650601. /30360826880.; C(4)=- 38430045./42505157632.
				C(5)=      1225./53368641.
			case default ! order > 10
				do m=1,s,1
					C(m) = 1./real(2*m-1)
					do n=1,s,1
						if (n .ne. m) then
						C(m) = C(m) * (real(2*n-1)**2.) &
						     / ( (real(2*n-1)**2.) - (real(2*m-1)**2.) )
						end if
					end do
				end do
		end select  
	end subroutine
!
!-----------------------------------------------------------------------	
	subroutine setup_bracket(stencil,nx,ny,delta,kx,ky)
! This subroutine setup the wave vectors "brackets" kx_br and ky_br corresponding
! a Finite difference centered scheme for discretizing a first order derivative :
! i kx_br^p*IDFT{B^i}^p = IDFT{(dB/dx)^i}^p where IDFT stands for Inverse Discrete 
! Fourier Transform and where the derivative in the real space reads : 
! (dB/dx)^i = (1/dx)*sum(C(m)*(B^{i+(m-1/2)}-B^{i-(m-1/2)}),m = 1..s) + O[dx^(2s)]
! with the coefficients C(m) given by the subroutine setup_coef_center_scheme
		implicit none
		integer, intent(in)                  :: stencil, nx, ny
		real, dimension(2), intent(in)       :: delta
		real, dimension(nx+2), intent(inout) :: kx
		real, dimension(ny+2), intent(inout) :: ky
! locals
		integer                  :: i,j,m
		real, dimension(stencil) :: C
		real, dimension(nx+2)    :: kx_br
		real, dimension(ny+2)    :: ky_br
!
		call setup_coef_center_scheme(stencil,C)
		do i=1,nx+2,1
			kx_br(i) = 0.
			do m=1,stencil,1
				kx_br(i) = kx_br(i) &
				         + ( ( C(m) * sin( real(2*m-1) * 0.5 * kx(i) * delta(1) ) ) &
				           / ( 0.5 * delta(1) ) )
			end do
		end do
		do j=1,ny+2,1
			ky_br(j) = 0.
			do m=1,stencil,1
				ky_br(j) = ky_br(j) &
				         + ( ( C(m) * sin( real(2*m-1) * 0.5 * ky(j) * delta(2) ) ) &
				           / ( 0.5 * delta(2) ) )
			end do
		end do
		kx = kx_br
		ky = ky_br
	end subroutine setup_bracket
!
!-----------------------------------------------------------------------                        
!      subroutine ompmove2(ppart,kpic,ncl,ihole,noff,nyp,tsort,tmov,kstrt&
!     &,nvp,nx,ny,mx,my,npbmx,nbmax,mx1,list,irc)                              ! M. Touati
      subroutine ompmove2(ppart,kpic,ncl,ihole,noff,nyp,tsort,tmov,kstrt&
     &,nvp,nx,ny,mx,my,npbmx,nbmax,mx1,list,x,y,irc)                           ! M. Touati
! reorder particles by tile with OpenMP and MPI
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: kstrt, nvp, nx, ny, mx, my, npbmx, nbmax
      integer, intent(in) :: mx1
      integer, intent(in) :: noff, nyp
      real, intent(inout) :: tsort, tmov
      integer, intent(inout) :: irc
      logical, intent(in) :: list
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      real, dimension(:), intent(in) :: x, y                                   ! M. Touati
! local data
      integer :: idimp, nppmx, mxyp1, ntmax
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mxyp1 = size(ppart,3)
      ntmax = size(ihole,2) - 1
! check if required size of buffer has increased
      if (szpbuf < idimp*npbmx*mxyp1) then
         if (szpbuf /= 0) deallocate(ppbuff)
! allocate new buffer
         allocate(ppbuff(idimp,npbmx,mxyp1))
         szpbuf = idimp*npbmx*mxyp1
      endif
! check if required size of buffers has increased
      if (szbufs < idimp*nbmax) then
         if (szbufs /= 0) deallocate(sbufl,sbufr,rbufl,rbufr)
! allocate new buffers
         allocate(sbufl(idimp,nbmax),sbufr(idimp,nbmax))
         allocate(rbufl(idimp,nbmax),rbufr(idimp,nbmax))
         szbufs = idimp*nbmax
      endif
! check if required size of buffers has increased
      if (sznbufs < 3*mx1) then
         if (sznbufs /= 0) deallocate(ncll,nclr,mcll,mclr)
! allocate new buffers
         allocate(ncll(3,mx1),nclr(3,mx1),mcll(3,mx1),mclr(3,mx1))
         sznbufs = 3*mx1
      endif
!
! first part of particle reorder on x and y cell with mx, my tiles:
! list of particles leaving tile already calculated by push
      if (list) then
! updates: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
         call mporderf2a(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,  &
     &tsort,irc)
         if (irc /= 0) then
            write (*,*) kstrt,'mporderf2a error: ntmax, irc=',ntmax,irc
         endif
! calculate list of particles leaving tile
      else
! updates ppart, ppbuff, sbufl, sbufr, ncl, ihole, ncll, nclr, irc
!         call mporder2a(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,   &
!     &nclr,noff,nyp,tsort,nx,ny,mx,my,irc)                                    ! M. Touati
         call mporder2a(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,   &
     &nclr,noff,nyp,tsort,nx,ny,mx,my,x,y,irc)                                 ! M. Touati
         if (irc /= 0) then
            write (*,*) kstrt,'mporder2a error: ntmax, irc=',ntmax,irc
         endif
      endif
      if (irc /= 0) then
         call PPABORT()
         stop
      endif
!
! move particles into appropriate spatial regions with MPI:
! updates rbufr, rbufl, mcll, mclr
      call mpmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,tmov,    &
     &kstrt,nvp)
!
! second part of particle reorder on x and y cell with mx, my tiles:
! updates ppart, kpic
      call mporder2b(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr, &
     &tsort,irc)
      if (irc /= 0) then
         write (*,*) kstrt,'mporder2b error: nppmx, irc=',nppmx,irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpfft2r(f,g,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
     &indy,kstrt,nvp,kyp,ny,mter,ierr)
! performs 2d real/complex FFT for scalar data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kyp, ny
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:), intent(inout) :: f
      complex, dimension(:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! inverse fourier transform: from real to complex
      if (isign < 0) then
! moves scalar grids from non-uniform to uniform partition
         call mpfmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,ierr&
     &)
! wrapper function for scalar 2d real/complex FFT
         call mpfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! forward fourier transform: from complex to real
      else if (isign > 0) then
! wrapper function for scalar 2d real/complex FFT
         call mpfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! moves scalar grids from uniform to non-uniform partition
         call mpfmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,ierr&
     &)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpfft2rn(f,g,noff,nyp,isign,mixup,sct,tfft,tfmov,indx,&
     &indy,kstrt,nvp,kyp,ny,mter,ierr)
! performs 2d real/complex FFT for n component vector data
! data in real space has a non-uniform partition,
! data in fourier space has a uniform partition
      implicit none
      integer, intent(in) :: isign, indx, indy, kstrt, nvp, kyp, ny
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: mter, ierr
      real, intent(inout) :: tfmov
      real, dimension(:,:,:), intent(inout) :: f
      complex, dimension(:,:,:), intent(inout) :: g
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, dimension(2), intent(inout) :: tfft
! inverse fourier transform
      if (isign < 0) then
! moves vector grids from non-uniform to uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
! wrapper function for n component vector 2d real/complex FFT
         call mpfft2rn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! forward fourier transform
      else if (isign > 0) then
! wrapper function for n component vector 2d real/complex FFT
         call mpfft2rn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,nvp,kyp)
! moves vector grids from uniform to non-uniform partition
         call mpfnmove2(f,noff,nyp,isign,tfmov,kyp,ny,kstrt,nvp,mter,   &
     &ierr)
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpcguard2(f,nyp,tguard,nx,kstrt,nvp)
! copy scalar guard cells to local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! copies data to guard cells in non-uniform partitions
      call mpcguard2(f,nyp,tguard,kstrt,nvp)
! replicates local periodic scalar field
      call mpdguard2x(f,nyp,tguard,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpncguard2(f,nyp,tguard,nx,kstrt,nvp)
! copy vector guard cells to local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! copies data to guard cells in non-uniform partitions
      call mpncguard2(f,nyp,tguard,kstrt,nvp)
! replicates local periodic vector field
      call mpcguard2x(f,nyp,tguard,nx)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpaguard2(f,nyp,tguard,nx,kstrt,nvp)
! add scalar guard cells from local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:), intent(inout) :: f
! accumulates local periodic scalar field
      call mpaguard2x(f,nyp,tguard,nx)
! adds scalar data from guard cells in non-uniform partitions
      call mpnaguard2(f,nyp,tguard,nx,kstrt,nvp)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmpnacguard2(f,nyp,tguard,nx,kstrt,nvp)
! add vector guard cells from local and remote partitions
      implicit none
      integer, intent(in) :: nyp, nx, kstrt, nvp
      real, intent(inout) :: tguard
      real, dimension(:,:,:), intent(inout) :: f
! accumulates local periodic vector field
      call mpacguard2x(f,nyp,tguard,nx)
! adds vector data from guard cells in non-uniform partitions
      call mpnacguard2(f,nyp,tguard,nx,kstrt,nvp)
      end subroutine
!
      end module
