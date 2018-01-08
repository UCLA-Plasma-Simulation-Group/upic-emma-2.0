!-----------------------------------------------------------------------
! Fortran library for spectral field solvers
! written by Michael Touati, UCLA, 2017
! 2D MPI/OpenMP PIC Code:

!-----------------------------------------------------------------------
	  subroutine MPPMAXWELB2SPECTRALPML(exy,bxy,cu,                     &
     &ffc,affp,ci,dt,wm,nx,ny,kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
! this subroutine solves 2d maxwell's equation in fourier space for
! transverse magnetic field for a half time step
! input: all, output: wm, bxy
! approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the magnetic field is first updated half a step using the equations:
! bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
! by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
! bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
! where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! j,k = fourier mode numbers, except for
! bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
! bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
! bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
! bxy(i,k,j) = i-th component of complex magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprical of velocity of light
! dt = time interval between successive calculations
! magnetic field energy is also calculated, using
! wm = (nx*ny)**2 * sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci, dt, ax, ay, wm
      complex exy, bxy, cu, ffc
      dimension exy(4,nyv,kxp), bxy(4,nyv,kxp), cu(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
      real kx, ky
      dimension kx(nx+2), ky(ny+2)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real cdth, anorm, dkx, dky
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7
      double precision wp, sum2
      complex bxy_1_1, bxy_j1_1, bxy_1_k1, bxy_j1_k1
      dimension bxy_1_1(4), bxy_j1_1(4), bxy_1_k1(4), bxy_j1_k1(4)
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      cdth = 0.5*dt/ci
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      sum2 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1,zt2,zt3,zt4,zt5,zt6,
!!$OMP& zt7,wp)
!!$OMP& REDUCTION(+:sum2)
      do 20 j = 1, kxps
      dkx = kx(j+joff+1)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = ky(k)
! update magnetic field half time step, ky > 0
! IDFT[Bx](kpx(j),ky(k)) = bxy(1,k,j) (kpx(j) = kx(j+joff+1))
         zt1 = exy(3,k,j) + exy(4,k,j)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,k,j)),real(exy(2,k,j)))
         zt3 = cmplx(-aimag(exy(1,k,j)),real(exy(1,k,j)))
         zt4 = bxy(1,k,j) - cdth*(dky*zt1)
         zt5 = bxy(2,k,j) + cdth*(dkx*zt1)
         zt6 = bxy(3,k,j) - cdth*(dkx*zt2)
         zt7 = bxy(4,k,j) + cdth*(dky*zt3)
         bxy(1,k,j) = zt4
         bxy(2,k,j) = zt5
         bxy(3,k,j) = zt6
         bxy(4,k,j) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  +(zt6+zt7)*conjg(zt6+zt7))
! update magnetic field half time step, ky < 0
! IDFT[Bx](kpx(j),-ky(k)) = bxy(1,k1,j) (kpx(j) = kx(j+joff+1))
         zt1 = exy(3,k1,j) + exy(4,k1,j)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,k1,j)),real(exy(2,k1,j)))
         zt3 = cmplx(-aimag(exy(1,k1,j)),real(exy(1,k1,j)))
         zt4 = bxy(1,k1,j) + cdth*(dky*zt1)
         zt5 = bxy(2,k1,j) + cdth*(dkx*zt1)
         zt6 = bxy(3,k1,j) - cdth*(dkx*zt2)
         zt7 = bxy(4,k1,j) - cdth*(dky*zt3)
         bxy(1,k1,j) = zt4
         bxy(2,k1,j) = zt5
         bxy(3,k1,j) = zt6
         bxy(4,k1,j) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + (zt6+zt7)*conjg(zt6+zt7))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
! ky = 0
! IDFT{Bx}(kpx(j),-ky(1)=0) = bxy(1,1,j) (kpx(j)=kx(j+joff+1))
		 zt1 = exy(3,1,j) + exy(4,1,j)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,1,j)),real(exy(2,1,j)))
         zt4 = zero
         zt5 = bxy(2,1,j) + cdth*(dkx*zt1)
         zt6 = bxy(3,1,j) - cdth*(dkx*zt2)
         zt7 = bxy(4,1,j)
         bxy(1,1,j) = zt4
         bxy(2,1,j) = zt5
         bxy(3,1,j) = zt6
         bxy(4,1,j) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  +(zt6+zt7)*conjg(zt6+zt7))
! ky = ny/2
! IDFT[Bx](kpx(j),-ky(nyh+1)) = bxy(1,k1,j)  (kpx(j) = kx(j+joff+1))
         bxy(1,k1,j) = zero
         bxy(2,k1,j) = zero
         bxy(3,k1,j) = zero
         bxy(4,k1,j) = zero
      endif
      sum2 = sum2 + wp
   20 continue
!!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         dkx = kx(nxh+1)
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = ky(k)
! update magnetic field half time step
! kx = 0
! IDFT[Bx](kx(1)=0, ky(k)) = bxy(1,k,1)
! IDFT[Bx](kx(1)=0,-ky(k)) = conjg(bxy(1,k,1))
		 zt1 = exy(3,k,1) + exy(4,k,1)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt3 = cmplx(-aimag(exy(1,k,1)),real(exy(1,k,1)))
         zt4 = bxy(1,k,1) - cdth*(dky*zt1)
         zt5 = zero
         zt6 = bxy(3,k,1)
         zt7 = bxy(4,k,1) + cdth*(dky*zt3)
         bxy(1,k,1) = zt4
         bxy(2,k,1) = zt5
         bxy(3,k,1) = zt6
         bxy(4,k,1) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  +(zt6+zt7)*conjg(zt6+zt7))
! kx = nx / 2
! IDFT{Bx}(kx(nxh+1),-ky(k))} = bxy(1,k1,1)
! IDFT{Bx}(kx(nxh+1), ky(k))  = conjg(bxy(1,k1,1))
         bxy(1,k1,1) = zero
         bxy(2,k1,1) = zero
         bxy(3,k1,1) = zero
         bxy(4,k1,1) = zero
   30    continue
         k1 = nyh + 1
         dky = ky(k1)
! kx = ky = 0
! IDFT[Bx](kx(1)=0  ,-ky(1)=0) = cmplx( real(bxy(1,1,1)),0.0)
         bxy_1_1(1) = zero
         bxy_1_1(2) = zero
         bxy_1_1(3) = zero
         bxy_1_1(4) = zero
! kx = nx/2 and ky = 0
! IDFT[Bx](kx(nxh+1),-ky(1)=0)   = cmplx(aimag(bxy(1,1,1)),0.0)
         bxy_j1_1(1) = zero
         bxy_j1_1(2) = zero
         bxy_j1_1(3) = zero
         bxy_j1_1(4) = zero
! packed :
         bxy(1,1,1)  = cmplx(real(bxy_1_1(1)),real(bxy_j1_1(1)))
     	 bxy(2,1,1)  = cmplx(real(bxy_1_1(2)),real(bxy_j1_1(2)))
         bxy(3,1,1)  = cmplx(real(bxy_1_1(3)),real(bxy_j1_1(3)))
         bxy(4,1,1)  = cmplx(real(bxy_1_1(4)),real(bxy_j1_1(4)))
! kx = 0 and ky = ny/2
! IDTF[Bx](kx(1)=0  ,-ky(nyh+1)) = cmplx( real(bxy(1,k1,1)),0.0)
         bxy_1_k1(1) = zero
         bxy_1_k1(2) = zero
         bxy_1_k1(3) = zero
         bxy_1_k1(4) = zero
! kx = nx/2 and ky = ny/2
! IDTF[Bx](kx(nxh+1),-ky(nyh+1)) = cmplx(aimag(bxy(1,k1,1)),0.0)
         bxy_j1_k1(1) = zero
         bxy_j1_k1(2) = zero
         bxy_j1_k1(3) = zero
         bxy_j1_k1(4) = zero
! packed :
         bxy(1,k1,1) = cmplx(real(bxy_1_k1(1)),real(bxy_j1_k1(1)))
         bxy(2,k1,1) = cmplx(real(bxy_1_k1(2)),real(bxy_j1_k1(2)))
         bxy(3,k1,1) = cmplx(real(bxy_1_k1(3)),real(bxy_j1_k1(3)))
         bxy(4,k1,1) = cmplx(real(bxy_1_k1(4)),real(bxy_j1_k1(4)))
      endif
      sum2 = sum2 + wp
   40 continue
      wm = real(nx)*real(ny)*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPMAXWELE2SPECTRALPML(exy,bxy,cu,ffc,                 &
     &affp,ci,dt,wf,nx,ny,kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
! this subroutine solves 2d maxwell's equation in fourier space for
! transverse electric for a whole time step
! input: all, output: wf, exy, bxy
! approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the electric field is updated a whole step using the equations:
! ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
!                       - affp*dt*cux(kx,ky)*s(kx,ky)
! ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
!                       - affp*dt*cuy(kx,ky)*s(kx,ky)
! ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
!                       - affp*dt*cuz(kx,ky)*s(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! j,k = fourier mode numbers, except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
! ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
! and similarly for bx, by, bz.
! cu(i,k,j) = i-th component of complex current density and
! exy(i,k,j) = i-th component of complex electric field,
! bxy(i,k,j) = i-th component of complex magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! aimag(ffc(k,j)) = finite-size particle shape factor s
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprical of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = (nx*ny)**2 * sum((1/affp)*|exyz(kx,ky,kz)|**2)
! magnetic field energy is also calculated, using
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci, dt, ax, ay, wf
      complex exy, bxy, cu, ffc
      dimension exy(4,nyv,kxp), bxy(4,nyv,kxp), cu(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
      real kx, ky
      dimension kx(nx+2), ky(ny+2)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real cdt, adt, anorm, dkx, dky, afdt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9, zt10
      double precision ws, sum1
      complex exy_1_1, exy_j1_1, exy_1_k1, exy_j1_k1
      dimension exy_1_1(4), exy_j1_1(4), exy_1_k1(4), exy_j1_k1(4)
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      cdt = dt/ci
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1,zt2,zt3,zt4,zt5,zt6,
!$OMP& zt7,zt8,zt9,zt10,ws)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxps
      dkx = kx(j+joff+1)
      ws = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = ky(k)
         afdt = adt*aimag(ffc(k,j))
! ky > 0
! IDFT[Ex](kpx(j),ky(k)) = exy(1,k,j) (kpx(j) = kx(j+joff+1))
         zt4  = bxy(1,k,j)
         zt5  = bxy(2,k,j)
         zt6  = bxy(3,k,j) + bxy(4,k,j)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k,j) + cdt*(dky*zt1) - afdt    *cu(1,k,j)
         zt8  = exy(2,k,j) - cdt*(dkx*zt1) - afdt    *cu(2,k,j)
         zt9  = exy(3,k,j) + cdt*(dkx*zt2) - afdt*0.5*cu(3,k,j)
         zt10 = exy(4,k,j) - cdt*(dky*zt3) - afdt*0.5*cu(3,k,j)
         exy(1,k,j) = zt7
         exy(2,k,j) = zt8
         exy(3,k,j) = zt9
         exy(4,k,j) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
! ky < 0
! IDFT[Ex](kpx(j),-ky(k)) = exy(1,k1,j) (kpx(j) = kx(j+joff+1))
         zt4  = bxy(1,k1,j)
         zt5  = bxy(2,k1,j)
         zt6  = bxy(3,k1,j) + bxy(4,k1,j)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k1,j) - cdt*(dky*zt1) - afdt    *cu(1,k1,j)
         zt8  = exy(2,k1,j) - cdt*(dkx*zt1) - afdt    *cu(2,k1,j)
         zt9  = exy(3,k1,j) + cdt*(dkx*zt2) - afdt*0.5*cu(3,k1,j)
         zt10 = exy(4,k1,j) + cdt*(dky*zt3) - afdt*0.5*cu(3,k1,j)
         exy(1,k1,j) = zt7
         exy(2,k1,j) = zt8
         exy(3,k1,j) = zt9
         exy(4,k1,j) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
   10    continue
! mode numbers ky = 0, ny/2
! ky = 0
! IDFT{Ex}(kpx(j),-ky(1)=0) = exy(1,1,j) (kpx(j)=kx(j+joff+1))
         afdt = adt*aimag(ffc(1,j))
         zt5  = bxy(2,1,j)
         zt6  = bxy(3,1,j) + bxy(4,1,j)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt7 = zero
         zt8  = exy(2,1,j) - cdt*(dkx*zt1) - afdt    *cu(2,1,j)
         zt9  = exy(3,1,j) + cdt*(dkx*zt2) - afdt*0.5*cu(3,1,j)
         zt10 = exy(4,1,j)                 - afdt*0.5*cu(3,1,j)
         exy(1,1,j) = zt7
         exy(2,1,j) = zt8
         exy(3,1,j) = zt9
         exy(4,1,j) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
! ky = ny/2
! IDFT[Ex](kpx(j),-ky(nyh+1)) = exy(1,k1,j)  (kpx(j) = kx(j+joff+1))
         k1 = nyh + 1
         exy(1,k1,j) = zero
         exy(2,k1,j) = zero
         exy(3,k1,j) = zero
         exy(4,k1,j) = zero
      endif
      sum1 = sum1 + ws
   20 continue
!$OMP END PARALLEL DO
      ws = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         dkx = kx(nxh+1)
         do 30 k = 2, nyh
! kx = 0
! IDFT[Ex](kx(1)=0, ky(k)) = exy(1,k,1)
! IDFT[Ex](kx(1)=0,-ky(k)) = conjg(exy(1,k,1))
         dky = ky(k)
         afdt = adt*aimag(ffc(k,1))
         zt4  = bxy(1,k,1)
         zt6  = bxy(3,k,1) + bxy(4,k,1)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k,1) + cdt*(dky*zt1) - afdt    *cu(1,k,1)
         zt8  = zero
         zt9  = exy(3,k,1)                 - afdt*0.5*cu(3,k,1)
         zt10 = exy(4,k,1) - cdt*(dky*zt3) - afdt*0.5*cu(3,k,1)
         exy(1,k,1) = zt7
         exy(2,k,1) = zt8
         exy(3,k,1) = zt9
         exy(4,k,1) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
! kx = nx / 2
! IDFT[Ex](kx(nxh+1),-ky(k)) = exy(1,k1,1)
! IDFT[Ex](kx(nxh+1), ky(k)) = conjg(exy(1,k1,1))
 	     k1 = ny2 - k
         exy(1,k1,1) = zero
         exy(2,k1,1) = zero
         exy(3,k1,1) = zero
         exy(4,k1,1) = zero
   30    continue
         k1 = nyh + 1
         dky = ky(k1)
! kx = ky = 0
! IDTF[Ex](kx(1)=0  ,-ky(1)=0) = cmplx(real(exy(1,1,1)),0.0)
         exy_1_1(1)  = zero
         exy_1_1(2)  = zero
         exy_1_1(3)  = zero
         exy_1_1(4)  = zero
! kx = nx/2 and ky = 0
! IDFT[Ex](kx(nxh+1),-ky(1)=0)   = cmplx(aimag(exy(1,1,1)),0.0)
         exy_j1_1(1)  = zero
         exy_j1_1(2)  = zero
         exy_j1_1(3)  = zero
         exy_j1_1(4)  = zero
! packed :
		 exy(1,1,1)  = cmplx(real(exy_1_1(1)),real(exy_j1_1(1)))
         exy(2,1,1)  = cmplx(real(exy_1_1(2)),real(exy_j1_1(2)))
         exy(3,1,1)  = cmplx(real(exy_1_1(3)),real(exy_j1_1(3)))
         exy(4,1,1)  = cmplx(real(exy_1_1(4)),real(exy_j1_1(4)))
! kx = 0 and ky = ny/2
! IDTF[Ex](kx(1)=0  ,-ky(nyh+1)) = cmplx( real(exy(1,k1,1)),0.0)
         exy_1_k1(1)  = zero
         exy_1_k1(2)  = zero
         exy_1_k1(3)  = zero
         exy_1_k1(4)  = zero
! kx = nx/2 and ky = ny/2
! IDTF[Ex](kx(nxh+1),-ky(nyh+1)) = cmplx(aimag(exy(1,k1,1)),0.0)
         exy_j1_k1(1)  = zero
         exy_j1_k1(2)  = zero
         exy_j1_k1(3)  = zero
         exy_j1_k1(4)  = zero
! packed :
         exy(1,k1,1) = cmplx(real(exy_1_k1(1)),real(exy_j1_k1(1)))
         exy(2,k1,1) = cmplx(real(exy_1_k1(2)),real(exy_j1_k1(2)))
         exy(3,k1,1) = cmplx(real(exy_1_k1(3)),real(exy_j1_k1(3)))
         exy(4,k1,1) = cmplx(real(exy_1_k1(4)),real(exy_j1_k1(4)))
      endif
      sum1 = sum1 + ws
   40 continue
      wf = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------
	  subroutine MPPMAXWELB2YEEPML(exy,bxy,exy_corr,bxy_corr,cu,ffc,    &
     &affp,ci,dt,wm,nx,ny,kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
! this subroutine solves 2d maxwell's equation in fourier space for
! transverse magnetic field for a half time step
! input: all, output: wm, bxy
! approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the magnetic field is first updated half a step using the equations:
! bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
! by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
! bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
! where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! j,k = fourier mode numbers, except for
! bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
! bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0,
! bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
! bxy(i,k,j) = i-th component of complex magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprical of velocity of light
! dt = time interval between successive calculations
! magnetic field energy is also calculated, using
! wm = (nx*ny)**2 * sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci, dt, ax, ay, wm
      complex exy, bxy, cu, ffc
      dimension exy(4,nyv,kxp), bxy(4,nyv,kxp), cu(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
      complex exy_corr, bxy_corr
      dimension exy_corr(4,4), bxy_corr(4,4)
      real kx, ky
      dimension kx(nx+2), ky(ny+2)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real cdth, anorm, dkx, dky
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7
      double precision wp, sum2
      complex bxy_1_1, bxy_j1_1, bxy_1_k1, bxy_j1_k1
      dimension bxy_1_1(4), bxy_j1_1(4), bxy_1_k1(4), bxy_j1_k1(4)
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      cdth = 0.5*dt/ci
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      sum2 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1,zt2,zt3,zt4,zt5,zt6,
!!$OMP& zt7,wp)
!!$OMP& REDUCTION(+:sum2)
      do 20 j = 1, kxps
      dkx = kx(j+joff+1)
      wp = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = ky(k)
! update magnetic field half time step, ky > 0
! IDFT[Bx](kpx(j),ky(k)) = bxy(1,k,j) (kpx(j) = kx(j+joff+1))
         zt1 = exy(3,k,j) + exy(4,k,j)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,k,j)),real(exy(2,k,j)))
         zt3 = cmplx(-aimag(exy(1,k,j)),real(exy(1,k,j)))
         zt4 = bxy(1,k,j) - cdth*(dky*zt1)
         zt5 = bxy(2,k,j) + cdth*(dkx*zt1)
         zt6 = bxy(3,k,j) - cdth*(dkx*zt2)
         zt7 = bxy(4,k,j) + cdth*(dky*zt3)
         bxy(1,k,j) = zt4
         bxy(2,k,j) = zt5
         bxy(3,k,j) = zt6
         bxy(4,k,j) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  +(zt6+zt7)*conjg(zt6+zt7))
! update magnetic field half time step, ky < 0
! IDFT[Bx](kpx(j),-ky(k)) = bxy(1,k1,j) (kpx(j) = kx(j+joff+1))
         zt1 = exy(3,k1,j) + exy(4,k1,j)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,k1,j)),real(exy(2,k1,j)))
         zt3 = cmplx(-aimag(exy(1,k1,j)),real(exy(1,k1,j)))
         zt4 = bxy(1,k1,j) + cdth*(dky*zt1)
         zt5 = bxy(2,k1,j) + cdth*(dkx*zt1)
         zt6 = bxy(3,k1,j) - cdth*(dkx*zt2)
         zt7 = bxy(4,k1,j) - cdth*(dky*zt3)
         bxy(1,k1,j) = zt4
         bxy(2,k1,j) = zt5
         bxy(3,k1,j) = zt6
         bxy(4,k1,j) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + (zt6+zt7)*conjg(zt6+zt7))
   10    continue
! mode numbers ky = 0, ny/2
         k1 = nyh + 1
! ky = 0
! IDFT{Bx}(kpx(j),-ky(1)=0) = bxy(1,1,j) (kpx(j)=kx(j+joff+1))
		 zt1 = exy(3,1,j) + exy(4,1,j)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,1,j)),real(exy(2,1,j)))
         zt4 = zero
         zt5 = bxy(2,1,j) + cdth*(dkx*zt1)
         zt6 = bxy(3,1,j) - cdth*(dkx*zt2)
         zt7 = bxy(4,1,j)
         bxy(1,1,j) = zt4
         bxy(2,1,j) = zt5
         bxy(3,1,j) = zt6
         bxy(4,1,j) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  +(zt6+zt7)*conjg(zt6+zt7))
! ky = ny/2
! IDFT[Bx](kpx(j),-ky(nyh+1)) = bxy(1,k1,j)  (kpx(j) = kx(j+joff+1))
         zt1 = exy(3,k1,j) + exy(4,k1,j)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,k1,j)),real(exy(2,k1,j)))
         zt3 = cmplx(-aimag(exy(1,k1,j)),real(exy(1,k1,j)))
         zt4 = bxy(1,k1,j) + cdth*(dky*zt1)
         zt5 = bxy(2,k1,j) + cdth*(dkx*zt1)
         zt6 = bxy(3,k1,j) - cdth*(dkx*zt2)
         zt7 = bxy(4,k1,j) - cdth*(dky*zt3)
         bxy(1,k1,j) = zt4
         bxy(2,k1,j) = zt5
         bxy(3,k1,j) = zt6
         bxy(4,k1,j) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  +(zt6+zt7)*conjg(zt6+zt7))
      endif
      sum2 = sum2 + wp
   20 continue
!!$OMP END PARALLEL DO
      wp = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         dkx = kx(nxh+1)
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = ky(k)
! update magnetic field half time step
! kx = 0
! IDFT[Bx](kx(1)=0, ky(k)) = bxy(1,k,1)
! IDFT[Bx](kx(1)=0,-ky(k)) = conjg(bxy(1,k,1))
		 zt1 = exy(3,k,1) + exy(4,k,1)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt3 = cmplx(-aimag(exy(1,k,1)),real(exy(1,k,1)))
         zt4 = bxy(1,k,1) - cdth*(dky*zt1)
         zt5 = zero
         zt6 = bxy(3,k,1)
         zt7 = bxy(4,k,1) + cdth*(dky*zt3)
         bxy(1,k,1) = zt4
         bxy(2,k,1) = zt5
         bxy(3,k,1) = zt6
         bxy(4,k,1) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  +(zt6+zt7)*conjg(zt6+zt7))
! kx = nx / 2
! IDFT{Bx}(kx(nxh+1),-ky(k))} = bxy(1,k1,1)
! IDFT{Bx}(kx(nxh+1), ky(k))  = conjg(bxy(1,k1,1))
         k1 = ny2 - k
         !--------------------------
         zt1 = exy(3,k1,1) + exy(4,k1,1)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy(2,k1,1)),real(exy(2,k1,1)))
         zt3 = cmplx(-aimag(exy(1,k1,1)),real(exy(1,k1,1)))
         zt4 = bxy(1,k1,1) + cdth*(dky*zt1)
         zt5 = bxy(2,k1,1) + cdth*(dkx*zt1)
         zt6 = bxy(3,k1,1) - cdth*(dkx*zt2)
         zt7 = bxy(4,k1,1) - cdth*(dky*zt3)
         bxy(1,k1,1) = zt4
         bxy(2,k1,1) = zt5
         bxy(3,k1,1) = zt6
         bxy(4,k1,1) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + (zt6+zt7)*conjg(zt6+zt7))
   30    continue
         k1 = nyh + 1
         dky = ky(k1)
! kx = ky = 0
! IDFT[Bx](kx(1)=0  ,-ky(1)=0) = cmplx( real(bxy(1,1,1)),0.0)
         bxy_corr(1,1) = bxy_corr(1,1)
         bxy_corr(1,2) = bxy_corr(1,2)
         bxy_corr(1,3) = bxy_corr(1,3)
         bxy_corr(1,4) = bxy_corr(1,4)
         zt4 = zero
         zt5 = zero
         zt6 = bxy_corr(1,3)
         zt7 = bxy_corr(1,4)
         bxy_1_1(1) = zt4
         bxy_1_1(2) = zt5
         bxy_1_1(3) = zt6
         bxy_1_1(4) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + (zt6+zt7)*conjg(zt6+zt7))
! kx = nx/2 and ky = 0
! IDFT[Bx](kx(nxh+1),-ky(1)=0)   = cmplx(aimag(bxy(1,1,1)),0.0)
         zt1 = exy_corr(2,3) + exy_corr(2,4)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy_corr(2,2)),real(exy_corr(2,2)))
         bxy_corr(2,1) = zero
         bxy_corr(2,2) = bxy_corr(2,2) + cdth*(dkx*zt1)
         bxy_corr(2,3) = bxy_corr(2,3) - cdth*(dkx*zt2)
         bxy_corr(2,4) = bxy_corr(2,4)
         zt4 = cmplx(real(bxy_corr(2,1)),0.0)
         zt5 = cmplx(real(bxy_corr(2,2)),0.0)
         zt6 = cmplx(real(bxy_corr(2,3)),0.0)
         zt7 = cmplx(real(bxy_corr(2,4)),0.0)
         bxy_j1_1(1) = zt4
         bxy_j1_1(2) = zt5
         bxy_j1_1(3) = zt6
         bxy_j1_1(4) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + (zt6+zt7)*conjg(zt6+zt7))
! packed :
         bxy(1,1,1)  = cmplx(real(bxy_1_1(1)),real(bxy_j1_1(1)))
     	 bxy(2,1,1)  = cmplx(real(bxy_1_1(2)),real(bxy_j1_1(2)))
         bxy(3,1,1)  = cmplx(real(bxy_1_1(3)),real(bxy_j1_1(3)))
         bxy(4,1,1)  = cmplx(real(bxy_1_1(4)),real(bxy_j1_1(4)))
! kx = 0 and ky = ny/2
! IDTF[Bx](kx(1)=0  ,-ky(nyh+1)) = cmplx( real(bxy(1,k1,1)),0.0)
         zt1 = exy_corr(3,3) + exy_corr(3,4)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt3 = cmplx(-aimag(exy_corr(3,1)),real(exy_corr(3,1)))
         bxy_corr(3,1) = bxy_corr(3,1) + cdth*(dky*zt1)
         bxy_corr(3,2) = zero
         bxy_corr(3,3) = bxy_corr(3,3)
         bxy_corr(3,4) = bxy_corr(3,4) - cdth*(dky*zt3)
         zt4 = cmplx(real(bxy_corr(3,1)),0.0)
         zt5 = cmplx(real(bxy_corr(3,2)),0.0)
         zt6 = cmplx(real(bxy_corr(3,3)),0.0)
         zt7 = cmplx(real(bxy_corr(3,4)),0.0)
         bxy_1_k1(1) = zt4
         bxy_1_k1(2) = zt5
         bxy_1_k1(3) = zt6
         bxy_1_k1(4) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + (zt6+zt7)*conjg(zt6+zt7))
! kx = nx/2 and ky = ny/2
! IDTF[Bx](kx(nxh+1),-ky(nyh+1)) = cmplx(aimag(bxy(1,k1,1)),0.0)
         !------------------------------------------
         zt1 = exy_corr(4,3) + exy_corr(4,4)
         zt1 = cmplx(-aimag(zt1),real(zt1))
         zt2 = cmplx(-aimag(exy_corr(4,2)),real(exy_corr(4,2)))
         zt3 = cmplx(-aimag(exy_corr(4,1)),real(exy_corr(4,1)))
         bxy_corr(4,1) = bxy_corr(4,1) + cdth*(dky*zt1)
         bxy_corr(4,2) = bxy_corr(4,2) + cdth*(dkx*zt1)
         bxy_corr(4,3) = bxy_corr(4,3) - cdth*(dkx*zt2)
         bxy_corr(4,4) = bxy_corr(4,4) - cdth*(dky*zt3)
         zt4 = cmplx(real(bxy_corr(4,1)),0.0)
         zt5 = cmplx(real(bxy_corr(4,2)),0.0)
         zt6 = cmplx(real(bxy_corr(4,3)),0.0)
         zt7 = cmplx(real(bxy_corr(4,4)),0.0)
         bxy_j1_k1(1) = zt4
         bxy_j1_k1(2) = zt5
         bxy_j1_k1(3) = zt6
         bxy_j1_k1(4) = zt7
         wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5)               &
     &                  + (zt6+zt7)*conjg(zt6+zt7))
! packed :
         bxy(1,k1,1) = cmplx(real(bxy_1_k1(1)),real(bxy_j1_k1(1)))
         bxy(2,k1,1) = cmplx(real(bxy_1_k1(2)),real(bxy_j1_k1(2)))
         bxy(3,k1,1) = cmplx(real(bxy_1_k1(3)),real(bxy_j1_k1(3)))
         bxy(4,k1,1) = cmplx(real(bxy_1_k1(4)),real(bxy_j1_k1(4)))
      endif
      sum2 = sum2 + wp
   40 continue
      wm = real(nx)*real(ny)*sum2
      return
      end
!-----------------------------------------------------------------------
      subroutine MPPMAXWELE2YEEPML(exy,bxy,exy_corr,bxy_corr,cu,ffc,    &
     &affp,ci,dt,wf,nx,ny,kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
! this subroutine solves 2d maxwell's equation in fourier space for
! transverse electric for a whole time step
! input: all, output: wf, exy, bxy
! approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
! where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
! the electric field is updated a whole step using the equations:
! ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
!                       - affp*dt*cux(kx,ky)*s(kx,ky)
! ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
!                       - affp*dt*cuy(kx,ky)*s(kx,ky)
! ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
!                       - affp*dt*cuz(kx,ky)*s(kx,ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! j,k = fourier mode numbers, except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0,
! ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
! and similarly for bx, by, bz.
! cu(i,k,j) = i-th component of complex current density and
! exy(i,k,j) = i-th component of complex electric field,
! bxy(i,k,j) = i-th component of complex magnetic field,
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! aimag(ffc(k,j)) = finite-size particle shape factor s
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! for fourier mode (jj-1,k-1), where jj = j + kxp*(kstrt - 1)
! affp = normalization constant = nx*ny/np, where np=number of particles
! ci = reciprical of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = (nx*ny)**2 * sum((1/affp)*|exyz(kx,ky,kz)|**2)
! magnetic field energy is also calculated, using
! nx/ny = system length in x/y direction
! kxp = number of data values per block
! kstrt = starting data block number
! nyv = first dimension of field arrays, must be >= ny
! nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, kstrt, nyv, kxp, nyhd
      real affp, ci, dt, ax, ay, wf
      complex exy, bxy, cu, ffc
      dimension exy(4,nyv,kxp), bxy(4,nyv,kxp), cu(3,nyv,kxp)
      dimension ffc(nyhd,kxp)
      complex exy_corr, bxy_corr
      dimension exy_corr(4,4), bxy_corr(4,4)
      real kx, ky
      dimension kx(nx+2), ky(ny+2)
! local data
      integer nxh, nyh, ny2, ks, joff, kxps, j, k, k1
      real cdt, adt, anorm, dkx, dky, afdt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9, zt10
      double precision ws, sum1
      complex exy_1_1, exy_j1_1, exy_1_k1, exy_j1_k1
      dimension exy_1_1(4), exy_j1_1(4), exy_1_k1(4), exy_j1_k1(4)
      complex cu_1_1, cu_j1_1, cu_1_k1, cu_j1_k1
      dimension cu_1_1(4), cu_j1_1(4), cu_1_k1(4), cu_j1_k1(4)
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      ks = kstrt - 1
      joff = kxp*ks
      kxps = min(kxp,max(0,nxh-joff))
      joff = joff - 1
      cdt = dt/ci
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      sum1 = 0.0d0
      if (kstrt.gt.nxh) go to 40
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!$OMP PARALLEL DO PRIVATE(j,k,k1,dkx,dky,zt1,zt2,zt3,zt4,zt5,zt6,
!$OMP& zt7,zt8,zt9,zt10,ws)
!$OMP& REDUCTION(+:sum1)
      do 20 j = 1, kxps
      dkx = kx(j+joff+1)
      ws = 0.0d0
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = ky(k)
         afdt = adt*aimag(ffc(k,j))
! ky > 0
! IDFT[Ex](kpx(j),ky(k)) = exy(1,k,j) (kpx(j) = kx(j+joff+1))
         zt4  = bxy(1,k,j)
         zt5  = bxy(2,k,j)
         zt6  = bxy(3,k,j) + bxy(4,k,j)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k,j) + cdt*(dky*zt1) - afdt    *cu(1,k,j)
         zt8  = exy(2,k,j) - cdt*(dkx*zt1) - afdt    *cu(2,k,j)
         zt9  = exy(3,k,j) + cdt*(dkx*zt2) - afdt*0.5*cu(3,k,j)
         zt10 = exy(4,k,j) - cdt*(dky*zt3) - afdt*0.5*cu(3,k,j)
         exy(1,k,j) = zt7
         exy(2,k,j) = zt8
         exy(3,k,j) = zt9
         exy(4,k,j) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
! ky < 0
! IDFT[Ex](kpx(j),-ky(k)) = exy(1,k1,j) (kpx(j) = kx(j+joff+1))
         zt4  = bxy(1,k1,j)
         zt5  = bxy(2,k1,j)
         zt6  = bxy(3,k1,j) + bxy(4,k1,j)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k1,j) - cdt*(dky*zt1) - afdt    *cu(1,k1,j)
         zt8  = exy(2,k1,j) - cdt*(dkx*zt1) - afdt    *cu(2,k1,j)
         zt9  = exy(3,k1,j) + cdt*(dkx*zt2) - afdt*0.5*cu(3,k1,j)
         zt10 = exy(4,k1,j) + cdt*(dky*zt3) - afdt*0.5*cu(3,k1,j)
         exy(1,k1,j) = zt7
         exy(2,k1,j) = zt8
         exy(3,k1,j) = zt9
         exy(4,k1,j) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
   10    continue
! mode numbers ky = 0, ny/2
! ky = 0
! IDFT{Ex}(kpx(j),-ky(1)=0) = exy(1,1,j) (kpx(j)=kx(j+joff+1))
         afdt = adt*aimag(ffc(1,j))
         zt5  = bxy(2,1,j)
         zt6  = bxy(3,1,j) + bxy(4,1,j)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt7 = zero
         zt8  = exy(2,1,j) - cdt*(dkx*zt1) - afdt    *cu(2,1,j)
         zt9  = exy(3,1,j) + cdt*(dkx*zt2) - afdt*0.5*cu(3,1,j)
         zt10 = exy(4,1,j)                 - afdt*0.5*cu(3,1,j)
         exy(1,1,j) = zt7
         exy(2,1,j) = zt8
         exy(3,1,j) = zt9
         exy(4,1,j) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
! ky = ny/2
! IDFT[Ex](kpx(j),-ky(nyh+1)) = exy(1,k1,j)  (kpx(j) = kx(j+joff+1))
         k1 = nyh + 1
         !-------------------
         dky = ky(nyh+1)
         afdt = adt*aimag(ffc(1,j))*exp(-0.5*(ay*dky)**2)
         zt4  = bxy(1,k1,j)
         zt5  = bxy(2,k1,j)
         zt6  = bxy(3,k1,j) + bxy(4,k1,j)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k1,j) - cdt*(dky*zt1) -afdt    *cu(1,k1,j)
         zt8  = exy(2,k1,j) - cdt*(dkx*zt1) -afdt    *cu(2,k1,j)
         zt9  = exy(3,k1,j) + cdt*(dkx*zt2) -afdt*0.5*cu(3,k1,j)
         zt10 = exy(4,k1,j) + cdt*(dky*zt3) -afdt*0.5*cu(3,k1,j)
         exy(1,k1,j) = zt7
         exy(2,k1,j) = zt8
         exy(3,k1,j) = zt9
         exy(4,k1,j) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
      endif
      sum1 = sum1 + ws
   20 continue
!$OMP END PARALLEL DO
      ws = 0.0d0
! mode numbers kx = 0, nx/2
      if (ks.eq.0) then
         dkx = kx(nxh+1)
         do 30 k = 2, nyh
! kx = 0
! IDFT[Ex](kx(1)=0, ky(k)) = exy(1,k,1)
! IDFT[Ex](kx(1)=0,-ky(k)) = conjg(exy(1,k,1))
         dky = ky(k)
         afdt = adt*aimag(ffc(k,1))
         zt4  = bxy(1,k,1)
         zt6  = bxy(3,k,1) + bxy(4,k,1)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k,1) + cdt*(dky*zt1) - afdt    *cu(1,k,1)
         zt8  = zero
         zt9  = exy(3,k,1)                 - afdt*0.5*cu(3,k,1)
         zt10 = exy(4,k,1) - cdt*(dky*zt3) - afdt*0.5*cu(3,k,1)
         exy(1,k,1) = zt7
         exy(2,k,1) = zt8
         exy(3,k,1) = zt9
         exy(4,k,1) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  +(zt9+zt10)*conjg(zt9+zt10))
! kx = nx / 2
! IDFT[Ex](kx(nxh+1),-ky(k)) = exy(1,k1,1)
! IDFT[Ex](kx(nxh+1), ky(k)) = conjg(exy(1,k1,1))
 	     k1 = ny2 - k
 	     !--------------------------
         afdt = adt*aimag(ffc(k,1))*exp(-0.5*(ax*dkx)**2)
         zt4  = bxy(1,k1,1)
         zt5  = bxy(2,k1,1)
         zt6  = bxy(3,k1,1) + bxy(4,k1,1)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         zt7  = exy(1,k1,1) - cdt*(dky*zt1) -afdt    *cu(1,k1,1)
         zt8  = exy(2,k1,1) - cdt*(dkx*zt1) -afdt    *cu(2,k1,1)
         zt9  = exy(3,k1,1) + cdt*(dkx*zt2) -afdt*0.5*cu(3,k1,1)
         zt10 = exy(4,k1,1) + cdt*(dky*zt3) -afdt*0.5*cu(3,k1,1)
         exy(1,k1,1) = zt7
         exy(2,k1,1) = zt8
         exy(3,k1,1) = zt9
         exy(4,k1,1) = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  + (zt9+zt10)*conjg(zt9+zt10))
   30    continue
         k1 = nyh + 1
         dky = ky(k1)
! kx = ky = 0
! IDTF[Ex](kx(1)=0  ,-ky(1)=0) = cmplx(real(exy(1,1,1)),0.0)
         afdt = adt
         cu_1_1(1)  = cmplx(real(cu(1,1,1)),0.0)
         cu_1_1(2)  = cmplx(real(cu(2,1,1)),0.0)
         cu_1_1(3)  = cmplx(real(cu(3,1,1)),0.0)
		 exy_corr(1,1) = exy_corr(1,1) - afdt*0.5*cu_1_1(1)
         exy_corr(1,2) = exy_corr(1,2) - afdt*0.5*cu_1_1(2)
         exy_corr(1,3) = exy_corr(1,3) - afdt*0.5*cu_1_1(3)
         exy_corr(1,4) = exy_corr(1,4) - afdt*0.5*cu_1_1(3)
         zt7  = zero
         zt8  = zero
         zt9  = exy_corr(1,3)
         zt10 = exy_corr(1,4)
         exy_1_1(1)  = zt7
         exy_1_1(2)  = zt8
         exy_1_1(3)  = zt9
         exy_1_1(4)  = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  + (zt9+zt10)*conjg(zt9+zt10))
! kx = nx/2 and ky = 0
! IDFT[Ex](kx(nxh+1),-ky(1)=0)   = cmplx(aimag(exy(1,1,1)),0.0)
 	     afdt = adt*exp(-0.5*((ax*dkx)**2.))
         cu_j1_1(1)  = cmplx(aimag(cu(1,1,1)),0.0)
         cu_j1_1(2)  = cmplx(aimag(cu(2,1,1)),0.0)
         cu_j1_1(3)  = cmplx(aimag(cu(3,1,1)),0.0)
         zt5  = bxy_corr(2,2)
         zt6  = bxy_corr(2,3) + bxy_corr(2,4)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt2  = cmplx(-aimag(zt5),real(zt5))
         exy_corr(2,1) = zero
         exy_corr(2,2) = exy_corr(2,2) - cdt*(dkx*zt1)                  &
     &                 - afdt    *cu_j1_1(2)
         exy_corr(2,3) = exy_corr(2,3) + cdt*(dkx*zt2)                  &
     &                 - afdt*0.5*cu_j1_1(3)
         exy_corr(2,4) = exy_corr(2,4)                                  &
     &                 - afdt*0.5*cu_j1_1(3)
         zt7   = cmplx(real(exy_corr(2,1)),0.0)
         zt8   = cmplx(real(exy_corr(2,2)),0.0)
         zt9   = cmplx(real(exy_corr(2,3)),0.0)
         zt10  = cmplx(real(exy_corr(2,4)),0.0)        
         exy_j1_1(1)  = zt7
         exy_j1_1(2)  = zt8
         exy_j1_1(3)  = zt9
         exy_j1_1(4)  = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  + (zt9+zt10)*conjg(zt9+zt10))
! packed :
		 exy(1,1,1)  = cmplx(real(exy_1_1(1)),real(exy_j1_1(1)))
         exy(2,1,1)  = cmplx(real(exy_1_1(2)),real(exy_j1_1(2)))
         exy(3,1,1)  = cmplx(real(exy_1_1(3)),real(exy_j1_1(3)))
         exy(4,1,1)  = cmplx(real(exy_1_1(4)),real(exy_j1_1(4)))
! kx = 0 and ky = ny/2
! IDTF[Ex](kx(1)=0  ,-ky(nyh+1)) = cmplx( real(exy(1,k1,1)),0.0)
         afdt = adt*exp(-0.5*((ay*dky)**2.))
         cu_1_k1(1)  = cmplx(real(cu(1,k1,1)),0.0)
         cu_1_k1(2)  = cmplx(real(cu(2,k1,1)),0.0)
         cu_1_k1(3)  = cmplx(real(cu(3,k1,1)),0.0)
         zt4  = bxy_corr(3,1)
         zt6  = bxy_corr(3,3) + bxy_corr(3,4)
         zt1  = cmplx(-aimag(zt6),real(zt6))
         zt3  = cmplx(-aimag(zt4),real(zt4))
         exy_corr(3,1) = exy_corr(3,1) - cdt*(dky*zt1)                  &
     &                 - afdt    *cu_1_k1(1)
         exy_corr(3,2) = zero
         exy_corr(3,3) = exy_corr(3,3)                                  &
     &                 - afdt*0.5*cu_1_k1(3)
         exy_corr(3,4) = exy_corr(3,4) + cdt*(dky*zt3)                  &
     &                 - afdt*0.5*cu_1_k1(3)
         zt7  = cmplx(real(exy_corr(3,1)),0.0)
         zt8  = cmplx(real(exy_corr(3,2)),0.0)
         zt9  = cmplx(real(exy_corr(3,3)),0.0)
         zt10 = cmplx(real(exy_corr(3,4)),0.0)
         exy_1_k1(1)  = zt7
         exy_1_k1(2)  = zt8
         exy_1_k1(3)  = zt9
         exy_1_k1(4)  = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  + (zt9+zt10)*conjg(zt9+zt10))
! kx = nx/2 and ky = ny/2
! IDTF[Ex](kx(nxh+1),-ky(nyh+1)) = cmplx(aimag(exy(1,k1,1)),0.0)
         afdt = adt*exp(-0.5*(((ax*dkx)**2)+((ay*dky)**2)))
         cu_j1_k1(1)  = cmplx(aimag(cu(1,k1,1)),0.0)
         cu_j1_k1(2)  = cmplx(aimag(cu(2,k1,1)),0.0)
         cu_j1_k1(3)  = cmplx(aimag(cu(3,k1,1)),0.0)
         zt1  = cmplx(-aimag(bxy_corr(4,3)),real(bxy_corr(4,3)))        &
     &        + cmplx(-aimag(bxy_corr(4,4)),real(bxy_corr(4,4)))
         zt2  = cmplx(-aimag(bxy_corr(4,2)),real(bxy_corr(4,2)))
         zt3  = cmplx(-aimag(bxy_corr(4,1)),real(bxy_corr(4,1)))
         exy_corr(4,1) = exy_corr(4,1) - cdt*(dky*zt1)                  &
     &                 - afdt    *cu_j1_k1(1)
         exy_corr(4,2) = exy_corr(4,2) - cdt*(dkx*zt1)                  &
     &                 - afdt    *cu_j1_k1(2)
         exy_corr(4,3) = exy_corr(4,3) + cdt*(dkx*zt2)                  &
     &                 - afdt*0.5*cu_j1_k1(3)
         exy_corr(4,4) = exy_corr(4,4) + cdt*(dky*zt3)                  &
     &                 - afdt*0.5*cu_j1_k1(3)
         zt7  = cmplx(real(exy_corr(4,1)),0.0)
         zt8  = cmplx(real(exy_corr(4,2)),0.0)
         zt9  = cmplx(real(exy_corr(4,3)),0.0)
         zt10 = cmplx(real(exy_corr(4,4)),0.0)
         exy_j1_k1(1)  = zt7
         exy_j1_k1(2)  = zt8
         exy_j1_k1(3)  = zt9
         exy_j1_k1(4)  = zt10
         ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8)               &
     &                  + (zt9+zt10)*conjg(zt9+zt10))
! packed :
         exy(1,k1,1) = cmplx(real(exy_1_k1(1)),real(exy_j1_k1(1)))
         exy(2,k1,1) = cmplx(real(exy_1_k1(2)),real(exy_j1_k1(2)))
         exy(3,k1,1) = cmplx(real(exy_1_k1(3)),real(exy_j1_k1(3)))
         exy(4,k1,1) = cmplx(real(exy_1_k1(4)),real(exy_j1_k1(4)))
      endif
      sum1 = sum1 + ws
   40 continue
      wf = real(nx)*real(ny)*sum1
      return
      end
!-----------------------------------------------------------------------