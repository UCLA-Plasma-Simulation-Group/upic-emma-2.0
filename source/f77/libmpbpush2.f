!-----------------------------------------------------------------------
! Fortran Library for pushing electromagnetic particles
! 2-1/2D MPI/OpenMP PIC Codes:
! PPPMOVIN2L sorts particles by x,y grid in tiles of mx, my and
!            copies to segmented array ppart
! PPPCOPYOUT2 copies segmented particle data ppart to the array part
! PPPCHECK2L performs a sanity check to make sure particles sorted
!            by x,y grid in tiles of mx, my, are all within bounds.
! PPGPPOST2L calculates particle charge density using linear
!            interpolation
! PPGBPPUSH23L updates magnetized particle co-ordinates and velocities
!              using leap-frog scheme in time and linear interpolation
!              in space with various particle boundary conditions
! PPGBPPUSHF23L updates magnetized particle co-ordinates and velocities
!               using leap-frog scheme in time and linear interpolation
!               in space with periodic particle boundary conditions,
!               determines list of particles which are leaving each tile
! PPGRBPPUSH23L updates relativistic magnetized particle co-ordinates
!               and momenta using leap-frog scheme in time and linear
!               interpolation in space with various particle boundary
!               conditions
! PPGRBPPUSHF23L updates relativistic magnetized particle co-ordinates
!                and momenta using leap-frog scheme in time and linear
!                interpolation in space with periodic particle boundary
!                conditions, determines list of particles which are
!                leaving each tile
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: january 18, 2016
!-----------------------------------------------------------------------
!      subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax, &
!     &mx,my,mx1,mxyp1,irc)                                                    ! M. Touati
      subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax, &
     &mx,my,mx1,mxyp1,delta,irc)                                               ! M. Touati
! this subroutine sorts particles by x,y grid in tiles of
! mx, my and copies to segmented array ppart
! linear interpolation, spatial decomposition in y direction
! input: all except ppart, kpic, output: ppart, kpic
! part/ppart = input/output particle arrays
! part(1,n) = position x of particle n in partition
! part(2,n) = position y of particle n in partition
! kpic = output number of particles per tile
! nppmx = maximum number of particles in tile
! npp = number of particles in partition
! noff = backmost global gridpoint in particle partition
! idimp = size of phase space = 4
! npmax = maximum number of particles in each partition
! mx/my = number of grids in sorting cell in x and y
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer nppmx, idimp, npmax, mx, my, mx1, mxyp1, irc
      integer kpic, npp, noff
      real part, ppart
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
      real delta                                                               ! M. Touati
      dimension delta(2)                                                       ! M. Touati
! local data
      integer i, j, k, n, m, mnoff, ip, ierr
      mnoff = noff
      ierr = 0
! clear counter array
      do 10 k = 1, mxyp1
      kpic(k) = 0
   10 continue
! find addresses of particles at each tile and reorder particles
      do 30 j = 1, npp
!      n = part(1,j)                                                           ! M. Touati
	  n = floor(part(1,j)/delta(1))                                            ! M. Touati
      n = n/mx + 1
!      m = part(2,j)                                                           ! M. Touati
	  m = floor(part(2,j)/delta(2))                                            ! M. Touati
      m = (m - mnoff)/my
      m = n + mx1*m
      ip = kpic(m) + 1
      if (ip.le.nppmx) then
         do 20 i = 1, idimp
         ppart(i,ip,m) = part(i,j)
   20    continue
      else
         ierr = max(ierr,ip-nppmx)
      endif
      kpic(m) = ip
   30 continue
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
      subroutine PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyp1&
     &,irc)
! for 2d code, this subroutine copies segmented particle data ppart to
! the array part with original tiled layout
! spatial decomposition in y direction
! input: all except part, npp, irc, output: part, npp, irc
! part(i,j) = i-th coordinate for particle j in partition
! ppart(i,j,k) = i-th coordinate for particle j in partition in tile k
! kpic = number of particles per tile
! npp = number of particles in partition
! npmax = maximum number of particles in each partition
! nppmx = maximum number of particles in tile
! idimp = size of phase space = 5
! mxyp1 = total number of tiles in partition
! irc = maximum overflow, returned only if error occurs, when irc > 0
      implicit none
      integer npp, npmax, nppmx, idimp, mxyp1, irc
      real part, ppart
      integer kpic
      dimension part(idimp,npmax), ppart(idimp,nppmx,mxyp1)
      dimension kpic(mxyp1)
! local data
      integer i, j, k, npoff, nppp, ne, ierr
      npoff = 0
      ierr = 0
! loop over tiles
      do 30 k = 1, mxyp1
      nppp = kpic(k)
      ne = nppp + npoff
      if (ne.gt.npmax) ierr = max(ierr,ne-npmax)
      if (ierr.gt.0) nppp = 0
! loop over particles in tile
      do 20 j = 1, nppp
      do 10 i = 1, idimp
      part(i,j+npoff) = ppart(i,j,k)
   10 continue
   20 continue
      npoff = npoff + nppp
   30 continue
      npp = npoff
      if (ierr.gt.0) irc = ierr
      return
      end
!-----------------------------------------------------------------------
!      subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1&
!     &,myp1,irc)                                                              ! M. Touati
      subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1&
     &,myp1,nypmx,x,yp,irc)                                                    ! M. Touati
! this subroutine performs a sanity check to make sure particles sorted
! by x,y grid in tiles of mx, my, are all within bounds.
! tiles are assumed to be arranged in 2D linear memory
! input: all except irc
! output: irc
! ppart(1,n,k) = position x of particle n in tile k
! ppart(2,n,k) = position y of particle n in tile k
! kpic(k) = number of reordered output particles in tile k
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! nx = system length in x direction
! mx/my = number of grids in sorting cell in x/y
! mx1 = (system length in x direction - 1)/mx + 1
! myp1 = (partition length in y direction - 1)/my + 1
! irc = particle error, returned only if error occurs, when irc > 0
      implicit none
      integer noff, nyp, idimp, nppmx, nx, mx, my, mx1, myp1, irc
      real ppart
      integer kpic
      dimension ppart(idimp,nppmx,mx1*myp1)
      dimension kpic(mx1*myp1)
      integer nypmx                                                            ! M. Touati
      real x, yp                                                               ! M. Touati
      dimension x(nx+2), yp(nypmx)                                             ! M. Touati
! local data
      integer mxyp1, noffp, moffp, nppp, j, k, ist, nn, mm
      real edgelx, edgely, edgerx, edgery, dx, dy
      mxyp1 = mx1*myp1
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,noffp,moffp,nppp,nn,mm,ist,edgelx,edgely,edgerx,
!$OMP& edgery,dx,dy)
      do 20 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
!      edgelx = noffp                                                          ! M. Touati
!      edgerx = noffp + nn                                                     ! M. Touati
!      edgely = noff + moffp                                                   ! M. Touati
!      edgery = noff + moffp + mm                                              ! M. Touati
      edgelx = x(noffp + 1)                                                    ! M. Touati
      edgerx = x(noffp + nn + 1)                                               ! M. Touati
      edgely = yp(moffp + 1)                                                   ! M. Touati
      edgery = yp(moffp + mm + 1)                                              ! M. Touati
! loop over particles in tile
      do 10 j = 1, nppp
      dx = ppart(1,j,k)
      dy = ppart(2,j,k)
! find particles going out of bounds
      ist = 0
      if (dx.lt.edgelx) ist = 1
      if (dx.ge.edgerx) ist = 2
      if (dy.lt.edgely) ist = ist + 3
      if (dy.ge.edgery) ist = ist + 6
      if (ist.gt.0) irc = k
   10 continue
   20 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
!      subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv, &
!     &nypmx,mx1,mxyp1)                                                        ! M. Touati
      subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv, &
     &nypmx,mx1,mxyp1,x,yp,delta)                                              ! M. Touati
! for 2d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! OpenMP version using guard cells, for distributed data
! data deposited in tiles
! particles stored segmented array
! 17 flops/particle, 6 loads, 4 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n,m)=qm*(1.-dx)*(1.-dy)
! q(n+1,m)=qm*dx*(1.-dy)
! q(n,m+1)=qm*(1.-dx)*dy
! q(n+1,m+1)=qm*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! q(j,k) = charge density at grid point (j,kk),
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! qm = charge on particle, in units of e
! idimp = size of phase space = 4
! nppmx = maximum number of particles in tile
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of charge array, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
      implicit none
      integer noff, idimp, nppmx, mx, my, nxv, nypmx, mx1, mxyp1
      real qm
      real ppart, q
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1), q(nxv,nypmx), kpic(mxyp1)
      real x, yp, delta                                                        ! M. Touati
      dimension x(nxv), yp(nypmx), delta(2)                                    ! M. Touati
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
!      real x, y, dxp, dyp, amx, amy                                           ! M. Touati
      real xx, yy, dxp, dyp, amx, amy                                          ! M. Touati
      real sq
!     dimension sq(MXV,MYV)
      dimension sq(mx+1,my+1)
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!!$OMP PARALLEL DO
!!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,x,y,dxp,dyp,amx,amy,
!!$OMP& sq)                                                                    ! M. Touati
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,mnoff,nn,mm,xx,yy,dxp,dyp,
!$OMP& amx,amy,sq)                                                             ! M. Touati
      do 80 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! zero out local accumulator
      do 20 j = 1, my+1
      do 10 i = 1, mx+1
      sq(i,j) = 0.0
   10 continue
   20 continue
! loop over particles in tile
      do 30 j = 1, nppp
! find interpolation weights
!      x = ppart(1,j,k)                                                        ! M. Touati
!      y = ppart(2,j,k)                                                        ! M. Touati
!      nn = x                                                                  ! M. Touati
!      mm = y                                                                  ! M. Touati
!      dxp = qm*(x - real(nn))                                                 ! M. Touati
!      dyp = y - real(mm)                                                      ! M. Touati
!      nn = nn - noffp + 1                                                     ! M. Touati
!      mm = mm - mnoff                                                         ! M. Touati
!      amx = qm - dxp                                                          ! M. Touati
!      amy = 1.0 - dyp                                                         ! M. Touati
      xx = ppart(1,j,k)                                                        ! M. Touati
      yy = ppart(2,j,k)                                                        ! M. Touati
      nn = 1 + floor(xx/delta(1))                                              ! M. Touati
      mm = 1 + floor(yy/delta(2))                                              ! M. Touati
      dxp = xx - x(nn)                                                         ! M. Touati
      dyp = yy - yp(mm-noff)                                                   ! M. Touati
      nn = nn - noffp                                                          ! M. Touati
      mm = mm - mnoff - 1                                                      ! M. Touati
      amx = (delta(1) - dxp) / (delta(1)**2.)                                  ! M. Touati
      dxp = dxp / (delta(1)**2.)                                               ! M. Touati
      amy = (delta(2) - dyp) / (delta(2)**2.)                                  ! M. Touati
      dyp = dyp / (delta(2)**2.)                                               ! M. Touati
! deposit charge within tile to local accumulator
!      x = sq(nn,mm) + amx*amy                                                 ! M. Touati
!      y = sq(nn+1,mm) + dxp*amy                                               ! M. Touati
!      sq(nn,mm) = x                                                           ! M. Touati
!      sq(nn+1,mm) = y                                                         ! M. Touati
!      x = sq(nn,mm+1) + amx*dyp                                               ! M. Touati
!      y = sq(nn+1,mm+1) + dxp*dyp                                             ! M. Touati
!      sq(nn,mm+1) = x                                                         ! M. Touati
!      sq(nn+1,mm+1) = y                                                       ! M. Touati
      xx = sq(nn,mm) + qm*amx*amy                                              ! M. Touati
      yy = sq(nn+1,mm) + qm*dxp*amy                                            ! M. Touati
      sq(nn,mm) = xx                                                           ! M. Touati
      sq(nn+1,mm) = yy                                                         ! M. Touati
      xx = sq(nn,mm+1) + qm*amx*dyp                                            ! M. Touati
      yy = sq(nn+1,mm+1) + qm*dxp*dyp                                          ! M. Touati
      sq(nn,mm+1) = xx                                                         ! M. Touati
      sq(nn+1,mm+1) = yy                                                       ! M. Touati
   30 continue
! deposit charge to interior points in global array
      nn = min(mx,nxv-noffp)
      mm = min(my,nypmx-moffp)
      do 50 j = 2, mm
      do 40 i = 2, nn
      q(i+noffp,j+moffp) = q(i+noffp,j+moffp) + sq(i,j)
   40 continue
   50 continue
! deposit charge to edge points in global array
      mm = min(my+1,nypmx-moffp)
      do 60 i = 2, nn
!$OMP ATOMIC
      q(i+noffp,1+moffp) = q(i+noffp,1+moffp) + sq(i,1)
      if (mm > my) then
!$OMP ATOMIC
         q(i+noffp,mm+moffp) = q(i+noffp,mm+moffp) + sq(i,mm)
      endif
   60 continue
      nn = min(mx+1,nxv-noffp)
      do 70 j = 1, mm
!$OMP ATOMIC
      q(1+noffp,j+moffp) = q(1+noffp,j+moffp) + sq(1,j)
      if (nn > mx) then
!$OMP ATOMIC
         q(nn+noffp,j+moffp) = q(nn+noffp,j+moffp) + sq(nn,j)
      endif
   70 continue
   80 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
!      subroutine PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,&
!     &idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                       ! M. Touati
      subroutine PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,ci,qbm,dt,dtc,&
     &ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)              ! M. Touati
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored segmented array
! 119 flops/particle, 1 divide, 29 loads, 5 stores
! input: all, output: ppart, ek
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
! omz = (q/m)*bz(x(t),y(t)).
! position equations used are:
! x(t+dt)=x(t) + vx(t+dt/2)*dt
! y(t+dt)=y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
!      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real ci, qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1)
      real x, y , delta                                                        ! M. Touati
      dimension x(nx+2), y(ny+2), delta(2)                                     ! M. Touati
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
!      real x, y                                                               ! M. Touati
      real xx, yy                                                              ! M. Touati
      real sfxy, sbxy
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      sum2 = 0.0d0
! set boundary values
!      edgelx = 0.0                                                            ! M. Touati
!      edgely = 1.0                                                            ! M. Touati
!      edgerx = real(nx)                                                       ! M. Touati
!      edgery = real(ny-1)                                                     ! M. Touati
      edgelx = x(1)                                                            ! M. Touati
      edgely = y(2)                                                            ! M. Touati
      edgerx = x(nx+1)                                                         ! M. Touati
      edgery = y(ny)                                                           ! M. Touati
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
!         edgelx = 1.0                                                         ! M. Touati
!         edgerx = real(nx-1)                                                  ! M. Touati
         edgelx = x(2)                                                         ! M. Touati
         edgerx = x(nx)                                                        ! M. Touati
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,  
!$OMP& dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,sum1,sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)
      do 60 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
!      x = ppart(1,j,k)                                                        ! M. Touati
!      y = ppart(2,j,k)                                                        ! M. Touati
!      nn = x                                                                  ! M. Touati
!      mm = y                                                                  ! M. Touati
!      dxp = x - real(nn)                                                      ! M. Touati
!      dyp = y - real(mm)                                                      ! M. Touati
!      nn = nn - noffp + 1                                                     ! M. Touati
!      mm = mm - mnoff                                                         ! M. Touati
!      amx = 1.0 - dxp                                                         ! M. Touati
!      amy = 1.0 - dyp                                                         ! M. Touati
      xx = ppart(1,j,k)                                                        ! M. Touati
      yy = ppart(2,j,k)                                                        ! M. Touati
      nn = 1 + floor(xx/delta(1))                                              ! M. Touati
      mm = 1 + floor(yy/delta(2))                                              ! M. Touati
      dxp = xx - x(nn)                                                         ! M. Touati
      dyp = yy - y(mm)                                                         ! M. Touati
      nn = nn - noffp                                                          ! M. Touati
      mm = mm - mnoff - 1                                                      ! M. Touati
      amx = (delta(1) - dxp) / delta(1)                                        ! M. Touati
      amy = (delta(2) - dyp) / delta(2)                                        ! M. Touati
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
! time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
! calculate cyclotron frequency
	  !omxt = qtmh*ox                                                          ! M. Touati
      !omyt = qtmh*oy                                                          ! M. Touati
      !omzt = qtmh*oz                                                          ! M. Touati
      omxt = ci*qtmh*ox                                                        ! M. Touati
      omyt = ci*qtmh*oy                                                        ! M. Touati
      omzt = ci*qtmh*oz                                                        ! M. Touati
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
! new position
!      dx = x + dx*dtc                                                         ! M. Touati
!      dy = y + dy*dtc                                                         ! M. Touati
      dx = xx + dx*dtc                                                         ! M. Touati
      dy = yy + dy*dtc                                                         ! M. Touati
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   50 continue
      sum2 = sum2 + sum1
   60 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
!      subroutine PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm&
!     &,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)       ! M. Touati
      subroutine PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,ci,&
     &qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,  &
     &x,y,delta,irc)                                                           ! M. Touati
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored segmented array
! 119 flops/particle, 1 divide, 29 loads, 5 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
! omz = (q/m)*bz(x(t),y(t)).
! position equations used are:
! x(t+dt)=x(t) + vx(t+dt/2)*dt
! y(t+dt)=y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = velocity vx of particle n in partition in tile m
! ppart(4,n,m) = velocity vy of particle n in partition in tile m
! ppart(5,n,m) = velocity vz of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
!      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real ci, qbm, dt, dtc, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
      real x, y, delta                                                         ! M. Touati
      dimension x(nx+2), y(ny+2), delta(2)                                     ! M. Touati
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtmh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
!      real x, y                                                               ! M. Touati
      real xx, yy                                                              ! M. Touati
      real sfxy, sbxy
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
!      anx = real(nx)                                                          ! M. Touati
!      any = real(ny)                                                          ! M. Touati
      anx = x(nx+1)                                                            ! M. Touati
      any = y(ny+1)                                                            ! M. Touati
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!!$OMP PARALLEL DO
!!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!!$OMP& amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1, 
!!$OMP& rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,
!!$OMP& edgery,sum1,sfxy,sbxy)
!!$OMP& REDUCTION(+:sum2)                                                      ! M. Touati
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,xx,yy,dxp,dyp,
!$OMP& amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,
!$OMP& rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,
!$OMP& edgerx,edgery,sum1,sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)                                                       ! M. Touati
      do 70 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
!      edgelx = noffp                                                          ! M. Touati
!      edgerx = noffp + nn                                                     ! M. Touati
!      edgely = noff + moffp                                                   ! M. Touati
!      edgery = noff + moffp + mm                                              ! M. Touati
      edgelx = x(noffp + 1)                                                    ! M. Touati
      edgerx = x(noffp + nn + 1)                                               ! M. Touati
      edgely = y(noff + moffp + 1)                                             ! M. Touati
      edgery = y(noff + moffp + mm + 1)                                        ! M. Touati
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
! load local fields from global arrays
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 60 j = 1, nppp
! find interpolation weights
!      x = ppart(1,j,k)                                                        ! M. Touati
!      y = ppart(2,j,k)                                                        ! M. Touati
!      nn = x                                                                  ! M. Touati
!      mm = y                                                                  ! M. Touati
!      dxp = x - real(nn)                                                      ! M. Touati
!      dyp = y - real(mm)                                                      ! M. Touati
!      nn = nn - noffp + 1                                                     ! M. Touati
!      mm = mm - mnoff                                                         ! M. Touati
!      amx = 1.0 - dxp                                                         ! M. Touati
!      amy = 1.0 - dyp                                                         ! M. Touati
      xx = ppart(1,j,k)                                                        ! M. Touati
      yy = ppart(2,j,k)                                                        ! M. Touati
      nn = 1 + floor(xx/delta(1))                                              ! M. Touati
      mm = 1 + floor(yy/delta(2))                                              ! M. Touati
      dxp = xx - x(nn)                                                         ! M. Touati
      dyp = yy - y(mm)                                                         ! M. Touati
      nn = nn - noffp                                                          ! M. Touati
      mm = mm - mnoff - 1                                                      ! M. Touati
      amx = (delta(1) - dxp) / delta(1)                                        ! M. Touati
      amy = (delta(2) - dyp) / delta(2)                                        ! M. Touati
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
! time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
! calculate cyclotron frequency
      !omxt = qtmh*ox                                                          ! M. Touati
      !omyt = qtmh*oy                                                          ! M. Touati
      !omzt = qtmh*oz                                                          ! M. Touati
      omxt = ci*qtmh*ox                                                        ! M. Touati
      omyt = ci*qtmh*oy                                                        ! M. Touati
      omzt = ci*qtmh*oz                                                        ! M. Touati
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
! new position
!      dx = x + dx*dtc                                                         ! M. Touati
!      dy = y + dy*dtc                                                         ! M. Touati
      dx = xx + dx*dtc                                                         ! M. Touati
      dy = yy + dy*dtc                                                         ! M. Touati
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   60 continue
      sum2 = sum2 + sum1
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   70 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + 0.5*sum2
      return
      end
!-----------------------------------------------------------------------
!      subroutine PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci&
!     &,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                   ! M. Touati
      subroutine PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci&
     &,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)          ! M. Touati
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored segmented array
! 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
! input: all, output: ppart, ek
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum vx of particle n in partition in tile m
! ppart(4,n,m) = momentum vy of particle n in partition in tile m
! ppart(5,n,m) = momentum vz of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! kpic = number of particles per tile
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ipbc
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1)
      real x, y, delta                                                         ! M. Touati
      dimension x(nx+2), y(ny+2), delta(2)                                     ! M. Touati
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, nn, mm
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
!      real x, y                                                               ! M. Touati
      real xx, yy                                                              ! M. Touati
      real sfxy, sbxy
      dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
!     dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum2 = 0.0d0
! set boundary values
!      edgelx = 0.0                                                            ! M. Touati
!      edgely = 1.0                                                            ! M. Touati
!      edgerx = real(nx)                                                       ! M. Touati
!      edgery = real(ny-1)                                                     ! M. Touati
      edgelx = x(1)                                                            ! M. Touati
      edgely = y(2)                                                            ! M. Touati
      edgerx = x(nx+1)                                                         ! M. Touati
      edgery = y(ny)                                                           ! M. Touati
      if ((ipbc.eq.2).or.(ipbc.eq.3)) then
!         edgelx = 1.0                                                         ! M. Touati
!         edgerx = real(nx-1)                                                  ! M. Touati
         edgelx = x(2)                                                         ! M. Touati
         edgerx = x(nx)                                                        ! M. Touati
      endif
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!!$OMP PARALLEL DO
!!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,x,y,dxp,dyp,amx,amy,
!!$OMP& dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,
!!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,sfxy,
!!$OMP& sbxy)
!!$OMP& REDUCTION(+:sum2)                                                      ! M. Touati
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,mnoff,xx,yy,dxp,dyp,amx,amy,
!$OMP& dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,rot2,
!$OMP& rot3,rot4,rot5,rot6,rot7,rot8,rot9,p2,gami,qtmg,dtg,sum1,sfxy,
!$OMP& sbxy)
!$OMP& REDUCTION(+:sum2)                                                      ! M. Touati
      do 60 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      mnoff = moffp + noff - 1
! load local fields from global arrays
      nn = min(mx,nx-noffp) + 1
      mm = min(my,nyp-moffp) + 1
      do 20 j = 1, mm
      do 10 i = 1, nn
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm
      do 30 i = 1, nn
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 50 j = 1, nppp
! find interpolation weights
!      x = ppart(1,j,k)                                                        ! M. Touati
!      y = ppart(2,j,k)                                                        ! M. Touati
!      nn = x                                                                  ! M. Touati
!      mm = y                                                                  ! M. Touati
!      dxp = x - real(nn)                                                      ! M. Touati
!      dyp = y - real(mm)                                                      ! M. Touati
!      nn = nn - noffp + 1                                                     ! M. Touati
!      mm = mm - mnoff                                                         ! M. Touati
!      amx = 1.0 - dxp                                                         ! M. Touati
!      amy = 1.0 - dyp                                                         ! M. Touati
      xx = ppart(1,j,k)                                                        ! M. Touati
      yy = ppart(2,j,k)                                                        ! M. Touati
      nn = 1 + floor(xx/delta(1))                                              ! M. Touati
      mm = 1 + floor(yy/delta(2))                                              ! M. Touati
      dxp = xx - x(nn)                                                         ! M. Touati
      dyp = yy - y(mm)                                                         ! M. Touati
      nn = nn - noffp                                                          ! M. Touati
      mm = mm - mnoff - 1                                                      ! M. Touati
      amx = (delta(1) - dxp) / delta(1)                                        ! M. Touati
      amy = (delta(2) - dyp) / delta(2)                                        ! M. Touati
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
      qtmg = qtmh*gami
! time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      !omxt = qtmg*ox                                                          ! M. Touati
      !omyt = qtmg*oy                                                          ! M. Touati
      !omzt = qtmg*oz                                                          ! M. Touati
      omxt = ci*qtmg*ox                                                        ! M. Touati
      omyt = ci*qtmg*oy                                                        ! M. Touati
      omzt = ci*qtmg*oz                                                        ! M. Touati
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
! update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
!      dx = x + dx*dtg                                                         ! M. Touati
!      dy = y + dy*dtg                                                         ! M. Touati
      dx = xx + dx*dtg                                                         ! M. Touati
      dy = yy + dy*dtg                                                         ! M. Touati
! reflecting boundary conditions
      if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = ppart(2,j,k)
            ppart(4,j,k) = -ppart(4,j,k)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = ppart(1,j,k)
            ppart(3,j,k) = -ppart(3,j,k)
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
   50 continue
      sum2 = sum2 + sum1
   60 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
!-----------------------------------------------------------------------
!      subroutine PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,  &
!     &qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax&
!     &,irc)                                                                   ! M. Touati
      subroutine PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,  &
     &qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax&
     &,x,y,delta,irc)                                                          ! M. Touati
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! with periodic boundary conditions.
! also determines list of particles which are leaving this tile
! OpenMP version using guard cells, for distributed data
! data read in tiles
! particles stored segmented array
! 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
! input: all except ncl, ihole, irc, output: ppart, ncl, ihole, ek, irc
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! ppart(1,n,m) = position x of particle n in partition in tile m
! ppart(2,n,m) = position y of particle n in partition in tile m
! ppart(3,n,m) = momentum vx of particle n in partition in tile m
! ppart(4,n,m) = momentum vy of particle n in partition in tile m
! ppart(5,n,m) = momentum vz of particle n in partition in tile m
! fxy(1,j,k) = x component of force/charge at grid (j,kk)
! fxy(2,j,k) = y component of force/charge at grid (j,kk)
! fxy(3,j,k) = z component of force/charge at grid (j,kk)
! that is, convolution of electric field over particle shape,
! where kk = k + noff - 1
! bxy(1,j,k) = x component of magnetic field at grid (j,kk)
! bxy(2,j,k) = y component of magnetic field at grid (j,kk)
! bxy(3,j,k) = z component of magnetic field at grid (j,kk)
! that is, the convolution of magnetic field over particle shape,
! where kk = k + noff - 1
! kpic(k) = number of particles in tile k
! ncl(i,k) = number of particles going to destination i, tile k
! ihole(1,:,k) = location of hole in array left by departing particle
! ihole(2,:,k) = destination of particle leaving hole
! ihole(1,1,k) = ih, number of holes left (error, if negative)
! noff = lowermost global gridpoint in particle partition.
! nyp = number of primary (complete) gridpoints in particle partition
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprical of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! nx/ny = system length in x/y direction
! mx/my = number of grids in sorting cell in x/y
! nxv = first dimension of field arrays, must be >= nx+1
! nypmx = maximum size of particle partition, including guard cells.
! mx1 = (system length in x direction - 1)/mx + 1
! mxyp1 = mx1*myp1, where myp1=(partition length in y direction-1)/my+1
! ntmax = size of hole array for particles leaving tiles
! irc = maximum overflow, returned only if error occurs, when irc > 0
! optimized version
      implicit none
      integer noff, nyp, idimp, nppmx, nx, ny, mx, my, nxv, nypmx
      integer mx1, mxyp1, ntmax, irc
      real qbm, dt, dtc, ci, ek
      real ppart, fxy, bxy
      integer kpic, ncl, ihole
      dimension ppart(idimp,nppmx,mxyp1)
      dimension fxy(3,nxv,nypmx), bxy(3,nxv,nypmx)
      dimension kpic(mxyp1), ncl(8,mxyp1)
      dimension ihole(2,ntmax+1,mxyp1)
      real x, y, delta                                                         ! M. Touati
      dimension x(nx+2), y(ny+2), delta(2)                                     ! M. Touati
! local data
      integer MXV, MYV
      parameter(MXV=33,MYV=33)
      integer noffp, moffp, nppp
      integer mnoff, i, j, k, ih, nh, nn, mm
      real qtmh, ci2, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anx, any, edgelx, edgely, edgerx, edgery
!      real x, y                                                               ! M. Touati
      real xx, yy                                                              ! M. Touati
      real sfxy, sbxy
!     dimension sfxy(3,MXV,MYV), sbxy(3,MXV,MYV)
      dimension sfxy(3,mx+1,my+1), sbxy(3,mx+1,my+1)
      double precision sum1, sum2
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
!      anx = real(nx)                                                          ! M. Touati
!      any = real(ny)                                                          ! M. Touati
      anx = x(nx+1)                                                            ! M. Touati
      any = y(ny+1)                                                            ! M. Touati
      sum2 = 0.0d0
! error if local array is too small
!     if ((mx.ge.MXV).or.(my.ge.MYV)) return
! loop over tiles
!!$OMP PARALLEL DO
!!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,x,y,dxp,dyp,amx,
!!$OMP& amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,rot1,
!!$OMP& rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,edgerx,
!!$OMP& edgery,p2,gami,qtmg,dtg,sum1,sfxy,sbxy)
!!$OMP& REDUCTION(+:sum2)                                                      ! M. Touati
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,noffp,moffp,nppp,nn,mm,ih,nh,mnoff,xx,yy,dxp,dyp,
!$OMP& amx,amy,dx,dy,dz,ox,oy,oz,acx,acy,acz,omxt,omyt,omzt,omt,anorm,
!$OMP& rot1,rot2,rot3,rot4,rot5,rot6,rot7,rot8,rot9,edgelx,edgely,
!$OMP& edgerx,edgery,p2,gami,qtmg,dtg,sum1,sfxy,sbxy)
!$OMP& REDUCTION(+:sum2)                                                      ! M. Touati
      do 70 k = 1, mxyp1
      noffp = (k - 1)/mx1
      moffp = my*noffp
      noffp = mx*(k - mx1*noffp - 1)
      nppp = kpic(k)
      nn = min(mx,nx-noffp)
      mm = min(my,nyp-moffp)
!      edgelx = noffp                                                          ! M. Touati
!      edgerx = noffp + nn                                                     ! M. Touati
!      edgely = noff + moffp                                                   ! M. Touati
!      edgery = noff + moffp + mm                                              ! M. Touati
      edgelx = x(noffp + 1)                                                    ! M. Touati
      edgerx = x(noffp + nn + 1)                                               ! M. Touati
      edgely = y(noff + moffp + 1)                                             ! M. Touati
      edgery = y(noff + moffp + mm + 1)                                        ! M. Touati
      ih = 0
      nh = 0
      mnoff = moffp + noff - 1
! load local fields from global arrays
      do 20 j = 1, mm+1
      do 10 i = 1, nn+1
      sfxy(1,i,j) = fxy(1,i+noffp,j+moffp)
      sfxy(2,i,j) = fxy(2,i+noffp,j+moffp)
      sfxy(3,i,j) = fxy(3,i+noffp,j+moffp)
   10 continue
   20 continue
      do 40 j = 1, mm+1
      do 30 i = 1, nn+1
      sbxy(1,i,j) = bxy(1,i+noffp,j+moffp)
      sbxy(2,i,j) = bxy(2,i+noffp,j+moffp)
      sbxy(3,i,j) = bxy(3,i+noffp,j+moffp)
   30 continue
   40 continue
! clear counters
      do 50 j = 1, 8
      ncl(j,k) = 0
   50 continue
      sum1 = 0.0d0
! loop over particles in tile
      do 60 j = 1, nppp
! find interpolation weights
!      x = ppart(1,j,k)                                                        ! M. Touati
!      y = ppart(2,j,k)                                                        ! M. Touati
!      nn = x                                                                  ! M. Touati
!      mm = y                                                                  ! M. Touati
!      dxp = x - real(nn)                                                      ! M. Touati
!      dyp = y - real(mm)                                                      ! M. Touati
!      nn = nn - noffp + 1                                                     ! M. Touati
!      mm = mm - mnoff                                                         ! M. Touati
!      amx = 1.0 - dxp                                                         ! M. Touati
!      amy = 1.0 - dyp                                                         ! M. Touati
      xx = ppart(1,j,k)                                                        ! M. Touati
      yy = ppart(2,j,k)                                                        ! M. Touati
      nn = 1 + floor(xx/delta(1))                                              ! M. Touati
      mm = 1 + floor(yy/delta(2))                                              ! M. Touati
      dxp = xx - x(nn)                                                         ! M. Touati
      dyp = yy - y(mm)                                                         ! M. Touati
      nn = nn - noffp                                                          ! M. Touati
      mm = mm - mnoff - 1                                                      ! M. Touati
      amx = (delta(1) - dxp) / delta(1)                                        ! M. Touati
      amy = (delta(2) - dyp) / delta(2)                                        ! M. Touati
! find electric field
      dx = amx*sfxy(1,nn,mm)
      dy = amx*sfxy(2,nn,mm)
      dz = amx*sfxy(3,nn,mm)
      dx = amy*(dxp*sfxy(1,nn+1,mm) + dx)
      dy = amy*(dxp*sfxy(2,nn+1,mm) + dy)
      dz = amy*(dxp*sfxy(3,nn+1,mm) + dz)
      acx = amx*sfxy(1,nn,mm+1)
      acy = amx*sfxy(2,nn,mm+1)
      acz = amx*sfxy(3,nn,mm+1)
      dx = dx + dyp*(dxp*sfxy(1,nn+1,mm+1) + acx) 
      dy = dy + dyp*(dxp*sfxy(2,nn+1,mm+1) + acy)
      dz = dz + dyp*(dxp*sfxy(3,nn+1,mm+1) + acz)
! find magnetic field
      ox = amx*sbxy(1,nn,mm)
      oy = amx*sbxy(2,nn,mm)
      oz = amx*sbxy(3,nn,mm)
      ox = amy*(dxp*sbxy(1,nn+1,mm) + ox)
      oy = amy*(dxp*sbxy(2,nn+1,mm) + oy)
      oz = amy*(dxp*sbxy(3,nn+1,mm) + oz)
      acx = amx*sbxy(1,nn,mm+1)
      acy = amx*sbxy(2,nn,mm+1)
      acz = amx*sbxy(3,nn,mm+1)
      ox = ox + dyp*(dxp*sbxy(1,nn+1,mm+1) + acx) 
      oy = oy + dyp*(dxp*sbxy(2,nn+1,mm+1) + acy)
      oz = oz + dyp*(dxp*sbxy(3,nn+1,mm+1) + acz)
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = ppart(3,j,k) + dx
      acy = ppart(4,j,k) + dy
      acz = ppart(5,j,k) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! renormalize magnetic field
      qtmg = qtmh*gami
! time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      !omxt = qtmg*ox                                                          ! M. Touati
      !omyt = qtmg*oy                                                          ! M. Touati
      !omzt = qtmg*oz                                                          ! M. Touati
      omxt = ci*qtmg*ox                                                        ! M. Touati
      omyt = ci*qtmg*oy                                                        ! M. Touati
      omzt = ci*qtmg*oz                                                        ! M. Touati
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      ppart(3,j,k) = dx
      ppart(4,j,k) = dy
      ppart(5,j,k) = dz
! update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
!      dx = x + dx*dtg                                                         ! M. Touati
!      dy = y + dy*dtg                                                         ! M. Touati
      dx = xx + dx*dtg                                                         ! M. Touati
      dy = yy + dy*dtg                                                         ! M. Touati
! find particles going out of bounds
      mm = 0
! count how many particles are going in each direction in ncl
! save their address and destination in ihole
! use periodic boundary conditions and check for roundoff error
! mm = direction particle is going
      if (dx.ge.edgerx) then
         if (dx.ge.anx) dx = dx - anx
         mm = 2
      else if (dx.lt.edgelx) then
         if (dx.lt.0.0) then
            dx = dx + anx
            if (dx.lt.anx) then
               mm = 1
            else
               dx = 0.0
            endif
         else
            mm = 1
         endif
      endif
      if (dy.ge.edgery) then
         if (dy.ge.any) dy = dy - any
         mm = mm + 6
      else if (dy.lt.edgely) then
         if (dy.lt.0.0) then
            dy = dy + any
            if (dy.lt.any) then
               mm = mm + 3
            else
               dy = 0.0
            endif
         else
            mm = mm + 3
         endif
      endif
! set new position
      ppart(1,j,k) = dx
      ppart(2,j,k) = dy
! increment counters
      if (mm.gt.0) then
         ncl(mm,k) = ncl(mm,k) + 1
         ih = ih + 1
         if (ih.le.ntmax) then
            ihole(1,ih+1,k) = j
            ihole(2,ih+1,k) = mm
         else
            nh = 1
         endif
      endif
   60 continue
      sum2 = sum2 + sum1
! set error and end of file flag
! ihole overflow
      if (nh.gt.0) then
         irc = ih
         ih = -ih
      endif
      ihole(1,1,k) = ih
   70 continue
!$OMP END PARALLEL DO
! normalize kinetic energy
      ek = ek + sum2
      return
      end
