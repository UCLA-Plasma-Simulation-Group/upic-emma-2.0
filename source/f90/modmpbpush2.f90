!-----------------------------------------------------------------------
!
      module modmpbpush2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpbpush2.f
! mpmovin2 reorder and copy particles to ordered array
!          calls PPPMOVIN2L
! mpcopyout2 copy ordered particles to unordered array
!            calls PPPCOPYOUT
! mpcheck2 verify particles are all in correct tiles
!          calls PPPCHECK2L
! mppost2 deposits charge density
!         calls PPGPPOST2L
! mpbpush2 push magnetized particles
!          calls PPGBPPUSH23L
! mpbpushf2 push magnetized particles and determine which particles are
!           leaving tile
!           calls PPGBPPUSHF23L
! mprbpush2 push relativistic, magnetized particles
!           calls PPGRBPPUSH23L
! mprbpushf2 push relativistic, magnetized particles and determine which
!            particles are leaving tile
!            calls PPGRBPPUSHF23L
! wmpbpush2 generic procedure to push magnetized particles
!           calls mprbpushf2, mpbpushf2, mprbpush2, or mpbpush2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 29, 2016
!
      use libmpbpush2_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
!      subroutine mpmovin2(part,ppart,kpic,npp,noff,mx,my,mx1,irc)             ! M. Touati
      subroutine mpmovin2(part,ppart,kpic,npp,noff,mx,my,mx1,delta,irc)        ! M. Touati
! order particles in part by tiles and copy to ppart
! store number of particles in each tile in kpic
      implicit none
      integer, intent(in) :: mx, my, mx1
      integer, intent(in) :: npp, noff
      integer, intent(inout) :: irc
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(inout) :: kpic
      real, dimension(:), intent(in) :: delta                                  ! M. Touati
! local data
      integer :: idimp, npmax, nppmx, mxyp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
! call low level procedure
!      call PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax,mx,my, &
!     &mx1,mxyp1,irc)                                                          ! M. Touati
      call PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,npmax,mx,my, &
     &mx1,mxyp1,delta,irc)                                                     ! M. Touati
! check for errors
      if (irc /= 0) then
         write (*,*) 'mpmovin2 overflow error, irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mpcopyout2(part,ppart,kpic,npp,irc)
! copy ordered particles in array ppart to unordered array part
! store total particle number in npp
      implicit none
      integer, intent(inout) :: npp, irc
      real, dimension(:,:), intent(inout) :: part
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
! local data
      integer :: idimp, npmax, nppmx, mxyp1
! extract dimensions
      idimp = size(part,1); npmax = size(part,2)
      nppmx = size(ppart,2)
      mxyp1 = size(kpic,1)
! call low level procedure
      call PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,mxyp1,irc)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,irc)               ! M. Touati
      subroutine mpcheck2(ppart,kpic,noff,nyp,nx,mx,my,mx1,nypmx,x,yp,irc)     ! M. Touati
! perform a sanity check to make sure particles ordered by tiles are all
! within bounds.
      implicit none
      integer, intent(in) :: nx, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer :: nypmx                                                         ! M. Touati
      real, dimension(:) :: x, yp                                              ! M. Touati
! local data
      integer :: idimp, nppmx, myp1
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      myp1 = size(kpic,1)/mx1
! call low level procedure
!      call PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1,myp1,&
!     &irc)                                                                    ! M. Touati
      call PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,mx1,myp1,&
     &nypmx,x,yp,irc)                                                          ! M. Touati
! check error
      if (irc /= 0) then
         write (*,*) 'mpcheck2 error: irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mppost2(ppart,q,kpic,noff,qm,tdpost,mx,my,mx1)               ! M. Touati
      subroutine mppost2(ppart,q,kpic,noff,qm,tdpost,mx,my,mx1,x,yp,delta)     ! M. Touati
! deposit charge
      implicit none
      integer, intent(in) :: mx, my, mx1
      integer, intent(in) :: noff
      real, intent(in) :: qm
      real, intent(inout) :: tdpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:), intent(inout) :: q
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, yp, delta                           ! M. Touati 
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(q,1); nypmx = size(q,2)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv,nypmx, &
!     &mx1,mxyp1)                                                              ! M. Touati
      call PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,nxv,nypmx, &
     &mx1,mxyp1,x,yp,delta)                                                    ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,    &
!     &tpush,nx,ny,mx,my,mx1,ipbc)                                             ! M. Touati
      subroutine mpbpush2(ppart,fxy,bxy,kpic,noff,nyp,ci,qbm,dt,dtc,ek,  &
     &tpush,nx,ny,mx,my,mx1,x,y,delta,ipbc)                                    ! M. Touati
! push magnetized particles with 2d electromagnetic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc
      integer, intent(in) :: noff, nyp
      real, intent(in) :: ci, qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,idimp,&
!     &nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                             ! M. Touati
      call PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,ci,qbm,dt,dtc,ek,   &
     &idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)                    ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
!     &dtc,ek,tpush,nx,ny,mx,my,mx1,irc)                                       ! M. Touati
      subroutine mpbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,ci,qbm,&
     &dt,dtc,ek,tpush,nx,ny,mx,my,mx1,x,y,delta,irc)                           ! M. Touati
! push magnetized particles with 2d electromagnetic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: ci, qbm, dt, dtc
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,  &
!     &dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)           ! M. Touati
      call PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,ci,qbm,  &
     &dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,x,y,     &
     &delta,irc)                                                               ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mprbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,&
!     &tpush,nx,ny,mx,my,mx1,ipbc)                                             ! M. Touati
      subroutine mprbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,&
     &tpush,nx,ny,mx,my,mx1,x,y,delta,ipbc)                                    ! M. Touati
! push relativistic, magnetized particles with 2d electromagnetic fields
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc
      integer, intent(in) :: noff, nyp
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,  &
!     &idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                       ! M. Touati
      call PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,  &
     &idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)              ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mprbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
!     &,dtc,ci,ek,tpush,nx,ny,mx,my,mx1,irc)                                   ! M. Touati
      subroutine mprbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
     &,dtc,ci,ek,tpush,nx,ny,mx,my,mx1,x,y,delta,irc)                          ! M. Touati
! push relativistic, magnetized particles with 2d electromagnetic fields
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(fxy,2); nypmx = size(fxy,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt, &
!     &dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)        ! M. Touati
      call PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt, &
     &dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,x,y,  &
     &delta,irc)                                                               ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine wmpbpush2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
!     &dtc,ci,ek,tpush,nx,ny,mx,my,mx1,ipbc,relativity,list,irc)               ! M. Touati
      subroutine wmpbpush2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
     &dtc,ci,ek,tpush,nx,ny,mx,my,mx1,ipbc,relativity,list,x,y,delta,irc)      ! M. Touati
! generic procedure to push magnetized particles
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, relativity
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      logical, intent(in) :: list
      real, intent(in) :: qbm, dt, dtc, ci
      real, intent(inout) :: ek, tpush
      real, dimension(:,:,:), intent(inout)  :: ppart
      real, dimension(:,:,:), intent(in) :: fxy, bxy
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout)  :: ncl
      integer, dimension(:,:,:), intent(inout)  :: ihole
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! also calculate list of particles leaving tile
      if (list) then
! updates ppart, wke, ncl, ihole, irc
         if (relativity==1) then
!            call mprbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
!     &,dtc,ci,ek,tpush,nx,ny,mx,my,mx1,irc)                                   ! M. Touati
            call mprbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt&
     &,dtc,ci,ek,tpush,nx,ny,mx,my,mx1,x,y,delta,irc)                          ! M. Touati
         else
!            call mpbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,qbm,dt,&
!     &dtc,ek,tpush,nx,ny,mx,my,mx1,irc)                                       ! M. Touati
            call mpbpushf2(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,ci,qbm,&
     &dt,dtc,ek,tpush,nx,ny,mx,my,mx1,x,y,delta,irc)                              ! M. Touati
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart and wke
         if (relativity==1) then
!            call mprbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,&
!     &tpush,nx,ny,mx,my,mx1,ipbc)                                             ! M. Touati
            call mprbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ci,ek,&
     &tpush,nx,ny,mx,my,mx1,x,y,delta,ipbc)                                    ! M. Touati
         else
!            call mpbpush2(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,ek,    &
!     &tpush,nx,ny,mx,my,mx1,ipbc)                                             ! M. Touati
            call mpbpush2(ppart,fxy,bxy,kpic,noff,nyp,ci,qbm,dt,dtc,ek, &
     &tpush,nx,ny,mx,my,mx1,x,y,delta,ipbc)                                    ! M. Touati
         endif
      endif
      end subroutine
!
      end module
