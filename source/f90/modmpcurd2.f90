!-----------------------------------------------------------------------
!
      module modmpcurd2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpcurd2.f
! mpdjpost2 deposit current and update particle positions
!           calls PPGJPPOST2L
! mpdjpostf2 deposit current and update particle positions
!            determine which particles are leaving tile
!            calls PPGJPPOSTF2L
! mprdjpost2 deposit current and update relativistic particle positions
!            calls PPGRJPPOST2L
! mprdjpostf2 deposit current and update relativistic particle positions
!             determine which particles are leaving tile
!             calls PPGRJPPOSTF2L
! mpgmjpost2 deposits momentum flux
!            calls PPGMJPPOST2L
! mpgrmjpost2 deposits relativistic momentum flux
!             calls PPGRMJPPOST2L
! wmpdjpost2 generic procedure to deposit current and update particle
!            positions
!            calls mprdjpostf2, mpdjpostf2, mprdjpost2, or mpdjpost2
! wmpgmjpost2 generic procedure to deposit momentum flux
!             calls mpgrmjpost2, or mpgmjpost2
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: february 1, 2016
!
      use libmpcurd2_h
      implicit none
!
      contains  
!
!-----------------------------------------------------------------------
!      subroutine mpdjpost2(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,mx,my,&
!     &mx1,ipbc)                                                               ! M. Touati
      subroutine mpdjpost2(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,mx,my,&
     &mx1,x,y,delta,ipbc)                                                      ! M. Touati
! deposit current and update particle positions
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc
      integer, intent(in) :: noff
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny,mx,my,&
!     &nxv,nypmx,mx1,mxyp1,ipbc)                                               ! M. Touati
      call PPGJPPOST2L(ppart,cu,kpic,noff,qm,dt,nppmx,idimp,nx,ny,mx,my,&
     &nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)                                      ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,     &
!     &tdjpost,nx,ny,mx,my,mx1,irc)                                            ! M. Touati
      subroutine mpdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,     &
     &tdjpost,nx,ny,mx,my,mx1,x,y,delta,irc)                                   ! M. Touati
! deposit current and update particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,nppmx,   &
!     &idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)                        ! M. Touati
      call PPGJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,nppmx,   &
     &idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,x,y,delta,irc)               ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
! check error
      if (irc /= 0) then
         write (*,*) 'mpdjpostf2 error: irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mprdjpost2(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,mx&
!     &,my,mx1,ipbc)                                                           ! M. Touati
      subroutine mprdjpost2(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,mx&
     &,my,mx1,x,y,delta,ipbc)                                                  ! M. Touati
! deposit current and update relativistic particle positions
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc
      integer, intent(in) :: noff
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx,ny,mx&
!     &,my,nxv,nypmx,mx1,mxyp1,ipbc)                                           ! M. Touati
      call PPGRJPPOST2L(ppart,cu,kpic,noff,qm,dt,ci,nppmx,idimp,nx,ny,mx&
     &,my,nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)                                  ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mprdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci, &
!     &tdjpost,nx,ny,mx,my,mx1,irc)                                            ! M. Touati
      subroutine mprdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci, &
     &tdjpost,nx,ny,mx,my,mx1,x,y,delta,irc)                                   ! M. Touati
! deposit current and update relativistic particle positions
! determine which particles are leaving tile
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! local data
      integer :: idimp, nppmx, nxv, nypmx, mxyp1, ntmax
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nxv = size(cu,2); nypmx = size(cu,3)
      mxyp1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci,nppmx&
!     &,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,irc)                       ! M. Touati
      call PPGRJPPOSTF2L(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci,nppmx&
     &,idimp,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,x,y,delta,irc)              ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
! check error
      if (irc /= 0) then
         write (*,*) 'mprdjpostf2 error: irc=', irc
         call PPABORT()
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpgmjpost2(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mx1)        ! M. Touati
      subroutine mpgmjpost2(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mx1,  &
     &x,yp,delta)                                                              ! M. Touati
! deposit momentum flux
      implicit none
      integer, intent(in) :: mx, my, mx1
      integer, intent(in) :: noff
      real, intent(in) ::  qm
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, yp, delta                           ! M. Touati
! local data
      integer :: idimp, nppmx, mdim, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mdim = size(amu,1); nxv = size(amu,2); nypmx = size(amu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(mdim)
      case (4)
!         call PPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my,nxv,&
!     &nypmx,mx1,mxyp1)                                                        ! M. Touati
         call PPGMJPPOST2L(ppart,amu,kpic,noff,qm,nppmx,idimp,mx,my,nxv,&
     &nypmx,mx1,mxyp1,x,yp,delta)                                              ! M. Touati
      case default
         write (*,*) 'mpgmjpost2: unsupported dimension mdim = ', mdim
      end select
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine mpgrmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,  &
!     &mx1)                                                                    ! M. Touati
      subroutine mpgrmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,  &
     &mx1,x,yp,delta)                                                          ! M. Touati
! deposit relativistic momentum flux
      implicit none
      integer, intent(in) :: mx, my, mx1
      integer, intent(in) :: noff
      real, intent(in) ::  qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, yp, delta                           ! M. Touati
! local data
      integer :: idimp, nppmx, mdim, nxv, nypmx, mxyp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mdim = size(amu,1); nxv = size(amu,2); nypmx = size(amu,3)
      mxyp1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(mdim)
      case (4)
!         call PPGRMJPPOST2L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx,my,&
!     &nxv,nypmx,mx1,mxyp1)                                                    ! M. Touati
         call PPGRMJPPOST2L(ppart,amu,kpic,noff,qm,ci,nppmx,idimp,mx,my,&
     &nxv,nypmx,mx1,mxyp1,x,yp,delta)                                          ! M. Touati
      case default
         write (*,*) 'mpgrmjpost2: unsupported dimension mdim = ', mdim
      end select
! record time
      call dtimer(dtime,itime,1)
      tdcjpost = tdcjpost + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine wmpdjpost2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci,  &
!     &tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,list,irc)                       ! M. Touati
      subroutine wmpdjpost2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci,  &
     &tdjpost,nx,ny,mx,my,mx1,ipbc,relativity,list,x,y,delta,irc)              ! M. Touati
! generic procedure to deposit current and update particle positions
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, ny, mx, my, mx1, ipbc, relativity
      integer, intent(in) :: noff, nyp
      integer, intent(inout) :: irc
      logical, intent(in) :: list
      real, intent(in) :: qm, dt, ci
      real, intent(inout) :: tdjpost
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: cu
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      real, dimension(:), intent(in) :: x, y, delta                            ! M. Touati
! also calculate list of particles leaving tile
      if (list) then
! updates ppart, cue, ncl, ihole, irc
         if (relativity==1) then
!            call mprdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci, &
!     &tdjpost,nx,ny,mx,my,mx1,irc)                                            ! M. Touati
            call mprdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,ci, &
     &tdjpost,nx,ny,mx,my,mx1,x,y,delta,irc)                                   ! M. Touati
         else
!            call mpdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,     &
!     &tdjpost,nx,ny,mx,my,mx1,irc)                                            ! M. Touati
            call mpdjpostf2(ppart,cu,kpic,ncl,ihole,noff,nyp,qm,dt,     &
     &tdjpost,nx,ny,mx,my,mx1,x,y,delta,irc)                                   ! M. Touati
         endif
! do not also calculate list of particles leaving tile
      else
! updates ppart, cue
         if (relativity==1) then
!            call mprdjpost2(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,mx&
!     &,my,mx1,ipbc)                                                           ! M. Touati
            call mprdjpost2(ppart,cu,kpic,noff,qm,dt,ci,tdjpost,nx,ny,mx&
     &,my,mx1,x,y,delta,ipbc)                                                  ! M. Touati
         else
!            call mpdjpost2(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,mx,my,&
!     &mx1,ipbc)                                                               ! M. Touati
            call mpdjpost2(ppart,cu,kpic,noff,qm,dt,tdjpost,nx,ny,mx,my,&
     &mx1,x,y,delta,ipbc)                                                      ! M. Touati
         endif
      endif
      end subroutine
!
!-----------------------------------------------------------------------
!      subroutine wmpgmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,  &
!     &mx1,relativity)                                                         ! M. Touati
      subroutine wmpgmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,  &
     &mx1,relativity,x,yp,delta)                                               ! M. Touati
! generic procedure to deposit momentum flux
      implicit none
      integer, intent(in) :: noff, mx, my, mx1, relativity
      real, intent(in) ::  qm, ci
      real, intent(inout) :: tdcjpost
      real, dimension(:,:,:), intent(in) :: ppart
      real, dimension(:,:,:), intent(inout) :: amu
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:), intent(in) :: x, yp, delta                           ! M. Touati
! updates amu
      if (relativity==1) then
!         call mpgrmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,mx1)       ! M. Touati
         call mpgrmjpost2(ppart,amu,kpic,noff,qm,ci,tdcjpost,mx,my,     &
     &mx1,x,yp,delta)                                                          ! M. Touati
      else
!         call mpgmjpost2(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mx1)           ! M. Touati
         call mpgmjpost2(ppart,amu,kpic,noff,qm,tdcjpost,mx,my,mx1,x,yp,&
     &delta)                                                                   ! M. Touati
      endif
      end subroutine
!
      end module
