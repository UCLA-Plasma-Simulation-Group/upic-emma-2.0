!-----------------------------------------------------------------------
!
      module modmpsort2
!
! Fortran90 wrappers to 2d MPI/OpenMP PIC library libmpsort2.f
! mporder2a creates list of particles which are leaving tile, and
!           buffers outgoing particles
!           calls PPPORDER2LA
! mporderf2a buffers outgoing particles
!            calls PPPORDERF2LA
! mporder2b copies buffers into particle array
!           calls PPPORDER2LB
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: january 14, 2016
!
      use libmpsort2_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
!      subroutine mporder2a(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,&
!     &nclr,noff,nyp,tsort,nx,ny,mx,my,irc)                                    ! M. Touati
      subroutine mporder2a(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,&
     &nclr,noff,nyp,tsort,nx,ny,mx,my,x,y,irc)                                 ! M. Touati
! performs first part of particle reordering into tiles
      implicit none
      integer, intent(in) :: nx, ny, mx, my
      integer, intent(in) :: noff, nyp
      real, intent(inout) :: tsort
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
      integer, dimension(:,:), intent(inout) :: ncll, nclr
      real, dimension(:), intent(in) :: x, y                                   ! M. Touati      
! local data
      integer :: idimp, nppmx, npbmx, nbmax, ntmax, mx1, myp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(sbufl,2)
      ntmax = size(ihole,2) - 1
      mx1 = size(ncll,2); myp1 = size(kpic,1)/mx1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
!      call PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr&
!     &,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,nbmax,irc)       ! M. Touati
      call PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr&
     &,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,nbmax,x,y,irc)    ! M. Touati
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporderf2a(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr&
     &,tsort,irc)
! performs first part of particle reordering into tiles,
! does not create list of particles which are leaving tile
      implicit none
      real, intent(inout) :: tsort
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      real, dimension(:,:), intent(inout) :: sbufl, sbufr
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
      integer, dimension(:,:), intent(inout) :: ncll, nclr
! local data
      integer :: idimp, nppmx, npbmx, nbmax, mx1, ntmax, myp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(sbufl,2)
      mx1 = size(ncll,2)
      ntmax = size(ihole,2) - 1; myp1 = size(ihole,3)/mx1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,   &
     &idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mporder2b(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,&
     &mclr,tsort,irc)
! performs second part of particle reordering into tiles
      implicit none
      real, intent(inout) :: tsort
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(in) :: ppbuff
      real, dimension(:,:), intent(in) :: rbufl, rbufr
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(in) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
      integer, dimension(:,:), intent(in) :: mcll, mclr
! local data
      integer :: idimp, nppmx, npbmx, nbmax, mx1, ntmax, myp1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2); nbmax = size(rbufl,2)
      mx1 = size(mcll,2)
      ntmax = size(ihole,2) - 1; myp1 = size(kpic,1)/mx1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr&
     &,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
      end module
