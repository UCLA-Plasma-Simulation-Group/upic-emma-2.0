!-----------------------------------------------------------------------
! Interface file for libmpbpush2.f
      module libmpbpush2_h
      implicit none
!
      interface
!         subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
!     &npmax,mx,my,mx1,mxyp1,irc)                                              ! M. Touati
         subroutine PPPMOVIN2L(part,ppart,kpic,npp,noff,nppmx,idimp,    &
     &npmax,mx,my,mx1,mxyp1,delta,irc)                                         ! M. Touati
         implicit none
         integer, intent(in) :: nppmx, idimp, npmax, mx, my, mx1, mxyp1
         integer, intent(in) :: npp, noff
         integer, intent(inout) :: irc
         real, dimension(idimp,npmax), intent(in) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         integer, dimension(mxyp1), intent(inout) :: kpic
         real, dimension(2), intent(in) :: delta
         end subroutine
      end interface
!
      interface
         subroutine PPPCOPYOUT2(part,ppart,kpic,npp,npmax,nppmx,idimp,  &
     &mxyp1,irc)
         implicit none
         integer, intent(in) :: npmax, nppmx, idimp, mxyp1
         integer, intent(inout) :: npp, irc
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         integer, dimension(mxyp1), intent(in) :: kpic
         end subroutine
      end interface
!
      interface
!         subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,&
!     &mx1,myp1,irc)                                                           ! M. Touati
         subroutine PPPCHECK2L(ppart,kpic,noff,nyp,idimp,nppmx,nx,mx,my,&
     &mx1,myp1,nypmx,x,yp,irc)                                                 ! M. Touati
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, my, mx1, myp1
         integer, intent(in) :: noff, nyp
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(in) :: ppart
         integer, dimension(mx1*myp1), intent(in) :: kpic
         integer, intent(in) :: nypmx                                          ! M. Touati 
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(nypmx), intent(in) :: yp                              ! M. Touati          
         end subroutine
      end interface
!
      interface
!         subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,  &
!     &nxv,nypmx,mx1,mxyp1)                                                    ! M. Touati
         subroutine PPGPPOST2L(ppart,q,kpic,noff,qm,idimp,nppmx,mx,my,  &
     &nxv,nypmx,mx1,mxyp1,x,yp,delta)                                                ! M. Touati
         implicit none
         integer, intent(in) :: idimp, nppmx, mx, my, nxv, nypmx
         integer, intent(in) :: mx1, mxyp1
         integer, intent(in) :: noff
         real, intent(in) :: qm
         real, dimension(idimp,nppmx,mxyp1), intent(in) :: ppart
         real, dimension(nxv,nypmx), intent(inout) :: q
         integer, dimension(mxyp1), intent(in) :: kpic
         real, dimension(nxv), intent(in) :: x
         real, dimension(nypmx), intent(in) :: yp
         real, dimension(2), intent(in) :: delta
         end subroutine
      end interface
!
      interface
!         subroutine PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc,&
!     &ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                    ! M. Touati
         subroutine PPGBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,ci,qbm,dt, &
     &dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)           ! M. Touati
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: ci, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(ny+2), intent(in) :: y                                ! M. Touati
         real, dimension(2), intent(in) :: delta                               ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,&
!     &qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax,  &
!     &irc)                                                                    ! M. Touati
         subroutine PPGBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp,&
     &ci,qbm,dt,dtc,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ntmax&
     &,x,y,delta,irc)                                                           ! M. Touati
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: ci, qbm, dt, dtc
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout)  :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout)  :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout)  :: ihole
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(ny+2), intent(in) :: y                                ! M. Touati
         real, dimension(2), intent(in) :: delta                               ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc&
!     &,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,ipbc)                ! M. Touati
         subroutine PPGRBPPUSH23L(ppart,fxy,bxy,kpic,noff,nyp,qbm,dt,dtc&
     &,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,x,y,delta,ipbc)       ! M. Touati
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ipbc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(ny+2), intent(in) :: y                                ! M. Touati
         real, dimension(2), intent(in) :: delta                               ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp&
!     &,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,    &
!     &ntmax,irc)                                                              ! M. Touati
         subroutine PPGRBPPUSHF23L(ppart,fxy,bxy,kpic,ncl,ihole,noff,nyp&
     &,qbm,dt,dtc,ci,ek,idimp,nppmx,nx,ny,mx,my,nxv,nypmx,mx1,mxyp1,    &
     &ntmax,x,y,delta,irc)                                                     ! M. Touati
         implicit none
         integer, intent(in) :: noff, nyp, idimp, nppmx, nx, ny, mx, my
         integer, intent(in) :: nxv, nypmx, mx1, mxyp1, ntmax
         integer, intent(inout) :: irc
         real, intent(in) :: qbm, dt, dtc, ci
         real, intent(inout) :: ek
         real, dimension(idimp,nppmx,mxyp1), intent(inout) :: ppart
         real, dimension(3,nxv,nypmx), intent(in) :: fxy, bxy
         integer, dimension(mxyp1), intent(in) :: kpic
         integer, dimension(8,mxyp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mxyp1), intent(inout) :: ihole
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(ny+2), intent(in) :: y                                ! M. Touati
         real, dimension(2), intent(in) :: delta                               ! M. Touati
         end subroutine
      end interface
!
      end module