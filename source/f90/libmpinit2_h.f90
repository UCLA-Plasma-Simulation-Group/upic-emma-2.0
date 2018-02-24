!-----------------------------------------------------------------------
! Interface file for libmpinit2.f
      module libmpinit2_h
      implicit none
!
      interface
         subroutine PDICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp,  &
     &idps)
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
         real, dimension(idps), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
!         subroutine PDNICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp, &
!     &idps)                                                                   ! M. Touati
		 subroutine PDNICOMP2L(edges,nyp,noff,nypmx,nypmn,ny,kstrt,nvp, &
      &idps,delta)                                                             ! M. Touati
         implicit none
         integer, intent(in) :: ny, kstrt, nvp, idps
         integer, intent(inout) :: nyp, noff, nypmx, nypmn
         real, dimension(2), intent(in) :: delta                               ! M. Touati 
         real, dimension(idps), intent(inout) :: edges
         end subroutine
      end interface
!
      interface
!         subroutine PDISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy, &
!     &nx,ny,idimp,npmax,idps,ipbc,ierr)                                       ! M. Touati
         subroutine PDISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy, &
     &nx,ny,idimp,npmax,idps,ipbc,x,y,ierr)                                    ! M. Touati
         implicit none
         integer, intent(in) :: nps, npx, npy, nx, ny, idimp, npmax
         integer, intent(in) :: idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(ny+2), intent(in) :: y                                ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
!     &npx,npy,nx,ny,idimp,npmax,idps,ipbc,ierr)                               ! M. Touati
         subroutine PDISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
     &npx,npy,nx,ny,idimp,npmax,idps,ipbc,x,y,ierr)                            ! M. Touati
         implicit none
         integer, intent(in) :: nps, npx, npy, nx, ny, idimp, npmax
         integer, intent(in) :: idps, ipbc
         integer, intent(inout) :: npp, ierr
         real, intent(in) :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax), intent(inout) :: part
         real, dimension(idps), intent(in) :: edges
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(ny+2), intent(in) :: y                                ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
!     &,mx1,mxyp1,irc)                                                         ! M. Touati
         subroutine PPDBLKP2L(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my&
     &,mx1,mxyp1,delta,irc)                                                    ! M. Touati
         implicit none
         integer, intent(in) :: idimp, npmax, mx, my, mx1, mxyp1, npp
         integer, intent(in) :: noff
         integer, intent(inout) :: nppmx, irc
         real, dimension(idimp,npmax), intent(in) :: part
         integer, dimension(mxyp1), intent(inout) :: kpic
         real, dimension(2), intent(in) :: delta                               ! M. Touati
         end subroutine
      end interface
!
      interface
         function ranorm()
         implicit none
         double precision :: ranorm
         end function
      end interface
!
      interface
         function randum()
         implicit none
         double precision :: randum
         end function
      end interface
!
      end module
