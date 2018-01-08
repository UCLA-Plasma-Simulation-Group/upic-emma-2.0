!-----------------------------------------------------------------------
! Interface file for libmpsort2.f
      module libmpsort2_h
      implicit none
!
      interface
!         subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,&  
!     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
!     &nbmax,irc)                                                              ! M. Touati
         subroutine PPPORDER2LA(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,&  
     &ncll,nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,myp1,npbmx,ntmax,  &
     &nbmax,x,y,irc)                                                           ! M. Touati
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, ny, mx, my, mx1, myp1
         integer, intent(in) :: npbmx, ntmax, nbmax, noff, nyp
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         integer, dimension(mx1*myp1), intent(in) :: kpic
         integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(inout) :: ihole
         integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         real, dimension(nx+2), intent(in) :: x                                ! M. Touati
         real, dimension(ny+2), intent(in) :: y                                ! M. Touati
         end subroutine
      end interface
!
      interface
         subroutine PPPORDERF2LA(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll&
     &,nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(inout) :: ppbuff
         real, dimension(idimp,nbmax), intent(inout) :: sbufl, sbufr
         integer, dimension(8,mx1*myp1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(3,mx1), intent(inout) :: ncll, nclr
         end subroutine
      end interface
!
      interface
         subroutine PPPORDER2LB(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,&
     &mcll,mclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, myp1, npbmx, ntmax
         integer, intent(in) :: nbmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1*myp1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1*myp1), intent(in) :: ppbuff
         real, dimension(idimp,nbmax), intent(in) :: rbufl, rbufr
         integer, dimension(mx1*myp1), intent(inout) :: kpic
         integer, dimension(8,mx1*myp1), intent(in) :: ncl
         integer, dimension(2,ntmax+1,mx1*myp1), intent(in) :: ihole
         integer, dimension(3,mx1), intent(in) :: mcll, mclr
         end subroutine
      end interface
!
      end module