!-----------------------------------------------------------------------
! Interface file for libmpfield2.f
      module libmpfield2_h
      implicit none
!
      interface
!         subroutine MPPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
!     &nyv,kxp,nyhd)                                                           ! M. Touati
         subroutine MPPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
     &nyv,kxp,nyhd,kx,ky)                                                      ! M. Touati
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         complex, dimension(nyv,kxp), intent(in)  :: q
         complex, dimension(2,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
!         subroutine MPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
!     &nyv,kxp,nyhd)                                                           ! M. Touati
          subroutine MPPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt, &
     &nyv,kxp,nyhd,kx,ky)                                                      ! M. Touati
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp
         real, intent(inout) :: we
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
!         subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)                        ! M. Touati
         subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp,kx,ky)                   ! M. Touati
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: cu
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine MPPCUPERP2(cu,nx,ny,kstrt,nyv,kxp)                        ! M. Touati
         subroutine MPPCUPERP2YEE(cu,nx,ny,kstrt,nyv,kxp,kx,ky)                ! M. Touati
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: cu
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine MIPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,  &
!     &nyhd)                                                                   ! M. Touati
         subroutine MIPPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,  &
     &nyhd,kx,ky)                                                              ! M. Touati
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         end subroutine
      end interface
!
      interface
!         subroutine MPPMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,   &
!     &kstrt,nyv,kxp,nyhd)
		 subroutine MPPMAXWEL2(exy,bxy,cu,ffc,&
     &affp,ci,dt,wf,wm,nx,ny,kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)                   ! M. Touati
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci, dt, ax, ay
         real, intent(inout) :: wf, wm
         complex, dimension(3,nyv,kxp), intent(inout) :: exy, bxy
         complex, dimension(3,nyv,kxp), intent(in)  :: cu
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         end subroutine
      end interface
!
	  interface
         subroutine MPPMAXWEL2YEE(exy,bxy,exy_corr,bxy_corr,cu,ffc,&
     &affp,ci,dt,wf,wm,nx,ny,kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)                   ! M. Touati
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci, dt, ax, ay
         real, intent(inout) :: wf, wm
         complex, dimension(3,nyv,kxp), intent(inout) :: exy, bxy
         complex, dimension(4,3), intent(inout)       :: exy_corr, bxy_corr    ! M. Touati
         complex, dimension(3,nyv,kxp), intent(in)  :: cu
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         end subroutine
      end interface
!
      interface
         subroutine MPPEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,  &
     &nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(3,nyv,kxp), intent(in) :: exy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPEMFIELD2YEE(fxy,exy,ffc,isign,kx,ky,ax,ay,nx,ny, &
     &kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(3,nyv,kxp), intent(in) :: exy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         real, intent(in)                  :: ax, ay
         real, dimension(nx+2), intent(in) :: kx
         real, dimension(ny+2), intent(in) :: ky
         end subroutine
      end interface
!
      interface
         subroutine MPPLSREMFIELD2(fxy,exy,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(in) :: fxy
         complex, dimension(3,nyv,kxp), intent(inout) :: exy
         end subroutine
      end interface
!
      interface
         subroutine PPBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx)
         implicit none
         integer, intent(in) :: nyp, nx, nxe, nypmx
         real, intent(in) :: omx, omy, omz
         real, dimension(3,nxe,nypmx), intent(inout) :: bxy
         end subroutine
      end interface
!
      interface
         subroutine PPADDVRFIELD2(a,b,c,ndim,nxe,nypmx)
         implicit none
         integer, intent(in) :: ndim, nxe, nypmx
         real, dimension(ndim,nxe,nypmx), intent(inout) :: a
         real, dimension(ndim,nxe,nypmx), intent(in) :: b, c
         end subroutine
      end interface
!
      interface
         subroutine MPPBBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,  &
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
!         subroutine MPPDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)                 ! M. Touati
         subroutine MPPDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp,kx,ky)            ! M. Touati
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: dcu
         complex, dimension(4,nyv,kxp), intent(in) :: amu
         real, dimension(nx+2), intent(in) :: kx                               ! M. Touati
         real, dimension(ny+2), intent(in) :: ky                               ! M. Touati
         end subroutine
      end interface
!
      interface
         subroutine MPPADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(inout) :: dcu
         complex, dimension(4,nyv,kxp), intent(in) :: amu
         end subroutine
      end interface
!
      interface
         subroutine MPPEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf, &
     &nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: isign, nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ax, ay, affp, wp0, ci
         real, intent(inout) :: wf
         complex, dimension(3,nyv,kxp), intent(in) :: dcu
         complex, dimension(3,nyv,kxp), intent(inout) :: exy
         complex, dimension(nyhd,kxp), intent(inout) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPOTP2(q,pot,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(nyv,kxp), intent(inout) :: pot
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPELFIELD22(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in)  :: q
         complex, dimension(2,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPELFIELD23(q,fxy,ffc,we,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(inout) :: we
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(3,nyv,kxp), intent(inout) :: fxy
         complex, dimension(nyhd,kxp), intent(inout) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(ndim,nyv,kxp), intent(in) :: f
         complex, dimension(nyv,kxp), intent(inout) :: df
         end subroutine
      end interface
!
      interface
         subroutine MPPGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, ndim, nyv, kxp
         complex, dimension(nyv,kxp), intent(in) :: df
         complex, dimension(ndim,nyv,kxp), intent(inout) :: f
         end subroutine
      end interface
!
      interface
         subroutine MPPCURLF2(f,g,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(in) :: f
         complex, dimension(3,nyv,kxp), intent(inout) :: g
         end subroutine
      end interface
!
      interface
         subroutine MPPAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp
         complex, dimension(3,nyv,kxp), intent(in) :: bxy
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         end subroutine
      end interface
!
      interface
         subroutine MPPAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,&
     &nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         complex, dimension(3,nyv,kxp), intent(in) :: bxy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPAPOTP23(cu,axy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,nyhd&
     &)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: ci
         real, intent(inout) :: wm
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: axy
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPETFIELD23(dcu,exy,ffe,affp,ci,wf,nx,ny,kstrt,nyv,&
     &kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         real, intent(in) :: affp, ci
         real, intent(inout) :: wf
         complex, dimension(3,nyv,kxp), intent(in) :: dcu
         complex, dimension(3,nyv,kxp), intent(inout) :: exy
         complex, dimension(nyhd,kxp), intent(in) :: ffe
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH2(q,qs,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(nyv,kxp), intent(in) :: q
         complex, dimension(nyv,kxp), intent(inout) :: qs
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MPPSMOOTH23(cu,cus,ffc,nx,ny,kstrt,nyv,kxp,nyhd)
         implicit none
         integer, intent(in) :: nx, ny, kstrt, nyv, kxp, nyhd
         complex, dimension(3,nyv,kxp), intent(in) :: cu
         complex, dimension(3,nyv,kxp), intent(inout) :: cus
         complex, dimension(nyhd,kxp), intent(in) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine PPRDMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,  &
     &kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, kstrt, nyv, kxp
         integer, intent(in) :: modesxpd, modesyd
         complex, dimension(nyv,kxp), intent(in) :: pot
         complex, dimension(modesyd,modesxpd), intent(inout) :: pott
         end subroutine
      end interface
!
      interface
         subroutine PPWRMODES2(pot,pott,nx,ny,modesx,modesy,kstrt,nyv,  &
     &kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, kstrt, nyv, kxp
         integer, intent(in) :: modesxpd, modesyd
         complex, dimension(nyv,kxp), intent(inout) :: pot
         complex, dimension(modesyd,modesxpd), intent(in) :: pott
         end subroutine
      end interface
!
      interface
         subroutine PPRDVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,    &
     &kstrt,nyv,kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, kstrt
         integer, intent(in) :: nyv, kxp, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp), intent(in) :: vpot
         complex, dimension(ndim,modesyd,modesxpd), intent(inout) ::     &
     &   vpott
         end subroutine
      end interface
!
      interface
         subroutine PPWRVMODES2(vpot,vpott,nx,ny,modesx,modesy,ndim,    &
     &kstrt,nyv,kxp,modesxpd,modesyd)
         implicit none
         integer, intent(in) :: nx, ny, modesx, modesy, ndim, kstrt
         integer, intent(in) :: nyv, kxp, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp), intent(inout) :: vpot
         complex, dimension(ndim,modesyd,modesxpd), intent(in) :: vpott
         end subroutine
      end interface
!
      end module
