!-----------------------------------------------------------------------
! Interface file for libpml2.f
!
module libmpfieldpml2_h
	implicit none
!
	interface
    	subroutine MPPMAXWELB2SPECTRALPML(exy,bxy,&
	                              cu,ffc,affp,ci,dt,wm,nx,ny,&
                                  kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
        	implicit none
         	integer, intent(in)                          :: nx, ny, kstrt, nyv, kxp, nyhd
         	real, intent(in)                             :: affp, ci, dt, ax, ay
         	real, intent(inout)                          :: wm
         	complex, dimension(4,nyv,kxp), intent(inout) :: exy, bxy
         	complex, dimension(3,nyv,kxp), intent(in)    :: cu
         	complex, dimension(nyhd,kxp), intent(in)     :: ffc
         	real, dimension(nx+2), intent(in)            :: kx
         	real, dimension(ny+2), intent(in)            :: ky
         end subroutine
      end interface
!-----------------------------------------------------------------------
!
	interface 
    	subroutine MPPMAXWELE2SPECTRALPML(exy,bxy,&
	                              cu,ffc,affp,ci,dt,wf,nx,ny,&
                                  kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
        	implicit none
         	integer, intent(in)                          :: nx, ny, kstrt, nyv, kxp, nyhd
         	real, intent(in)                             :: affp, ci, dt, ax, ay
         	real, intent(inout)                          :: wf
         	complex, dimension(4,nyv,kxp), intent(inout) :: exy, bxy
         	complex, dimension(3,nyv,kxp), intent(in)    :: cu
         	complex, dimension(nyhd,kxp), intent(in)     :: ffc
         	real, dimension(nx+2), intent(in)            :: kx
         	real, dimension(ny+2), intent(in)            :: ky
         end subroutine
    end interface	
!
!-----------------------------------------------------------------------
!
	interface
    	subroutine MPPMAXWELB2YEEPML(exy,bxy,exy_corr,bxy_corr,&
	                              cu,ffc,affp,ci,dt,wm,nx,ny,&
                                  kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
        	implicit none
         	integer, intent(in)                          :: nx, ny, kstrt, nyv, kxp, nyhd
         	real, intent(in)                             :: affp, ci, dt, ax, ay
         	real, intent(inout)                          :: wm
         	complex, dimension(4,nyv,kxp), intent(inout) :: exy, bxy
         	complex, dimension(4,4), intent(inout)       :: exy_corr, bxy_corr
         	complex, dimension(3,nyv,kxp), intent(in)    :: cu
         	complex, dimension(nyhd,kxp), intent(in)     :: ffc
         	real, dimension(nx+2), intent(in)            :: kx
         	real, dimension(ny+2), intent(in)            :: ky
         end subroutine
      end interface
!-----------------------------------------------------------------------
!
	interface 
    	subroutine MPPMAXWELE2YEEPML(exy,bxy,exy_corr,bxy_corr,&
	                              cu,ffc,affp,ci,dt,wf,nx,ny,&
                                  kstrt,nyv,kxp,nyhd,kx,ky,ax,ay)
        	implicit none
         	integer, intent(in)                          :: nx, ny, kstrt, nyv, kxp, nyhd
         	real, intent(in)                             :: affp, ci, dt, ax, ay
         	real, intent(inout)                          :: wf
         	complex, dimension(4,nyv,kxp), intent(inout) :: exy, bxy
         	complex, dimension(4,4), intent(inout)       :: exy_corr, bxy_corr
         	complex, dimension(3,nyv,kxp), intent(in)    :: cu
         	complex, dimension(nyhd,kxp), intent(in)     :: ffc
         	real, dimension(nx+2), intent(in)            :: kx
         	real, dimension(ny+2), intent(in)            :: ky
         end subroutine
    end interface	
!
!-----------------------------------------------------------------------
end module