module antenna
!
	implicit none
	contains
!
	subroutine LaunchLaser(fxyze,bxyze,wl,x,yp,affp,delta,ci,dt,phtime,L_PML,theta,polardir,tlaunch,FWHMt,&
	                      FWHMs,xfocal,yfocal,omega0,a0,BCx,BCy,nx,nyp,nxe,nypmx,ntime,propdir,shape,laserpulse)
		implicit none
		integer,intent(in)                          :: BCx, BCy, nx, nyp, nxe, nypmx
		integer, intent(in)                         :: ntime, propdir, shape
		real, intent(in)                            :: theta, polardir
		real, intent(in)                            :: tlaunch, FWHMt
		real, intent(in)                            :: FWHMs, xfocal, yfocal
		real, intent(in)                            :: omega0, a0
		real, intent(in)                            :: affp,ci, dt, phtime, L_PML
		real, dimension(nxe), intent(in)            :: x
		real, dimension(nypmx), intent(in)          :: yp
		real, dimension(2), intent(in)              :: delta
		real, dimension(3,nxe,nypmx), intent(inout) :: fxyze, bxyze
		real, intent(inout)                         :: wl
		logical, intent(inout)                      :: laserpulse
! locals
		integer                      :: ntlas, ifant
		integer                      :: antlength
!
		!ntlas = 1
		!ifant = ntime - ntlas*(ntime/ntlas) ! use antenna every ntlas steps
		!if (ifant == 0) then
			if (shape == 0) then 
! plane antenna
				call LaunchPlaneAntenna_fst(fxyze,bxyze,wl,x,yp,affp,delta,polardir,tlaunch,&
				                        FWHMt,omega0,a0,ci,dt,phtime,L_PML,&
				                        ntime,propdir,BCx,BCy,nx,nyp,nxe,nypmx,laserpulse)
			else 
! Gaussian antenna
				! call LaunchGaussianAntenna(fxyze,bxyze,x,yp,delta,polardir,tlaunch,&
! 				                           FWHMt,omega0,a0,ci,dt,phtime,L_PML,&
! 				                           ntime,propdir,BCx,BCy,nxe,nypmx)
			endif
		!endif
	end subroutine LaunchLaser
!
	subroutine LaunchPlaneAntenna(fxyze,bxyze,wl,x,yp,affp,delta,polardir,tlaunch,&
				                  FWHMt,omega0,a0,ci,dt,phtime,L_PML,&
				                  ntime,propdir,BCx,BCy,nx,nyp,nxe,nypmx,laserpulse)
		implicit none
		integer,intent(in)                          :: BCx, BCy, nx, nyp, nxe, nypmx
		integer, intent(in)                         :: ntime, propdir
		real, intent(in)                            :: polardir
		real, intent(in)                            :: tlaunch, FWHMt
		real, intent(in)                            :: omega0, a0
		real, intent(in)                            :: affp, ci, dt, phtime, L_PML
		real, dimension(nxe), intent(in)            :: x
		real, dimension(nypmx), intent(in)          :: yp
		real, dimension(2), intent(in)              :: delta
		real, dimension(3,nxe,nypmx), intent(inout) :: fxyze, bxyze
		real, intent(inout)                         :: wl
		logical, intent(inout)                      :: laserpulse
! locals
		integer  :: ii, jj, dim, dim1, dim2, dim3, Npt
		logical  :: antenna
		real     :: a, pi, Ft, xi, x0, x1, x2, omega, k, phi0, phi, k0, alpha, beta, anorm, fl1, bl1, fl2, bl2
!
	    k0 = omega0 * ci
	    omega = sin(omega0*dt)/dt
		k = omega * ci
		pi = 3.14159265359
		anorm = 0.5 *delta(1) * delta(2) / affp
		if (propdir == 1) then
			if (BCx == 0) then
				x1 = L_PML
				x2 = L_PML + (2.*pi/k0)
! 				x2 = real(nx)*delta(1)
			else
! 			    x1 = delta(1)
! 				x2 = delta(1) + 2.*pi/k0
! 				x1 =   delta(1)
! 				x2 = 2*delta(1)
				x1 = floor(2.*pi/(k0*delta(1)))*delta(1)
				x2 = x1 + floor(2.*pi/(k0*delta(1)))*delta(1)
				Npt = 1 + floor((x2-x1)/delta(1))
				a = a0 / Npt
			end if
			dim1 = 2
			dim2 = 3
			dim3 = 1
		else if (propdir == 2) then
			if (BCy == 0) then
				x1 = L_PML
				x2 = L_PML + (2.*pi/k0)
! 				x2 = real(nx)*delta(1)
			else
			    x1 = delta(2)
				x2 = delta(2) + 2.*pi/k0
			end if
			dim1 = 3
			dim2 = 1
			dim3 = 2
		end if
		wl = 0.
		!$OMP  PARALLEL DO PRIVATE(ii,jj,antenna,xi,phi0,Ft,phi,alpha,beta) &
		!$OMP& REDUCTION(+:wl)
		do ii = 1,nx
			do jj = 1,nyp
				if (propdir == 1) then
					antenna = ( x(ii)  >= x1) .and. ( x(ii)  <= x2)
					xi = x(ii)
				else if (propdir == 2) then
					antenna = (yp(jj)  >= x1) .and. (yp(jj)  <= x2)
					xi = yp(jj)
				else
					antenna = .false.
				end if
				if (antenna) then
					phi0 = 0.5*pi*((phtime-tlaunch)-((xi-x1)*ci)-FWHMt)/FWHMt
					Ft  = cos(phi0) ** 2.
! 					Ft = exp(-4.*log(2.)*(((phtime-tlaunch)-((xi-x1)*ci)-2.5*FWHMt)/FWHMt)**2.)
!                     Ft = 1.
! 					phi = omega0 * phtime - k0 * xi
					phi = omega0*(phtime-tlaunch) - (k0*(xi-x1))
					alpha = sqrt(     abs(polardir) )
					beta  = sqrt(1. - abs(polardir) )
					if (polardir >= 0.) then
						fxyze(dim1,ii,jj) =   a0 * omega0 * alpha * sin(phi) * Ft 
						fxyze(dim2,ii,jj) = - a0 * omega0 * beta  * cos(phi) * Ft
						fxyze(dim3,ii,jj) =   0.
						bxyze(dim1,ii,jj) =   a0 * omega0 * beta  * cos(phi) * Ft
						bxyze(dim2,ii,jj) =   a0 * omega0 * alpha * sin(phi) * Ft
						bxyze(dim3,ii,jj) =   0.
					else
						print*,'not yet supported'
						stop
! 						fxyze(dim1,ii,jj) =   a0 * omega0 * alpha * sin(phi) * Ft 
! 						bxyze(dim2,ii,jj) = - a0 * omega0 * beta  * cos(phi) * Ft
					endif
						wl = wl + anorm * (  (fxyze(dim1,ii,jj)**2.) & 
						                   + (fxyze(dim2,ii,jj)**2.) & 
						                   + (bxyze(dim1,ii,jj)**2.) & 
						                   + (bxyze(dim2,ii,jj)**2.)  )
				else
					fxyze(dim1,ii,jj) = 0.
					fxyze(dim2,ii,jj) = 0.
					fxyze(dim3,ii,jj) = 0.
					bxyze(dim1,ii,jj) = 0.
					bxyze(dim2,ii,jj) = 0.
					bxyze(dim3,ii,jj) = 0.
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!$OMP PARALLEL DO PRIVATE(jj)
		do jj = 1,nyp,1
			do dim=1,3,1
				fxyze(dim,nx+1,jj) = fxyze(dim,1,jj)
				bxyze(dim,nx+1,jj) = bxyze(dim,1,jj)
				fxyze(dim,nx+2,jj) = fxyze(dim,1,jj)
				bxyze(dim,nx+2,jj) = bxyze(dim,1,jj)
			end do
		end do
		!$OMP END PARALLEL DO
		if (phtime >= (tlaunch + 2.*FWHMt)) then
			laserpulse = .false.
			wl = 0.
		endif
	end subroutine LaunchPlaneAntenna


	subroutine LaunchPlaneAntenna_fst(fxyze,bxyze,wl,x,yp,affp,delta,polardir,tlaunch,&
				                  FWHMt,omega0,a0,ci,dt,phtime,L_PML,&
				                  ntime,propdir,BCx,BCy,nx,nyp,nxe,nypmx,laserpulse)
		implicit none
		integer,intent(in)                          :: BCx, BCy, nx, nyp, nxe, nypmx
		integer, intent(in)                         :: ntime, propdir
		real, intent(in)                            :: polardir
		real, intent(in)                            :: tlaunch, FWHMt
		real, intent(in)                            :: omega0, a0
		real, intent(in)                            :: affp, ci, dt, phtime, L_PML
		real, dimension(nxe), intent(in)            :: x
		real, dimension(nypmx), intent(in)          :: yp
		real, dimension(2), intent(in)              :: delta
		real, dimension(3,nxe,nypmx), intent(inout) :: fxyze, bxyze
		real, intent(inout)                         :: wl
		logical, intent(inout)                      :: laserpulse
! locals
		integer  :: ii, jj, dim, dim1, dim2, dim3, Npt
		logical  :: antenna
		real     :: a, pi, Ft, xi, x0, x1, x2, omega, k, phi0, phi, k0, alpha, beta, anorm, fl1, bl1, fl2, bl2
		integer  :: ix, iy
		real     :: phase, cos_phase,sin_phase
!
	    k0 = omega0 * ci
	    omega = sin(omega0*dt)/dt
		k = omega * ci
		pi = 3.14159265359
		anorm = 0.5 *delta(1) * delta(2) / affp
		
		alpha = sqrt(     abs(polardir) )
		beta  = sqrt(1. - abs(polardir) )
        phase = omega0*(phtime-tlaunch)
        cos_phase = cos(phase)
        sin_phase = sin(phase)
		if (propdir == 1) then
			if (BCx == 0) then
				ix = L_PML/delta(1)+3
				
			else
                ix = 2
			end if
			dim1 = 2
			dim2 = 3
			dim3 = 1
		    do jj = 1,nyp
		    	fxyze(dim1,ix,jj) =    fxyze(dim1,ix,jj) + a0 * omega0 * omega0 * alpha * cos_phase
			fxyze(dim2,ix,jj) =    fxyze(dim2,ix,jj) + a0 * omega0 * omega0 * beta  * cos_phase 
			! fxyze(dim3,ii,jj) =   0.
		    end do
		else if (propdir == 2) then
			if (BCy == 0) then
				
				iy=L_PML/delta(2)+1
! 				x2 = real(nx)*delta(1)
			else
				iy = 2
			     
			end if
			dim1 = 3
			dim2 = 1
			dim3 = 2
			do ii=1,nx
			    fxyze(dim1,ii,iy) =    a0 * omega0 * alpha * sin(phi) 
				fxyze(dim2,ii,iy) =    a0 * omega0 * beta  * sin(phi)
			end do
		end if

		do jj = 1,nyp,1
			do dim=1,3,1
				fxyze(dim,nx+1,jj) = fxyze(dim,1,jj)
				bxyze(dim,nx+1,jj) = bxyze(dim,1,jj)
				fxyze(dim,nx+2,jj) = fxyze(dim,1,jj)
				bxyze(dim,nx+2,jj) = bxyze(dim,1,jj)
			end do
		end do

		if (phtime >= (tlaunch + 2.*FWHMt)) then
			laserpulse = .false.
			wl = 0.
		endif
	end subroutine LaunchPlaneAntenna_fst
!
! 	subroutine LaunchGaussianAntenna(fxyze,bxyze,lasxstart,lasxend,lask,lase0, &
! 	                                &inorder,ci,laspol,ntime,dt,launchlas,rayl,&
! 	                                &gamma,lasw0,lasfocus,antlength)
! 		implicit none
! 		real, dimension(:,:,:,:), pointer :: fxyze,bxyze
! 		real :: ci, lask, lase0, dt,rayl,gamma,lasw0
! 		integer :: laspol, lasxstart,lasxend, lasfocus,antlength
! 		integer :: inorder, ntime, launchlas,LINEAR
! ! local data
! 		integer :: xstart,xend,ystart,yend, ant_length,xoffset,yoffset
! 		real :: laserk,omegaenv,a0,beta,xx,xx2,yy,y0,wz,curv,RL,ttovci
! 		integer :: ii, jj, yboxsize,s1,s2,s3
! 		real :: per_env, lon_env, tt, xfocus,xcenter, yj
! 		real :: pi = 3.1415926536
! ! common block for parallel processing
! 		integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! 		parameter(lstat=10)
! 		common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble,      &
! 		&lworld
! ! locals
! 		integer istatus, nvp, idproc, np, ioff, id, nrec0, i, j, k, l
! 		integer ierr
! 		dimension istatus(lstat)
! 		call MPI_COMM_RANK(lworld,idproc,ierr)
! 		call MPI_COMM_SIZE(lworld,nvp,ierr)
! 
! 		s1 = size(fxyze,1); s2 = size(fxyze,2); s3 = size(fxyze,3);
! 
! 		yend = s3 - 2; ystart = 2; yboxsize = s3-3
! 		xoffset = 1; yoffset = 1
! 		if (inorder==LINEAR) then
! 			yend = s3 - 1; ystart = 1; yboxsize = s3-1
! 			xoffset = 0; yoffset = 0
! 		endif
! 		if (nvp == 1) then
! 			y0 = yboxsize/2 + 0.5
! 		else
! 			y0 = nvp/2*yboxsize + 0.5
! 		endif
! 		if (idproc == 1) then
! 			print *, ' - Launching Gaussian Antenna.'
! 		endif
! 
! 		tt = real(ntime) * dt; ant_length = antlength
! 		beta = sqrt(1.0 - 1.0/gamma/gamma)
! 		ttovci = tt/ci; 
! 		! for ntime == 0, length of antenna is ant_length
! 		! for ntime > 0, length of antenna is ant_length/2
! 		if (ntime == 0) then
! 			if (ttovci - int(ttovci) > 0.999) then
! 				xend = lasxend+xoffset - ( int(ttovci)+1 )
! 			else 
! 				xend = lasxend+xoffset - int(ttovci)
! 			endif
! 			xstart = xend - ant_length
! 		else
! 			if (ttovci - int(ttovci) > 0.999) then
! 				xend = lasxend+xoffset - ( int(ttovci)+1 ) - ant_length/2
! 			else 
! 				xend = lasxend+xoffset - int(ttovci) - ant_length/2
! 			endif
! 			xstart = xend - ant_length/2
! 		endif
! 		xcenter = real(lasxstart+lasxend) * 0.5 + xoffset + ttovci
! 		laserk = lask; omegaenv = pi / real(lasxend-lasxstart); a0 = lase0; RL = rayl
! 		xfocus = real(lasfocus+xoffset) - ttovci * beta
! !
! 		do ii = xstart,xend
! 			do jj = ystart, yend
! 				xx = real(ii)-xfocus; xx2 = real(ii)-xcenter
! ! yj is the y-coordinate in global coordinate without guard cell
! 				yj = real(idproc) * real(yboxsize) + real(jj - yoffset)
! 				yy = (y0 - yj)
! 				wz = lasw0 * sqrt(1+xx*xx/RL/RL); 
! 				if ( abs(real(ii)-xfocus) .le. 0.00001) then
! 					curv = 0.0
! 				else
! 					curv = 0.5 * yy*yy/ (xx+RL*RL/xx) * (1+beta)
! 				endif
! 				lon_env = cos(omegaenv*xx2) ** 2
! 				per_env = sqrt(lasw0/wz)*exp(-yy*yy/wz/wz) * &
! 				& cos(laserk*xx2+laserk*curv-atan(xx/RL))
! 				if (laspol == 2) then ! e2 polarized
! 					fxyze(2,ii,jj,1) = a0 * per_env * lon_env
! 					bxyze(3,ii,jj,1) = ci * a0 * per_env * lon_env
! 				else if (laspol == 3) then ! e3 polarized
! 					fxyze(3,ii,jj,1) = a0 * per_env * lon_env
! 					bxyze(2,ii,jj,1) = -ci * a0 * per_env * lon_env
! 				endif
! 			enddo
! 		enddo
! 		if ((xcenter-real(xstart))*2 .ge. real(lasxend-lasxstart)) then
! 			launchlas = 0 ! stop using antenna
! 		endif
! 	end subroutine LaunchGaussianAntenna
	
end module
