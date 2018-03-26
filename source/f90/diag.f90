module diag

! DIAG_INIT open files that will be written by each MPI node to store 
!           the simulation results in the folder /results/
!
! 
	implicit none
	contains
!
!-----------------------------------------------------------------------
      subroutine DIAG_INIT(movion, kstrt, N_rho, N_rho_ele, N_rho_ion,&
           				   N_jx, N_jy, N_jz,&  
           				   N_jex, N_jey, N_jez,&
           				   N_jix, N_jiy, N_jiz,&
           				   N_rho_Fourier, N_rho_ele_Fourier, N_rho_ion_Fourier,&
           				   N_jx_Fourier, N_jy_Fourier,N_jz_Fourier,&
           				   N_histogram_ele, N_histogram_ion,N_jez_Fourier,&
           			       N_jix_Fourier, N_jiy_Fourier,N_jiz_Fourier,&
           				   N_rho_Fourier_arg, N_rho_ele_Fourier_arg, N_rho_ion_Fourier_arg,&
				           N_jx_Fourier_arg, N_jy_Fourier_arg, N_jz_Fourier_arg,&
				           N_jex_Fourier_arg, N_jey_Fourier_arg, N_jez_Fourier_arg,&
           				   N_jix_Fourier_arg, N_jiy_Fourier_arg, N_jiz_Fourier_arg,&
           				   N_Ex_Fourier, N_Ey_Fourier,N_Ez_Fourier,&
           				   N_Bx_Fourier, N_By_Fourier, N_Bz_Fourier,&
           				   N_Ex_Fourier_arg, N_Ey_Fourier_arg, N_Ez_Fourier_arg,&
           				   N_Bx_Fourier_arg, N_By_Fourier_arg, N_Bz_Fourier_arg,&
           				   N_Ex, N_Ey, N_Ez, N_Bx, N_By, N_Bz, N_Pix, N_Piy, N_Piz,&
           				   N_ele_px_vs_xy, N_ele_py_vs_xy, N_ele_pz_vs_xy,&
           				   N_ion_px_vs_xy, N_ion_py_vs_xy, N_ion_pz_vs_xy,&
           				   N_es_fields, N_el_fields, N_ma_fields,&
           				   N_kinetic_energy, N_ele_kinetic_energy, N_ion_kinetic_energy,&
           				   N_el_fields_pml, N_ma_fields_pml, N_fields_src,&
           				   N_fields_esc, N_el_dumped, N_ma_dumped)
! This subroutine open files that will be written by each MPI node to store 
! the simulation results in the folder /results/
! input    : kstrt
! inoutput : N_Q for each quantity Q to be stored 
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      logical, intent(in)    :: movion
      integer, intent(in)    :: kstrt
      integer, intent(inout) :: N_rho, N_rho_ele, N_rho_ion
      integer, intent(inout) :: N_jx, N_jy, N_jz
      integer, intent(inout) :: N_jex, N_jey, N_jez
      integer, intent(inout) :: N_jix, N_jiy, N_jiz
      integer, intent(inout) :: N_rho_Fourier, N_rho_ele_Fourier, N_rho_ion_Fourier
      integer, intent(inout) :: N_jx_Fourier, N_jy_Fourier,N_jz_Fourier
      integer, intent(inout) :: N_histogram_ele, N_histogram_ion,N_jez_Fourier
      integer, intent(inout) :: N_jix_Fourier, N_jiy_Fourier,N_jiz_Fourier
      integer, intent(inout) :: N_rho_Fourier_arg, N_rho_ele_Fourier_arg, N_rho_ion_Fourier_arg
      integer, intent(inout) :: N_jx_Fourier_arg, N_jy_Fourier_arg, N_jz_Fourier_arg
      integer, intent(inout) :: N_jex_Fourier_arg, N_jey_Fourier_arg, N_jez_Fourier_arg
      integer, intent(inout) :: N_jix_Fourier_arg, N_jiy_Fourier_arg, N_jiz_Fourier_arg
      integer, intent(inout) :: N_Ex_Fourier, N_Ey_Fourier,N_Ez_Fourier
      integer, intent(inout) :: N_Bx_Fourier, N_By_Fourier, N_Bz_Fourier
      integer, intent(inout) :: N_Ex_Fourier_arg, N_Ey_Fourier_arg, N_Ez_Fourier_arg
      integer, intent(inout) :: N_Bx_Fourier_arg, N_By_Fourier_arg, N_Bz_Fourier_arg
      integer, intent(inout) :: N_Ex, N_Ey, N_Ez, N_Bx, N_By, N_Bz, N_Pix, N_Piy, N_Piz
      integer, intent(inout) :: N_ele_px_vs_xy, N_ele_py_vs_xy, N_ele_pz_vs_xy
      integer, intent(inout) :: N_ion_px_vs_xy, N_ion_py_vs_xy, N_ion_pz_vs_xy
      integer, intent(inout) :: N_es_fields, N_el_fields, N_ma_fields
      integer, intent(inout) :: N_kinetic_energy, N_ele_kinetic_energy, N_ion_kinetic_energy
      integer, intent(inout) :: N_el_fields_pml, N_ma_fields_pml, N_fields_src
      integer, intent(inout) :: N_fields_esc, N_el_dumped, N_ma_dumped
! local data
      character(len=20) :: str
! create a folder /results      
      call system('mkdir -p results')
! open files with the unit N_Q entitled dig/"Q""kstrt".dat where "kstrt" is the processor 
! reference number (<<MPI) and "Q" is the physical quantity to be stored 
	  write (str, *) kstrt
	  str = adjustl(str)
! histogram of initialized momentum distributions in selected spatial cells
	  N_histogram_ele = N_histogram_ele * kstrt
 	  open(unit=N_histogram_ele,file ='results/histogram_ele'//trim(str)//'.dat', &
 		   form='formatted',status='unknown')
 	  N_histogram_ion = N_histogram_ion * kstrt
 	  open(unit=N_histogram_ion,file ='results/histogram_ion'//trim(str)//'.dat', &
 		   form='formatted',status='unknown')
! densities in the real space :
	  N_rho = N_rho * kstrt
	  open(unit=N_rho,file ='results/rho'//trim(str)//'.dat', &
	       form='formatted',status='unknown')
	  N_rho_ele = N_rho_ele * kstrt
	  if (movion) open(unit=N_rho_ele,file ='results/rho_ele'//trim(str)//'.dat', &
	       form='formatted',status='unknown')
	  N_rho_ion = N_rho_ion * kstrt
	  if (movion) open(unit=N_rho_ion,file ='results/rho_ion'//trim(str)//'.dat', &
	       form='formatted',status='unknown')
! current densities in the real space
	  N_jx = N_jx * kstrt
	  open(unit=N_jx,file ='results/jx'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jy = N_jy * kstrt
	  open(unit=N_jy,file ='results/jy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jz = N_jz * kstrt
	  open(unit=N_jz,file ='results/jz'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jex = N_jex * kstrt
	  open(unit=N_jex,file ='results/jex'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jey = N_jey * kstrt
	  open(unit=N_jey,file ='results/jey'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jez = N_jez * kstrt
	  open(unit=N_jez,file ='results/jez'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jix = N_jix * kstrt
	  open(unit=N_jix,file ='results/jix'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jiy = N_jiy * kstrt
	  open(unit=N_jiy,file ='results/jiy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jiz = N_jiz * kstrt
	  open(unit=N_jiz,file ='results/jiz'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
! densities Fourier transforms modulus and arguments
	  N_rho_Fourier = N_rho_Fourier * kstrt
	  open(unit=N_rho_Fourier,file ='results/rho_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_rho_Fourier_arg = N_rho_Fourier_arg * kstrt
	  open(unit=N_rho_Fourier_arg,file ='results/rho_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_rho_ele_Fourier = N_rho_ele_Fourier * kstrt
	  !open(unit=N_rho_ele_Fourier,file ='results/rho_ele_Fourier'//trim(str)//'.dat', &
	  !    form='formatted',status='unknown')
	  N_rho_ele_Fourier_arg = N_rho_ele_Fourier_arg * kstrt
	  !open(unit=N_rho_ele_Fourier_arg,file ='results/rho_ele_Fourier_arg'//trim(str)//'.dat', &
	  !	   form='formatted',status='unknown')
	  N_rho_ion_Fourier = N_rho_ion_Fourier * kstrt
	  !open(unit=N_rho_ion_Fourier,file ='results/rho_ion_Fourier'//trim(str)//'.dat', &
	  !	   form='formatted',status='unknown')
	  N_rho_ion_Fourier_arg = N_rho_ion_Fourier_arg * kstrt
	  !open(unit=N_rho_ion_Fourier_arg,file ='results/rho_ion_Fourier_arg'//trim(str)//'.dat', &
	  !	   form='formatted',status='unknown')
! densities Fourier transforms modulus and arguments		   
	  N_jx_Fourier = N_jx_Fourier * kstrt
	  open(unit=N_jx_Fourier,file ='results/jx_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jy_Fourier = N_jy_Fourier * kstrt
	  open(unit=N_jy_Fourier,file ='results/jy_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jz_Fourier = N_jz_Fourier * kstrt
	  open(unit=N_jz_Fourier,file ='results/jz_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jx_Fourier_arg = N_jx_Fourier_arg * kstrt
	  open(unit=N_jx_Fourier_arg,file ='results/jx_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jy_Fourier_arg = N_jy_Fourier_arg * kstrt
	  open(unit=N_jy_Fourier_arg,file ='results/jy_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_jz_Fourier_arg = N_jz_Fourier_arg * kstrt
	  open(unit=N_jz_Fourier_arg,file ='results/jz_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
 	  N_jez_Fourier = N_jez_Fourier * kstrt
! 	  open(unit=N_jez_Fourier,file ='results/jez_Fourier'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jex_Fourier_arg = N_jex_Fourier_arg * kstrt
! 	  open(unit=N_jex_Fourier_arg,file ='results/jex_Fourier_arg'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jey_Fourier_arg = N_jey_Fourier_arg * kstrt
! 	  open(unit=N_jey_Fourier_arg,file ='results/jey_Fourier_arg'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jez_Fourier_arg = N_jez_Fourier_arg * kstrt
! 	  open(unit=N_jez_Fourier_arg,file ='results/jez_Fourier_arg'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jix_Fourier = N_jix_Fourier * kstrt
! 	  open(unit=N_jix_Fourier,file ='results/jix_Fourier'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jiy_Fourier = N_jiy_Fourier * kstrt
! 	  open(unit=N_jiy_Fourier,file ='results/jiy_Fourier'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jiz_Fourier = N_jiz_Fourier * kstrt
! 	  open(unit=N_jiz_Fourier,file ='results/jiz_Fourier'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jix_Fourier_arg = N_jix_Fourier_arg * kstrt
! 	  open(unit=N_jix_Fourier_arg,file ='results/jix_Fourier_arg'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jiy_Fourier_arg = N_jiy_Fourier_arg * kstrt
! 	  open(unit=N_jiy_Fourier_arg,file ='results/jiy_Fourier_arg'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
 	  N_jiz_Fourier_arg = N_jiz_Fourier_arg * kstrt
! 	  open(unit=N_jiz_Fourier_arg,file ='results/jiz_Fourier_arg'//trim(str)//'.dat', &
! 		   form='formatted',status='unknown')
! electromagnetic fields Fourier transforms modulus and arguments		   
	  N_Ex_Fourier = N_Ex_Fourier * kstrt
	  open(unit=N_Ex_Fourier,file ='results/Ex_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Ey_Fourier = N_Ey_Fourier * kstrt
	  open(unit=N_Ey_Fourier,file ='results/Ey_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Ez_Fourier = N_Ez_Fourier * kstrt
	  open(unit=N_Ez_Fourier,file ='results/Ez_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Ex_Fourier_arg = N_Ex_Fourier_arg * kstrt
	  open(unit=N_Ex_Fourier_arg,file ='results/Ex_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Ey_Fourier_arg = N_Ey_Fourier_arg * kstrt
	  open(unit=N_Ey_Fourier_arg,file ='results/Ey_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Ez_Fourier_arg = N_Ez_Fourier_arg * kstrt
	  open(unit=N_Ez_Fourier_arg,file ='results/Ez_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Bx_Fourier = N_Bx_Fourier * kstrt
	  open(unit=N_Bx_Fourier,file ='results/Bx_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_By_Fourier = N_By_Fourier * kstrt
	  open(unit=N_By_Fourier,file ='results/By_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Bz_Fourier = N_Bz_Fourier * kstrt
	  open(unit=N_Bz_Fourier,file ='results/Bz_Fourier'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Bx_Fourier_arg = N_Bx_Fourier_arg * kstrt
	  open(unit=N_Bx_Fourier_arg,file ='results/Bx_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_By_Fourier_arg = N_By_Fourier_arg * kstrt
	  open(unit=N_By_Fourier_arg,file ='results/By_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Bz_Fourier_arg = N_Bz_Fourier_arg * kstrt
	  open(unit=N_Bz_Fourier_arg,file ='results/Bz_Fourier_arg'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
! electromagnetic fields
	  N_Ex = N_Ex * kstrt
	  open(unit=N_Ex,file ='results/Ex'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Ey = N_Ey * kstrt
	  open(unit=N_Ey,file ='results/Ey'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Ez = N_Ez * kstrt
	  open(unit=N_Ez,file ='results/Ez'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Bx = N_Bx * kstrt
	  open(unit=N_Bx,file ='results/Bx'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_By = N_By * kstrt
	  open(unit=N_By,file ='results/By'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Bz = N_Bz * kstrt
	  open(unit=N_Bz,file ='results/Bz'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Pix = N_Pix * kstrt
	  open(unit=N_Pix,file ='results/Pix'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Piy = N_Piy * kstrt
	  open(unit=N_Piy,file ='results/Piy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_Piz = N_Piz * kstrt
	  open(unit=N_Piz,file ='results/Piz'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
! macro electrons phase-space
	  N_ele_px_vs_xy = N_ele_px_vs_xy * kstrt
	  open(unit=N_ele_px_vs_xy,file ='results/ele_px_vs_xy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_ele_py_vs_xy = N_ele_py_vs_xy * kstrt
	  open(unit=N_ele_py_vs_xy,file ='results/ele_py_vs_xy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_ele_pz_vs_xy = N_ele_pz_vs_xy * kstrt
	  open(unit=N_ele_pz_vs_xy,file ='results/ele_pz_vs_xy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
! macro ions phase-space
	  N_ion_px_vs_xy = N_ion_px_vs_xy * kstrt
	  open(unit=N_ion_px_vs_xy,file ='results/ion_px_vs_xy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_ion_py_vs_xy = N_ion_py_vs_xy * kstrt
	  open(unit=N_ion_py_vs_xy,file ='results/ion_py_vs_xy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
	  N_ion_pz_vs_xy = N_ion_pz_vs_xy * kstrt
	  open(unit=N_ion_pz_vs_xy,file ='results/ion_pz_vs_xy'//trim(str)//'.dat', &
		   form='formatted',status='unknown')
! energies
	  if (kstrt == 1) then
	    open(unit=N_es_fields,file ='results/Ues.dat', &
		     form='formatted',status='unknown')
		open(unit=N_el_fields,file ='results/Uel.dat', &
		     form='formatted',status='unknown')
	    open(unit=N_ma_fields,file ='results/Uma.dat', &
		     form='formatted',status='unknown')
		open(unit=N_kinetic_energy,file ='results/Uk.dat', &
		     form='formatted',status='unknown')
		open(unit=N_ele_kinetic_energy,file ='results/Uke.dat', &
		     form='formatted',status='unknown')
		open(unit=N_ion_kinetic_energy,file ='results/Uki.dat', &
		     form='formatted',status='unknown')
		open(unit=N_el_fields_pml,file ='results/Uel_pml.dat', &
		     form='formatted',status='unknown')
		open(unit=N_ma_fields_pml,file ='results/Uma_pml.dat', &
		     form='formatted',status='unknown')
		open(unit=N_fields_src,file ='results/Usrc.dat', &
		     form='formatted',status='unknown')
		open(unit=N_fields_esc,file ='results/Uesc.dat', &
		     form='formatted',status='unknown')
		open(unit=N_el_dumped,file ='results/Uel_dumped.dat', &
		     form='formatted',status='unknown')
		open(unit=N_ma_dumped,file ='results/Uma_dumped.dat', &
		     form='formatted',status='unknown')   
      end if
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine DIAG_REAL_MOMENTS(N_rho, N_jx, N_jy , N_jz, nyp, nx, nxe, nypmx,&
                                   de, phtime, yp, x, q, cu, tdiag)
! This subroutine stores the charge and current densities (real space) in the 
! corresponding files
! input    : N_rho, N_jx, N_jy , N_jz, nyp, nx, nxe, nypmx, yp, x, qe, cue
! inoutput : tdiag
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      integer, intent(in)                      :: N_rho, N_jx, N_jy , N_jz
      integer, intent(in)                      :: nyp, nx, nxe, nypmx
      real, intent(in)                         :: de, phtime
      real, dimension(nypmx), intent(in)       :: yp
      real, dimension(nxe), intent(in)         :: x
      real, dimension(nxe,nypmx), intent(in)   :: q
      real, dimension(3,nxe,nypmx), intent(in) :: cu
      real, intent(inout)                      :: tdiag
! local data
      integer                                  :: i, j
      integer, dimension(4)                    :: itime
      double precision                         :: dtime
      real                                     :: zero
!
      zero = 1.e-15
      call dtimer(dtime,itime,-1)
      do j=1,nyp,1
      	do i=1,nx,1
      		if (abs(q(i,j)) .gt. zero) then 
      			write (N_rho,'(4E23.15)') phtime, yp(j), x(i), q(i,j)
      		else
      			write (N_rho,'(4E23.15)') phtime, yp(j), x(i), zero
      		end if
      		if (abs(cu(1,i,j)) .gt. zero) then 
	    		write (N_jx,'(4E23.15)')  phtime, yp(j), x(i), cu(1,i,j)
	    	else
      			write (N_jx,'(4E23.15)')  phtime, yp(j), x(i), zero
      		end if
	    	if (abs(cu(2,i,j)) .gt. zero) then 
	    		write (N_jy,'(4E23.15)')  phtime, yp(j), x(i), cu(2,i,j)
	    	else
      			write (N_jy,'(4E23.15)')  phtime, yp(j), x(i), zero
      		end if
	    	if (abs(cu(3,i,j)) .gt. zero) then 
	    		write (N_jz,'(4E23.15)')  phtime, yp(j), x(i), cu(3,i,j)
	    	else
      			write (N_jz,'(4E23.15)')  phtime, yp(j), x(i), zero
      		end if
      	end do
      end do
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine DIAG_FOURIER_MOMENTS(N_rho_Fourier, N_rho_Fourier_arg,&
                                      N_jx_Fourier, N_jx_Fourier_arg,&
                                      N_jy_Fourier, N_jy_Fourier_arg,&
                                      N_jz_Fourier, N_jz_Fourier_arg,&
                                      modesyd, modesxpd, kxp, nye,&
                                      de, phtime, kpx, ky, qf, cuf, tdiag)
! This subroutine stores the charge and current densities (Fourier space) in the 
! corresponding files
! input    : N_rho_Fourier, N_rho_Fourier_arg, N_jx_Fourier, N_jx_Fourier_arg,
!            N_jy_Fourier, N_jy_Fourier_arg, N_jz_Fourier, N_jz_Fourier_arg,
!            modesxpd, modesyd, kxp, nye, phtime, kpx, ky, qf, cuf
! inoutput : tdiag
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      integer, intent(in)                                :: N_rho_Fourier, N_rho_Fourier_arg
      integer, intent(in)                                :: N_jx_Fourier, N_jx_Fourier_arg
      integer, intent(in)                                :: N_jy_Fourier, N_jy_Fourier_arg
      integer, intent(in)                                :: N_jz_Fourier, N_jz_Fourier_arg
      integer, intent(in)                                :: modesxpd, modesyd, kxp, nye
      real, intent(in)                                   :: de, phtime
      real, dimension(kxp), intent(in)                   :: kpx
      real, dimension(nye), intent(in)                   :: ky
      complex, dimension(modesyd,modesxpd), intent(in)   :: qf
      complex, dimension(3,modesyd,modesxpd), intent(in) :: cuf
      real, intent(inout)                                :: tdiag
! local data
      integer                                            :: i, j
      integer, dimension(4)                              :: itime
      double precision                                   :: dtime
      real                                               :: zero
      call dtimer(dtime,itime,-1)     
      zero = 1.e-15
      do i=1,modesxpd,1
	    do j=1+int(modesyd/2),1,-1
	    	if (abs(qf(2*j-1,i)) .gt. zero) then
	    		write (N_rho_Fourier,'(4E23.15)') phtime, kpx(i),&
       										      -ky(j), abs(qf(2*j-1,i))
       		else
       			write (N_rho_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(cuf(1,2*j-1,i)) .gt. zero) then
       	    	write (N_jx_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), abs(cuf(1,2*j-1,i))
       		else
       			write (N_jx_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(cuf(2,2*j-1,i)) .gt. zero) then
       			write (N_jy_Fourier,'(4E23.15)') phtime, kpx(i),&
       											 -ky(j), abs(cuf(2,2*j-1,i))
       		else
       			write (N_jy_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(cuf(3,2*j-1,i)) .gt. zero) then
	       		write (N_jz_Fourier,'(4E23.15)') phtime, kpx(i),&
    	   										 -ky(j), abs(cuf(3,2*j-1,i))
    	   	else
       			write (N_jz_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		write (N_rho_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		       -ky(j), atan2(aimag(qf(2*j-1,i)),real(qf(2*j-1,i)))
       		write (N_jx_Fourier_arg,'(4E23.15)')  phtime, kpx(i),&
       		       -ky(j), atan2(aimag(cuf(1,2*j-1,i)),real(cuf(1,2*j-1,i)))
       		write (N_jy_Fourier_arg,'(4E23.15)')  phtime, kpx(i),&
       		       -ky(j), atan2(aimag(cuf(2,2*j-1,i)),real(cuf(2,2*j-1,i)))
       		write (N_jz_Fourier_arg,'(4E23.15)')  phtime, kpx(i),&
       		       -ky(j), atan2(aimag(cuf(3,2*j-1,i)),real(cuf(3,2*j-1,i)))
       	end do
       	do j=2,1+int(modesyd/2),1
       		if (abs(qf(2*j,i)) .gt. zero) then
	    		write (N_rho_Fourier,'(4E23.15)') phtime, kpx(i),&
       		    								 ky(j), abs(qf(2*j,i))
       		else
       			write (N_rho_Fourier,'(4E23.15)') phtime, kpx(i),&
       		    								 ky(j), zero
       		end if
       		if (abs(cuf(1,2*j,i)) .gt. zero) then
       	    	write (N_jx_Fourier,'(4E23.15)')  phtime, kpx(i),&
       											 ky(j), abs(cuf(1,2*j,i))
       		else
       			write (N_jx_Fourier,'(4E23.15)')  phtime, kpx(i),&
       		    								 ky(j), zero
       		end if
       		if (abs(cuf(2,2*j,i)) .gt. zero) then
       	    	write (N_jy_Fourier,'(4E23.15)')  phtime, kpx(i),&
       										     ky(j), abs(cuf(2,2*j,i))
       		else
       			write (N_jy_Fourier,'(4E23.15)')  phtime, kpx(i),&
       		    								 ky(j), zero
       		end if
       		if (abs(cuf(3,2*j,i)) .gt. zero) then
       	    	write (N_jz_Fourier,'(4E23.15)')  phtime, kpx(i),&
       										     ky(j), abs(cuf(3,2*j,i))
       		else
       			write (N_jz_Fourier,'(4E23.15)')  phtime, kpx(i),&
       		    								 ky(j), zero
       		end if
       		write (N_rho_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		       ky(j), atan2(aimag(qf(2*j,i)),real(qf(2*j,i)))
       		write (N_jx_Fourier_arg,'(4E23.15)')  phtime, kpx(i),&
       		       ky(j), atan2(aimag(cuf(1,2*j,i)),real(cuf(1,2*j,i)))
       		write (N_jy_Fourier_arg,'(4E23.15)')  phtime, kpx(i),&
       		       ky(j), atan2(aimag(cuf(2,2*j,i)),real(cuf(2,2*j,i)))
       		write (N_jz_Fourier_arg,'(4E23.15)')  phtime, kpx(i),&
       		       ky(j), atan2(aimag(cuf(3,2*j,i)),real(cuf(3,2*j,i)))
       	end do
      end do
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!------------------------------------------------------------------------------
	  subroutine UPDATE_POYNTING(ci,affp,dt,L_PML,delta,x,yp,fxyze,bxyze,&
	                             BCx,BCy,nx,ny,nyp,nxe,nypmx,pixyze,&
	                             wpesc,wpescx,wpescy,tdiag)
		implicit none
		integer, intent(in)                       :: BCx,BCy, nx, ny, nyp, nxe,nypmx
		real, intent(in)                          :: ci, affp, dt, L_PML
		real, dimension(2), intent(in)            :: delta
		real, dimension(nxe), intent(in)          :: x
		real, dimension(nypmx), intent(in)        :: yp
		real, dimension(3,nxe,nypmx), intent(in)  :: fxyze, bxyze
		real, dimension(3,nxe,nypmx), intent(out) :: pixyze
		real, intent(inout)                       :: wpesc, wpescx, wpescy, tdiag
! local
		integer               :: i, j
      	integer, dimension(4) :: itime
      	double precision      :: dtime
		real :: x1, x2, y1, y2, anorm, anormx, anormy
!
		call dtimer(dtime,itime,-1)
		x1 = L_PML
		x2 = real(nx)*delta(1)-L_PML
		y1 = L_PML
		y2 = real(ny)*delta(2)-L_PML
      	anorm = 1. / (ci*affp)
      	anormx = delta(1)*dt
      	anormy = delta(2)*dt
      	!$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:wpescx,wpescy)
      	do j=1,nyp,1
      		do i=1,nx,1
      			pixyze(1,i,j) = anorm * ( (fxyze(2,i,j)*bxyze(3,i,j)) - (fxyze(3,i,j)*bxyze(2,i,j)) )
      			pixyze(2,i,j) = anorm * ( (fxyze(3,i,j)*bxyze(1,i,j)) - (fxyze(1,i,j)*bxyze(3,i,j)) )
      			pixyze(3,i,j) = anorm * ( (fxyze(1,i,j)*bxyze(2,i,j)) - (fxyze(2,i,j)*bxyze(1,i,j)) )
      			if ( (BCx == 0) .and. (BCy == 1) )then
      				if      (x(i) == x1) then
      					wpescx = wpescx - ( anormy * pixyze(1,i,j) )
      				else if (x(i) == x2) then
      					wpescx = wpescx + ( anormy * pixyze(1,i,j) )
      				end if
      			end if
      			if ( (BCx == 1) .and. (BCy == 0) )then
      				if      (yp(j) == y1) then
      					wpescy = wpescy - ( anormx * pixyze(2,i,j) )
      				else if (yp(j) == y2) then
      					wpescy = wpescy + ( anormx * pixyze(2,i,j) )
      				end  if   			
      			end if
!      			if ( (BCx == 0) .and. (BCy == 0) )then
      				if      ( (x(i) == x1)  .and. (yp(j).ge.y1) .and. (yp(j).le.y2) ) then
      					wpescx = wpescx - ( anormy * pixyze(1,i,j) )
      				else if ( (x(i) == x2)  .and. (yp(j).ge.y1) .and. (yp(j).le.y2) ) then
      					wpescx = wpescx + ( anormy * pixyze(1,i,j) )
      				else if ( (yp(j) == y1) .and. (x(i).ge.x1)  .and. (x(i).le.x2)  ) then
      					wpescy = wpescy - ( anormx * pixyze(2,i,j) )
      				else if ( (yp(j) == y2) .and. (x(i).ge.x1)  .and. (x(i).le.x2)  ) then
      					wpescy = wpescy + ( anormx * pixyze(2,i,j) )
 !     				end if    			
      			end if
      		end do
      	end do
      	!$OMP END PARALLEL DO
      	wpesc = wpescx + wpescy
      	call dtimer(dtime,itime,1)
      	tdiag = tdiag + real(dtime)
	  end subroutine 
!
!-----------------------------------------------------------------------
      subroutine DIAG_REAL_FIELD(N_Ex, N_Ey, N_Ez, N_Bx, N_By, N_Bz,&
                             nxe, nypmx, nyp, nx,&
                             de, phtime, yp, x, fxyze, bxyze, tdiag)
! This subroutine stores the electromagnetic field (real space) in the 
! corresponding files
! input    : N_Ex, N_Ey, N_Ez, N_Bx, N_By, N_Bz,
!            nxe, nypmx, nyp, nx, phtime, yp, x, fxyze, bxyze
! inoutput : tdiag
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      integer, intent(in)                      :: N_Ex, N_Ey, N_Ez, N_Bx, N_By, N_Bz
      integer, intent(in)                      :: nxe, nypmx, nyp, nx
      real, intent(in)                         :: de, phtime
      real, dimension(nypmx), intent(in)       :: yp
      real, dimension(nxe), intent(in)         :: x
      real, dimension(3,nxe,nypmx), intent(in) :: fxyze, bxyze
      real, intent(inout)                      :: tdiag
! local data
      integer                                  :: i, j
      integer, dimension(4)                    :: itime
      double precision                         :: dtime
      real                                     :: zero
      call dtimer(dtime,itime,-1)     
      zero = 1.e-15
      do j=1,nyp
      	do i=1,nx
      		if (abs(fxyze(1,i,j)) .gt. zero) then
		    	write (N_Ex,'(4E23.15)') phtime, yp(j), x(i), fxyze(1,i,j)
		    else
		    	write (N_Ex,'(4E23.15)') phtime, yp(j), x(i), zero
		    end if
	    	if (abs(fxyze(2,i,j)) .gt. zero) then
	    		write (N_Ey,'(4E23.15)') phtime, yp(j), x(i), fxyze(2,i,j)
	    	else
		    	write (N_Ey,'(4E23.15)') phtime, yp(j), x(i), zero
		    end if
	    	if (abs(fxyze(3,i,j)) .gt. zero) then
	    		write (N_Ez,'(4E23.15)') phtime, yp(j), x(i), fxyze(3,i,j)
	    	else
		    	write (N_Ez,'(4E23.15)') phtime, yp(j), x(i), zero
		    end if
	    	if (abs(bxyze(1,i,j)) .gt. zero) then
	    		write (N_Bx,'(4E23.15)') phtime, yp(j), x(i), bxyze(1,i,j)
	    	else
	    		write (N_Bx,'(4E23.15)') phtime, yp(j), x(i), zero
	    	end if
	    	if (abs(bxyze(2,i,j)) .gt. zero) then
	    		write (N_By,'(4E23.15)') phtime, yp(j), x(i), bxyze(2,i,j)
	    	else
	    		write (N_By,'(4E23.15)') phtime, yp(j), x(i), zero
	    	end if
	    	if (abs(bxyze(3,i,j)) .gt. zero) then
	    		write (N_Bz,'(4E23.15)') phtime, yp(j), x(i), bxyze(3,i,j)
	    	else
	    		write (N_Bz,'(4E23.15)') phtime, yp(j), x(i), zero
	    	end if
    	end do
      end do
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine DIAG_POYNTING(N_Pix,N_Piy,N_Piz,&
                          	   nxe, nypmx, nyp, nx,&
                          	   de, phtime, yp, x, pixyze,tdiag)
! This subroutine stores the Poynting vector (real space) in the 
! corresponding files
! input    : N_Pix, N_Piy, N_Piz
!            nxe, nypmx, nyp, nx, phtime, yp, x, fxyze, bxyze
! inoutput : tdiag
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      integer, intent(in)                      :: N_Pix, N_Piy, N_Piz
      integer, intent(in)                      :: nxe, nypmx, nyp, nx
      real, intent(in)                         :: de, phtime
      real, dimension(nypmx), intent(in)       :: yp
      real, dimension(nxe), intent(in)         :: x
      real, dimension(3,nxe,nypmx), intent(in) :: pixyze
      real, intent(inout)                      :: tdiag
! local data
      integer                                  :: i, j
      integer, dimension(4)                    :: itime
      double precision                         :: dtime
      real                                     :: zero
      call dtimer(dtime,itime,-1)     
      zero = 1.e-15
      do j=1,nyp
      	do i=1,nx
      		if (abs(pixyze(1,i,j)) .gt. zero) then
		    	write (N_Pix,'(4E23.15)') phtime, yp(j), x(i), pixyze(1,i,j)
		    else
		    	write (N_Pix,'(4E23.15)') phtime, yp(j), x(i), zero
		    end if
	    	if (abs(pixyze(2,i,j)) .gt. zero) then
	    		write (N_Piy,'(4E23.15)') phtime, yp(j), x(i), pixyze(2,i,j)
	    	else
		    	write (N_Piy,'(4E23.15)') phtime, yp(j), x(i), zero
		    end if
	    	if (abs(pixyze(3,i,j)) .gt. zero) then
	    		write (N_Piz,'(4E23.15)') phtime, yp(j), x(i), pixyze(3,i,j)
	    	else
		    	write (N_Piz,'(4E23.15)') phtime, yp(j), x(i), zero
		    end if
    	end do
      end do
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------     
      subroutine DIAG_FOURIER_FIELDS(N_Ex_Fourier, N_Ey_Fourier, N_Ez_Fourier,&
                                N_Bx_Fourier, N_By_Fourier, N_Bz_Fourier,&
                                N_Ex_Fourier_arg, N_Ey_Fourier_arg, N_Ez_Fourier_arg,&
                                N_Bx_Fourier_arg, N_By_Fourier_arg, N_Bz_Fourier_arg,&
                                modesxpd, modesyd, kxp, nye,&
                                de, phtime, kpx, ky, exyzf, bxyzf, tdiag)
! This subroutine stores the electromagnetic field (Fourier space) in the 
! corresponding files
! input    : N_Ex_Fourier, N_Ey_Fourier, N_Ez_Fourier,
!            N_Bx_Fourier, N_By_Fourier, N_Bz_Fourier,
!            N_Ex_Fourier_arg, N_Ey_Fourier_arg, N_Ez_Fourier_arg
!            N_Bx_Fourier_arg, N_By_Fourier_arg, N_Bz_Fourier_arg
!            modesxpd, modesyd, kxp, nye, de, phtime, kpx, ky, exyzf, bxyzf
! inoutput : tdiag
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      integer, intent(in)                                :: N_Ex_Fourier, N_Ey_Fourier
      integer, intent(in)                                :: N_Ez_Fourier, N_Bx_Fourier
      integer, intent(in)                                :: N_By_Fourier, N_Bz_Fourier
      integer, intent(in)                                :: N_Ex_Fourier_arg
      integer, intent(in)                                :: N_Ey_Fourier_arg
      integer, intent(in)                                :: N_Ez_Fourier_arg
      integer, intent(in)                                :: N_Bx_Fourier_arg
      integer, intent(in)                                :: N_By_Fourier_arg
      integer, intent(in)                                :: N_Bz_Fourier_arg
      integer, intent(in)                                :: modesxpd, modesyd, kxp, nye
      real, intent(in)                                   :: de, phtime
      real, dimension(kxp), intent(in)                   :: kpx
      real, dimension(nye), intent(in)                   :: ky
      complex, dimension(3,modesyd,modesxpd), intent(in) :: exyzf, bxyzf
      real, intent(inout)                                :: tdiag
! local data
      integer                                            :: i, j
      integer, dimension(4)                              :: itime
      double precision                                   :: dtime
      real                                               :: zero
      call dtimer(dtime,itime,-1)     
      zero = 1.e-15
      do i=1,modesxpd,1
	    do j=1+int(modesyd/2),1,-1
	    	if (abs(exyzf(1,2*j-1,i)) .gt. zero ) then
       	    	write (N_Ex_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), abs(exyzf(1,2*j-1,i))
       		else
       			write (N_Ex_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(exyzf(2,2*j-1,i)) .gt. zero ) then
       	    	write (N_Ey_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), abs(exyzf(2,2*j-1,i))
       		else
       			write (N_Ey_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(exyzf(3,2*j-1,i)) .gt. zero ) then
       	    	write (N_Ez_Fourier,'(4E23.15)') phtime, kpx(i),&
       										 -ky(j), abs(exyzf(3,2*j-1,i))
       		else
       			write (N_Ez_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(bxyzf(1,2*j-1,i)) .gt. zero ) then
       	    	write (N_Bx_Fourier,'(4E23.15)') phtime, kpx(i),&
       										 -ky(j), abs(bxyzf(1,2*j-1,i))
       		else
       			write (N_Bx_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(bxyzf(2,2*j-1,i)) .gt. zero ) then
       	    	write (N_By_Fourier,'(4E23.15)') phtime, kpx(i),&
       										 -ky(j), abs(bxyzf(2,2*j-1,i))
       		else
       			write (N_By_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(bxyzf(3,2*j-1,i)) .gt. zero ) then
       	    	write (N_Bz_Fourier,'(4E23.15)') phtime, kpx(i),&
       										 -ky(j), abs(bxyzf(3,2*j-1,i))
       		else
       			write (N_Bz_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		write (N_Ex_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		-ky(j), atan2(aimag(exyzf(1,2*j-1,i)),real(exyzf(1,2*j-1,i)))
       		write (N_Ey_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		-ky(j), atan2(aimag(exyzf(2,2*j-1,i)),real(exyzf(2,2*j-1,i)))
       		write (N_Ez_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		-ky(j), atan2(aimag(exyzf(3,2*j-1,i)),real(exyzf(3,2*j-1,i)))
       		write (N_Bx_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		-ky(j), atan2(aimag(bxyzf(1,2*j-1,i)),real(bxyzf(1,2*j-1,i)))
       		write (N_By_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		-ky(j), atan2(aimag(bxyzf(2,2*j-1,i)),real(bxyzf(2,2*j-1,i)))
       		write (N_Bz_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		-ky(j), atan2(aimag(bxyzf(3,2*j-1,i)),real(bxyzf(3,2*j-1,i)))
       	end do
       	do j=2,1+int(modesyd/2),1
       	    if (abs(exyzf(1,2*j,i)) .gt. zero ) then
       	    	write (N_Ex_Fourier,'(4E23.15)') phtime, kpx(i),&
       										   ky(j), abs(exyzf(1,2*j,i))
       		else
       			write (N_Ex_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(exyzf(2,2*j,i)) .gt. zero ) then
       	    	write (N_Ey_Fourier,'(4E23.15)') phtime, kpx(i),&
       										   ky(j), abs(exyzf(2,2*j,i))
       		else
       			write (N_Ey_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(exyzf(3,2*j,i)) .gt. zero ) then
       	    	write (N_Ez_Fourier,'(4E23.15)') phtime, kpx(i),&
       										   ky(j), abs(exyzf(3,2*j,i))
       		else
       			write (N_Ez_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(bxyzf(1,2*j,i)) .gt. zero ) then
       	    	write (N_Bx_Fourier,'(4E23.15)') phtime, kpx(i),&
       										   ky(j), abs(bxyzf(1,2*j,i))
       		else
       			write (N_Bx_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(bxyzf(2,2*j,i)) .gt. zero ) then
       	    	write (N_By_Fourier,'(4E23.15)') phtime, kpx(i),&
       										   ky(j), abs(bxyzf(2,2*j,i))
       		else
       			write (N_By_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		if (abs(bxyzf(3,2*j,i)) .gt. zero ) then
       	    	write (N_Bz_Fourier,'(4E23.15)') phtime, kpx(i),&
       										   ky(j), abs(bxyzf(3,2*j,i))
       		else
       			write (N_Bz_Fourier,'(4E23.15)') phtime, kpx(i),&
       										     -ky(j), zero
       		end if
       		write (N_Ex_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		ky(j), atan2(aimag(exyzf(1,2*j,i)),real(exyzf(1,2*j,i)))
       		write (N_Ey_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		ky(j), atan2(aimag(exyzf(2,2*j,i)),real(exyzf(2,2*j,i)))
       		write (N_Ez_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		ky(j), atan2(aimag(exyzf(3,2*j,i)),real(exyzf(3,2*j,i)))
       		write (N_Bx_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		ky(j), atan2(aimag(bxyzf(1,2*j,i)),real(bxyzf(1,2*j,i)))
       		write (N_By_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		ky(j), atan2(aimag(bxyzf(2,2*j,i)),real(bxyzf(2,2*j,i)))
       		write (N_Bz_Fourier_arg,'(4E23.15)') phtime, kpx(i),&
       		ky(j), atan2(aimag(bxyzf(3,2*j,i)),real(bxyzf(3,2*j,i)))
       	end do
      end do
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------     
      subroutine DIAG_PART(N_px_vs_xy, N_py_vs_xy, N_pz_vs_xy, &
                           npp, idimp, nppmx0, mxyp1, Ndiag, de, phtime, ppart, tdiag)
! This subroutine stores in the corresponding files the velocity components in the 
! x and y-direction of macro-particles versus their position
! input    : N_px_vs_xy, N_py_vs_xy, N_pz_vs_xy
!            npp, idimp, nppmx0, mxyp1, Ndiag, phtime, ppart
! inoutput : tdiag
! written by Michael J. Touati, UCLA, 2016, March
      implicit none
      integer, intent(in)                             :: N_px_vs_xy, N_py_vs_xy, N_pz_vs_xy
      integer, intent(in)                             :: npp, idimp, nppmx0, mxyp1, Ndiag
      real, intent(in)                                :: de, phtime
      real, dimension(idimp,nppmx0,mxyp1), intent(in) :: ppart
      real, intent(inout)                             :: tdiag
! local
      integer, dimension(4)                           :: itime
      double precision                                :: dtime
      integer                                         :: l, tile
      call dtimer(dtime,itime,-1)
      do tile=1,mxyp1
        do l=1,nppmx0,Ndiag
        	if ((ppart(1,l,tile).ne.0.).and.(ppart(2,l,tile).ne.0.)&
        	    .and.(ppart(3,l,tile).ne.0.).and.(ppart(4,l,tile).ne.0.)&
        	    .and.(ppart(5,l,tile).ne.0.)) then
     	  		write (N_px_vs_xy,'(4E23.15)') phtime, ppart(1,l,tile), ppart(2,l,tile),&
     	        	                           ppart(3,l,tile)
     	  		write (N_py_vs_xy,'(4E23.15)') phtime, ppart(1,l,tile), ppart(2,l,tile),&
     	        	                           ppart(4,l,tile)
     	  		write (N_pz_vs_xy,'(4E23.15)') phtime, ppart(1,l,tile), ppart(2,l,tile),&
     	                                       ppart(5,l,tile)
     	    end if
        end do
      end do
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!----------------------------------------------------------------------- 
      subroutine DIAG_HISTOGRAM(ci,vtx,vty,vtz,N_histogram,npp,idimp,npmax,nxe,nye,&
                                part,x,y,delta, tdiag)
      	implicit none
      	real, intent(in)                         :: ci,vtx,vty,vtz
      	integer, intent(in)                      :: N_histogram
      	integer, intent(in)                      :: npp, idimp, npmax, nxe, nye
      	real, dimension(idimp,npmax), intent(in) :: part
      	real, dimension(nxe), intent(in)         :: x
      	real, dimension(nye), intent(in)         :: y
      	real, dimension(2), intent(in)           :: delta
      	real, intent(inout) :: tdiag
! local
		integer :: N, l, i, j, k, N_part
		real(8) :: mu, vT, v_max, w_max, dv
		real, dimension(:,:), allocatable :: part_cell,histogram
		logical :: in_cell
		integer, dimension(4)                    :: itime
      	double precision                         :: dtime
      	call dtimer(dtime,itime,-1)
		vT = sqrt( ( (vtx**2.)+(vty**2.)+(vtz**2.) ) / 3. )
		mu = ci*ci * vT*vT
		mu = 1. / mu
		v_max = 4.e0 * sqrt(2.e0) * sqrt( 1.e0 + sqrt( 1.e0 + (mu**2.) ) ) / ( mu * ci )
		! cell (3,3) of this MPI tile for example :
		i = 3
		j = 3
		N_part = 0
		do l=1,npp,1
			in_cell = (part(1,l) .gt. x(i)) 
			in_cell = ( (in_cell) .and. (part(1,l) .le. x(i+1)) )
			in_cell = ( (in_cell) .and. (part(2,l) .gt. y(j)) )
			in_cell = ( (in_cell) .and. (part(2,l) .le. y(j+1)) )
			if (in_cell) then
				N_part = N_part + 1
			end if
		end do
		N  = 1000
		allocate(part_cell(3,N_part),histogram(2,N))
		part_cell = 0.
		k    = 0
		do l=1,npp,1
			in_cell = ( (part(1,l) .gt. x(i)) .and. (part(1,l) .le. x(i+1)) )
			in_cell = ( in_cell .and. (part(2,l) .gt. y(j)) .and. (part(2,l) .le. y(j+1)) )
			if (in_cell) then
				k = k + 1
				part_cell(1:3,k) = part(3:5,l)
			end if
		end do
		histogram(1:2,1:N) = 0.
		dv        = v_max / N
		do i = 1,N-1
			histogram(1,i+1) = histogram(1,i) + dv
			do l=1,N_part
				w_max = sqrt((part_cell(1,l)**2.)+(part_cell(2,l)**2.)+(part_cell(3,l)**2.))
				if ( (w_max .gt. histogram(1,i)) .and. (w_max .le. histogram(1,i+1)) ) then
					histogram(2,i+1) = histogram(2,i+1) + 1.
				end if
			end do
			histogram(2,i+1) = histogram(2,i+1) / (N*dv)
			write (N_histogram,'(2E23.15)') histogram(1,i),histogram(2,i)
		end do
		write (N_histogram,'(2E23.15)') histogram(1,N),histogram(2,N)
		deallocate(part_cell,histogram)
		call dtimer(dtime,itime,1)
		tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------      
      subroutine DIAG_ENERGY(BCx,BCy,movion,phtime,tdiag, &
                             N_es_fields, N_el_fields, N_ma_fields, &
                             N_kinetic_energy,N_ele_kinetic_energy, N_ion_kinetic_energy, &
                             N_el_fields_pml, N_ma_fields_pml, N_fields_src, &
           				     N_fields_esc, N_el_dumped, N_ma_dumped, &
                             we, wf, wm, wk, wke, wki, wsrc, wfsim, wmsim, &
                             wesc,wfpml,wmpml,wfdump,wmdump)
! This subroutine stores in the corresponding files the velocity components in the 
! x and y-direction of macro-particles versus their position
! input    : N_es_fields, N_el_fields, N_ma_fields, N_kinetic_energy,
!            N_ele_kinetic_energy, N_ion_kinetic_energy, N_el_fields_pml,
!            N_ma_fields_pml, N_fields_src,
!            phtime, we, wf, wm, wk, wke, wki, wth, wthe, wthi
! inoutput : tdiag
! written by Michael J. Touati, UCLA, 2016, March
		  implicit none
		  integer, intent(in) :: BCx, BCy
		  logical, intent(in) :: movion
		  integer, intent(in) :: N_es_fields, N_el_fields, N_ma_fields
		  integer, intent(in) :: N_kinetic_energy, N_ele_kinetic_energy, N_ion_kinetic_energy
		  integer, intent(in) :: N_el_fields_pml, N_ma_fields_pml, N_fields_src
		  integer, intent(in) :: N_fields_esc, N_el_dumped, N_ma_dumped
		  real, intent(in)    :: phtime, we, wf, wm, wk, wke, wki
		  real, intent(in)    :: wsrc, wfsim, wmsim, wesc
		  real, intent(in)    :: wfpml, wmpml, wfdump, wmdump
		  real, intent(inout) :: tdiag
! local data
		  integer, dimension(4)                    :: itime
		  double precision                         :: dtime
		  call dtimer(dtime,itime,-1)
		  write (N_es_fields,'(2E23.15)')                      phtime, we
		  if ( (BCx == 1) .and. (BCy == 1) ) then
			  write (N_el_fields,'(2E23.15)')                  phtime, wf
			  write (N_ma_fields,'(2E23.15)')                  phtime, wm
			  write (N_fields_src,'(2E23.15)')                 phtime, wsrc
			  write (N_fields_esc,'(2E23.15)')                 phtime, wesc
		  else
			  write (N_el_fields,'(2E23.15)')                  phtime, wfsim
			  write (N_ma_fields,'(2E23.15)')                  phtime, wmsim
			  write (N_fields_src,'(2E23.15)')                 phtime, wsrc
			  write (N_fields_esc,'(2E23.15)')                 phtime, wesc
			  write (N_el_fields_pml,'(2E23.15)')              phtime, wfpml
			  write (N_ma_fields_pml,'(2E23.15)')              phtime, wmpml
			  write (N_el_dumped,'(2E23.15)')                  phtime, wfdump
			  write (N_ma_dumped,'(2E23.15)')                  phtime, wmdump
		  end if
		  write (N_kinetic_energy,'(2E23.15)')                 phtime, wk
		  if (movion) write (N_ele_kinetic_energy,'(2E23.15)') phtime, wke
		  if (movion) write (N_ion_kinetic_energy,'(2E23.15)') phtime, wki
		  call dtimer(dtime,itime,1)
		  tdiag = tdiag + real(dtime)
      end subroutine

end module
