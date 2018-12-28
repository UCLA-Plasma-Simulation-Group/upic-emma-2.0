
      module upic_m_parameters

!
! File numbers for diagnostics 
!
! N_rho, N_jx, ... = indices for diagnostic file units
!                    (each MPI node store results in a different file)
!
      integer :: N_rho = 10, N_rho_ele = 11, N_rho_ion = 12 
      integer  :: N_jx=13, N_jy=14, N_jz=15
      integer  :: N_jex=16, N_jey=17, N_jez=18
      integer  :: N_jix=19, N_jiy=20, N_jiz=21
      integer  :: N_rho_Fourier=22, N_rho_ele_Fourier=23, N_rho_ion_Fourier=24
      integer  :: N_jx_Fourier=25, N_jy_Fourier=26, N_jz_Fourier=27
      integer  :: N_histogram_ele=28, N_histogram_ion=29, N_jez_Fourier=30
      integer  :: N_jix_Fourier=31, N_jiy_Fourier=32, N_jiz_Fourier=33
      integer  :: N_rho_Fourier_arg=34, N_rho_ele_Fourier_arg=35, N_rho_ion_Fourier_arg=36
      integer  :: N_jx_Fourier_arg=37, N_jy_Fourier_arg=38, N_jz_Fourier_arg=39
      integer  :: N_jex_Fourier_arg=40, N_jey_Fourier_arg=41, N_jez_Fourier_arg=42
      integer  :: N_jix_Fourier_arg=43, N_jiy_Fourier_arg=44, N_jiz_Fourier_arg=45
      integer  :: N_Ex_Fourier=46, N_Ey_Fourier=47, N_Ez_Fourier=48
      integer  :: N_Bx_Fourier=49, N_By_Fourier=50, N_Bz_Fourier=51
      integer  :: N_Ex_Fourier_arg=52, N_Ey_Fourier_arg=53, N_Ez_Fourier_arg=54
      integer  :: N_Bx_Fourier_arg=55, N_By_Fourier_arg=56, N_Bz_Fourier_arg=57
      integer  :: N_Ex = 58, N_Ey = 59, N_Ez = 60
      integer  :: N_Bx=61, N_By = 62, N_Bz = 63
      integer  :: N_Pix=64, N_Piy = 65, N_Piz = 66
      integer  :: N_ele_px_vs_xy = 67, N_ele_py_vs_xy = 68, N_ele_pz_vs_xy = 69
      integer  :: N_ion_px_vs_xy = 70, N_ion_py_vs_xy = 71, N_ion_pz_vs_xy = 72
      integer  :: N_es_fields = 73, N_el_fields = 74, N_ma_fields = 75
      integer  :: N_kinetic_energy = 76, N_ele_kinetic_energy = 77, N_ion_kinetic_energy = 78
      integer  :: N_el_fields_pml = 79, N_ma_fields_pml = 80, N_fields_src = 81
      integer  :: N_fields_esc = 82, N_el_dumped = 83, N_ma_dumped = 84

      real, parameter :: part_factor = 2.5

      integer, parameter :: p_x_dim = 2
      integer, parameter :: p_p_dim = 3

      end module upic_m_parameters


