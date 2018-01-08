! OpenMP utility library
! written by Viktor K. Decyk, UCLA
! copyright 2016, regents of the university of california
! update: january 12, 2016
      module omplib
      use omp_lib
      implicit none
!
! common data for parallel processing
! nthreads = number of threads used
      integer :: nthreads
      save
!
      private
      public :: INIT_OMP, SETNTHSIZE, GETNTHSIZE
!
      contains
!
!-----------------------------------------------------------------------
!      subroutine INIT_OMP(nth)                                 ! M. Touati
      subroutine INIT_OMP(kstrt,nvp, nth)                       ! M. Touati
! initialize openmp library
! use nth threads if nth > 0; otherwise, use the number found
!             and if nccpus/nvp >= nth
      implicit none
!      integer, intent(in) :: nth                               ! M. Touati
	  integer, intent(in) :: kstrt, nvp, nth                    ! M. Touati
! local data
      integer :: ncpus
! determine how many processors are available
      ncpus = omp_get_num_procs()
!      write (*,*) 'number of cpus found = ', ncpus             ! M. Touati
      nthreads = omp_get_max_threads()
      nthreads = ncpus / nvp                                    ! M. Touati
!      write (*,*) 'maximum number of threads = ', nthreads     ! M. Touati
!	  if (nth > 0) nthreads = nth                               ! M. Touati
      if ( (nth > 0) .and. (nth .lt. nthreads) ) nthreads = nth ! M. Touati
      call omp_set_num_threads(nthreads)
      write (*,*) 'MPI node', kstrt, 'is using ', nthreads, 'OpenMP thread(s)'
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine SETNTHSIZE(nth)
! set number of threads
      implicit none
      integer, intent(in) :: nth
      if (nth > 0) nthreads = nth
      call omp_set_num_threads(nthreads)
      end subroutine
!
!-----------------------------------------------------------------------
      function GETNTHSIZE()
! get number of threads
      implicit none
      integer :: GETNTHSIZE
      GETNTHSIZE = nthreads
      end function
!
      end module omplib
