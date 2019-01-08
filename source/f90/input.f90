module input

	implicit none
	public :: read_input_deck
	
!=====================================================================================
!                      Public variables known in the whole code 
!            (see the input deck file 'input_deck' for the definitions)
!=====================================================================================
! Simulation options :
		integer, public          :: N_threads
		integer, public          :: FDTD
		real, public             :: cfl
		logical, public          :: relativistic
		logical, public          :: moving_ions
		logical, public          :: laserpulse
		logical, public          :: plasma
		real, public             :: delta_x
		real, public             :: delta_y
      	integer, public          :: indx
      	integer, public          :: indy
		real, public             :: t_end
		character(len=3), public :: units
! Plasma properties  :
		real, public             :: atomic_weight
		real, public             :: atomic_number
		real, public             :: ionization_state
		real, public             :: param
		real, public             :: ax, ay
		real, public             :: den_me
		real, public             :: Npic
		real, public             :: ptx, pty, ptz
		real, public             :: px0, py0, pz0
! Laser pulse properties
		integer, public          :: propdir
		real, public             :: theta
        real, public             :: polardir
		integer, public          :: shape
		real, public             :: tlaunch
		real, public             :: FWHMt, FWHMs
		real, public             :: xfocal, yfocal
		real, public             :: omega0, a0
! Boundary conditions
		integer, public          :: BCx, BCy
		integer, public          :: PML_scheme
		real, public             :: n_cond_PML
		real, public             :: L_PML
! Diagnostics :
		real, public             :: Delta_t_diag		
! density
        character(len=2048), public :: density_x, density_y

	contains

	subroutine read_input_deck()
		implicit none
		Character(len=80)         :: str
		Integer                   :: i, istr
! read the file 'init_parameters'
		open (unit=1, file='input_deck',&
			& position='rewind', access='sequential',&
			& form='formatted', action='read',status='old')
		read: do
		call get_str(str)
		if (str == 'end') exit read
		if (str(1:1) /= '#') cycle read
		istr = 80
		do i = 1, len(str)
			if (str(i:i) == ' ') then
			istr = i
			exit
			end if
		end do
		select case (str(1:istr))
! Simulation options :
			case ('#N_threads')
				N_threads = get_integer(str(istr+1:))
			case ('#FDTD')
				FDTD = get_integer(str(istr+1:))
			case ('#cfl')
				cfl = get_real(str(istr+1:))
			case ('#laserpulse')  
				laserpulse = get_logical(str(istr+1:))
			case ('#plasma')
				plasma = get_logical(str(istr+1:))
			case ('#relativistic')  
				relativistic = get_logical(str(istr+1:))
			case ('#moving_ions')
				moving_ions = get_logical(str(istr+1:))
			case ('#Delta_x')
				delta_x = get_real(str(istr+1:))
			case ('#Delta_y')
				delta_y = get_real(str(istr+1:))
			case ('#indx')
				indx = get_real(str(istr+1:))
			case ('#indy')
				indy = get_real(str(istr+1:))
			case ('#t_end')
				t_end = get_real(str(istr+1:))
			case ('#units')
				units = get_char(str(istr+1:))
! Plasma properties  :
			case ('#atomic_weight')
				Atomic_weight = get_real(str(istr+1:))
			case ('#atomic_number')
				Atomic_number = get_real(str(istr+1:))
			case ('#ionization_state')
				ionization_state = get_real(str(istr+1:))
			case ('#param')
				param = get_real(str(istr+1:))
			case ('#den_me')
				den_me = get_real(str(istr+1:))
			case ('#ax')
				ax = get_real(str(istr+1:))
			case ('#ay')
				ay = get_real(str(istr+1:))
			case ('#Npic')
				Npic = get_real(str(istr+1:))
			case ('#ptx')
				ptx = get_real(str(istr+1:))
			case ('#pty')
				pty = get_real(str(istr+1:))
			case ('#ptz')
				ptz = get_real(str(istr+1:))
			case ('#px0')
				px0 = get_real(str(istr+1:))
			case ('#py0')
				py0 = get_real(str(istr+1:))
			case ('#pz0')
				pz0 = get_real(str(istr+1:))
! Laser pulse properties  :
			case ('#propdir')
				propdir = get_integer(str(istr+1:))
			case ('#theta')
				theta = get_real(str(istr+1:))
        	case ('#polardir') 
        		polardir = get_real(str(istr+1:))
			case ('#shape') 
				shape = get_integer(str(istr+1:))
			case ('#tlaunch')
				tlaunch = get_real(str(istr+1:))
			case ('#FWHMt') 
				FWHMt = get_real(str(istr+1:))
			case ('#FWHMs') 
				FWHMs = get_real(str(istr+1:))
			case ('#xfocal')
				xfocal = get_real(str(istr+1:))
		 	case ('#yfocal')
				yfocal = get_real(str(istr+1:))
			case ('#omega0')
				omega0 = get_real(str(istr+1:)) 
			case ('#a0')
				a0 = get_real(str(istr+1:))
! Boundary conditions
			case ('#BCx')
				BCx = get_integer(str(istr+1:))
			case ('#BCy')
				BCy = get_integer(str(istr+1:))
			case ('#PML_scheme')
				PML_scheme = get_integer(str(istr+1:))
			case ('#L_PML')
				L_PML = get_real(str(istr+1:))
			case ('#n_cond_PML')
				n_cond_PML = get_real(str(istr+1:))
! Diagnostics :
			case ('#Delta_t_diag')
				Delta_t_diag = get_real(str(istr+1:))
! Density :
			case ('#density_x')
				density_x = get_line(str(istr+1:))
			case ('#density_y')
				density_y = get_line(str(istr+1:))
		end select
		end do read
		close(1)
	end subroutine read_input_deck

	subroutine get_str(str)
		character(len=*), intent(inout) :: str
		read (1,'(A)',end=10) str
		str = adjustl(str)
		return
		continue ; 10 str = 'end'
	end subroutine get_str

	elemental function get_logical(str)
		character(len=*), intent(in) :: str
		logical :: get_logical
		read (str, *) get_logical
		return
	end function get_logical

	elemental function get_integer(str)
		character(len=*), intent(in) :: str
		integer :: get_integer
		read (str, *) get_integer
		return
	end function get_integer

	elemental function get_real(str)
		character(len=*), intent(in) :: str
		real :: get_real
		read (str, *) get_real
		return
	end function get_real

	elemental function get_char(str)
		character(len=*), intent(in) :: str
		Character(len=2) :: get_char
		read (str, *) get_char
		return
		end function get_char

        function get_line(str)
                character(len=*), intent(in) :: str
                Character(len=2048) :: get_line
                read (str, '(A/)') get_line
                return
                end function get_line
end module input
