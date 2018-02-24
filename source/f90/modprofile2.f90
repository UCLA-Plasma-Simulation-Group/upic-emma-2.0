!-----------------------------------------------------------------------
module libprofile_h
use m_rvm
implicit none

public :: getdensity1d
contains

subroutine getdensity1d(expr, varname, edges, dx, nx, npx, sr, x)
    character (len=*), intent(in) :: expr, varname
    double precision, intent(in) :: edges
    double precision, dimension(:), intent(inout) :: x
    double precision, intent(in) :: dx
    integer, intent(in) :: npx, nx, sr

    ! local
    integer :: i, ierr
    integer :: part_nx, gppe
    double precision, parameter :: epslon = 1D-14
    double precision, dimension(:), allocatable :: pcharge
    double precision :: qlow, den_ltot, pvol, qx, dxp_2

    type(t_parser) :: dens

    call setup(dens, expr, (/varname/), ierr)

    ! get half distance between sampling grid points
    dxp_2 = 0.5/sr
    !enddo

    ! find total number sampling grid points per cell
    gppe = nx * sr + 1

    ! initialize temp buffers
    allocate( pcharge(gppe) )
    den_ltot = 0.0
    ! Get position of sampling points for the full box
    ! Inject particles

    ! evaluate the density expression
    do i = 1, gppe
        ! overwrite in one direction
        pcharge(i) = eval(dens, (/ edges + (2 * i - 2) * dxp_2 * dx /))
    enddo
    ! numerically integrate the density expression
    den_ltot = 0
    do i = 1, gppe-1
        den_ltot = den_ltot + 0.5 * (pcharge(i)+pcharge(i+1))
    enddo
    qx = den_ltot / npx

    part_nx = 0
    qlow = qx * 0.5
    i = 1
    do while (.true.)
        if (qlow < qx) then
            qlow = qlow + 0.5 * (pcharge(i) + pcharge(i+1))
            if ( i > gppe) exit
            i = i + 1
        else
            part_nx = part_nx + 1
            qlow = qlow - qx
            ! here pvol is used as tmp storing the density gradient.
            pvol = pcharge(i)-pcharge(i-1) 
            if (pvol <= epslon) then  ! piece-wise constant
                x(part_nx) = dxp_2 * 2.0 * (i - qlow/pcharge(i) - 1)
            else  ! assuming linear gradient within this piece
                x(part_nx) = dxp_2 * 2.0 * (i - (sqrt(pcharge(i)**2 + &
                        & 2*qlow*pvol) - pcharge(i))/pvol - 1)
            endif
        endif
    enddo
    
    call delete( dens )
    deallocate( pcharge )


end subroutine getdensity1d

end module

