MODULE trajectory
  use globals
  use constants
  implicit none
  private
  public :: trajectory_calc

  contains

  !============================================================================
  ! trajectory_calc (START)
  !============================================================================
  subroutine trajectory_calc(jet)
  !============================================================================
  ! This subroutine calculates a vertical jet trajectory. 
  ! The input is the '' (jet) or 'CJ' (counter-jet). 
  !============================================================================
    character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet
    real(dp) :: dir
    integer :: i

    ! The jet goes upwards, the CJ goes downwards (negative z)
    if (jet == '') then
      dir = 1.d0
    elseif (jet == 'CJ') then
      dir = -1.d0
    else
      print *, 'Invalid jet type'
      stop
    endif

    z(1) = dir*z0               ! jet height
    rj_out(1) = rout_0
    rj_in(1) = rin_0
    do i = 2, Ns
      z(i) = z(i-1) * dlogz
      rj_out(i) = r_fun_out(z(i))
      rj_in(i) = r_fun_in(z(i))
    end do

    open(10, file=trim(opath('trajectory'//jet//'.dat')))
    do i = 1, Ns
      write(10,*) z(i), rj_out(i), rj_in(i)
    enddo
    close(10)

  end subroutine trajectory_calc
  !============================================================================
  ! trajectory_calc (END)
  !============================================================================


END MODULE trajectory
