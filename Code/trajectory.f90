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
    rj(1) = r0
    do i = 1, Ns-1
      z(i+1) = z(i) * dlogZ
      rj(i+1) = r_fun(z(i+1))
    enddo

    open(10, file='trajectory'//jet//'.dat')
    do i = 1, Ns
      write(10,*) z(i), rj(i)
    enddo
    close(10)

  end subroutine trajectory_calc
  !============================================================================
  ! trajectory_calc (END)
  !============================================================================


END MODULE trajectory
