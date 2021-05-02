MODULE trajectory
  use globals
  implicit none
  private
  public :: trajectory_calc

  contains

  !============================================================================
  ! trajectory_calc (START)
  !============================================================================
  ! This subroutine calculates a vertical jet trajectory. 
  ! The input is the '' (jet) or 'cj' (counter-jet). 
  !============================================================================

  subroutine trajectory_calc(jet)

    character(LEN=*), intent(in) :: jet ! 'cj' or 'jet'
    real(dp) :: dir
    allocate( z(Ns), rj(Ns) )

    ! The jet goes upwards, the CJ goes downwards (negative z)
    if (jet == '') then
      write(*,*) '\n Calculating jet trajectory'
    elseif (jet == 'cj') then
      write(*,*) '\n Calculating counter-jet trajectory'
    endif

    !System-observer vector (unitary)
    SO = (/sin(jo*pi/180d0), sin(jo*pi/180), cos(jo*pi/180d0)/)

    ! The jet goes upwards, the CJ goes downwards (negative z)
    if (jet == '') then
      dir = 1.d0
    elseif (jet == 'cj') then
      dir = -1.d0
    else
      write(*,*) 'Invalid jet type'
      stop
    endif

    z(1) = dir*z0
    do i = 1, Ns
      rj(i) = z(i)*tan(thetaj)           !Jet radius     
      z(i+1) = z(i) + dir*vj*dt
    enddo

    open(10, file='trajectory'//jet//'.dat')
    do i = 1, Ns
      write(10,*) z(i), rj(i)
    enddo
    close(10)

    deallocate( z, rj )
    return
  end subroutine trajectory_calc


  !============================================================================
  ! trajectory_calc (END)
  !============================================================================


END MODULE trajectory
