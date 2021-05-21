PROGRAM main
  implicit none
  real(kind(1.d0)) :: B

  interface
    subroutine jet_calc(B_value, verbose) 
      real(kind(1.d0)), intent(in) :: B_value
      logical, intent(in) :: verbose
    end subroutine jet_calc
  end interface

  B = 1.d-2
  call jet_calc(B, .true.)

END PROGRAM main

  !==================================================================

  subroutine jet_calc(B_value, verbose)
    use global
    use controles
    use dist, only: diste_run
    use radsyn, only: radsyn_run, radsyn_data_run

    implicit none
    real(dp), intent(in) :: B_value
    logical, intent(in) :: verbose
    real(dp) :: L_sy_tot
    !integer, parameter :: m_data=16 ! Number of data points
    integer, parameter :: m_data=3 ! Number of data points
    real(dp) :: E_data(m_data), S_data(m_data)

    !================================================================


    if (verbose) write(*,*) '\nJET_CALC started'
    if (dist .and. verbose) write(*,*) '\n Calculate NT electron distribution'
    if (dist) call diste_run 
 
    if (radsyn .and. verbose) write(*,*) 'Calculate synchrotron SED'
    if (radsyn) call radsyn_run(B_value, L_sy_tot)! Give explicit value of B
    if (radsyn .and. verbose) write(*,'(A10,X,Es10.2)') 'L_sy =', L_sy_tot

    if (radsyn_data) then

      !E_data = (/4.13d-8, 2.37d-7, 1.36d-6, 7.8d-6, 4.5d-5, 2.56d-4, 1.47d-3, 8.4d-3, 4.8d-2, 2.8d-1, 1.6d0, 9.1d0, 52.3d0, 299.7d0, 1718.7d0, 9855.2d0/)
      E_data = (/4.13d-8, 8.4d-3, 9.1d0/)

      if (verbose) write(*,*) 'Calculate synchrotron SED at specific energies'
      call radsyn_data_run(B_value, E_data, S_data)! Give explicit value of B
      if (verbose) write(*,*) 'S_data =', S_data
    endif

    !================================================================

    if (verbose) then
      write(*,*) 'B [uG] = ', B_value*1.d6
      write(*,*) 'E_max [TeV] = ', Emax_e/eV/1.d12
      write(*,*) '\nJET_CALC finished\n'
    endif
  end subroutine

  !===================================================================


   
    
