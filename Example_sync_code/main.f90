PROGRAM main
  implicit none
  real(kind(1.d0)) :: B, E_data(3), S_data(3)

  interface
    subroutine jet_calc(B_value, verbose) 
      real(kind(1.d0)), intent(in) :: B_value
      logical, intent(in) :: verbose
    end subroutine jet_calc
  end interface

  B = 1.d-2
  call jet_calc(B, .true.)

  E_data = (/1d-19, 1d-14, 1.d-10/)

  call sed_model(B, E_data, S_data) ! Energies in cgs

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

    !================================================================


    if (verbose) write(*,*) '\nJET_CALC started'
    if (dist .and. verbose) write(*,*) '\n Calculate NT electron distribution'
    if (dist) call diste_run 
 
    if (radsyn .and. verbose) write(*,*) 'Calculate synchrotron SED'
    if (radsyn) call radsyn_run(B_value, L_sy_tot)! Give explicit value of B
    if (radsyn .and. verbose) write(*,'(A10,X,Es10.2)') 'L_sy =', L_sy_tot

    !================================================================

    if (verbose) then
      write(*,*) 'B [uG] = ', B_value*1.d6
      write(*,*) 'E_max [TeV] = ', Emax_e/eV/1.d12
      write(*,*) '\nJET_CALC finished\n'
    endif
  end subroutine jet_calc


  !==================================================================


  subroutine sed_model(B_value, E_data, S_data)
    !===================================================================
    ! Subroutine to calculate synchrotron SED at specific energies
    ! for the given input parameters
    !===================================================================
    use global
    use controles
    use dist, only: diste_run
    use radsyn, only: radsyn_run, radsyn_data_run

    implicit none
    integer, parameter :: m_data=3 ! Number of data points
    real(dp), intent(in) :: B_value
    real(dp), intent(in) :: E_data(m_data)
    real(dp), intent(out) :: S_data(m_data)

    call diste_run 

    call radsyn_data_run(B_value, E_data, S_data)! Give explicit value of B

    !write(*,*) 'E_data[eV] =', E_data/eV
    !write(*,*) 'S_data[cgs] =', S_data

  end subroutine sed_model

  !===================================================================



   
    
