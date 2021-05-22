PROGRAM main
  use global
  implicit none
  real(dp) :: B_value
  real(dp), allocatable :: E_data(:)
  real(dp) :: S_data(m_data_max)

  interface
    subroutine jet_calc(B_value, verbose) 
      real(kind(1.d0)), intent(in) :: B_value
      logical, intent(in) :: verbose
    end subroutine jet_calc
  end interface

  interface
    subroutine sed_model(B_value, E_vector, S_vector) 
      use global
      real(dp), intent(in) :: B_value
      real(kind(1.d0)), intent(in) :: E_vector(:)
      !real(kind(1.d0)), intent(out) :: S_vector(:)
      real(dp), intent(out) :: S_vector(m_data_max)
    end subroutine sed_model
  end interface 

  !-------------------------------------------------------------

  B_value = 1.d-2
  call jet_calc(B_value, .true.)

  allocate( E_data(4) )
  E_data = (/1d-19, 1d-14, 1.d-10, 2.d-8/)

  call sed_model(B_value, E_data, S_data) ! Energies in cgs

  !deallocate(E_data,S_data)

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

    if (verbose) write(*,*) '\nJET_CALC started'
    if (dist .and. verbose) write(*,*) '\n Calculate NT electron distribution'
    if (dist) call diste_run 
 
    if (radsyn .and. verbose) write(*,*) 'Calculate synchrotron SED'
    if (radsyn) call radsyn_run(B_value, L_sy_tot)
    if (radsyn .and. verbose) write(*,'(A10,X,Es10.2)') 'L_sy =', L_sy_tot

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
    real(dp), intent(in) :: B_value
    real(dp), intent(in) :: E_data(:)
    real(dp), intent(out) :: S_data(m_data_max)

    call diste_run 

    call radsyn_data_run(B_value, E_data, S_data)

    write(*,*) 'E_data[eV] =', E_data/eV
    write(*,*) 'S_data[cgs] =', S_data

  end subroutine sed_model

  !===================================================================



   
    
