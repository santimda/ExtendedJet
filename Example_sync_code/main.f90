PROGRAM main
  use global
  use constants
  implicit none
  real(dp) :: B_value, U_value
  real(dp), allocatable :: E_data(:)
  real(dp) :: S_data(m_data_max)

  interface
    subroutine jet_calc(B_value, U_value, verbose) 
      real(kind(1.d0)), intent(in) :: B_value, U_value
      logical, intent(in) :: verbose
    end subroutine jet_calc
  end interface

  interface
    subroutine sed_model(B_value, U_value, E_vector, S_vector) 
      use global
      real(dp), intent(in) :: B_value, U_value
      real(kind(1.d0)), intent(in) :: E_vector(:)
      real(dp), intent(out) :: S_vector(m_data_max)
    end subroutine sed_model
  end interface 

  !-------------------------------------------------------------

  B_value = 1.d-2
  U_value = 1.d6*eV
  call jet_calc(B_value, U_value, .true.)

  allocate( E_data(7) )
  E_data = (/1d-19, 1d-14, 1.d-10, 2.d-8, 1.d-5, 1.d-3, 1.d0/)
  call sed_model(B_value, U_value, E_data, S_data) ! Energies in cgs
  deallocate( E_data )

END PROGRAM main

  !==================================================================

  subroutine jet_calc(B_value, U_value, verbose)
    use global
    use controles
    use dist, only: diste_run
    use radsyn, only: radsyn_run, radsyn_data_run
    use radssc, only: radssc_run, radssc_data_run

    implicit none
    real(dp), intent(in) :: B_value, U_value
    logical, intent(in) :: verbose

    allocate( N_e(mE_e), E_e(mE_e) )

    if (verbose) write(*,*) '\nJET_CALC started'
    if (dist_calc .and. verbose) write(*,*) '\n Calculate the NT electron distribution'
    if (dist_calc) call diste_run(U_value)
 
    if (radsyn_calc .and. verbose) write(*,*) 'Calculate the synchrotron SED'
    if (radsyn_calc .or. radssc_calc) call radsyn_run(B_value, E_sy, L_sy)
    if (radsyn_calc .and. verbose) write(*,'(A10,X,Es10.2)') 'L_sy =', L_sy_tot

    if (radssc_calc .and. verbose) write(*,*) 'Calculate the SSC SED'
    if (radssc_calc) call radssc_run(E_sy, L_sy, E_ic, L_ic)
    if (radssc_calc .and. verbose) write(*,'(A10,X,Es10.2)') 'L_SSC =', L_ic_tot

    if (verbose) then
      write(*,*) 'B [uG] = ', B_value*1.d6
      write(*,*) 'U [eV/cm3] = ', U_value/eV
      write(*,*) 'E_max [TeV] = ', Emax_e/eV/1.d12
      write(*,*) '\nJET_CALC finished\n'
    endif

    deallocate( N_e, E_e )

  end subroutine jet_calc


  !==================================================================

  subroutine sed_model(B_value, U_value, E_data, S_data)
    !===================================================================
    ! Subroutine to calculate synchrotron SED at specific energies
    ! for the given input parameters
    !===================================================================
    use global
    use controles
    use constants
    use dist, only: diste_run
    use radsyn, only: radsyn_run, radsyn_data_run
    use radssc, only: radssc_run, radssc_data_run

    implicit none
    real(dp), intent(in) :: B_value, U_value
    real(dp), intent(in) :: E_data(:)
    real(dp), intent(out) :: S_data(m_data_max)
    real(dp) :: S_data_syn(m_data_max), S_data_ssc(m_data_max)
    logical, parameter :: verbose=.False.

    allocate( N_e(mE_e), E_e(mE_e) )

    call diste_run(U_value)

    call radsyn_data_run(B_value, E_data, S_data_syn)
    call radssc_data_run(B_value, E_data, S_data_ssc)
    S_data = S_data_syn + S_data_ssc

    if(verbose) write(*,*) 'E_data[eV] =', E_data/eV
    if(verbose) write(*,*) 'S_data[cgs] =', S_data

    deallocate( N_e, E_e ) 

  end subroutine sed_model

  !===================================================================



   
    
