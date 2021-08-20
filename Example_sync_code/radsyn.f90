module radsyn
  use constants
  use global
  use controles
  
  implicit none
  private
  public :: radsyn_run, radsyn_data_run

contains
  
  !=============================================================================
  ! radsyn_run
  !=============================================================================
  ! This routine calculates the intrinsic synchrotron emission in the fluid frame
  !-----------------------------------------------------------------------------
  subroutine radsyn_run(B_new, E_sy, L_sy)
    !$ use omp_lib
    implicit none
    real(dp), intent(in) :: B_new
    real(dp), intent(out) :: E_sy(mE_sy), L_sy(mE_sy)
    real(dp) :: E_sy_min, E_sy_max, dE_sy, E_sy_int, gamma_e(mE_e), x
    real(dp) :: nuc, nuc_max, sum, cte_P_sy, P_sy, a, B_par
    logical :: write_ecsv = .False.

    B = B_new
    B_par = B*sqrt(2./3.) ! Assuming isotropy

    !------------------------------------------------------------------------
    !  Reading information here  

    if (.not. dist_calc ) then
      ! N_e(i) = number of e- with energy E_e(i) per unit energy
      open (14,file='dist_e.dat')
      do i = 1, mE_e
        read(14,*) E_e(i), N_e(i)
        E_e(i) = E_e(i)*eV    !   [erg]
      enddo
      close(14) 
    endif
    Eint_e = E_e(2)/E_e(1)
    gamma_e = E_e/mec2

    !------------------------------------------------------------------------
    !  Radiative calculations start here  
    
    ! Define the limits of the photon energies and create an energy vector
    nuc_max = 4.22d6*B_par*gamma_e(mE_e)**2
    E_sy_max = 10.d0*h*nuc_max  ! 
    E_sy_min = 10.d6*h          ! 10 MHz
    call vector_log(E_sy_min, E_sy_max, E_sy_int, E_sy)

    ! Calculate the L_sy at each photon energy E_sy(i) 
    cte_P_sy = 1.85*sqrt(3.)*qe**3*B_par/h/mec2 !4.39d-22*B_par/h
    L_sy_tot = 0.d0               ! Total luminosity: L_sy = int P_sy(E_sy) dE_sy
    do i = 1, mE_sy              ! E_sy(i), or nu(i)
       dE_sy = E_sy(i)*(E_sy_int-1.d0)
       sum = 0.d0
       
       ! L_sy(E_sy) = integral_E_e(min)^E_e(max) N_e(E_e) P_sy(E_e,E_sy) dE_e [erg/s/erg]
       do j = 1, mE_e             ! N_e(j), nuc(j), gamma_e(j)
          nuc = 4.22d6*B_par*gamma_e(j)**2 ! synchrotron characteristic frequency from Eq. (5.38)
          x = E_sy(i)/(h*nuc)     ! 0.1 < x = nu/nu_c < 10 --> Use approximation from Eq. (5.41)
          if (x <= 12.d0) then
            P_sy = cte_P_sy * x**(1.d0/3.d0) * exp(-x) ! Eq. (5.37)
            dE_e = E_e(j) * (Eint_e - 1.d0)
            sum = N_e(j) * P_sy * dE_e + sum
          endif
       enddo                ! N_e 
       L_sy(i) = sum         ! [erg/s/erg]
       L_sy_tot = L_sy_tot + L_sy(i)*dE_sy
       
    enddo  ! E_sy(i),nu(i)

    if (write_files) then
      open(17,file='radsyn.dat')
      do i = 1, mE_sy
        write(17,*) E_sy(i)/eV, L_sy(i) ! L(E) for the SED
      enddo
      close(17)
    endif

    ! Make an ecsv file containing a header, E, F(E), F(E)_error
    if (write_ecsv) then
      a = 4.d0*pi*d**2
      open(18,file='data2.ecsv')
      write(18, *) '# %ECSV 0.9'
      write(18, *) '# ---'
      write(18, *) '# datatype:'
      write(18, *) '# - {name: energy, unit: eV, datatype: float64}'
      write(18, *) '# - {name: flux, unit: 1 / (cm2 s erg), datatype: float64}'
      write(18, *) '# - {name: flux_error, unit: 1 / (cm2 s erg), datatype: float64}'
      write(18, *) '# meta: !!omap'
      write(18, *) '# - keywords: !!omap'
      write(18, *) '#   - cl: {value: 0.64}'
      write(18, *) 'energy flux flux_error '
      do i = 1, mE_sy, int(mE_sy/16) ! Just some values in the SED, assume errors of ~30%
        write(18,*) E_sy(i)/eV, L_sy(i)/a/E_sy(i), 0.3*L_sy(i)/a/E_sy(i) 
        !write(18,*) E_sy(i)/eV, L_sy(i)/a, 0.1*L_sy(i)/a 
      enddo
      close(17)
    endif
    
  END subroutine radsyn_run
   


  !=============================================================================
  ! radsyn_data_run
  !=============================================================================
  ! This routine calculates the intrinsic synchrotron emission in the fluid frame
  ! at the energies E_data, and returns the flux at each of those energies
  !-----------------------------------------------------------------------------
  subroutine radsyn_data_run(B_new, E_data, S_data)
    !$ use omp_lib
    implicit none
    real(dp), intent(in) :: B_new
    real(dp), intent(in) :: E_data(:)
    real(dp), intent(out) :: S_data(m_data_max)
    real(dp) :: gamma_e(mE_e), x
    real(dp) :: nuc, sum, P_sy, cte_P_sy, a, B_par
    integer :: m_data

    B = B_new
    B_par = B*sqrt(2./3.) ! Assuming isotropy

    !------------------------------------------------------------------------
    !  Reading information here  

    if (.not. dist_calc ) then
      ! N_e(i) = number of e- with energy E_e(i) per unit energy
      open (14,file='dist_e.dat')
      do i = 1, mE_e
        read(14,*) E_e(i), N_e(i)
        E_e(i) = E_e(i)*eV    !   [erg]
      enddo
      close(14) 
    endif
    Eint_e = E_e(2)/E_e(1)
    gamma_e = E_e/mec2

    !------------------------------------------------------------------------
    !  Radiative calculations start here  
    
    m_data = SIZE(E_data)
    if (m_data > m_data_max) then
      print*, 'Too many points to fit, make m_data_max larger in global.f90!'
      stop
    endif
    allocate( L_sy_data(m_data) )

    cte_P_sy = 1.85*sqrt(3.)*qe**3*B_par/h/mec2 !4.39d-22*B_par/h

    ! Calculate the L_sy at each photon energy E_sy(i)     
    !-----------------------------------------------------------------------
    !!!!!!!!!   HACER UNA PRUEBA CON OPENMP ACA   !!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, m_data           ! E_sy(i), or nu(i)
       sum = 0.d0
       
       do j = 1, mE_e             ! N_e(j), nuc(j), gamma_e(j)
          nuc = 4.22d6*B_par*gamma_e(j)**2 ! synchrotron characteristic frequency from Eq. (5.38)
          x = E_data(i)/(h*nuc)       ! 0.1 < x = nu/nu_c < 10 --> Use approximation from Eq. (5.41)
          if (x <= 12.d0) then
            P_sy = cte_P_sy * x**(1.d0/3.d0) * exp(-x) ! Eq. (5.37)
            dE_e = E_e(j) * (Eint_e - 1.d0)
            sum = N_e(j) * P_sy * dE_e + sum
          endif
       enddo                ! N_e 
       L_sy_data(i) = sum         ! [erg/s/erg]
       
    enddo  ! E_sy(i),nu(i)
    !!!!!!!!!   HACER UNA PRUEBA CON OPENMP ACA   !!!!!!!!!!!!!!!!!!!!!!!!!!
    !-----------------------------------------------------------------------

    ! Define the output vector for the first m_data positions
    S_data = 0.d0
    a = 4.d0*pi*d**2
    do i = 1, m_data
      S_data(i) = L_sy_data(i)/a
    enddo

    deallocate( L_sy_data )
    
  END subroutine radsyn_data_run
   
  
!==============================
  
END module radsyn



