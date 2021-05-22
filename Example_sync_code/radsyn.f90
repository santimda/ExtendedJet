module radsyn

  use global
  
  implicit none
  private
  public :: radsyn_run, radsyn_data_run

contains
  
  !=============================================================================
  ! radsyn_run
  !=============================================================================
  ! This routine calculates the intrinsic synchrotron emission in the fluid frame
  !-----------------------------------------------------------------------------
  subroutine radsyn_run(B_new, L_sy_tot)
    use global
    !$ use omp_lib
    implicit none
    ! mEe -> electrons, m_nu -> emitted photons
    real(dp), intent(in) :: B_new
    real(dp), intent(out) :: L_sy_tot
    real(dp) :: eps(m_nu), dEps, eps_int, gamma_e(mE_e), x
    real(dp) :: nuc, nuc_max, nu_min, nu_max, sum, sum_sy, cte_Psy, a, B_par
    allocate( N_e(mE_e), E_e(mE_e) )
    allocate( nu(m_nu), Lsy(m_nu) )

    B = B_new
    B_par = B*sqrt(2./3.) ! Assuming isotropy

    !------------------------------------------------------------------------
    !  Reading information here  

    ! N_e(i) = number of e- with energy E_e(i) per unit energy
    open (14,file='dist_e.dat')
    do i = 1, mE_e
       read(14,*) E_e(i), N_e(i)
       E_e(i) = E_e(i)*eV    !   [erg]
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma_e = E_e/mec2

    !------------------------------------------------------------------------
    !  Radiative calculations start here  
    
    ! Define nu_min and nu_max:
    nuc_max = 4.22d6*B_par*gamma_e(mE_e)**2
    nu_max = 10.d0*nuc_max ! 
    nu_min = 10.d6          ! 10 MHz
    ! Define a vector with the photon energies epsilon: eps(i) = h*nu(i)
    call vector_log(h*nu_min, h*nu_max, eps_int, eps)
    nu = eps/h             ! [nu] = Hz 

    cte_Psy = 1.85*sqrt(3.)*qe**3*B_par/h/mec2 !4.39d-22*B_par/h
    sum_sy = 0.d0 ! To estimate the total luminosity: L_sy = int P_sy(eps) dEps
    ! We calculate the luminosity at each value of the photon energy eps(i) between eps_min and eps_max
    do i = 1, m_nu           ! eps(i), or nu(i)
       dEps = eps(i)*(eps_int-1.d0)
       sum = 0.d0
       
       ! P(eps) = integral_E_e(min)^E_e(max) N_e(E_e) P(E_e,eps) dE_e [erg/s/erg]
       ! Define: N_e(E_e)*P(E_e,eps)*dE_e = dP_eps , so that P(eps) = integral_E_e dP_eps
       do j = 1, mE_e             ! N_e(j), nuc(j), gamma_e(j)
          nuc = 4.22d6*B_par*gamma_e(j)**2 ! synchrotron characteristic frequency from Eq. (5.38)
          x = nu(i)/nuc       ! 0.1 < x = nu/nu_c < 10 --> Use approximation from Eq. (5.41)
          if (x <= 12.d0) then
             ! Calculate P_eps [erg/s/Hz] and transform to erg/s/erg dividing by h
             P_eps = cte_Psy * x**(1.d0/3.d0) * exp(-x) ! Eq. (5.37)
             dE_e = E_e(j)*(Eint_e - 1.d0)
             dP_eps = N_e(j)*P_eps*dE_e
             sum = dP_Eps + sum
          endif
       enddo                ! N_e 
       Lsy(i) = sum         ! [erg/s/erg]
       sum_sy = sum_sy + Lsy(i)*deps
       
    enddo  ! eps(i),nu(i)
    L_sy_tot = sum_sy

    open(17,file='radsyn.dat')
    a = 4.d0*pi*d**2
    do i = 1, m_nu
       write(17,*) eps(i)/eV, Lsy(i)/a ! F(E) for the SED
    enddo
    close(17)

    ! Make an ecsv file containing a header, E, F(E), F(E)_error
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
    do i = 1, m_nu, int(m_nu/16) ! Just some values in the SED, assume errors of ~30%
       write(18,*) eps(i)/eV, Lsy(i)/a/eps(i), 0.3*Lsy(i)/a/eps(i) 
       !write(18,*) eps(i)/eV, Lsy(i)/a, 0.1*Lsy(i)/a 
    enddo
    close(17)

    deallocate( nu,Lsy )
    deallocate( N_e,E_e )
    
  END subroutine radsyn_run
   
  
!==============================

  !=============================================================================
  ! radsyn_data_run
  !=============================================================================
  ! This routine calculates the intrinsic synchrotron emission in the fluid frame
  ! at the energies E_data, and returns the flux at each of those energies
  !-----------------------------------------------------------------------------
  subroutine radsyn_data_run(B_new, E_data, S_data)
    use global
    !$ use omp_lib
    implicit none
    ! mEe -> electrons, m_nu -> emitted photons
    real(dp), intent(in) :: B_new
    real(dp), intent(in) :: E_data(:)
    real(dp) :: gamma_e(mE_e), x
    real(dp) :: nuc, sum, cte_Psy, a, B_par
    real(dp), intent(out) :: S_data(m_data_max)
    integer :: m_data
    allocate( N_e(mE_e), E_e(mE_e) )

    B = B_new
    B_par = B*sqrt(2./3.) ! Assuming isotropy

    !------------------------------------------------------------------------
    !  Reading information here  

    ! N_e(i) = number of e- with energy E_e(i) per unit energy
    open (14,file='dist_e.dat')
    do i = 1, mE_e
       read(14,*) E_e(i), N_e(i)
       E_e(i) = E_e(i)*eV    !   [erg]
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma_e = E_e/mec2

    !------------------------------------------------------------------------
    !  Radiative calculations start here  
    
    m_data = SIZE(E_data)
    if (m_data > m_data_max) then
      print*, 'Too many points to fit, make m_data_max larger in global.f90!'
      stop
    endif
    allocate( nu(m_data), Lsy(m_data) )
    nu = E_data/h             ! [nu] = Hz 

    cte_Psy = 1.85*sqrt(3.)*qe**3*B_par/h/mec2 !4.39d-22*B_par/h
    ! Calculate the luminosity at each value of the photon energy E_data(i)
    do i = 1, m_data           ! eps(i), or nu(i)
       sum = 0.d0
       
       do j = 1, mE_e             ! N_e(j), nuc(j), gamma_e(j)
          nuc = 4.22d6*B_par*gamma_e(j)**2 ! synchrotron characteristic frequency from Eq. (5.38)
          x = nu(i)/nuc       ! 0.1 < x = nu/nu_c < 10 --> Use approximation from Eq. (5.41)
          if (x <= 12.d0) then
             ! Calculate P_eps [erg/s/Hz] and transform to erg/s/erg dividing by h
             P_eps = cte_Psy * x**(1.d0/3.d0) * exp(-x) ! Eq. (5.37)
             dE_e = E_e(j)*(Eint_e - 1.d0)
             dP_eps = N_e(j)*P_eps*dE_e
             sum = dP_eps + sum
          endif
       enddo                ! N_e 
       Lsy(i) = sum         ! [erg/s/erg]
       
    enddo  ! eps(i),nu(i)

    ! Define the output vector for the first m_data positions
    S_data = 0.d0
    a = 4.d0*pi*d**2
    do i = 1, m_data
      S_data(i) = Lsy(i)/a
    enddo

    deallocate( nu, Lsy )
    deallocate( N_e, E_e )
    
  END subroutine radsyn_data_run
   
  
!==============================
  
END module radsyn



