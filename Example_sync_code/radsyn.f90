module radsyn
  ! This module calculates the intrinsic synchrotron emission in the fluid frame
  use global
  
  implicit none
  private
  public :: radsyn_run

contains
  
  !============================================================================
  ! radsyn_run
  !============================================================================
  
  subroutine radsyn_run  
    use global
    !$ use omp_lib
    implicit none
    ! mEe -> electrons, m_nu -> emitted photons
    real(dp) :: eps(m_nu), dEps, eps_int, gamma_e(mE_e), x
    real(dp) :: nuc, nuc_max, nu_min, nu_max, sum, sum_sy, cte_P
    allocate( N_e(mE_e), E_e(mE_e) )
    allocate( nu(m_nu), Lsy(m_nu) )

80  format(A10,X,Es10.2)

    !------------------------------------------------------------------------
    !
    !  Reading information here  
    !

    ! N_e(i) = number of e- with energy E_e(i) per unit energy
    open (14,file='dist_e.dat')
    do i = 1, mE_e
       read(14,*) E_e(i),N_e(i)
       E_e(i) = E_e(i)*eV    !   [erg]
    enddo
    close(14) 
    Eint_e = E_e(2)/E_e(1)
    gamma_e = E_e/mec2
    
    !------------------------------------------------------------------------
    !
    !  Radiative calculations start here  
    !
    
    ! Define nu_min and nu_max:
    nuc_max = 4.22d6*B_par*gamma_e(mE_e)**2
    nu_max = 20.d0*nuc_max ! 
    nu_min = 1.d6          ! 1 MHz
    ! Define a vector with the photon energies epsilon: eps(i) = h*nu(i)
    call vector_log(h*nu_min, h*nu_max, eps_int, eps)
    nu = eps/h             ! [nu] = Hz 

    cte_P = sqrt(3.)*qe**3/(h*mec2)
    sum_sy = 0.d0 ! To estimate the total luminosity: L_sy = int P_sy(eps) dEps
    ! We calculate the luminosity at each value of the photon energy eps(i) between eps_min and eps_max
    do i = 1,m_nu           ! eps(i), or nu(i)
       dEps = eps(i)*(eps_int-1.d0)
       sum = 0.d0
       
       ! P(eps) = integral_E_e(min)^E_e(max) N_e(E_e) P(E_e,eps) dE_e [erg/s/erg]
       ! Define: N_e(E_e)*P(E_e,eps)*dE_e = dP_eps , so that P(eps) = integral_E_e dP_eps
       do j = 1,mE_e             ! N_e(j), nuc(j), gamma_e(j)
          nuc = 4.22d6*B_par*gamma_e(j)**2 ! synchrotron characteristic frequency from Eq. (5.38)
          x = nu(i)/nuc       ! 0.1 < x = nu/nu_c < 10 --> Use approximation from Eq. (5.41)
          if (x > 10.d0) then
             dP_eps = 0.d0
          else
             ! Calculate P_eps [erg/s/Hz] and transform to erg/s/erg dividing by h
             P_eps = cte_P * B_par * ( 1.85*x**(1./3.)*exp(-x) )! Eq. (5.37)
             dE_e = E_e(j)*(Eint_e - 1.d0)
             dP_eps = N_e(j)*P_eps*dE_e
          endif
          sum = dP_Eps + sum
       enddo                ! N_e 
       Lsy(i) = sum         ! [erg/s/erg]
       sum_sy = sum_sy + Lsy(i)*deps
       
    enddo  ! eps(i),nu(i)
    write(*,80) 'L_sy =',sum_sy

    ! Ya tengo la emisión total no absorbida:
    open(17,file='radsyn.dat')
    do i = 1, m_nu
       write(17,*) eps(i)/eV, eps(i)*Lsy(i) ! Eps*L(Eps) for the SED
    enddo
    close(17)

    deallocate( nu,Lsy )
    deallocate( N_e,E_e )
    
  END subroutine radsyn_run
   
  
!==============================

  
END module radsyn



