module dist

  use global

  implicit none
  private
  public :: diste_run

  contains

  !============================================================================
  ! diste_run
  !============================================================================

  subroutine diste_run  
    use global
    implicit none
    real(dp) :: integral, Ka
    allocate( N_e(mE_e),E_e(mE_e) )

    write(*,*) 'Calculo la distribución de electrones NT'

    ! Define a logarithmically spaced energy vector
    ! The subroutine vector_log is defined in global.f90
    call vector_log(Emin_e, 10.d0*Emax_e, Eint_e, E_e)

    ! Obtain the normalization constant:
    ! N = Ka*E**(-alpha), Ka from: int N(E')*E'dE' = E_NT = U*Vol     
    integral = 0.d0
    do i = 1, mE_e
      N_e(i) = E_e(i)**(-p) * (exp(-E_e(i)/Emax_e))
      dE_e = E_e(i)*(Eint_e-1.d0)
      integral = integral + E_e(i)*N_e(i)*dE_e
    enddo
    Ka = U*Vol/integral   
    N_e = Ka*N_e       

    ! Save the particle energy distribution
    open(79,file='dist_e.dat')
    do i = 1, mE_e
      write(79,*) E_e(i)/eV, N_e(i)
    enddo
    close(79)

    deallocate( N_e,E_e ) 

  END subroutine diste_run

END module dist
