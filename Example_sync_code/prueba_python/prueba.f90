subroutine prueba
  !============================================================================
  ! calculate electron distribution radiation
  !============================================================================
  !use global
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: mE_e = 100
  real(dp), parameter :: Emin_e=1.d3, Emax_e=1.d6, p=2.d0
  real(dp) :: integral, N_e(mE_e), E_e(mE_e), dE_e, Eint_e
  integer :: i
  write(*,*) 'Calculo la distribución de electrones NT'

  ! Obtain the normalization constant:
  ! N = Ka*E**(-alpha), Ka from: int N(E')*E'dE' = E_NT = U*Vol     
  integral = 0.d0
  E_e(1) = 1.d0
  Eint_e = 2.d0
  do i = 1, mE_e
    N_e(i) = E_e(i)**(-p) 
    dE_e = E_e(i)*(Eint_e-1.d0)
    integral = integral + E_e(i)*N_e(i)*dE_e
  enddo    


end subroutine prueba