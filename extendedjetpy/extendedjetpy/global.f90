module global

use constants
implicit none
public

!===================================================

real(dp), parameter :: U = 100.d0*eV      ! [U_e_NT] = erg/cm^3
real(dp) :: B = 100.d-6        ! B = 100 uG
!real(dp) :: B_par = B*sqrt(2./3.) ! Assuming isotropy
real(dp), parameter :: p = 2.d0           ! Injection spectral index
real(dp), parameter :: R = 0.1*pc         ! Hot-spot size
real(dp), parameter :: Vol = (4.d0/3.d0)*pi*R**3 ! Assume spherical region 
real(dp), dimension(:), allocatable :: E_e, N_e
real(dp), parameter :: Emin_e = 1.d8*eV   ! 100 MeV
real(dp) :: Emax_e = 1.d13*eV  ! 10 TeV
real(dp) :: Eint_e,dE_e
integer, parameter :: mE_e = 400          ! Energy resolution for electrons
integer :: i,j,l
real(dp), parameter :: d = 1.d3*pc      ! 1 kpc


!===================================================
! emission

real(dp) :: P_Eps, dP_Eps
integer, parameter :: m_nu = 800  ! frequency resolution for synchrotron photons
real(dp), dimension(:), allocatable :: nu, Lsy

integer, parameter :: m_data_max = 50 ! maximum number of points to fit

!===================================================

contains

!============================================================================
! vector_log  --> Sets a logarithmic step 
!============================================================================
  
  SUBROUTINE vector_log(Emin,Emax,Eint,E_vec)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: Emin,Emax
    integer :: N_puntos, i_vec
    real(dp), intent(out) :: Eint
    real(dp), intent(out) :: E_vec(:)

    N_puntos = SIZE(E_vec)
    Eint = (Emax/Emin)**(1.d0/dfloat(N_puntos-1))

    E_vec(1) = Emin
    do i_vec=2,N_puntos
       E_vec(i_vec) = E_vec(i_vec-1)*Eint    
    enddo

    RETURN

  END SUBROUTINE vector_log


!============================================================================
  

END module global
