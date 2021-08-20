module global

use constants
implicit none
public

!===================================================
! free parameters

real(dp) :: U = 1.d6*eV      ! [U_e_NT] = erg/cm^3
real(dp) :: B = 100.d-6      ! B = 100 uG

!==================================================
! parameters 

real(dp), parameter :: p = 2.d0           ! Injection spectral index
real(dp), parameter :: R = 0.01*pc        ! Hot-spot size
real(dp), parameter :: Vol = (4.d0/3.d0)*pi*R**3 ! Assume spherical region 
real(dp), parameter :: Emin_e = 1.d8*eV   ! 100 MeV
real(dp) :: Emax_e = 1.d13*eV             ! 10 TeV
real(dp) :: Eint_e,dE_e
integer, parameter :: mE_e = 400          ! Energy resolution for electrons
integer :: i,j,l
real(dp), parameter :: d = 1.d3*pc        ! 1 kpc

real(dp), dimension(:), allocatable :: E_e, N_e

!===================================================
! emission

integer, parameter :: mE_sy = 600       ! energy resolution for synchrotron photons
integer, parameter :: mE_ph = 250       ! energy resolution for SSC photons
integer, parameter :: m_data_max = 100  ! maximum number of points to fit

real(dp) :: E_sy(mE_sy), L_sy(mE_sy), L_sy_tot
real(dp) :: E_ic(mE_ph), L_ic(mE_ph), L_ic_tot
real(dp), dimension(:), allocatable :: L_sy_data
real(dp), dimension(:), allocatable :: L_ic_data

!===================================================

real(dp) :: dummy ! dummy variable

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
