module global

implicit none
public

!===================================================

integer, parameter :: dp = kind(1.d0)

real(dp), parameter :: pi = acos(-1.d0)
real(dp), parameter :: eV = 1.602176d-12
real(dp), parameter :: c = 2.99792458d10
real(dp), parameter :: k = 1.380658d-16      ! erg/K , boltzmann's constant
real(dp), parameter :: qe = 4.8032068d-10
real(dp), parameter :: mec2 = 0.511d6*eV
real(dp), parameter :: re = 2.818d-13        ! electron classical radius 
real(dp), parameter :: sigma = 5.6705d-5     ! stephan-boltzmann's constant
real(dp), parameter :: sigma_T = 6.6524d-25  ! Thomson cross-section
real(dp), parameter :: h = 6.6260755d-27     ! planck's constant (erg*s)
real(dp), parameter :: h_bar = h/(2.*pi)
real(dp), parameter :: alfa = 7.29735308d-3 
real(dp), parameter :: Msolyr = 1.989d33/31556926.d0 ! (m_sol/yr) in CGS
real(dp), parameter :: Rsol = 6.96d10        ! solar radius (cm)
real(dp), parameter :: Lsol = 4.d33          ! solar luminosity (CGS)
real(dp), parameter :: G_N = 6.67259d-8      ! Gravitational constant
real(dp), parameter :: UA = 1.496d13         ! 1 UA in cm
real(dp), parameter :: pc = 3.0856d18        ! 1 pc = 206265.d0*UA
real(dp), parameter :: mp = 1.673d-24        ! proton mass in CGS
real(dp), parameter :: ga = 4.d0/3.d0        ! the adiabatic coefficient 
real(dp), parameter :: m_pic2 = 134.98d6*eV  ! neutral pion rest-mass energy
real(dp), parameter :: mpc2 = 1836.d0*mec2   ! proton rest-mass energy

!===================================================

real(dp), parameter :: U = 100.d0*eV      ! [U_e_NT] = erg/cm^3
real(dp) :: B = 100.d-6        ! B = 100 uG
!real(dp) :: B_par = B*sqrt(2./3.) ! Assuming isotropy
real(dp), parameter :: p = 2.d0           ! Injection spectral index
real(dp), parameter :: R = 0.1*pc         ! Hot-spot size
real(dp), parameter :: Vol = (4.d0/3.d0)*pi*R**3 ! Assume spherical region 
real(dp), dimension(:), allocatable :: E_e, N_e, Q_e
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
       E_vec(i_vec)=E_vec(i_vec-1)*Eint    
    enddo

    RETURN

  END SUBROUTINE vector_log


!============================================================================
  

END module global
