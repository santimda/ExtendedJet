MODULE globals
  use constants
  use system_parameters
  implicit none
  public

  !====================================================================
  ! This module contains constants (from constants.f90), variables and 
  ! parameters that all other modules use. 
  ! It also has defined some subroutines that many modules call.
  !====================================================================

  integer, parameter :: nth = 4   ! number of threads (parallelization)

  !===========================================================
  ! Simulation grid/vectors size
  ! Ns, numbins and dttFF_min are now runtime variables set by the &run
  ! namelist (see system_parameters.f90); their defaults live there.
  integer, parameter :: nIn = 150      !Number of points in the electron energy distribution (default: 150)
  integer, parameter :: numdt = 100000  !Number of segments for integrating the energy distribution (default: 100000)
  integer, parameter :: m_data_max = 300 !Maximum number of points to fit
  integer, parameter :: NSED_sy = 150    !Number of points for the synchrotron SED
  integer, parameter :: NSED_IC = 100    !Number of points for the IC SED

  real(dp), dimension(:), allocatable :: deltaZ, deltatFF

  !System-observer vector (unitary)
  real(dp) :: SO(3) 

  !===SEDs======================================================================================================

  real(dp), parameter :: Emin_sy = 1d-8*eV  !Minimum energy for the synchrotron SED in the OF
  real(dp), parameter :: Emax_sy = 2d10*eV  !Maximum energy for the synchrotron SED in the OF

  real(dp), parameter :: Emin_IC = 1d0*eV   !Minimum energy for the IC SED in the OF
  real(dp), parameter :: Emax_IC = 1d14*eV  !Maximum energy for the IC SED in the OF

  !=============================================================================================================

  real(dp), dimension(:), allocatable :: z, rj_out, rj_in      !Jet position, radius

  real(dp), dimension(:), allocatable :: delta_obs  !Doppler factor w.r.t. observer
  real(dp), dimension(:), allocatable :: boost      !Doppler boosting (FF -> OF)

  real(dp), dimension(:), allocatable :: gj, bj     !Jet Lorentz factor and velocity (beta)
  real(dp), dimension(:), allocatable :: BFF, vnorm, VFF !Magnetic field, jet expansion velocity, segment volume
   
  real(dp), dimension(:), allocatable :: Kappa, EphMFF, nfFF !Dilution factor, photon energy and number density in the FF
  real(dp), dimension(:), allocatable :: joFF, sjoFF      !Jet-observer angle and sin(jet-observer angle) in the FF

  logical, parameter :: ani = .True., iso = .False. ! To define isotropic or anisotropic IC
  real(dp), parameter :: cte_IC_K = 2d0*re**2*mec2**3/(pi*h_bar**3*c**2) !IC constant in Khangulyan 2014
  integer, parameter :: N_bands = 5                   !Number of flux bands
  integer :: Ninj

  real(dp) :: dummy 

  !=============================================================================================================

  contains

  !=============================================================================================================
  ! Energy vector
  !=============================================================================================================

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
    do i_vec = 2, N_puntos
      E_vec(i_vec) = E_vec(i_vec-1)*Eint    
    enddo
    RETURN

  END SUBROUTINE vector_log

  !=============================================================================================================
  ! Integration
  !=============================================================================================================

  SUBROUTINE integrate(N,H,FI,S)
    ! Subroutine for integration over f(x) with the Simpson rule.
    ! FI: integrand f(x); H: interval; S: integral.
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: N
    real(dp), intent(in) :: H,FI(N)
    real(dp), intent(out) :: S
    real(dp) :: S0,S1,S2
    integer :: I

    S  = 0.d0
    S0 = 0.d0
    S1 = 0.d0
    S2 = 0.d0
    DO I = 2, N-1, 2
      S1 = S1 + FI(I-1)
      S0 = S0 + FI(I)
      S2 = S2 + FI(I+1)
    END DO
    S = H*(S1 + 4.*S0 + S2)/3.
    !If N is even, add the last slice separately
    IF(MOD(N,2).eq.0) S = S + H*(5.d0*FI(N) + 8.d0*FI(N-1) - FI(N-2))/12.d0
    RETURN

  End SUBROUTINE integrate

  !=============================================================================================================
  ! IC radiation (Khangulyan 2014)
  !=============================================================================================================
  !----------------------------------------------------------------------------
  ! For the energy losses (dE/dt) 
  !----------------------------------------------------------------------------
  real(dp) function F_ani(u)  !Anisotropic IC (eq. 35)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: u
  real(dp), parameter :: c_ani = 6.13d0

  F_ani = c_ani*u*log(1d0+2.16d0*u/c_ani)/(1d0+c_ani*u/0.822d0)

  return
  end function F_ani
  !----------------------------------------------------------------------------
  real(dp) function F_iso(u)  !Isotropic IC (eqs. 38,20)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: u
  real(dp) :: g_u
  real(dp), parameter :: c_iso = 5.68d0

  g_u = (1d0 + (-0.362d0*u**0.682d0)/(1d0 + 0.826d0*u**1.281d0))**(-1d0)
  F_iso = c_iso*u*log(1d0+0.722d0*u/c_iso)/(1d0+c_iso*u/0.822d0) * g_u

  return
  end function F_iso

  !----------------------------------------------------------------------------
  ! For the interaction rate per energy unit (dN/dEdt) 
  !----------------------------------------------------------------------------   
  real(dp) function F_1(x0)  !Anisotropic IC (eqs. 18-20)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: x0
  real(dp) :: g

  g = (1d0 + (0.153d0*x0**0.857d0)/(1d0 + 0.254d0*x0**1.84d0))**(-1d0)
  F_1 = (pi**2/6d0 + x0)*exp(-x0) * g

  return
  end function F_1
  !----------------------------------------------------------------------------
  real(dp) function F_2(x0)  !Anisotropic IC (eqs. 18-20)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: x0
  real(dp) :: g

  g = (1d0 + (1.33d0*x0**0.691d0)/(1d0 + 0.534d0*x0**1.668d0))**(-1d0)
  F_2 = (pi**2/6d0 + x0)*exp(-x0) * g

  return
  end function F_2
  !----------------------------------------------------------------------------  
  real(dp) function F_3(x0)  !Isotropic IC (eqs. 24,25,20)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: x0
  real(dp) :: g
  real(dp), parameter :: c_3 = 0.319d0

  g = (1d0 + (0.443d0*x0**0.606d0)/(1d0 + 0.540d0*x0**1.481d0))
  F_3 = pi**2/6d0*(1d0+c_3*x0)/(1d0+pi**2*c_3*x0/6d0)*exp(-x0) * g

  return
  end function F_3
  !----------------------------------------------------------------------------
  real(dp) function F_4(x0)  !Isotropic IC (eqs. 24,25,20)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: x0
  real(dp) :: g
  real(dp), parameter :: c_4 = 6.62d0

  g = (1d0 + (0.726d0*x0**0.461d0)/(1d0 + 0.382d0*x0**1.457d0))**(-1d0)
  F_4 = pi**2/6d0*(1d0+c_4*x0)/(1d0+pi**2*c_4*x0/6d0)*exp(-x0) * g

  return
  end function F_4

  !======================================================================================

  !real(dp) function nBB(Epho,T)
  !use constants
  !real(dp), intent(in) :: Epho,T
  !nBB = 2d0*pi/((h*c)**3)*(Epho**2/(exp(Epho/(kB*T))-1d0))
  !return
  !end function nBB

  !====================================================================================

  SUBROUTINE versorresta(V1,V2,N,versor)
    ! subroutine to calculate the versor: (V1-V2)/|V1-V2|
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: N
    integer :: i,j
    real(dp), intent(in) :: V1,V2
    real(dp), intent(out) :: versor
    real(dp) :: resta,suma
    DIMENSION V1(N),V2(N),versor(N),resta(N)

    suma = 0.d0
    do i = 1, N
      resta(i) = V1(i) - V2(i)
      suma = suma + resta(i)**2
    enddo
    suma = sqrt(suma)          ! suma = || versor_resta ||
    do j = 1,N
      versor(j) = resta(j)/suma
    enddo

    RETURN
  END SUBROUTINE versorresta

  !============================================================================
  ! p-p cross section

  Function sigma_pp(L1,Ep)
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: L1, Ep
    real(dp) :: sigma_pp
    if (Ep < E_th) then !E_th= 1.22 GeV
      sigma_pp = 0.d0
    else           
      sigma_pp = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27 !cm²
      sigma_pp = sigma_pp*(1.d0-(E_th/Ep)**4)**2 ! Add factor for better fit at low energies 
    endif
    return
  end Function sigma_pp


  !============================================================================
  ! Logarithmic interpolation
  integer function locate(xx, x)
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: x, xx(:)
    ! Given an array xx(1:n) , and given a value x , returns a value j such that x is between
    ! xx(j) and xx(j+1) . xx(1:n) must be monotonic, either increasing or decreasing. j=0
    ! or j=n is returned to indicate that x is out of range.
    integer :: jl, jm, ju
    integer :: j, n

    n = SIZE(xx) ! x and v have to have the same size
    jl = 0                            ! Initialize lower
    ju = n+1                          ! and upper limits.
    
    do while (ju-jl > 1)
      jm = (ju + jl)/2
      if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
        jl = jm
      else
        ju = jm
      endif
    enddo
    
    if (x == xx(1)) then
      j = 1
    else if (x == xx(n)) then
      j = n-1
    else
      j = jl
    endif
    locate = j
    end function locate

  !============================================================================
  ! Logarithmic interpolation
  real(dp) function Interpol1d(xp, x, v)
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: xp, v(:), x(:)
  integer :: j_x, n
  real(dp) :: slope
  n = SIZE(v) ! x and v have to have the same size
  j_x = locate(x, xp)
  if (j_x == 0) then
    slope = (v(2) - v(1)) / (x(2) - x(1))
    Interpol1d = max(v(1) + slope * (xp - x(1)), 0d0)
  else if (j_x == n) then
    slope = (v(n) - v(n-1)) / (x(n) - x(n-1))
    Interpol1d = max(v(n) + slope * (xp - x(n)), 0d0)
  else
    if (v(j_x) > 0d0 .and. v(j_x + 1) > 0d0) then
      slope = log(v(j_x+1)/v(j_x))/log(x(j_x+1)/x(j_x))
      Interpol1d = v(j_x) * (xp/x(j_x))**slope
    else
      slope = (v(j_x+1) - v(j_x)) / (x(j_x+1) - x(j_x))
      Interpol1d = v(j_x) + slope * (xp - x(j_x))
    endif
  endif

  return
  end function Interpol1d

  !============================================================================
  
END MODULE globals
