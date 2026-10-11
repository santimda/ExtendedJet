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

END MODULE globals
