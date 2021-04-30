MODULE globals
  implicit none
  
  integer :: i,i2,j,j2,l,n

!===CONSTANTS IN CGS UNITS====================================================================================
  
  real*8, parameter :: pi = 4d0*atan(1d0)		!Number pi
  real*8, parameter :: c = 2.99792d10			!Speed of light in vacuum
  real*8, parameter :: eV = 1.60218d-12			!Electronvolt
  real*8, parameter :: me = 9.10938d-28			!Electron mass
  real*8, parameter :: mec2 = 5.10999d5*eV      !Electron rest energy
  real*8, parameter :: re = 2.81794d-13         !Electron classical radius
  real*8, parameter :: mp = 1.67262d-24			!Proton mass
  real*8, parameter :: mpc2 = 9.38272d8*eV		!Proton rest energy
  real*8, parameter :: qe = 4.80320d-10			!Elemental charge
  real*8, parameter :: Msun = 1.9885d33			!Solar mass
  real*8, parameter :: Msun_yr = 6.30118d25		!Msun/yr
  real*8, parameter :: Rsun = 6.96392d10		!Solar radius 
  real*8, parameter :: Lsun = 3.828d33			!Solar luminosity
  real*8, parameter :: AU = 1.49598d13			!Astronomical Unit
  real*8, parameter :: pc = 3.0857d18			!Parsec
  real*8, parameter :: kpc = 3.0857d21			!Kiloparsec
  real*8, parameter :: kB = 1.38065d-16			!Boltzmann constant
  real*8, parameter :: sigma = 5.67037d-5		!Stefan-Boltzamnn constant
  real*8, parameter :: h = 6.6260755d-27		!Planck's constant
  real*8, parameter :: hbar = h/(2d0*pi)		!Reduced Planck constant
  real*8, parameter :: Ry = 2.17987d-11         !Rydberg unit = 2.18e-18 J
  
!===JET PROPERTIES============================================================================================
  
  real*8, parameter :: gj = 3.d0						!Lorentz factor
  real*8, parameter :: vj = sqrt(gj**2-1d0)*c/gj        !Fluid bulk speed
  real*8, parameter :: Bj = vj/c						!Beta
  real*8, parameter :: thetaj = 0.1d0					!Half opening angle (rad)  
  real*8, parameter :: Lj = 5d36						!Jet luminosity (erg/s)
  real*8, parameter :: fnT = 1d-1						!Nonthermal energy fraction
  real*8, parameter :: LnT = Lj*fnT						!Available nonthermal luminosity
  real*8, parameter :: LnTFF = LnT/gj**2				!Available nonthermal luminosity in the FF
  real*8, parameter :: Pdotj = Lj*gj*Bj/(c*(gj-1d0))    !Jet lateral momentum flux
  
!===STAR PROPERTIES===========================================================================================
  
  real*8, parameter :: Lstar = 1d39										!Luminosity (erg/s)
  real*8, parameter :: Tstar = 4d4										!Temperature (K)
  real*8, parameter :: vW_inf = 2d8										!Terminal wind speed (cm/s)
  real*8, parameter :: EphM = 2.7d0*kB*Tstar							!Mean photon energy
  real*8, parameter :: Mdot = 1d-6*Msun_yr								!Mass loss rate
  real*8, parameter :: Rstar = sqrt(Lstar/(pi*sigma))/(2d0*Tstar**2)	!Radius
  
!===SYSTEM PARAMETERS=========================================================================================
  
  real*8 :: orbPhase										!Orbital phase	
  real*8, parameter :: dist = 3d0*kpc						!Distance to the system
  real*8, parameter :: a = 3d12								!Star-CO distance (cm)
  real*8, parameter :: T = 4d0*24d0*3600d0					!Period (s)
  real*8, parameter :: inc = 30d0							!Observer's angle with the z axis (0-90º)	
  real*8, parameter :: azi = 90d0							!Observer's angle with the x axis (0-360º)
  real*8, parameter :: dl = 1d11							!Length of jet segment/blob (cm)
  real*8, parameter :: lf = 50d0*a 							!Final jet length
  integer, parameter :: Ns = int(lf/dl) 					!Number of jet segments
  real*8, parameter :: dt = dl/vj							!Segment advection time for constant velocity
  real*8, parameter :: dtFF = dt/gj							!Advection time in the FF
  real*8, parameter :: PdotW_CO = Mdot*vW_inf*(1d0-Rstar/a)**0.8d0	!Wind thrust at the CO postion
  real*8, parameter :: &											!Height of the recol. shock (B-R 2016)
  &	zrec = a*sqrt((PdotW_CO*gj**2-8d0*Pdotj-sqrt(PdotW_CO*gj**2*(PdotW_CO*gj**2-16d0*Pdotj)))/(8d0*Pdotj))
  real*8, parameter :: nf0 = Lstar/(4d0*pi*(a**2+zrec**2)*c*EphM)	!Target photon density at the jet base
  real*8, parameter :: fB = 1.0d0									!Fraction of magnetic energy
  real*8, parameter :: B0 = sqrt(fB*8d0*pi*EphM*nf0)				!Magnetic field at the jet base
  !real*8, parameter :: B0 = sqrt(fB*8d0*Lj/(vj*(thetaj*zrec)**2))	!Magnetic field at the jet base
  real*8, parameter :: Xj = thetaj*PdotW_CO/(4d0*pi*Pdotj)			!Bosch-Ramon 2016
  
!===NONTHERMAL PARTICLES======================================================================================
  
  !Particle injection
  real*8, parameter :: facc = 1d-1			!Acceleration efficiency
  real*8, parameter :: EminFF = 1d6*eV		!Minimum energy in the FF
  real*8, parameter :: EmaxFF = 1d14*eV		!Maximum energy in the FF
  real*8, parameter :: p = 2d0				!Power law exponent
  
  !Particle distribution
  real*8, parameter :: E1 = EminFF	!Minimum energy considered in the FF
  real*8, parameter :: E2 = EmaxFF	!Maximum energy considered in the FF
  integer, parameter :: nNE = 50	!Energy steps (for plotting N(E))
  integer, parameter :: nIn = 100	!Energy steps (for the energy integrals)
  
!===SEDs======================================================================================================
  
  real*8, parameter :: Emin_sy = 1d-6*eV	!Minimum energy for the synchrotron SED in the OF
  real*8, parameter :: Emax_sy = 1d10*eV	!Maximum energy for the synchrotron SED in the OF
  integer, parameter :: NSED_sy = 100		!Number of points for the synchrotron SED
  
  real*8, parameter :: Emin_IC = 1d0*eV		!Minimum energy for the IC SED in the OF
  real*8, parameter :: Emax_IC = 1d14*eV	!Maximum energy for the IC SED in the OF
  integer, parameter :: NSED_IC = 100		!Number of points for the IC SED
  
!=============================================================================================================
  
  real*8, dimension(:), allocatable :: x,y,z,rj,lenj,d		!Jet position, radius, length and distance to star
  real*8, dimension(:), allocatable :: jw,wo,jo,cwoFF,	&	!Jet-wind, wind-observer and jet-observer angles,
  &									   sjoFF				!and cosine of wo in the FF
  logical, dimension(:), allocatable :: isOccult			!Is the jet occulted by the star?
  
  real*8, dimension(:), allocatable :: delta_star,delta_obs 	!Doppler factors with star and observer
  real*8, dimension(:), allocatable :: BFF,nfFF,vnorm,EphMFF, &	!Magnetic and photon fields, jet expansion
  &									   VFF,Kappa				!velocity, mean energy of the stellar photons,
  																!segment volume, dilution factor
  
  real*8, parameter :: cte_IC_K = 2d0*re**2*mec2**3/(pi*hbar**3*c**2)	!IC constant in Khangulyan 2014
  integer, parameter :: N_bands = 5										!Number of flux bands
  
!===ABSORTION PARAMETERS======================================================================================

  !Gamma-gamma
  integer, parameter :: mE_gg = NSED_IC, m_ang_gg = 400
  real*8, dimension(:,:), allocatable :: tau_gg
  real*8, parameter :: Ef_min_gg = Emin_IC
  real*8, parameter :: Ef_max_gg = Emax_IC
  real*8 :: Ef_gg(mE_gg),Ef_int_gg

  !Free-free
  integer, parameter :: mE_ff = 200, m_ang_ff = 400
  real*8, dimension(:,:), allocatable :: tau_ff
  real*8, parameter :: Ef_min_ff = Emin_sy  
  real*8, parameter :: Ef_max_ff = h*100d12	! 100 THz 
  real*8, parameter :: mu_w = 1.3d0			! Leitherer 1995 ?
  real*8, parameter :: Zq_w = 1.0d0			! Z = rms ionic charge ?
  real*8, parameter :: gamma_w = 1.0d0		! Leitherer 1995 ?
  real*8, parameter :: filling = 2d-1		! ..... 20% ?

!=============================================================================================================
  
contains
  
!=============================================================================================================
! Energy vector
!=============================================================================================================
  
  SUBROUTINE vector_log(Emin,Emax,Eint,E_vec)
    real*8, intent(in) :: Emin,Emax
    integer :: N_puntos, i_vec
    real*8, intent(out) :: Eint
    real*8, intent(out) :: E_vec(:)

    N_puntos = SIZE(E_vec)
    Eint = (Emax/Emin)**(1d0/dfloat(N_puntos-1))

    E_vec(1) = Emin
    do i_vec=2,N_puntos
       E_vec(i_vec)=E_vec(i_vec-1)*Eint    
    enddo

    RETURN

  END SUBROUTINE vector_log

!=============================================================================================================
! Integration
!=============================================================================================================

  SUBROUTINE integrate(N,H,FI,S)
    ! Subroutine for integration over f(x) with the Simpson rule.
    ! FI: integrand f(x); H: interval; S: integral.
    integer, intent(in) :: N
    real*8, intent(in) :: H,FI(N)
    real*8, intent(out) :: S
    real*8 :: S0,S1,S2

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
!	IC radiation (Khangulyan 2014)
!=============================================================================================================
!----------------------------------------------------------------------------
!	For the energy losses (dE/dt) 
!----------------------------------------------------------------------------
  	real*8 function F_ani(u)  !Anisotropic IC (eq. 35)
    real*8, intent(in) :: u
    real*8, parameter :: c_ani = 6.13d0

    F_ani = c_ani*u*log(1d0+2.16d0*u/c_ani)/(1d0+c_ani*u/0.822d0)

	return
  	end function F_ani
!----------------------------------------------------------------------------
  	real*8 function F_iso(u)  !Isotropic IC (eqs. 38,20)
  	real*8, intent(in) :: u
    real*8 :: g_u
    real*8, parameter :: c_iso = 5.68d0
    
    g_u = (1d0 + (-0.362d0*u**0.682d0)/(1d0 + 0.826d0*u**1.281d0))**(-1d0)
    F_iso = c_iso*u*log(1d0+0.722d0*u/c_iso)/(1d0+c_iso*u/0.822d0) * g_u
    
    return
  	end function F_iso
  	
!----------------------------------------------------------------------------
!	For the interaction rate per energy unit (dN/dEdt) 
!----------------------------------------------------------------------------   
  	real*8 function F_1(x0)  !Anisotropic IC (eqs. 18-20)
    real*8, intent(in) :: x0
    real*8 :: g
    
    g = (1d0 + (0.153d0*x0**0.857d0)/(1d0 + 0.254d0*x0**1.84d0))**(-1d0)
    F_1 = (pi**2/6d0 + x0)*exp(-x0) * g
    
    return
  	end function F_1
!----------------------------------------------------------------------------
	real*8 function F_2(x0)  !Anisotropic IC (eqs. 18-20)
    real*8, intent(in) :: x0
    real*8 :: g
    
    g = (1d0 + (1.33d0*x0**0.691d0)/(1d0 + 0.534d0*x0**1.668d0))**(-1d0)
    F_2 = (pi**2/6d0 + x0)*exp(-x0) * g
    
    return
    end function F_2
!----------------------------------------------------------------------------  
  	real*8 function F_3(x0)  !Isotropic IC (eqs. 24,25,20)
    real*8, intent(in) :: x0
    real*8 :: g
    real*8, parameter :: c_3 = 0.319d0
    
    g = (1d0 + (0.443d0*x0**0.606d0)/(1d0 + 0.540d0*x0**1.481d0))
    F_3 = pi**2/6d0*(1d0+c_3*x0)/(1d0+pi**2*c_3*x0/6d0)*exp(-x0) * g
    
    return
	end function F_3
!----------------------------------------------------------------------------
  	real*8 function F_4(x0)  !Isotropic IC (eqs. 24,25,20)
    real*8, intent(in) :: x0
    real*8 :: g
    real*8, parameter :: c_4 = 6.62d0
    
    g = (1d0 + (0.726d0*x0**0.461d0)/(1d0 + 0.382d0*x0**1.457d0))**(-1d0)
    F_4 = pi**2/6d0*(1d0+c_4*x0)/(1d0+pi**2*c_4*x0/6d0)*exp(-x0) * g
    
	return
	end function F_4
  
!=============================================================================================================
     
  real*8 function nBB(Epho,T)
    real*8, intent(in) :: Epho,T
    nBB = 2d0*pi/((h*c)**3)*(Epho**2/(exp(Epho/(kB*T))-1d0))
    return
  end function nBB

!=============================================================================================================
	
  
END MODULE globals
