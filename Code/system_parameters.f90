module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   !===SYSTEM PARAMETERS=========================================================================================

   real(dp), parameter :: dist = 16.8*Mpc          !Distance to the system
   real(dp), parameter :: jo = 17.d0*pi/180.d0     !Observer's angle with the jet axis (rad)
   real(dp), parameter :: M_BH = 6.5d9             !Black hole mass [Msun]
   real(dp), parameter :: L_Edd = 1.d38*M_BH       !Eddington luminosity
   real(dp), parameter :: R_g = 1.5d13*(M_BH/1.d8) !Gravitational radius (M_BH [Msun])
   real(dp), parameter :: z0 = 1.d1*R_g            !Jet launching height
   real(dp), parameter :: zmax = 1.d3*R_g          !Jet maximum height
   real(dp), parameter :: EphM = kB*T_CMB          !Target photon energy in the obs frame
   real(dp), parameter :: nf0 = u_CMB/EphM         !Target photon number density in the obs frame

   !===JET PROPERTIES============================================================================================

   real(dp), parameter :: gj = 10.d0              !Lorentz factor
   real(dp), parameter :: vj = sqrt(gj**2-1d0)*c/gj        !Fluid bulk speed
   real(dp), parameter :: Bj = vj/c              !Beta
   real(dp), parameter :: thetaj = 0.1d0         !Half opening angle (rad)  
   real(dp), parameter :: Lj = L_Edd*5.d-3       !Jet luminosity (erg/s)
   real(dp), parameter :: fnT = 1d-1             !Nonthermal energy fraction
   real(dp), parameter :: LnT = Lj*fnT           !Available nonthermal luminosity
   real(dp), parameter :: LnTFF = LnT/gj**2      !Available nonthermal luminosity in the FF
   real(dp), parameter :: Pdotj = Lj*gj*Bj/(c*(gj-1d0))    !Jet lateral momentum flux
   real(dp), parameter :: fB = 1.0d-2              !Fraction of magnetic energy
   real(dp), parameter :: B0 = sqrt(fB*8d0*Lj/(vj*(thetaj*z0)**2))  !Magnetic field at the jet base
   real(dp), parameter :: exp_B = 1.d0             ! B = B0*(z/z0)^exp_B

   !===NONTHERMAL PARTICLES======================================================================================

   !Particle injection
   real(dp), parameter :: facc = 1d-1        !Acceleration efficiency
   real(dp), parameter :: EminFF = 1d6*eV    !Minimum energy in the FF
   real(dp), parameter :: EmaxFF = 2d14*eV   !Maximum energy in the FF
   real(dp), parameter :: p = 2d0            !Power law exponent

   !Particle distribution
   real(dp), parameter :: E1 = EminFF  !Minimum energy considered in the FF
   real(dp), parameter :: E2 = EmaxFF  !Maximum energy considered in the FF

   !============================================================================


END module system_parameters
