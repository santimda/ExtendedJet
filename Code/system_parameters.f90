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
   real(dp), parameter :: L_Edd = 1.25d38*M_BH     !Eddington luminosity
   real(dp), parameter :: R_g = 1.477d5*M_BH       !Gravitational radius (M_BH [Msun])
   real(dp), parameter :: z0 = 4*R_g               !Jet launching height
   real(dp), parameter :: zmax = 6.5d2*R_g         !Jet maximum height
   real(dp), parameter :: zmax_inj = 5.d1*R_g       !Jet maximum cell of injection
   real(dp), parameter :: EphM = kB*T_CMB          !Target photon energy in the obs frame
   real(dp), parameter :: nf0 = u_CMB/EphM         !Target photon number density in the obs frame

   !===JET PROPERTIES============================================================================================

   real(dp), parameter :: gj = 1.5d0              !Lorentz factor
   real(dp), parameter :: vj = sqrt(gj**2-1d0)*c/gj        !Fluid bulk speed
   real(dp), parameter :: Bj = vj/c              !Beta
   real(dp), parameter :: thetaj = 0.1d0         !Half opening angle (rad)  
   real(dp), parameter :: r0 = tan(thetaj) * z0  !Jet launching radius
   real(dp), parameter :: geomIndex = 0.6d0        !Jet opening parameter
   real(dp), parameter :: Lj = L_Edd*1.d-5       !Jet luminosity (erg/s)
   real(dp), parameter :: fnT = 3d-2             !Nonthermal energy fraction
   real(dp), parameter :: LnT = Lj*fnT           !Available nonthermal luminosity
   real(dp), parameter :: LnTFF = LnT!/gj**2      !Available nonthermal luminosity in the FF
   real(dp), parameter :: Pdotj = Lj*gj*Bj/(c*(gj-1d0))    !Jet lateral momentum flux
   real(dp), parameter :: fB = 2.5d-3              !Fraction of magnetic energy
   real(dp), parameter :: B0 = sqrt(fB*8d0*Lj/(vj*(thetaj*z0)**2))  !Magnetic field at the jet base
   real(dp), parameter :: exp_B = 1.d0             ! B = B0*(z/z0)^exp_B

   !===NONTHERMAL PARTICLES======================================================================================

   !Particle injection
   real(dp), parameter :: facc = 1d-1        !Acceleration efficiency
   real(dp), parameter :: EminFF = 1d6*eV    !Minimum energy in the FF
   real(dp), parameter :: EmaxFF = 2d14*eV   !Maximum energy in the FF
   real(dp), parameter :: p = 2.13d0            !Power law exponent

   !Particle distribution
   real(dp), parameter :: E1 = EminFF  !Minimum energy considered in the FF
   real(dp), parameter :: E2 = EmaxFF  !Maximum energy considered in the FF

   !============================================================================

   contains
      !*************************************************************************************
      !   Jet radius as a function of z
      !*************************************************************************************
      real(dp) function r_fun(zz)
         real(dp), intent(in) :: zz
         r_fun = r0 * (zz/z0)**geomIndex
         return
      end function
      !*************************************************************************************
      !   Jet FF magnetic field at height z
      !*************************************************************************************
      real(dp) function BFF_fun(zz)
         real(dp), intent(in) :: zz
         real(dp) :: rr
         rr = r_fun(zz)
         BFF_fun = (B0/gj) * (r0/rr)**exp_B
         return
      end function
      !*************************************************************************************
      !   Jet lateral velocity at height z ( = dr/dz * vj)
      !*************************************************************************************
      real(dp) function v_lat(zz)
         real(dp), intent(in) :: zz
         v_lat = vj * tan(thetaj) * (z0/zz)**(1d0-geomIndex)
         return
      end function

END module system_parameters
