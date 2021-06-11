module system_parameters
   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   !===SYSTEM PARAMETERS=========================================================================================

   real(dp), parameter :: dist = 16.8*Mpc          !Distance to the system (EHT+2018?)
   real(dp), parameter :: jo = 17.d0*pi/180.d0     !Observer's angle with the jet axis (rad) (EHT+2018?)
   real(dp), parameter :: M_BH = 6.5d9             !Black hole mass [Msun] (EHT+2018?)
   real(dp), parameter :: L_Edd = 1.25d38*M_BH     !Eddington luminosity
   real(dp), parameter :: R_g = 1.477d5*M_BH       !Gravitational radius (M_BH [Msun])
   real(dp), parameter :: EphM = kB*T_CMB          !Target photon energy in the obs frame
   real(dp), parameter :: nf0 = u_CMB/EphM         !Target photon number density in the obs frame
   
   !===JET PROPERTIES============================================================================================
   
   real(dp), parameter :: gj = 1.5d0               !Lorentz factor  (Algaba+2021: Gamma ~ 1-3)
   real(dp), parameter :: z0 = gj**2*R_g           !Jet launching height
   real(dp), parameter :: zmax = 3d5*R_g           !Jet maximum height
   real(dp), parameter :: zmax_inj = 5.d4*R_g      !Jet maximum cell of injection
   real(dp), parameter :: vj = sqrt(gj**2-1d0)*c/gj!Fluid bulk speed
   real(dp), parameter :: Bj = vj/c                !Beta
   real(dp), parameter :: thetaj = 0.6d0           !Half opening angle (rad)  
   real(dp), parameter :: r0 = 5.d0*R_g            !Jet launching radius (Kim+2018)
   real(dp), parameter :: geomIndex = 0.5d0        !Jet opening parameter (Kim+2018); for r > 10^3 R_S see Mertens+2016
   real(dp), parameter :: exp_B = 1.d0             !B = B0*(r/r0)^exp_B

   !===FREE PARAMETERS=========================================================================================

   real(dp) :: Lj = L_Edd*4d-5            !Jet luminosity (erg/s) (Algaba+2021? arbitrary)
   real(dp) :: fnT != 1.d-2               !Nonthermal energy fraction
   real(dp) :: LnT != Lj*fnT              !Available nonthermal luminosity
   real(dp) :: LnTFF != LnT!/gj**2        !Available nonthermal luminosity in the FF

   real(dp) :: fB = 1d-1                                 !Fraction of magnetic energy in the jet
   real(dp) :: B0 != sqrt(fB*8d0*Lj/(vj*(thetaj*z0)**2))  !Magnetic field at the jet base

   !===NONTHERMAL PARTICLES======================================================================================

   !Particle injection
   real(dp), parameter :: facc = 1d-2           !Acceleration efficiency
   real(dp), parameter :: EminFF = 1d6*eV       !Minimum electron energy in the FF
   real(dp), parameter :: Emin_injFF = 1d6*eV   !Minimum electron energy of injected particles in the FF
   real(dp), parameter :: EmaxFF = 1d15*eV      !Maximum electron energy in the FF
   real(dp), parameter :: EpminFF = 1d9*eV      !Minimum proton energy in the FF
   real(dp), parameter :: Epmin_injFF = 1d9*eV  !Minimum proton energy of injected particles in the FF
   real(dp), parameter :: EpmaxFF = 1d17*eV     !Maximum proton energy in the FF
   real(dp), parameter :: p = 2.12d0            !Power law exponent for Q(E) = cte * E**-p
   real(dp), parameter :: a_inj = 1.00001d0     !Power law exponent for L(z) = cte_inj * z**-a_inj, a_inj =/ 1.d0 (else error!)

   !============================================================================

   contains

      !*************************************************************************************
      !   Jet radius as a function of z
      !*************************************************************************************
      real(dp) function r_fun(zz)
         real(dp), intent(in) :: zz
         real(dp) :: zz_
         zz_ = abs(zz)
         r_fun = r0 * (zz_/z0)**geomIndex
         return
      end function
      !*************************************************************************************
      !   Jet FF magnetic field at height z
      !*************************************************************************************
      real(dp) function BFF_fun(zz)
         real(dp), intent(in) :: zz
         real(dp) :: rr, zz_
         !print*,B0
         zz_ = abs(zz)
         rr = r_fun(zz_)
         BFF_fun = (B0/gj) * (r0/rr)**exp_B
         return
      end function
      !*************************************************************************************
      !   Jet lateral velocity at height z ( = dr/dz * vj)
      !*************************************************************************************
      real(dp) function v_lat(zz)
         real(dp), intent(in) :: zz
         real(dp) :: zz_
         zz_ = abs(zz)
         v_lat = vj * (r0/z0) * geomIndex * (z0/zz_)**(1d0-geomIndex)
         return
      end function

   !============================================================================

END module system_parameters
