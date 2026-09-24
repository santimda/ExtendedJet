module system_parameters
   use constants
   use controls
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled (CenA)
   !====================================================================

   !===SYSTEM PARAMETERS=========================================================================================

   real(dp), parameter :: dist = 16.8*Mpc !16.8*Mpc          !Distance to the system (EHT+2018?)
   real(dp), parameter :: jo = 17.d0*pi/180.d0 !17.d0*pi/180.d0     !Observer's angle with the jet axis (rad) (EHT+2018?)
   real(dp), parameter :: M_BH = 6.5d9!6.5d9             !Black hole mass [Msun] (EHT+2018?)
   real(dp), parameter :: L_Edd = 1.25d38*M_BH     !Eddington luminosity
   real(dp), parameter :: R_g = 1.477d5*M_BH       !Gravitational radius (M_BH [Msun])
   real(dp), parameter :: EphM = kB*T_CMB          !Target photon energy in the obs frame
   real(dp), parameter :: nf0 = u_CMB/EphM         !Target photon number density in the obs frame
   
   !===JET PROPERTIES============================================================================================
   
   real(dp), parameter :: gj0 = 1.01d0             !Lorentz factor at the base of the jet  (Algaba+2021: Gamma ~ 1-3)
   real(dp), parameter :: vj0  = sqrt(gj0**2-1d0)*c/gj0 !Fluid bulk speed at the base of the jet
   real(dp) :: zmax_inj = 1e4*R_g                  !Jet maximum cell of injection
   real(dp), parameter :: r_br = 23.0 * 2*R_g
   real(dp), parameter :: z_br = 70 * 2*R_g
   real(dp) :: z0 = 30d0 * R_g
   real(dp), parameter :: z_acc = 2d2 * R_g        !Acceleration region starts at z_acc (Algaba+2021)
   !real(dp), parameter :: geomIndex_out = 0.6d0        !Jet opening parameter (Kim+2018); for r > 10^3 R_S see Mertens+2016
   real(dp) :: geomIndex_out = 0.56d0                    !Jet opening parameter (Kim+2018); for r > 10^3 R_S see Mertens+2016
   real(dp) :: rout_0                                    !Jet launching radius (Kim+2018) - computed at runtime in initialize_allocate
   real(dp), parameter :: rin_0 = 0.0*R_g             
   real(dp), parameter :: geomIndex_in = 0.6d0        

   !===NONTHERMAL PARTICLES======================================================================================

   !Particle injection
   real(dp), parameter :: LnTpFF = 1d43
   real(dp) :: facc = 1d-2                      !Acceleration efficiency
   real(dp), parameter :: EminFF = 1d6*eV       !Minimum electron energy in the FF  !~~~~
   real(dp), parameter :: Emin_injFF = 2d6*eV   !Minimum electron energy of injected particles in the FF
                                                !~~~~
   real(dp), parameter :: EmaxFF = 5d15*eV      !Maximum electron energy in the FF
                                                !~~~~
   real(dp) :: fnT = 1.d-1                      !Nonthermal energy fraction
   real(dp), parameter :: EpminFF = 1d9*eV      !Minimum proton energy in the FF
   real(dp), parameter :: Epmin_injFF = 2d9*eV  !Minimum proton energy of injected particles in the FF
   real(dp), parameter :: EpmaxFF = 1d18*eV     !Maximum proton energy in the FF
   real(dp) :: p = 2.05d0                        !Power law exponent for Q(E) = cte * E**-p

   !===FREE PARAMETERS=========================================================================================

   real(dp) :: Lj = 1d43                        !Jet luminosity (erg/s) (Algaba+2021? arbitrary)
   real(dp) :: sigma_0 = 100d0                    !Initial Magnetization, sigma = B'^2/ (4pi rho' c^2 h')
   real(dp) :: mu  != gj0 * (10d0 + 1d0)!= gj0 * (sigma_0 + 1d0)    !!  !Total energy flux per unit of rest-mass energy flux
   real(dp) :: Mdot_j  != Lj / (c**2 * gj0 * (sigma_0 + 1d0)) !! !Mass flux
   real(dp) :: a_inj = 0.01d0                    !Power law exponent for dLe_dz(z) ~ z**(-(1 + a_inj)),
   real(dp) :: zmax = 1.2d6*R_g           !Jet maximum height !~~~~
                                                !a_inj =/ 0 (else error!)
   real(dp) :: gammaIndex = 0.0016              !Power law exponent for Gamma(z) = (z/z0)**gammaIndex

   !===GLOBAL VARIABLES (TIED TO FREE PARAMETERS)================================================================

   !real(dp) :: B0 != sqrt(fB*8d0*Lj/(vj*(thetaj*z0)**2))  !Magnetic field at the jet base
   real(dp) :: LnT != Lj*fnT              !Available nonthermal luminosity
   real(dp) :: LnTFF != LnT!/gj**2        !Available nonthermal luminosity in the FF
   real(dp) :: max_enthalpy
   !max_enthalpy = 1d0 + LnT / (mu/1.12d0 * Mdot_j * c**2)
   real(dp) :: dlogz != (zmax/z0)**(1d0/(Ns-1d0))
   ! L(z) = cte_inj * z**-a_inj --> LnTFF = integral L(z) dz = log(zmax_inj/z0) iif a_inj == 1.d0, else: 
   real(dp) :: cte_inj != LnTFF/( (zmax_inj**(-a_inj+1.) - z0**(-a_inj+1.))/(-a_inj+1.) )
   real(dp) :: dt != dl/vj           !Segment advection time for constant velocity
   real(dp) :: dtFF != dt/gj         !Advection time in the FF

   !==============================================================================================================

   contains

      !*************************************************************************************
      !   Jet outer radius as a function of z
      !*************************************************************************************
      real(dp) function r_fun_out(zz)
         real(dp), intent(in) :: zz
         real(dp) :: r_fun
         !r_fun1 = rout_0 * (abs(zz)/z0)**geomIndex_out
         !r_fun2 = 0.0 !40 * (abs(zz)/z0)**0.15
          !r_fun1 = r_br * (abs(zz)/z_br)**0.76
          r_fun = rout_0 * (abs(zz)/z0)**geomIndex_out
          r_fun_out = r_fun
          !r_fun_out = min(r_fun1, r_fun2)
      end function

      !*************************************************************************************
      !   Jet inner radius as a function of z
      !*************************************************************************************
      real(dp) function r_fun_in(zz)
         real(dp), intent(in) :: zz
         real(dp) :: r_fun1, r_fun2
         if (inner_radius_calc) then
            r_fun1 = rin_0 * (abs(zz)/z0)**geomIndex_in
         else
            r_fun1 = 0.0
         end if
         !r_fun1 = rin_0 * (abs(zz)/z0)**geomIndex_in
         r_fun2 = 0.0 !40 * (abs(zz)/z0)**0.15
         r_fun_in = max(r_fun1, r_fun2)
      end function


      !*************************************************************************************
      !   Jet Lorentz factor as a function of z
      !*************************************************************************************
      real(dp) function lorentzFactor(zz)
         real(dp), intent(in) :: zz
         !lorentzFactor = sqrt(1d0 + betaGamma(abs(zz))**2)
         !lorentzFactor = max(gj0, (abs(zz)/(2e2*R_g))**gammaIndex)
         lorentzFactor = (gj0 + (abs(zz)/z_acc))**gammaIndex
         !lorentzFactor = min(lorentzFactor, 0.99*mu/max_enthalpy) 
         ! The factor 0.99 is included to avoid numerical errors if B = 0
      end function

      !*************************************************************************************
      !   Jet normalized velocity as a function of z
      !*************************************************************************************
      real(dp) function betaJet(zz)
         real(dp), intent(in) :: zz
         !betaJet = betaGamma(zz) / lorentzFactor(zz)
         betaJet = sqrt(1d0 - 1d0/lorentzFactor(abs(zz))**2)
      end function

      !*************************************************************************************
      !   Jet Beta times Lorentz factor as a function of z
      !*************************************************************************************
      real(dp) function betaGamma(zz)
         real(dp), intent(in) :: zz
        ! betaGamma = 0.37 * (abs(zz)/R_g)**gammaIndex
         betaGamma = betaJet(abs(zz)) * lorentzFactor(abs(zz))
      end function

      !*************************************************************************************
      ! Power injected in electrons per unit zz
      !*************************************************************************************
      real(dp) function dLe_dz(zz)
         real(dp), intent(in) :: zz
         real(dp) L0
         L0 = LnT * a_inj / (1d0 - (z0/zmax_inj)**a_inj)
         dLe_dz = L0/z0 * (abs(zz)/z0)**(-(1d0+a_inj))
      end function

      !*************************************************************************************
      ! Power injected in electrons up to a height zz
      !*************************************************************************************
      real(dp) function Le_z(zz)
         real(dp), intent(in) :: zz
         if (zz < zmax_inj) then
            Le_z = LnT * (1d0 - (z0/abs(zz))**a_inj)/(1d0 - (z0/zmax_inj)**a_inj)
         else
            Le_z = LnT
         end if
      end function

      !*************************************************************************************
      ! Power in cold matter
      !*************************************************************************************
      real(dp) function Lp(zz)
         real(dp), intent(in) :: zz
         Lp = lorentzFactor(abs(zz)) * Mdot_j * c**2
      end function

      !*************************************************************************************
      ! Specific enthalpy
      !*************************************************************************************
      real(dp) function enthalpy(zz)
         real(dp), intent(in) :: zz
         enthalpy = 1d0 + Le_z(zz)/Lp(zz)
      end function

      !*************************************************************************************
      ! Magnetization
      !*************************************************************************************
      real(dp) function sigma_B(zz)
         real(dp), intent(in) :: zz
         !write(*,*) "ent_ratio = ", enthalpy(abs(zz))/max_enthalpy
         sigma_B = max(mu / (lorentzFactor(abs(zz)) * enthalpy(abs(zz))) - 1d0, 1d-10)
      end function

      !*************************************************************************************
      !   Jet FF magnetic field at height z
      !*************************************************************************************
      real(dp) function BFF_fun(zz)
         real(dp), intent(in) :: zz
         real(dp) :: rr, Gamma, sB
         rr = sqrt(r_fun_out(zz)**2-r_fun_in(zz)**2)
         Gamma = lorentzFactor(abs(zz))
         sB = sigma_B(abs(zz))
         BFF_fun = 2d0/(Gamma*rr) * sqrt( Lj/(betaJet(abs(zz))*c) * sB/(sB+1d0))
      end function

      !*************************************************************************************
      !  Jet magnetic power at height z
      !*************************************************************************************
      real(dp) function L_B(zz)
          real(dp), intent(in) :: zz
          real(dp) :: rr, Bff
          rr = sqrt(r_fun_out(zz)**2-r_fun_in(zz)**2)
          Bff = BFF_fun(zz)
          L_B = (Bff**2/(4d0*pi)) * (pi*rr**2) * betaJet(abs(zz)) * c * lorentzFactor(abs(zz))**2
        end function
      
      !*************************************************************************************
      !   Jet lateral velocity at height z ( = dr/dz * vj)
      !*************************************************************************************
      real(dp) function v_lat(zz)
         real(dp), intent(in) :: zz
         v_lat = betaJet(zz) * c * (rout_0/z0) * geomIndex_out * (z0/abs(zz))**(1d0-geomIndex_out)
      end function

   !============================================================================

END module system_parameters
