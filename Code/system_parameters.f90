module system_parameters
   use constants
   use controls
   implicit none
   public

   !====================================================================
   ! This module contains every model parameter (jet + RIAF-adjacent
   ! constants) as well as the run-time controls formerly hard-coded in
   ! globals.f90/controls.f90. All physical parameters are plain module
   ! variables (not PARAMETER) so that read_namelist() can override them;
   ! their defaults below are the M87 Jet+RIAF values (A's
   ! system_parameters_M87.f90). Quantities that are derived from other
   ! parameters are recomputed by derive_parameters(), never read directly
   ! from a namelist.
   !====================================================================

   !===SYSTEM PARAMETERS=========================================================================================

   real(dp) :: dist_mpc = 16.8d0                   !Distance to the system [Mpc]
   real(dp) :: incl_deg = 17.d0                     !Observer's angle with the jet axis [deg]
   real(dp) :: M_BH = 6.5d9                         !Black hole mass [Msun]

   real(dp) :: dist                                 !Distance to the system [cm] (derived)
   real(dp) :: jo                                   !Observer's angle with the jet axis [rad] (derived)
   real(dp) :: L_Edd                                !Eddington luminosity (derived)
   real(dp) :: R_g                                  !Gravitational radius [cm] (derived)
   real(dp), parameter :: EphM = kB*T_CMB           !Target photon energy in the obs frame
   real(dp), parameter :: nf0 = u_CMB/EphM          !Target photon number density in the obs frame

   !===JET PROPERTIES============================================================================================

   real(dp) :: gj0 = 1.01d0                         !Lorentz factor at the base of the jet  (Algaba+2021: Gamma ~ 1-3)
   real(dp) :: vj0                                  !Fluid bulk speed at the base of the jet (derived)

   ! Lengths in the namelist are given in units of R_g; the physical (cm)
   ! values below are computed by derive_parameters()/read_namelist().
   real(dp) :: z0_rg = 50.d0                        !Jet launching height [R_g]
   real(dp) :: zmax_rg = 1.2d6                      !Jet maximum height [R_g]
   real(dp) :: zmax_inj_rg = 1.d5                   !Jet maximum cell of injection [R_g]
   real(dp) :: r_br_rg = 46.d0                      !=23*2, jet break radius [R_g]
   real(dp) :: z_br_rg = 140.d0                     !=70*2, jet break height [R_g]
   real(dp) :: z_acc_rg = 2.d2                      !Acceleration region starts at z_acc [R_g] (Algaba+2021)
   real(dp) :: rin_0_rg = 0.d0                      !Jet inner launching radius [R_g]

   real(dp) :: z0                                   !Jet launching height [cm] (derived)
   real(dp) :: zmax                                 !Jet maximum height [cm] (derived)
   real(dp) :: zmax_inj                             !Jet maximum cell of injection [cm] (derived)
   real(dp) :: r_br                                 !Jet break radius [cm] (derived)
   real(dp) :: z_br                                 !Jet break height [cm] (derived)
   real(dp) :: z_acc                                !Acceleration region start [cm] (derived)
   real(dp) :: rin_0                                !Jet inner launching radius [cm] (derived)

   real(dp) :: geomIndex_out = 0.56d0               !Jet opening parameter (Kim+2018); for r > 10^3 R_S see Mertens+2016
   real(dp) :: rout_0                               !Jet launching (outer) radius [cm] (derived)
   real(dp) :: geomIndex_in = 0.6d0

   !===NONTHERMAL PARTICLES======================================================================================

   !Particle injection
   real(dp) :: LnTpFF = 1d43
   real(dp) :: facc = 1d-2                          !Acceleration efficiency
   real(dp) :: EminFF_eV = 1d6                      !Minimum electron energy in the FF [eV]
   real(dp) :: Emin_injFF_eV = 2d6                  !Minimum electron energy of injected particles in the FF [eV]
   real(dp) :: EmaxFF_eV = 5d15                     !Maximum electron energy in the FF [eV]
   real(dp) :: fnT = 1.d-1                          !Nonthermal energy fraction
   real(dp) :: EpminFF_eV = 1d9                     !Minimum proton energy in the FF [eV]
   real(dp) :: Epmin_injFF_eV = 2d9                 !Minimum proton energy of injected particles in the FF [eV]
   real(dp) :: EpmaxFF_eV = 1d18                    !Maximum proton energy in the FF [eV]
   real(dp) :: p = 2.05d0                            !Power law exponent for Q(E) = cte * E**-p

   real(dp) :: EminFF, Emin_injFF, EmaxFF           !Electron energy bounds in the FF [erg] (derived)
   real(dp) :: EpminFF, Epmin_injFF, EpmaxFF        !Proton energy bounds in the FF [erg] (derived)

   !===FREE PARAMETERS=========================================================================================

   real(dp) :: Lj = 1.2d43                          !Jet luminosity (erg/s) (Algaba+2021? arbitrary)
   real(dp) :: sigma_0 = 100d0                       !Initial Magnetization, sigma = B'^2/ (4pi rho' c^2 h')
   real(dp) :: mu                                    !Total energy flux per unit of rest-mass energy flux (derived)
   real(dp) :: Mdot_j                                !Mass flux (derived)
   real(dp) :: a_inj = 0.01d0                        !Power law exponent for dLe_dz(z) ~ z**(-(1 + a_inj))
                                                     !a_inj =/ 0 (else error!)
   real(dp) :: gammaIndex = 0.16d0                   !Power law exponent for Gamma(z) = (z/z0)**gammaIndex

   !===GLOBAL VARIABLES (TIED TO FREE PARAMETERS)================================================================
   ! All recomputed by derive_parameters(); never set directly.

   real(dp) :: LnT                                   !Available nonthermal luminosity
   real(dp) :: LnTFF                                 !Available nonthermal luminosity in the FF
   real(dp) :: max_enthalpy
   real(dp) :: dlogz
   real(dp) :: cte_inj
   real(dp) :: taper_norm = 1d0                    !Rescaling of the injection so that the exp(-z/zmax_inj)-tapered power integrates to LnT (fix_norm_taper)
   real(dp) :: dt
   real(dp) :: dtFF

   !===RUN-TIME CONTROLS (namelist &run)==========================================================================
   ! These replace the old compile-time PARAMETERs in globals.f90.

   integer :: Ns = 200                               !Number of jet segments, log-spaced from z0 to zmax
   integer :: numbins = 500                          !Number of integration bins for the energy distribution
   real(dp) :: dttFF_min = 1.d-6                     !Precision in timestep for N(E)

   integer, parameter :: max_apertures = 10
   real(dp) :: apertures_rg(max_apertures) = (/ 7.3d0, 26.d0, 65.d0, 195.d0, 650.d0, 3.5d5, &
                                                 0.d0, 0.d0, 0.d0, 0.d0 /)
   integer :: n_ap = 6                               !Number of apertures actually used (<= max_apertures)

   character(len=8) :: sync_approx = 'full'          !'full' (A's jsync/alpha_sync integrals) or 'delta' (B's override)

   character(len=512) :: outdir = '.'                !All output files are written under this directory
   character(len=512) :: scriptdir = '.'             !Directory holding SortNe.sh etc. (used if sort_outputs)

   !==============================================================================================================

   contains

      !*************************************************************************************
      !   Read the &model and &run namelist groups from a file, then derive everything else.
      !*************************************************************************************
      subroutine read_namelist(path)
         character(len=*), intent(in) :: path
         integer :: lun, ios
         character(len=512) :: iomsg_

         namelist /model/ M_BH, dist_mpc, incl_deg, gj0, &
                           z0_rg, zmax_rg, zmax_inj_rg, r_br_rg, z_br_rg, z_acc_rg, rin_0_rg, &
                           geomIndex_out, geomIndex_in, &
                           LnTpFF, facc, EminFF_eV, Emin_injFF_eV, EmaxFF_eV, fnT, &
                           EpminFF_eV, Epmin_injFF_eV, EpmaxFF_eV, p, &
                           Lj, sigma_0, a_inj, gammaIndex

         namelist /run/ outdir, scriptdir, Ns, numbins, dttFF_min, apertures_rg, n_ap, sync_approx, &
                         JET_, CJ_, NE_calc, Np_calc, radsyn_calc, radpsyn_calc, radic_calc, &
                         print_tlosses, sort_outputs, fix_norm_taper, fix_inj_emin

         open(newunit=lun, file=trim(path), status='old', action='read', iostat=ios, iomsg=iomsg_)
         if (ios /= 0) then
            print '(A,A,A,A)', 'FATAL ERROR: cannot open namelist file "', trim(path), '": ', trim(iomsg_)
            stop 1
         end if

         read(lun, nml=model, iostat=ios, iomsg=iomsg_)
         if (ios /= 0) then
            print '(A,A,A,A)', 'FATAL ERROR: malformed or unknown entry in &model of "', trim(path), '": ', trim(iomsg_)
            stop 1
         end if

         rewind(lun)
         read(lun, nml=run, iostat=ios, iomsg=iomsg_)
         if (ios /= 0) then
            print '(A,A,A,A)', 'FATAL ERROR: malformed or unknown entry in &run of "', trim(path), '": ', trim(iomsg_)
            stop 1
         end if
         close(lun)

         if (n_ap > max_apertures) then
            print '(A,I0)', 'FATAL ERROR: n_ap exceeds max_apertures = ', max_apertures
            stop 1
         end if

         call execute_command_line('mkdir -p '//trim(outdir))

         call derive_parameters()

      end subroutine read_namelist

      !*************************************************************************************
      !   Set (and create) the output directory. Exposed to f2py so that several Python
      !   processes sharing the same in-process module state can each redirect outputs.
      !*************************************************************************************
      subroutine set_outdir(path)
         character(len=*), intent(in) :: path
         outdir = trim(path)
         call execute_command_line('mkdir -p '//trim(outdir))
      end subroutine set_outdir

      !*************************************************************************************
      !   Prefix a bare filename with the current output directory.
      !*************************************************************************************
      function opath(fname) result(full)
         character(len=*), intent(in) :: fname
         character(len=1024) :: full
         full = trim(outdir)//'/'//trim(fname)
      end function opath

      !*************************************************************************************
      !   Recompute every quantity that is tied to the free parameters above. Called after
      !   the namelist is read, and again from initialize_allocate/syn_model/power_model
      !   whenever those routines change a parameter directly (f2py entry points keep their
      !   original argument-based interface, so they update the physical (cm/erg) variables
      !   in place rather than through the *_rg / *_eV namelist mirrors).
      !*************************************************************************************
      subroutine derive_parameters()
         L_Edd = 1.25d38*M_BH
         R_g = 1.477d5*M_BH
         dist = dist_mpc*Mpc
         jo = incl_deg*pi/180.d0

         z0 = z0_rg*R_g
         zmax = zmax_rg*R_g
         zmax_inj = zmax_inj_rg*R_g
         r_br = r_br_rg*R_g
         z_br = z_br_rg*R_g
         z_acc = z_acc_rg*R_g
         rin_0 = rin_0_rg*R_g

         EminFF = EminFF_eV*eV
         Emin_injFF = Emin_injFF_eV*eV
         EmaxFF = EmaxFF_eV*eV
         EpminFF = EpminFF_eV*eV
         Epmin_injFF = Epmin_injFF_eV*eV
         EpmaxFF = EpmaxFF_eV*eV

         vj0 = sqrt(gj0**2-1d0)*c/gj0
         rout_0 = r_br * (z0/z_br)**geomIndex_out

         LnT = Lj*fnT
         LnTFF = LnT

         mu = gj0 * (sigma_0 + 1d0)
         Mdot_j = Lj / (c**2 * gj0 * (sigma_0 + 1d0))
         max_enthalpy = 1d0 + LnT / (lorentzFactor(zmax_inj)*Mdot_j * c**2)

         cte_inj = LnTFF/( (zmax_inj**(-a_inj+1.) - z0**(-a_inj+1.))/(-a_inj+1.) )

         dlogz = (zmax/z0)**(1d0/(Ns-1d0))

         taper_norm = 1d0
         if (fix_norm_taper) taper_norm = taper_correction()

      end subroutine derive_parameters

      !*************************************************************************************
      !   Injection is multiplied by exp(-z/zmax_inj) (edistribution: const_zjz) but dLe_dz is
      !   normalised without it, so the injected power falls short of LnT. This returns the
      !   factor LnT / int_{z0}^{zmax} dLe_dz(z) exp(-z/zmax_inj) dz that restores it.
      !*************************************************************************************
      real(dp) function taper_correction()
         integer, parameter :: n = 4000
         real(dp) :: lz, dlz, z1, z2, f1, f2, tot
         integer :: i
         dlz = log(zmax/z0)/n
         tot = 0d0
         z1 = z0
         f1 = dLe_dz(z1)*z1*exp(-z1/zmax_inj)
         do i = 1, n
            z2 = z0*exp(i*dlz)
            f2 = dLe_dz(z2)*z2*exp(-z2/zmax_inj)
            tot = tot + 0.5d0*(f1+f2)*dlz
            f1 = f2
         end do
         taper_correction = LnT/tot
      end function taper_correction

      !*************************************************************************************
      !   Jet outer radius as a function of z
      !*************************************************************************************
      real(dp) function r_fun_out(zz)
         real(dp), intent(in) :: zz
         real(dp) :: r_fun
         r_fun = rout_0 * (abs(zz)/z0)**geomIndex_out
         r_fun_out = r_fun
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
         r_fun2 = 0.0
         r_fun_in = max(r_fun1, r_fun2)
      end function


      !*************************************************************************************
      !   Jet Lorentz factor as a function of z
      !*************************************************************************************
      real(dp) function lorentzFactor(zz)
         real(dp), intent(in) :: zz
         lorentzFactor = (gj0 + (abs(zz)/z_acc))**gammaIndex
      end function

      !*************************************************************************************
      !   Jet normalized velocity as a function of z
      !*************************************************************************************
      real(dp) function betaJet(zz)
         real(dp), intent(in) :: zz
         betaJet = sqrt(1d0 - 1d0/lorentzFactor(abs(zz))**2)
      end function

      !*************************************************************************************
      !   Jet Beta times Lorentz factor as a function of z
      !*************************************************************************************
      real(dp) function betaGamma(zz)
         real(dp), intent(in) :: zz
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
