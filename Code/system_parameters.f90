module system_parameters

   use constants
   implicit none
   public

   !====================================================================
   ! This module contains the parameters of the specific system modelled
   !====================================================================

   !===JET PROPERTIES============================================================================================

   real(dp), parameter :: gj = 3.d0             !Lorentz factor
   real(dp), parameter :: vj = sqrt(gj**2-1d0)*c/gj        !Fluid bulk speed
   real(dp), parameter :: Bj = vj/c              !Beta
   real(dp), parameter :: thetaj = 0.1d0         !Half opening angle (rad)  
   real(dp), parameter :: Lj = 5d36              !Jet luminosity (erg/s)
   real(dp), parameter :: fnT = 1d-1             !Nonthermal energy fraction
   real(dp), parameter :: LnT = Lj*fnT           !Available nonthermal luminosity
   real(dp), parameter :: LnTFF = LnT/gj**2      !Available nonthermal luminosity in the FF
   real(dp), parameter :: Pdotj = Lj*gj*Bj/(c*(gj-1d0))    !Jet lateral momentum flux

   !===SYSTEM PARAMETERS=========================================================================================

   real(dp), parameter :: dist = 3d0*kpc           !Distance to the system
   real(dp), parameter :: inc = 30d0               !Observer's angle with the z axis (0-90º) 
   real(dp), parameter :: azi = 90d0               !Observer's angle with the x axis (0-360º)
   real(dp), parameter :: dl = 1d11                !Length of jet segment/blob (cm)
   real(dp), parameter :: lf = 50d0*a              !Final jet length
   integer, parameter :: Ns = int(lf/dl)           !Number of jet segments
   real(dp), parameter :: dt = dl/vj               !Segment advection time for constant velocity
   real(dp), parameter :: dtFF = dt/gj             !Advection time in the FF
   real(dp), parameter :: &                        !Height of the recol. shock (B-R 2016)
   & zrec = a*sqrt((PdotW_CO*gj**2-8d0*Pdotj-sqrt(PdotW_CO*gj**2*(PdotW_CO*gj**2-16d0*Pdotj)))/(8d0*Pdotj))
   real(dp), parameter :: nf0 = Lstar/(4d0*pi*(a**2+zrec**2)*c*EphM) !Target photon density at the jet base
   real(dp), parameter :: fB = 1.0d0                 !Fraction of magnetic energy
   real(dp), parameter :: B0 = sqrt(fB*8d0*pi*EphM*nf0)        !Magnetic field at the jet base
   !real(dp), parameter :: B0 = sqrt(fB*8d0*Lj/(vj*(thetaj*zrec)**2))  !Magnetic field at the jet base
   real(dp), parameter :: Xj = thetaj*PdotW_CO/(4d0*pi*Pdotj)      !Bosch-Ramon 2016

   !===NONTHERMAL PARTICLES======================================================================================

   !Particle injection
   real(dp), parameter :: facc = 1d-1      !Acceleration efficiency
   real(dp), parameter :: EminFF = 1d6*eV    !Minimum energy in the FF
   real(dp), parameter :: EmaxFF = 1d14*eV   !Maximum energy in the FF
   real(dp), parameter :: p = 2d0        !Power law exponent

   !Particle distribution
   real(dp), parameter :: E1 = EminFF  !Minimum energy considered in the FF
   real(dp), parameter :: E2 = EmaxFF  !Maximum energy considered in the FF
   integer, parameter :: nNE = 50  !Energy steps (for plotting N(E))
   integer, parameter :: nIn = 100 !Energy steps (for the energy integrals)

   !============================================================================


END module system_parameters
