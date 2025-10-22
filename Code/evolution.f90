MODULE evolution
   use constants
   use globals
   implicit none

   contains

   !============================================================================
   ! Quantities (START)
   !============================================================================
   subroutine Quantities(jet)
   !============================================================================
   ! This subroutine calculates the evolution of different quantities 
   ! from cell to cell. 
   ! The input is the '' (jet) or 'CJ' (counter-jet). 
   !============================================================================
      character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet
      real(dp) :: dir, cjo
      allocate( joFF(Ns) )

      open(30,file='quant_FF'//jet//'.dat')
      open(31,file='quant_FF_norm'//jet//'.dat')

      ! The jet goes upwards, the CJ goes downwards (negative v_z)
      if (jet == '') then
         dir = 1.d0
      elseif (jet == 'CJ') then
         dir = -1.d0
      else
         print *, 'Invalid jet type'
         stop
      endif

      cjo = cos(jo)
      do i = 1, Ns

         if (abs(z(i)) < zmax_inj) Ninj = i

         gj(i) = lorentzFactor(z(i))
         bj(i) = betaJet(z(i))
         deltaZ(i) = abs( z(i)) * (sqrt(dlogZ) - 1d0/sqrt(dlogZ) ) !abs value for CJ
         deltatFF(i) = deltaZ(i) / ( c*bj(i)*gj(i) )
         !print*, (deltaZ(i)/gj)/rj(i), z(i)/rj(i)/gj
         !Doppler factor between the jet and the observer
         delta_obs(i) = 1d0 / ( gj(i) * (1d0 - dir*bj(i)*cjo) )

         !Magnetic field
         BFF(i) = BFF_fun(z(i))

         !Stellar photon density and energy 
         nfFF(i) = nf0/delta_obs(i)  !!! CHECK!!!!!
         EphMFF(i) = EphM/delta_obs(i)

         !Dilution factor in the FF (Khangulyan 2014)
         Kappa(i) = delta_obs(i)**2 ! CHECK, FOR CMB

         !Velocity perpendicular to the jet propagation
         vnorm(i) = v_lat(z(i))

         !Jet observer angle in the FF
         joFF(i) = acos( (cjo - bj(i)) / (1d0 - bj(i)*cjo) )
         sjoFF(i) = sin(joFF(i))        !Sinus of joFF

         !Segment volume in the FF
         VFF(i) = pi*(rj_out(i)**2-rj_in(i)**2)*(deltaZ(i)/gj(i))

         !Doppler boosting factor (used to transform FF -> OF); the expression assumes that delta z is defined in the LF
         boost(i) = delta_obs(i)**2 

         write(30,*) z(i), BFF(i), rj_out(i), rj_in(i), gj(i), delta_obs(i), sjoFF(i), VFF(i), EphMFF(i), nfFF(i), deltaZ(i), deltatFF(i), sigma_B(z(i)), Le_z(z(i))
         write(31,*) z(i)/z(1), BFF(i)/BFF(1), vnorm(i)/vnorm(1), rj_out(i)/rj_out(1), rj_in(i)/rj_in(1), VFF(i)/VFF(1), gj(i)/gj(1), delta_obs(i)/delta_obs(1)

      end do
      close(30)
      close(31)

      deallocate( joFF )

   end subroutine Quantities
   !============================================================================
   ! Quantities (END)
   !============================================================================

END MODULE evolution
