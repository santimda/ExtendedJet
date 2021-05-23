MODULE evolution
   use globals
   implicit none

   contains

   !============================================================================
   ! Quantities (START)
   !============================================================================
   subroutine Quantities(jet)
   !============================================================================
   ! This subroutine calculates the evolution of different quantities 
   ! from cell to cell
   !============================================================================
      character(LEN=*), intent(in) :: jet ! 'cj' for counterjet or '' for jet
      real(dp) :: B, dir, cjo, boost
      allocate( z(Ns), rj(Ns) )
      allocate( delta_obs(Ns) )
      allocate( deltaZ(Ns) )
      allocate( deltatFF(Ns) )
      allocate( BFF(Ns), vnorm(Ns), VFF(Ns), nfFF(Ns) )
      allocate( EphMFF(Ns), Kappa(Ns), joFF(Ns), sjoFF(Ns) )

      open(10,file='trajectory'//jet//'.dat')
      open(20,file='deltas'//jet//'.dat')
      open(30,file='quant_FF'//jet//'.dat')
      open(31,file='quant_FF_norm'//jet//'.dat')

      ! The jet goes upwards, the CJ goes downwards (negative v_z)
      if (jet == '') then
         dir = 1.d0
      elseif (jet == 'cj') then
         dir = -1.d0
      else
         write(*,*) 'Invalid jet type'
         stop
      endif

      cjo = cos(jo)
      do i = 1, Ns

         read(10,*) z(i), rj(i)

         if (z(i) < zmax_inj) Ninj = i

         deltaZ(i) = z(i) * (sqrt(dlogZ)-1d0/sqrt(dlogZ))
         deltatFF(i) = deltaZ(i) / (vj*gj)

         !Doppler factor between the jet and the observer
         delta_obs(i) = 1d0/(gj*(1d0-dir*Bj*cjo))

         !Magnetic field
         B = B0*(rj(1)/rj(i))**exp_B
         BFF(i) = B/gj

         !Stellar photon density and energy 
         nfFF(i) = nf0/delta_obs(i)  !!! CHECK!!!!!
         EphMFF(i) = EphM/delta_obs(i)

         !Dilution factor in the FF (Khangulyan 2014)
         Kappa(i) = delta_obs(i)**2 ! CHECK, FOR CMB

         !Velocity perpendicular to the jet propagation
         if (i == 1) then
            vnorm(i) = vj*tan(thetaj)   !Conical expansion???
         else
            vnorm(i) = (rj(i)-rj(i-1))/(deltatFF(i)*gj)
         endif

         !Jet observer angle in the FF
         joFF(i) = acos( (cos(jo)-Bj)/(1d0-Bj*cos(jo)) )
         sjoFF(i) = sin(joFF(i))        !Sinus of joFF

         !Segment volume in the FF
         VFF(i) = pi*rj(i)**2*deltatFF(i)*vj

         !boost = delta_obs(i)**2/gj  ! CHECK THIS!
         boost = delta_obs(i)**2  ! CHECK THIS! (delta z defined in LF)

         write(20,*) delta_obs(i), boost
         write(30,*) z(i), BFF(i), rj(i), vnorm(i), Kappa(i), sjoFF(i), VFF(i), EphMFF(i), nfFF(i)
         write(31,*) z(i)/R_g, (BFF(i)/BFF(1))**2, (vnorm(i)/rj(i))/(vnorm(1)/rj(1)), VFF(i)/VFF(1)

      end do

      close(10)
      close(20)
      close(30)
      close(31)

      deallocate( z, rj )
      deallocate( delta_obs )
      deallocate( BFF, vnorm, VFF, nfFF ) 
      deallocate( EphMFF, Kappa, joFF, sjoFF )

   end subroutine Quantities
   !============================================================================
   ! Quantities (END)
   !============================================================================

END MODULE evolution
