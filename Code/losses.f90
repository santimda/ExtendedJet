MODULE losses
    use constants
    use globals
    implicit none

contains

    !============================================================================
    !   Total losses
    !============================================================================
    real(dp) function TLosses(B_FF, z_, Eph_MFF, Kapp, EFF, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: B_FF, Eph_MFF, Kapp, EFF, z_

    TLosses = SyncLoss(B_FF, EFF, particle) + AdLoss(z_, EFF)! + ICLoss(Eph_MFF,Kapp,EFF)

    end function TLosses

    !============================================================================
    !   Total losses (modified)
    !============================================================================
    real(dp) function TLosses_mod(B_FF, z_, Eph_MFF, Kapp, EFF, Ng_, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: B_FF, Eph_MFF, Kapp, EFF, z_, Ng_

    TLosses_mod = SyncLoss_mod(B_FF, EFF, z_, Ng_, particle) + AdLoss(z_, EFF)! + ICLoss(Eph_MFF,Kapp,EFF)

    end function TLosses_mod

    !============================================================================
    !   Synchrotron losses
    !============================================================================
    real(dp) function SyncLoss(B_FF, EFF, particle)
    real(dp), intent(in) :: B_FF, EFF
    character(1), intent(in) :: particle
    real(dp) :: Lgamma, factor

    Lgamma = EFF/mec2    !Lorentz factor of the particles in the fluid frame
    factor = 1d0
    if (particle == 'p') then
        factor = 2.966d-7    ! (me/mp)^2
        Lgamma = EFF/mpc2
    endif
    SyncLoss = - 1.05804d-15 * factor * B_FF**2 * (Lgamma**2 - 1d0)   !4/9*c*re^2 * B^2 * gamma^2 [erg/s]

    end function SyncLoss

    !============================================================================
    !   Synchrotron losses (modified)
    !============================================================================
    real(dp) function SyncLoss_mod(B_FF, EFF, zz_, Ng_, particle)
    real(dp), intent(in) :: B_FF, EFF, zz_, Ng_
    character(1), intent(in) :: particle
    real(dp) :: Lgamma, factor, alpha, B_cr, eps, tau

    Lgamma = EFF/mec2    !Lorentz factor of the particles in the fluid frame
    factor = 1d0
    if (particle == 'p') then
        factor = 2.966d-7    ! (me/mp)^2
        Lgamma = EFF/mpc2
    endif
    B_cr = 4.4e13
    eps = B_FF / B_cr * Lgamma**2
    alpha = 137d0 / 2d0 * pi * sigma_T * B_FF/B_cr / eps**2 * Ng_
    tau = alpha * r_fun_out(zz_)
    SyncLoss_mod = -1.05804d-15 * factor * B_FF**2 * (Lgamma**2-1d0) / (1d0 + tau)  !4/9*c*re^2 * B^2 * gamma^2 [erg/s]
 
    end function SyncLoss_mod
    
    !============================================================================
    !   Inverse Compton losses (Eq. (36) in Khangulyan 2014)                    
    !============================================================================
    real(dp) function ICLoss(Eph_MFF, Kapp, EFF)
    real(dp), intent(in) :: Eph_MFF, Kapp, EFF
    real(dp) :: TFF, t_iso, Eiso_ct

    TFF = Eph_MFF/(2.7d0*mec2)   !Dimensionless stellar temperature in the FF (EphMFF = 2.7*kB*TFF)
    t_iso = 4d0*TFF*(EFF/mec2)
    Eiso_ct = cte_IC_K*Kapp*TFF**2
    ICLoss = - Eiso_ct*F_iso(t_iso)*mec2    ![erg/s]

    end function ICLoss

    !============================================================================
    !   Adiabatic losses
    !============================================================================
    real(dp) function AdLoss(z_, EFF)
    real(dp), intent(in) :: z_, EFF
    real(dp) :: gg,zabs

    gg = EFF - mec2**2/EFF
    zabs = abs(z_)
    AdLoss = - 2d0/3d0 * geomIndex_out * lorentzFactor(zabs) * c * betaJet(zabs) / zabs * gg
    
    end function AdLoss
        
    !============================================================================
    !   Cell convection time
    !============================================================================
    real(dp) function t_cell(segN)
    integer, intent(in) :: segN
    if (segN==Ns) then
        t_cell = abs(z(segN) - z(segN-1))/(c*bj(segN))/gj(segN)
    else   
        t_cell = abs(z(segN+1) - z(segN))/(c*bj(segN))/gj(segN)
    endif
    
    end function t_cell
        
    !============================================================================
    !   Save to a file the different cooling times for a given segment
    !============================================================================
    subroutine Print_tloss(segN, particle)
    integer, intent(in) :: segN
    character(1), intent(in) :: particle
    real(dp) :: EFF, Eint, tSyn, tIC, tAd, tTot, tacc, tdiff, tcell, tconv
    character(3) :: num_seg
    
    !Initialization
    EFF = EminFF
    Eint = (EmaxFF/EminFF)**(1d0/(nIn-1))
    
    write(num_seg,'(I3)') segN
    open(20,file='losses_'//trim(adjustl(num_seg))//'.dat')  
    tcell = t_cell(segN)
    tconv = abs(z(segN))/(c*bj(segN))/gj(segN)

    do i = 1, nIn
        !Cooling time for each process for segment segN
        tSyn = EFF/abs(SyncLoss(BFF(segN),EFF,particle))    
        tIC = EFF/abs(ICLoss(EphMFF(segN),Kappa(segN),EFF))
        tAd = EFF/abs(AdLoss(z(segN),EFF))
        tTot = EFF/abs(TLosses(BFF(segN),z(segN),EphMFF(segN),Kappa(segN),EFF,particle))
        tacc = EFF/(facc*BFF(segN)*qe*c)                    !Acceleration timescale
        tdiff = 3d0*qe*BFF(segN)*rj_out(segN)**2/(2d0*c*EFF)    !Diffusion timescale
        write(20,*) EFF/eV, tSyn, tIC, tAd, tTot, tcell, tacc, tdiff, tconv
                    
        EFF = EFF*Eint  !Energy increase   
    enddo
    close(20)

    end subroutine Print_tloss
  
    !============================================================================

END MODULE losses
