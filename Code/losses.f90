MODULE losses
    use constants
    use globals
    implicit none

contains

    !**********************************************************************************************
    !   Total losses
    !**********************************************************************************************
    real(dp) function TLosses(B_FF, rj, v_norm, Eph_MFF, Kapp, EFF, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: B_FF, rj, v_norm, Eph_MFF, Kapp, EFF

    TLosses = SyncLoss(B_FF,EFF,particle) + AdLoss(rj, v_norm, EFF)! + ICLoss(Eph_MFF,Kapp,EFF)

    return
    end function TLosses

    !**********************************************************************************************
    !   Synchrotron losses
    !**********************************************************************************************
    real(dp) function SyncLoss(B_FF, EFF, particle)
    real(dp), intent(in) :: B_FF, EFF
    character(1), intent(in) :: particle
    real(dp) :: Lgamma, factor

    !Lorentz factor of the particles in the fluid frame
    Lgamma = EFF/mec2
    factor = 1d0
    if (particle == 'p') then
        factor = 2.966d-7    ! (me/mp)^2
        Lgamma = EFF/mpc2
    endif
    !Synchrotron losses in erg/s
    SyncLoss = - 1.05804d-15 * factor * B_FF**2 * Lgamma**2   !4/9*c*re^2 * B^2 * gamma^2

    return
    end function SyncLoss
    
    !**********************************************************************************************
    !   Inverse Compton losses (Eq. (36) in Khangulyan 2014)                    
    !**********************************************************************************************
    real(dp) function ICLoss(Eph_MFF, Kapp, EFF)
    real(dp), intent(in) :: Eph_MFF, Kapp, EFF
    real(dp) :: TFF, t_iso, Eiso_ct

    TFF = Eph_MFF/(2.7d0*mec2)   !Dimensionless stellar temperature in the FF (EphMFF = 2.7*kB*TFF)
    t_iso = 4d0*TFF*(EFF/mec2)
    Eiso_ct = cte_IC_K*Kapp*TFF**2
    
    !IC losses (in erg/s due to the mec2 multiplying)
    ICLoss = - Eiso_ct*F_iso(t_iso)*mec2

    return
    end function ICLoss

    !**********************************************************************************************
    !   Adiabatic losses
    !**********************************************************************************************
    real(dp) function AdLoss(r_j, v_norm, EFF)
    real(dp), intent(in) :: r_j, v_norm, EFF

    !Adiabatic losses
    AdLoss = - 2d0/3d0 * v_norm/r_j * EFF * gj    !E/t = E'/t' = ct

    end function AdLoss
        
    !**********************************************************************************************
    !   Write the different cooling times for a given segment
    !**********************************************************************************************
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
    if (segN==Ns) then
        tcell = (z(segN) - z(segN-1))/vj/gj
    else   
        tcell = (z(segN+1) - z(segN))/vj/gj
    endif
    tconv = z(segN)/vj/gj

    do i = 1, nIn
        !Cooling time for each process for segment segN
        tSyn = EFF/abs(SyncLoss(BFF(segN),EFF, 'e'))
        tIC = EFF/abs(ICLoss(EphMFF(segN),Kappa(segN),EFF))
        tAd = EFF/abs(AdLoss(rj(segN),vnorm(segN),EFF))
        tTot = EFF/abs(TLosses(BFF(segN),rj(segN),vnorm(segN),EphMFF(segN),Kappa(segN),EFF, particle))
        tacc = EFF/(facc*BFF(segN)*qe*c)                    !Acceleration timescale
        tdiff = 3d0*qe*BFF(segN)*rj(segN)**2/(2d0*c*EFF)    !Diffusion timescale
        write(20,*) EFF/eV, tSyn, tIC, tAd, tTot, tcell, tacc, tdiff, tconv
                    
        EFF = EFF*Eint  !Energy increase   
    enddo
    close(20)

    end subroutine Print_tloss

END MODULE losses
