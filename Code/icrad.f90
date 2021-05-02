MODULE icrad
  use globals
  implicit none

contains

  !============================================================================
  ! computation of the IC SED
  !============================================================================

  subroutine SEDIC ()
    real(dp) :: Eg_FF(NSED_IC),Eg_int_FF,Eg(NSED_IC),Eg_int
    real(dp) :: E_e(nIn),NE_e(nIn,Ns),NEcut_e(nIn)
    real(dp) :: LIC_FF(NSED_IC),LIC_cut(NSED_IC),LIC1(NSED_IC),LIC2(NSED_IC),LIC(NSED_IC)
    real(dp) :: LIC_bol_FF,LIC_bol_cut_FF
    real(dp) :: dummy, boost(Ns)
    integer :: j_s
    allocate( z(Ns), nfFF(Ns), EphMFF(Ns), Kappa(Ns), delta_obs(Ns) )

    write(*,*) '\n Calculating the IC SED along the jet'
        
    !Quantities at each segment
    open(10,file='quant_FF.dat')
    open(20,file='deltas.dat')
    do i = 1, Ns
      read(10,*) z(i), dummy, dummy, dummy, Kappa(i), dummy, dummy, EphMFF(i), nfFF(i)
      read(20,*) delta_obs(i), boost(i)
    enddo
    close(10)
    close(20)


    !N(E) of each cut
    open(30,file='NE_intcuts.dat')
    do i = 1,nIn
       do j = 1,Ns
          read(30,*) E_e(i),NE_e(i,j)
       enddo
    enddo
    close(30)

    !-------------------------------------------------------------------!
    !                     Inverse Compton                               !
    !-------------------------------------------------------------------!

    open(40,file='SEDIC_cuts.dat')
    open(50,file='SEDIC_jet.dat')

    write(*,'(1X,A,I3,A)') 'IC SED (',NSED_IC,')'

    call vector_log( Emin_ic,Emax_ic,Eg_int,Eg )

    LIC = 0.d0
    do j_s = 1,Ns !Loop for the segments

       !Photon energies in the FF such that they all match the same energies in the OF
       call vector_log( Emin_ic/delta_obs(j_s), Emax_ic/delta_obs(j_s), Eg_int_FF, Eg_FF )
       !!!!!!!!!!!!!!!!!!!!!!Por qué no Eg_FF = Eg/delta y ya está?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       do i = 1,nIn
          NEcut_e(i) = NE_e(i,j_s)
       enddo
              
       !Unabsorbed IC luminosity for each segment in LIC_FF
       call IC_Khang( E_e, NEcut_e, Eg_FF, 0.5d0, EphMFF(j_s), Kappa(j_s), LIC_FF, LIC_bol_cut_FF )
       ! IMPORTANT!!! CHANGE "0.5d0" FOR SOMETHING COMPATIBLE WITH ISOTROPIC PHOTON FIELD
       LIC_bol_FF = LIC_bol_FF + LIC_bol_cut_FF

       !Unabsorbed boosted luminosity for Eg (DIFFERENT for jet and c-jet)
       LIC_cut = delta_obs(j_s)**2/gj*LIC_FF ! These are vectors evaluated in E_g(i)
       LIC = LIC + LIC_cut

       do i = 1,NSED_IC
          write(40,*) Eg(i)/eV,Eg(i)*LIC_cut(i),j_s
       enddo
       write(40,*) 
       
       !Jet contribution of the segments up to given lengths
       if (z(j) <= 0.01d0*zmax) LIC1 = LIC
       if (z(j) <= 0.1d0*zmax) LIC2 = LIC
       
    enddo
    !write(*,'(3X,A,ES9.2E2,1X,A)') 'Jet unabsorbed IC luminosity (FF):',LIC_bol_FF,'erg/s'
    
    do i = 1,NSED_ic
       write(50,*) Eg(i)/eV,Eg(i)*LIC1(i),Eg(i)*LIC2(i),Eg(i)*LIC(i)
    enddo
    
    close(40)
    close(50)
    
    deallocate( z, EphMFF, nfFF, Kappa, delta_obs )

  END subroutine SEDIC
  

!/////////////////////////////////////////////////////////////////////////////
!                                                                           //
!                              SUBRUTINAS                                   //
!                                                                           //
!/////////////////////////////////////////////////////////////////////////////

! IMPORTANT: CHANGE THIS TO ISOTROPIC IC

subroutine IC_Khang(Ee,Ne,Ef,cang,Eph0,Kapp,L_IC,sum_IC)
  Real*8, intent(in) :: Ee(nIn),Ne(nIn),Ef(NSED_IC),cang,Eph0,Kapp
  Real*8, intent(out) :: L_IC(NSED_IC),sum_IC
  Real*8 :: P_eps,dP_eps,dEe,Eint_e
  Real*8 :: gamma_e(nIn),Eps(NSED_IC),Eps_int,dEps
  Real*8 :: T_ad,T_ad_ang,cte_Ndot,Ndot,x0,z_IC,t_K,sum

  T_ad = Eph0/(2.7*mec2)     ! adimensional temperature, Eph0 = 2.7*kT
  T_ad_ang = T_ad*(1.-cang)  ! anisotropic IC
  
  Eps = Ef/mec2              ! adimensional gamma photon energy
  Eps_int = Eps(2)/Eps(1)
  gamma_e = Ee/mec2
  Eint_e = Ee(2)/Ee(1)
  sum_ic = 0.d0

  cte_Ndot = cte_IC_K*Kapp*T_ad**2
  
  do i = 1,NSED_IC                      ! gamma-ray photon energies
     dEps = Eps(i)*(Eps_int-1.d0)
     sum = 0.d0
     do j = 1,nIn                       ! electron energies
        z_IC = Eps(i)/gamma_e(j)
        if (z_IC < 1d0) then    !
           t_K = 2d0*gamma_e(j)*T_ad_ang  !anisotropic
           x0 = z_IC/( (1.-z_IC)*t_K )
           Ndot = (cte_Ndot/gamma_e(j)**2) * ( (z_IC**2/(2.*(1.-z_IC)))*F_1(x0) + F_2(x0) ) !anisotropic
           P_eps = Eps(i)*Ndot           ! N^dot = nf*c*dsigmaIC
           dEe = Ee(j)*(Eint_e-1.d0)
           dP_eps = Ne(j)*P_eps*dEe
        else
           dP_eps = 0.d0
        endif
        sum = dP_eps + sum
     enddo
     L_ic(i) = sum
     !   P=[m_e c**2/s], lo pasamos a [erg/s]. Eps lo pasamos a eV.
     sum_ic = sum_ic + L_ic(i)*dEps*mec2
  enddo

end subroutine IC_Khang

!=====================================================================

END MODULE icrad