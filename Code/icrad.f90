MODULE icrad
	use globals
	implicit none

contains

  !============================================================================
  ! computation of the IC SED
  !============================================================================

  subroutine SEDIC ()
    real*8 :: Eg_FF(NSED_IC),Eg_int_FF,Eg(NSED_IC),Eg_int
    real*8 :: E_e(nIn),NE_e(nIn,Ns),NEcut_e(nIn)
    real*8 :: LIC_FF(NSED_IC),LIC_cut(NSED_IC),LIC1(NSED_IC),LIC2(NSED_IC),LIC(NSED_IC)
    real*8 :: LIC_bol_FF,LIC_bol_cut_FF
    real*8 :: dummy
    integer :: j_s
    allocate (lenj(Ns),nfFF(Ns),EphMFF(Ns),cwoFF(Ns),Kappa(Ns),delta_star(Ns),delta_obs(Ns))

    !Quantities at each segment
    open(10,file='Quant_FF.dat')
    do i = 1,Ns
       read(10,*) lenj(i),dummy,nfFF(i),dummy,dummy,EphMFF(i),cwoFF(i),Kappa(i)
    enddo
    close(10)

    !Quantities at each segment
    open(20,file='Deltas.dat')
    do i = 1,Ns
       read(20,*) delta_star(i),delta_obs(i)
    enddo
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
    do j_s = 1,Ns	!Loop for the segments

       !Photon energies in the FF such that they all match the same energies in the OF
       call vector_log( Emin_ic/delta_obs(j_s),Emax_ic/delta_obs(j_s),Eg_int_FF,Eg_FF )
       !!!!!!!!!!!!!!!!!!!!!!Por qué no Eg_FF = Eg/delta y ya está?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       do i = 1,nIn
          NEcut_e(i) = NE_e(i,j_s)
       enddo
              
       !Unabsorbed IC luminosity for each segment in LIC_FF
       call IC_Khang(E_e,NEcut_e,Eg_FF,cwoFF(j_s),EphMFF(j_s),Kappa(j_s),LIC_FF,LIC_bol_cut_FF)
       !call IC_OLD(E_e,NEcut_e,Eg_FF,cwoFF(j_s),EphMFF(j_s),nfFF(j_s),LIC_FF,LIC_bol_cut_FF)
       LIC_bol_FF = LIC_bol_FF + LIC_bol_cut_FF

       !Unabsorbed boosted luminosity for Eg (DIFFERENT for jet and c-jet)
       LIC_cut = delta_obs(j_s)**2/gj*LIC_FF ! These are vectors evaluated in E_g(i)
       LIC = LIC + LIC_cut

       do i = 1,NSED_IC
          write(40,*) Eg(i)/eV,Eg(i)*LIC_cut(i),j_s
       enddo
       write(40,*) 
       
       !Jet contribution of the segments up to given lengths
       if (lenj(j_s)-lenj(1) <= 0.1d0*a) LIC1 = LIC
       if (lenj(j_s)-lenj(1) <= a) LIC2 = LIC
       
    enddo
    !write(*,'(3X,A,ES9.2E2,1X,A)') 'Jet unabsorbed IC luminosity (FF):',LIC_bol_FF,'erg/s'
    
    do i = 1,NSED_ic
       write(50,*) Eg(i)/eV,Eg(i)*LIC1(i),Eg(i)*LIC2(i),Eg(i)*LIC(i)
    enddo
    
    close(40)
    close(50)
    
    deallocate (lenj,EphMFF,nfFF,cwoFF,Kappa,delta_star,delta_obs)

  END subroutine SEDIC
  
  !============================================================================
  ! computation of the IC SED
  !============================================================================

  subroutine SEDIC_CJ ()
    real*8 :: Eg_FF(NSED_IC),Eg_int_FF,Eg(NSED_IC),Eg_int
    real*8 :: E_e(nIn),NE_e(nIn,Ns),NEcut_e(nIn)
    real*8 :: LIC_FF(NSED_IC),LIC_cut(NSED_IC),LIC1(NSED_IC),LIC2(NSED_IC),LIC(NSED_IC)
    real*8 :: LIC_bol_FF,LIC_bol_cut_FF
    real*8 :: dummy
    integer :: j_s
    allocate (lenj(Ns),nfFF(Ns),EphMFF(Ns),cwoFF(Ns),Kappa(Ns),delta_star(Ns),delta_obs(Ns))

    !Quantities at each segment
    open(10,file='Quant_FFCJ.dat')
    do i = 1,Ns
       read(10,*) lenj(i),dummy,nfFF(i),dummy,dummy,EphMFF(i),cwoFF(i),Kappa(i)
    enddo
    close(10)

    !Quantities at each segment
    open(20,file='DeltasCJ.dat')
    do i = 1,Ns
       read(20,*) delta_star(i),delta_obs(i)
    enddo
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

    open(40,file='SEDIC_cutsCJ.dat')
    open(50,file='SEDIC_CJ.dat')

    call vector_log( Emin_ic,Emax_ic,Eg_int,Eg )

    LIC = 0.d0
    do j_s = 1,Ns !Loop for the segments

       !Photon energies in the FF such that they all match the same energies in the OF
       call vector_log( Emin_ic/delta_obs(j_s),Emax_ic/delta_obs(j_s),Eg_int_FF,Eg_FF )
       !!!!!!!!!!!!!!!!!!!!!!Por que no Eg_FF = Eg/delta y ya esta?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       do i = 1,nIn
          NEcut_e(i) = NE_e(i,j_s)
       enddo
              
       !Unabsorbed IC luminosity for each segment in LIC_FF
       call IC_Khang(E_e,NEcut_e,Eg_FF,cwoFF(j_s),EphMFF(j_s),Kappa(j_s),LIC_FF,LIC_bol_cut_FF)
       !call IC_OLD(E_e,NEcut_e,Eg_FF,cwoFF(j_s),EphMFF(j_s),nfFF(j_s),LIC_FF,LIC_bol_cut_FF)
       LIC_bol_FF = LIC_bol_FF + LIC_bol_cut_FF

       !Unabsorbed boosted luminosity for Eg (DIFFERENT for jet and c-jet)
       LIC_cut = delta_obs(j_s)**2/gj*LIC_FF ! These are vectors evaluated in E_g(i)
       LIC = LIC + LIC_cut

       do i = 1,NSED_IC
          write(40,*) Eg(i)/eV,Eg(i)*LIC_cut(i),j_s
       enddo
       write(40,*) 
       
       !Jet contribution of the segments up to given lengths
       if (lenj(j_s)-lenj(1) <= 0.1d0*a) LIC1 = LIC
       if (lenj(j_s)-lenj(1) <= a) LIC2 = LIC
       
    enddo
    !write(*,'(3X,A,ES9.2E2,1X,A)') 'C-jet unabsorbed IC luminosity (FF):',LIC_bol_FF,'erg/s'
    
    do i = 1,NSED_ic
       write(50,*) Eg(i)/eV,Eg(i)*LIC1(i),Eg(i)*LIC2(i),Eg(i)*LIC(i)
    enddo
    
    close(40)
    close(50)
    
    deallocate (lenj,EphMFF,nfFF,cwoFF,Kappa,delta_star,delta_obs)

  END subroutine SEDIC_CJ

  !==============================

!/////////////////////////////////////////////////////////////////////////////
!                                                                           //
!                              SUBRUTINAS                                   //
!                                                                           //
!/////////////////////////////////////////////////////////////////////////////

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
        if (z_IC < 1d0) then 		!
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

!--------------------------------------------------------------------------------

SUBROUTINE IC_OLD(Ee,Ne,Ef,cang,Eph0,nf,L_IC,sum_IC)
  integer, parameter :: dp = kind(1.d0)
  Real(dp),intent(in) :: Ee(nIn),Ne(nIn),Eph0,Ef(NSED_IC),cang,nf
  Real(dp),intent(out) :: L_IC(NSED_IC),sum_IC
  Real(dp) :: Eeint,Eps,Epsint
  Real(dp) :: dEps,dEe,dP_IC,sum,gamma,PepsIC,cs,sn
  Real(dp) :: u2,beta,s,ee0,pk,gc,wc,uc,vc,ec,wmax,wmin,csg,csa,ai1,ai2
  Real(dp) :: dsigmaIC,divis,cte,factor_ang
  cs = cang
  sn = sin(acos(cang))

  Eeint = Ee(2)/Ee(1)
  cte = pi*2.8d-13**2.
  factor_ang = c*(1.d0-cs)

  !     Calculo IC para distintas energías de Ef

  Epsint = Ef(2)/Ef(1)
  !     Eph será la energía del fotón incidente, Eps la del fotón emitido y 
  !     gamma la energía del electrón, todas se usan en unidades de mec2
  sum_IC = 0.d0
  do i=1,NSED_IC                 ! Ef
     Eps = Ef(i)/mec2
     dEps = Eps*(Epsint-1.d0)
     sum = 0.d0
     do j = 1,nIn            ! E
        gamma = Ee(j)/mec2  ! gm
        PepsIC = 0.d0
        u2 = sqrt(1.d0-1./gamma**2.) ! u
        beta=(1.d0/gamma**2.+u2**2.*sn**2.)/(1.+u2*cs) ! los usa en sigmaIC
        s=0.d0           ! acá arranca la subrutina para calcular sigmaIC
        ee0=Eph0/mec2
        pk=gamma*ee0*beta
        divis = sqrt(1.d0+2.d0*pk)
        gc=(gamma+ee0)/divis
        wc=pk/divis
        uc=pk/(1.d0+pk)
        vc=sqrt(ee0**2.+(gamma*u2)**2.+2.d0*gamma*ee0*u2*cs)/(gamma+ee0)
        ec=(1.d0+pk)/divis
        wmax=gc*wc*(1.+vc)
        wmin=wc/gc/(1.+vc)
        if((Eps > wmax).or.(Eps < wmin)) goto 7777
        csg=-(1.-gamma/gc/ec)/uc/vc
        csa= (Eps/gc/wc-1.)/vc
        ai1=1./sqrt((1+uc-2.d0*uc*csa*csg)/(1.+uc)/ec**2.+&
             uc**2.*(csa-csg)**2.)
        ai2=(1.-uc*csa*csg)*ai1**3.
        s= ( (1.d0/pk**2.-2.*ai1/pk/ec/wc+ai2/(ec*wc)**2.) + & 
             (2.d0*(1./pk-ai1/ec/wc)+(pk*ai1/ec/wc+ec*wc*(1.-uc*csa*csg)/pk)) )* &
             (cte/pk**2.*wc/gc/vc)
7777    continue
        dsigmaIC=s/mec2

        PepsIC = Ef(i)*nf*factor_ang*dsigmaIC
        dEe = Ee(j)*(Eeint-1.d0)
        dP_IC = Ne(j)*PepsIC*dEe
        sum = dP_IC + sum
     enddo
     L_IC(i) = sum     ! Es la P , luego grafico E(j)*LIC(Ej)
     sum_IC = sum_IC + sum*dEps*mec2
  enddo

  RETURN

end SUBROUTINE IC_OLD

END MODULE icrad
