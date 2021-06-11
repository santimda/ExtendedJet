MODULE icrad
  use constants
  use globals
  implicit none

contains

  !============================================================================
  ! Computation of the IC SED
  !============================================================================
  ! It calculates anisotropic IC if aniso=.True., else isotropic. 
  !============================================================================

  subroutine SEDIC(aniso)
    real(dp) :: Eg_FF(NSED_IC), Eg_int_FF, Eg(NSED_IC), Eg_int
    real(dp) :: E_e(nIn), NE_e(nIn,Ns), NEcut_e(nIn)
    real(dp) :: LIC_FF(NSED_IC), LIC_cut(NSED_IC), LIC1(NSED_IC), LIC2(NSED_IC), LIC(NSED_IC)
    real(dp) :: LIC_bol_FF, LIC_bol_cut_FF
    real(dp) :: dummy, boost(Ns), cang
    integer :: j_s
    logical, intent(in) :: aniso 
    allocate( delta_obs(Ns) )

    write(*,*) '\n Calculating the IC SED along the jet'
        
    !Quantities at each segment
    open(20,file='deltas.dat')
    do i = 1, Ns
      read(20,*) delta_obs(i), boost(i)
    enddo
    close(10)
    close(20)


    !N(E) of each cut
    open(30,file='NE_cuts.dat')
    do i = 1, nIn
       do j = 1, Ns
          read(30,*) E_e(i), NE_e(i,j), dummy
       enddo
    enddo
    close(30)

    cang = cos(jo) ! cosine of interaction angle, almost dummy variable for now

    !-------------------------------------------------------------------!
    !                     Inverse Compton                               !
    !-------------------------------------------------------------------!

    open(40,file='SEDIC_cuts.dat')
    open(50,file='SEDIC_jet.dat')

    write(*,'(1X,A,I3,A)') 'IC SED (',NSED_IC,')'

    call vector_log( Emin_ic,Emax_ic,Eg_int,Eg )

    LIC = 0.d0
    do j_s = 1, Ns !Loop for the segments

       !Photon energies in the FF such that they all match the same energies in the OF
       call vector_log( Emin_ic/delta_obs(j_s), Emax_ic/delta_obs(j_s), Eg_int_FF, Eg_FF )
       !!!!!!!!!!!!!!!!!!!!!!Por qu no Eg_FF = Eg/delta y ya esta?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       do i = 1,nIn
          NEcut_e(i) = NE_e(i,j_s)
       enddo
              
       !Unabsorbed IC luminosity for each segment in LIC_FF
       call IC_Khang( E_e, NEcut_e, Eg_FF, aniso, cang, EphMFF(j_s), Kappa(j_s), LIC_FF, LIC_bol_cut_FF )
       LIC_bol_FF = LIC_bol_FF + LIC_bol_cut_FF

       !Unabsorbed boosted luminosity for Eg (DIFFERENT for jet and c-jet)
       LIC_cut = boost(j_s)*LIC_FF ! These are vectors evaluated in E_g(i)
       LIC = LIC + LIC_cut

       do i = 1,NSED_IC
          write(40,*) Eg(i)/eV,Eg(i)*LIC_cut(i),j_s
       enddo
       write(40,*) 
       
       !Jet contribution of the segments up to given lengths
       if (z(j_s) <= 0.01d0*zmax) LIC1 = LIC
       if (z(j_s) <= 0.1d0*zmax) LIC2 = LIC
       
    enddo
    !write(*,'(3X,A,ES9.2E2,1X,A)') 'Jet unabsorbed IC luminosity (FF):',LIC_bol_FF,'erg/s'
    
    do i = 1,NSED_ic
       write(50,*) Eg(i)/eV,Eg(i)*LIC1(i),Eg(i)*LIC2(i),Eg(i)*LIC(i)
    enddo
    
    close(40)
    close(50)
    
    deallocate( delta_obs )

  END subroutine SEDIC
  

!/////////////////////////////////////////////////////////////////////////////
!                                                                           //
!                              SUBRUTINAS                                   //
!                                                                           //
!/////////////////////////////////////////////////////////////////////////////


!============================================================================
! This subroutine calculates the intrinsic IC emission using Khangulyan+2014
! formulae for black-body-like photon fields. 
! It calculates anisotropic IC if aniso=.True., else isotropic. 
!============================================================================

subroutine IC_Khang(Ee, Ne, Ef, aniso, cang, Eph0, Kapp, L_IC, sum_IC)
  real(dp), intent(in) :: Ee(nIn), Ne(nIn), Ef(NSED_IC), cang, Eph0, Kapp
  real(dp), intent(out) :: L_IC(NSED_IC), sum_IC
  logical, intent(in) :: aniso
  real(dp) :: P_eps, dP_eps, dEe, Eint_e
  real(dp) :: gamma_e(nIn), Eps(NSED_IC), Eps_int, dEps
  real(dp) :: T_ad, T_ad_ang, cte_Ndot, Ndot, x0, z_IC, t_K, sum
  integer, parameter :: nIn = 80       ! Number of points in the electron energy distribution

  T_ad = Eph0/(2.7*mec2)     ! adimensional temperature, Eph0 = 2.7*kT
  
  T_ad_ang = T_ad  !  dummy definition so compiler doesn't warn possible uninitialization
  if (aniso) T_ad_ang = T_ad*(1.-cang) ! anisotropic IC
  
  Eps = Ef/mec2              ! adimensional gamma photon energy
  Eps_int = Eps(2)/Eps(1)
  gamma_e = Ee/mec2
  Eint_e = Ee(2)/Ee(1)
  sum_ic = 0.d0

  cte_Ndot = cte_IC_K*Kapp*T_ad**2
  
  do i = 1,NSED_IC                      ! gamma-ray photon energies
     dEps = Eps(i)*(Eps_int-1.d0)
     sum = 0.d0
     do j = 1, nIn                       ! electron energies
        z_IC = Eps(i)/gamma_e(j)
        if (z_IC < 1d0) then    
          if (aniso) then                !anisotropic IC
            t_K = 2d0*gamma_e(j)*T_ad_ang
            x0 = z_IC/( (1.-z_IC)*t_K )
            Ndot = (cte_Ndot/gamma_e(j)**2) * ( (z_IC**2/(2.*(1.-z_IC)))*F_1(x0) + F_2(x0) ) 
          else
            t_K = 4.*gamma_e(j)*T_ad      !isotropic IC
            x0 = z_IC/( (1.-z_IC)*t_K )
            Ndot = (cte_Ndot/gamma_e(j)**2) * ( (z_IC**2/(2.*(1.-z_IC)))*F_3(x0) + F_4(x0) ) 
          endif
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
