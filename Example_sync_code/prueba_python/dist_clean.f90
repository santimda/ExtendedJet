subroutine diste_run  
  !============================================================================
  ! calculate electron distribution radiation
  !============================================================================
  use global
  implicit none
  real(dp) :: integral, Ka
  allocate( N_e(mE_e),E_e(mE_e) )

  call vector_log(Emin_e, 10.d0*Emax_e, Eint_e, E_e)
   
  integral = 0.d0
  do i = 1, mE_e
    N_e(i) = E_e(i)**(-p) * (exp(-E_e(i)/Emax_e))
    dE_e = E_e(i)*(Eint_e-1.d0)
    integral = integral + E_e(i)*N_e(i)*dE_e
  enddo
  Ka = U*Vol/integral   
  N_e = Ka*N_e       

  open(79,file='dist_e.dat')
  do i = 1, mE_e
    write(79,*) E_e(i)/eV, N_e(i)
  enddo
  close(79)

  deallocate( N_e,E_e ) 

end subroutine diste_run