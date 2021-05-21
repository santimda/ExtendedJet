PROGRAM main
  implicit none
  real(kind(1.d0)) :: B

  interface
    subroutine jet_calc(B_value) 
      real(kind(1.d0)), intent(in) :: B_value
    end subroutine jet_calc
  end interface

  B = 1.d-2
  call jet_calc(B)

  B = 1.d0
  call jet_calc(B)

END PROGRAM main

  !==================================================================

  subroutine jet_calc(B_value)
    use global
    use controles
    use dist, only: diste_run
    use radsyn, only: radsyn_run

    implicit none
    real(dp), intent(in) :: B_value

    !================================================================


    write(*,*) '\n ELECTRONES:'
          
    if (dist) call diste_run 

    if (radsyn) call radsyn_run(B_value)! Give explicit value of B

    !================================================================

    write(*,*) 'B [uG] = ', B_value*1.d6
    write(*,*) 'E_max [TeV] = ', Emax_e/eV/1.d12

    write(*,*) 'Terminado todo'
  end subroutine

  !===================================================================


   
    
