PROGRAM main

  use global
  use controles
  use dist, only: diste_run
  use radsyn, only: radsyn_run

  implicit none

  !================================================================


  write(*,*) 'ELECTRONES:'
        
  if (dist) call diste_run 

  if (radsyn) call radsyn_run       ! Use default value of B in system_parameters
  !if (radsyn) call radsyn_run(1.d-1)! Give explicit value of B

  !================================================================

  write(*,*) 'B [uG] = ', B*1.d6
  write(*,*) 'E_max [TeV] = ', Emax_e/eV/1.d12

  write(*,*) 'Terminado todo'
  
END PROGRAM main

   
    
