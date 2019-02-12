program main

use workproc
implicit none

!Local variables
integer      i

!Initialize MPI
call MPI_INIT(Glob_MPIErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  write (*,*) 'Program Expilitly Correlated Real Gaussians has started'
  write (*,*) 'Number of parallel processes running ',Glob_NumOfProcs
  write (*,*)
endif

call ReadIOFile()
call ProgramDataInit()

allocate(Glob_ExpVals(1:Glob_NumOfDRMCSteps))
do i=1,Glob_NumOfDRMCSteps

  Glob_CurrDRMCStep=i
  select case (Glob_DRMC(i)%Action)
  
     case('OP_DIPOLE')
        call ComputeExpValL0L1()
	if (Glob_ProcID==0) then
           write(*,*) 'Glob_n:         ',Glob_n
           write(*,*) 'NumOfDRMCSteps: ',Glob_NumOfDRMCSteps
           write(*,*) ' '
           write(*,*) 'Symmetry L=0: ',Glob_YOperatorString0
           write(*,*) 'Symmetry L=1: ',Glob_YOperatorString1
           write(*,*) ' '
           write(*,*) 'NonlinParam(:,2) L=0'
           write(*,*) Glob_NonlinParam0(:,2)
           write(*,*) 'NonlinParam(:,2) L=1'
           write(*,*) Glob_NonlinParam0(:,2)
           write(*,*) ' '
           write(*,*) 'FileName1: ',Glob_DRMC(i)%FileName1
	   write(*,*) 'FileName2: ',Glob_DRMC(i)%FileName2
	   write(*,*) 'FileName3: ',Glob_DRMC(i)%FileName3
	   write(*,*) 'FileName4: ',Glob_DRMC(i)%FileName4
           write(*,*) ' '
           write(*,*) 'S0(2,1:5):'
           write(*,*) Glob_S0(2,1:5)
           write(*,*) ' '
           write(*,*) 'S1(2,1:5):'
           write(*,*) Glob_S1(2,1:5)
           write(*,*) ' '
           write(*,*) 'c0'
           write(*,*) Glob_c0
           write(*,*) ' '
           write(*,*) 'c1'
           write(*,*) Glob_c1
           write(*,*) ' '
           write(*,*) 'Glob_YCoeff0: ',Glob_YCoeff0
           write(*,*) ' '
           write(*,*) 'Glob_YMatr0'
	   write(*,*) 'Glob_YMatr0(1): ',Glob_YMatr0(1:Glob_n,1:Glob_n,1)
	   write(*,*) 'Glob_YMatr0(2): ',Glob_YMatr0(1:Glob_n,1:Glob_n,2)
           write(*,*) ' '
           write(*,*) 'Glob_YCoeff1: ',Glob_YCoeff1
           write(*,*) ' '
           write(*,*) 'Glob_YMatr1'
	   write(*,*) 'Glob_YMatr1(1): ',Glob_YMatr1(1:Glob_n,1:Glob_n,1)
	   write(*,*) 'Glob_YMatr1(2): ',Glob_YMatr1(1:Glob_n,1:Glob_n,2)
	   write(*,*) ' '
	   write(*,*) 'Expectation Value (DIPOLE): ',Glob_ExpVals(Glob_CurrDRMCStep)
	endif
	
  endselect
enddo

if (Glob_ProcID==0) then
write(*,*) ' '
write(*,*) 'Data Reader and Matrix Calculator Program is completed'
write(*,*) 'Program has stopped'
end if

call MPI_FINALIZE(Glob_MPIErrCode)

end program main
