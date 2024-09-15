program main

use workproc
implicit none

!Local variables
integer        :: i,j
character(70)  :: ReadChar
real(dprec)    :: LineStrength
real(dprec)    :: OscillatorStrength
real(dprec)    :: DeltaE
real(dprec)    :: DeltaE_cm


!Initialize MPI
call MPI_INIT(Glob_MPIErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)

if (Glob_ProcID==0) then
	write (*,*) 'Program Expilitly Correlated Real Gaussians has started'
	write (*,*) 'Number of parallel processes running ',Glob_NumOfProcs
	write (*,*)
endif



call Readwf0wf1()
call ProgramDataInit()
if (Glob_ProcID==0) then
	write(*,*) ' '
	write(*,*)'Young string L=1: ',Glob_YOperatorString0
	write(*,*)'Young string L=2: ',Glob_YOperatorString1
endif


call ComputeSpinDep()

!call ComputeScalar()

call MPI_FINALIZE(Glob_MPIErrCode)


end program main


