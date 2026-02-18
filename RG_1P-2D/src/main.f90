program main

use workproc
implicit none

!Local variables
integer        :: i,j
character(70)  :: ReadChar
real(dprec)    :: LineStrength
real(dprec)    :: OscillatorStrength
real(dprec)    :: DeltaE


!Initialize MPI
call MPI_INIT(Glob_MPIErrCode)
call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)


if (Glob_ProcID==0) then
	write (*,*)
	write (*,*) ' Program Expilitly Correlated Real Gaussians has started'
	write (*,*) ' Number of parallel processes running ',Glob_NumOfProcs
	write (*,*)
    WRITE(*,*) ' Transition dipole moment integral calculations are starting for:'
    WRITE(*,*) '   L = 1, M_L = 0  ==>  L = 2, M_L = 0'
    WRITE(*,*) ' transition.'
endif


call Readwf0wf1()
call ProgramDataInit()


if (Glob_ProcID==0) then
	write(*,*)
	write(*,*)'Symmetry of the initial state, L=1: ', Glob_YOperatorString1
	write(*,*)'Symmetry of the  final  state, L=2: ', Glob_YOperatorString2
    write(*,*) 
 endif



DeltaE = Glob_CurrEnergy1-Glob_CurrEnergy2


call ComputeTranDipoL1L2()

	
if (Glob_ProcID==0) then

        write (*,*)
        write(*,'(a30)',advance='no') ' E{initial} - E{final}  =     '
        call writereal(6,abs(DeltaE))
        write(*,'(a10)') '   Hartree'
		

        write(*,*) 
		write(*,*) 
		write(*,*) 'Transition Dipole Moment Integral in Length Gauge:  '
		write(*,*) 
        write(*,'(a35)') ' <L=0|x|L=1>    =       0.0'
        write(*,'(a35)') ' <L=0|y|L=1>    =       0.0'
        write(*,'(a30)',advance='no') ' <L=0|z|L=1>    =     '
        call writerealadv(6,Glob_ExpVals1)

        write(*,*) 
		write(*,*) 
		write(*,*) 'Transition Dipole Moment Integral in Velocity Gauge:  '
		write(*,*) '  '
        write(*,'(a30)') ' <L=0|P(x)|L=1>    =     0.0'
        write(*,'(a30)') ' <L=0|P(y)|L=1>    =     0.0'		
        write(*,'(a30)',advance='no') ' <L=0|P(z)|L=1>    = -i *'
        call writerealadv(6,Glob_ExpVals2)

        write(*,*) 
        write(*,*) 
        write(*,*) 
        write(*,*) 'Data Reader and Matrix Calculator Program is completed'
        write(*,*) 'Program has stopped'
        write(*,*) 
endif



call MPI_FINALIZE(Glob_MPIErrCode)


end program main


