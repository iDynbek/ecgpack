
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
    write(*,*)
    write(*,'(A)') '  +----------------------------------------------------------------+'
    write(*,'(A)') '  |                                                                |'
    write(*,'(A)') '  |     Explicitly correlated real gaussians code has started      |'
    write(*,'(A)') '  |                                                                |'
    write(*,'(A)') '  |                calculating the x component of                  |'
    write(*,'(A)') '  |               transition dipole moment integral                |'
    write(*,'(A)') '  |                                                                |'
    write(*,'(A)') '  |           P^o, y component  -->  P^e, z component              |'
    write(*,'(A)') '  +----------------------------------------------------------------+'
    write(*,*)
    write(*,'(4X,A,I4)') 'Number of parallel processes:  ', Glob_NumOfProcs
    write(*,*)
endif


call READwf1wf2()
call ProgramDataInit()


DeltaE = Glob_CurrEnergy1-Glob_CurrEnergy2
call ComputeTranDipoL1L2()
	

if (Glob_ProcID==0) then
    write(*,*)
    write(*,*)
    write(*,'(A)') '   Energy difference:     '
    write(*,'(A)') '  ------------------------'
    write(*,*)
    write(*,'(4X,A)',advance='no') '    E(P) - E(D)  = '
    call writereal(6,abs(DeltaE))
    write(*,'(A)') '   Hartree'

    write(*,*)
    write(*,*)
    write(*,'(A)') '   Length Gauge:          '
    write(*,'(A)') '  ------------------------'
    write(*,*)    
    write(*,'(4X,A)',advance='no') ' < Po | x | Pe >   = '
    call writerealadv(6,Glob_ExpVals1)                   
    write(*,'(4X,A)') ' < Po | y | Pe >   =   0.0'
    write(*,'(4X,A)') ' < Po | z | Pe >   =   0.0'


    write(*,*)
    write(*,*)
    write(*,'(A)') '   Velocity Gauge:        '
    write(*,'(A)') '  ------------------------'
    write(*,*)
    write(*,'(4X,A)',advance='no') '< Po | P(z) | Pe >  ='
    call writerealadv(6,Glob_ExpVals2)
    write(*,'(4X,A)') '< Po | P(y) | Pe >  =  0.0'
    write(*,'(4X,A)') '< Po | P(x) | Pe >  =  0.0'


    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,'(A)') ' Program completed successfully.      '
    write(*,*)

endif


call MPI_FINALIZE(Glob_MPIErrCode)



end program main
