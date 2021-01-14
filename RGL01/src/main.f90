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
	write(*,*)'Symmetry L=0: ',Glob_YOperatorString0
	write(*,*)'Symmetry L=1: ',Glob_YOperatorString1
 endif


DeltaE = Glob_CurrEnergy0-Glob_CurrEnergy1
DeltaE_cm = DeltaE * 219474.6				! hartree to cm-1
call ComputeTranDipoL0L1()
LineStrength=Glob_ExpVals**TWO

!OscillatorStrength=TWO/THREE*Glob_ExpVals**TWO*abs(Glob_CurrEnergy0-Glob_CurrEnergy1)
OscillatorStrength = TWO/THREE*LineStrength*abs(DeltaE)


	
	
if (Glob_ProcID==0) then
        write(*,*) ' '
        write(*,'(a30)',advance='no') ' E{initial} - E{final}  =     '
        call writereal(6,abs(DeltaE))
        write(*,'(a10)') '   Hartree'
        write(*,'(a30)',advance='no') '                              '
        call writereal(6,abs(DeltaE_cm))
        write(*,'(a10)') 'cm-1'
        write(*,'(a35)') ' <L=0|Dipole{x}|L=1>    =       0.0'
        write(*,'(a35)') ' <L=0|Dipole{y}|L=1>    =       0.0'
        write(*,'(a30)',advance='no') ' <L=0|Dipole{z}|L=1>    =     '
        call writerealadv(6,Glob_ExpVals)
        write(*,'(a30)',advance='no') ' <L=0|Dipole{z}|L=1>**2 =     '
        call writerealadv(6,Glob_ExpVals**2)
        write(*,*) '  '
        write(*,*) ' To compute  Oscillator Strength :'
        write(*,*) '                           2'
        write(*,*) '    OscillatorStrength= -------- * S(multiple)* Delta_E'
        write(*,*) '                         3*g{i}'
        write(*,*) '    g{i} = 2J{i}+1'
        write(*,*) '    S(multiple) = Sum(Ji,Jk) S(line)'
        write(*,*) '    S(line)= <L=0|Dipole{z}|L=1>**2'
        write(*,*) '  '
        write(*,*) '  '
        write(*,*) 'Data Reader and Matrix Calculator Program is completed'
        write(*,*) 'Program has stopped'
        write(*,*) ' '
endif



call MPI_FINALIZE(Glob_MPIErrCode)


end program main


