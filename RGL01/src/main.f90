program main

use workproc
implicit none

!Local variables
integer        :: i,j
character(70)  :: ReadChar
real(dprec)    :: OscillatorStrength

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

if (Glob_ProcID==0) then
   write(*,*) 'Glob_n:         ',Glob_n
   write(*,*) 'NumOfDRMCSteps: ',Glob_NumOfDRMCSteps
   write(*,*) ' '
   write(*,*) 'Symmetry L=0: ',Glob_YOperatorString0
   write(*,*) 'Symmetry L=1: ',Glob_YOperatorString1
   write(*,*) ' '
   
   open(1,file=Glob_DataFileName,status='old',action='readwrite')
   do j=1,7+Glob_NumOfDRMCSteps
      read(1,*) ReadChar
   enddo
   write(1,*) 'Basis Size | Transition Dipole Matrix | Oscillator Strength'
   close(1)
endif

allocate(Glob_ExpVals(1:Glob_NumOfDRMCSteps))
do i=1,Glob_NumOfDRMCSteps

  Glob_CurrDRMCStep=i
  if (Glob_ProcID==0) then
     write(*,*) '================================================'
     write(*,*) ' '
     write(*,*) 'Current DRMC Step No.',Glob_CurrDRMCStep
     write(*,*) ' '
  endif
  call ReadIOFileForDRMCAction()
  select case (Glob_DRMC(i)%Action)
  
     case('OP_DIPOLE')
        call ComputeExpValL0L1()
	OscillatorStrength=TWO/THREE*Glob_ExpVals(Glob_CurrDRMCStep)**TWO*&
	                    abs(Glob_E0(Glob_CurrDRMCStep)-Glob_E1(Glob_CurrDRMCStep))
	if (Glob_ProcID==0) then
	   write(*,*) 'Dipole:'
	   write(*,*) '<|DipVec_x|>=<|DipVec_y|>=    0'
	   write(*,*) '<|DipVec_z|>=              ',Glob_ExpVals(Glob_CurrDRMCStep)
	   write(*,*) 'Transition Dipole Moment:  ',abs(Glob_ExpVals(Glob_CurrDRMCStep))
	   write(*,*) 'Oscillator Strength:       ',OscillatorStrength
	   open(1,file=Glob_DataFileName,status='old',action='readwrite')
	   do j=1,7+Glob_NumOfDRMCSteps+Glob_CurrDRMCStep
	      read(1,*) ReadChar
	   enddo
	   write(1,*) Glob_CurrBasisSize,Glob_ExpVals(Glob_CurrDRMCStep),OscillatorStrength
	   close(1)
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
