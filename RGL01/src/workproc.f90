module workproc
!This module contains basic work subroutines
use matelem
implicit none


contains


subroutine ReadIOFile()
!Subroutine ReadIOFile reads data (nonlinear variational
!parameters and other information) from the input/output 
!file whose name is specified by global variable 
!Glob_DataFileName and other files. If there is no such a file in the 
!current directory then the program stops.

!Local variables:
integer        :: OpenFileErr
real(dprec)    :: ReadRealA
integer        :: ReadInt,ReadErr
integer        :: WorkInt(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
integer        :: WorkInt0(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
integer        :: WorkInt1(max(max(Glob_YOperatorStringLength,20),Glob_FileNameLength))
integer        :: i,j,k,l,Line,j1,j2,j3,j4
character(70)  :: ReadChar
logical        :: ErrorInDataFile,IsDRMCStep

ErrorInDataFile=.false.

if (Glob_ProcID==0) then
  open(1,file=Glob_DataFileName,status='old',iostat=OpenFileErr)
  if (OpenFileErr/=0) then
    write (*,*) 'Error in DataFileInit: data file not found - ',Glob_DataFileName
    ErrorInDataFile=.true.
  endif
endif

call MPI_BCAST(ErrorInDataFile,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (ErrorInDataFile) stop

!Reading information
if (Glob_ProcID==0) Line=0
if (Glob_ProcID==0) then
  write(*,*) 'Reading initial conditions from data file ',Glob_DataFileName
  read(1,*) ReadChar(1:9),ReadInt
  write(*,'(1x,a9,1x,i6)') ReadChar(1:9),ReadInt
  Line=Line+1
  Glob_n=ReadInt-1 !Glob_n is the number of pseudoparticles
  if ((Glob_n<1).or.(ReadChar(1:9)/='PARTICLES')) then
    write(*,*) 'Error in data file, line ',Line   
    ErrorInDataFile=.true.
  endif
endif 
if (Glob_n>Glob_MaxAllowedNumOfPseudoParticles) then
  if (Glob_ProcID==0) then
    write (*,*) 'The version of the code you are running was compiled for the case'
    write (*,*) 'when the number of particles in the system is smaller or equal to', &
                 Glob_MaxAllowedNumOfParticles
    write (*,*) 'while the number of particles specified in the input file is',Glob_n+1
    write (*,*) 'Please make appropriate changes. Program will now stop.'
  endif
  ErrorInDataFile=.true.
endif
if (ErrorInDataFile) stop
call MPI_BCAST(Glob_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
Glob_np=Glob_n*(Glob_n+1)/2
Glob_npt=Glob_np
Glob_2raised3n2=TWO**((3*Glob_n)/TWO)

allocate(Glob_Mass(Glob_n+1))
if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:6),Glob_Mass(1:Glob_n+1)
  write(*,'(1x,a6)',advance='no') ReadChar(1:6)
  call writerealarradv(6,Glob_Mass,Glob_n+1)
endif
call MPI_BCAST(Glob_Mass,Glob_n+1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

allocate(Glob_PseudoCharge(Glob_n))
if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:7),Glob_PseudoCharge0,Glob_PseudoCharge(1:Glob_n)
  write(*,'(1x,a7)',advance='no') ReadChar(1:7)
  call writereal(6,Glob_PseudoCharge0)
  call writerealarradv(6,Glob_PseudoCharge,Glob_n)
endif
call MPI_BCAST(Glob_PseudoCharge0,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_PseudoCharge,Glob_n,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  Glob_RepulsionScalingParam=1.0_dprec
  Glob_RepScalParamSupplied=.false.
  Glob_RepulsionScalingParamPlus=1.0_dprec
  Glob_RepScalParamPlusSupplied=.false.  
  Glob_RepulsionScalingParamMinus=1.0_dprec
  Glob_RepScalParamMinusSupplied=.false.  
  do i=1,3
    read(1,*,iostat=ReadErr) ReadChar(1:29),ReadRealA  
    if ((ReadErr/=0).or.(ReadChar(1:23)/='REPULSION_SCALING_PARAM')) then
      backspace 1
    else
      if (ReadChar(1:28)=='REPULSION_SCALING_PARAM_PLUS') then
        Glob_RepulsionScalingParamPlus=ReadRealA
        Glob_RepScalParamPlusSupplied=.true.
        write(*,'(1x,a28)',advance='no') ReadChar(1:28)
        call writerealadv(6,Glob_RepulsionScalingParamPlus)      
      elseif (ReadChar(1:29)=='REPULSION_SCALING_PARAM_MINUS') then
        Glob_RepulsionScalingParamMinus=ReadRealA
        Glob_RepScalParamMinusSupplied=.true.
        write(*,'(1x,a28)',advance='no') ReadChar(1:29)
        call writerealadv(6,Glob_RepulsionScalingParamMinus)        
      else
        Glob_RepulsionScalingParam=ReadRealA
        Glob_RepScalParamSupplied=.true.
        write(*,'(1x,a23)',advance='no') ReadChar(1:23)
        call writerealadv(6,Glob_RepulsionScalingParam)      
      endif
    endif 
  enddo 
endif
call MPI_BCAST(Glob_RepScalParamSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_RepulsionScalingParam,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_RepScalParamPlusSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_RepulsionScalingParamPlus,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_RepScalParamMinusSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_RepulsionScalingParamMinus,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*,iostat=ReadErr) ReadChar(1:24),Glob_AttractionScalingParam
  if ((ReadErr/=0).or.(ReadChar(1:24)/='ATTRACTION_SCALING_PARAM')) then
    Glob_AttractionScalingParam=1.0_dprec
    Glob_AttrScalParamSupplied=.false.
    backspace 1
  else
    Glob_AttrScalParamSupplied=.true.
    write(*,'(1x,a24)',advance='no') ReadChar(1:24)
    call writerealadv(6,Glob_AttractionScalingParam)
  endif
endif
call MPI_BCAST(Glob_AttrScalParamSupplied,1,MPI_LOGICAL,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_AttractionScalingParam,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)

if (Glob_ProcID==0) then
  read(1,*) ReadChar(1:10),Glob_CurrBasisSize
  write(*,'(1x,a10,1x,i6)')  ReadChar(1:10),Glob_CurrBasisSize
endif
call MPI_BCAST(Glob_CurrBasisSize,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
if (Glob_ProcID==0) read(1,*) ReadChar

!Reading Data Reader and Matrix Calculator Program
ReadChar(1:70)=' '
if (Glob_ProcID==0) then
  Glob_NumOfDRMCSteps=0
  IsDRMCStep=.true.
  do while (IsDRMCStep)
    read(1,*) ReadChar(1:9)
    if (ReadChar(1:9)=='OP_DIPOLE') then
	Glob_NumOfDRMCSteps=Glob_NumOfDRMCSteps+1
    else
      IsDRMCStep=.false.
    endif
  enddo
  do i=1,Glob_NumOfDRMCSteps+1 
    backspace 1
  enddo
endif
call MPI_BCAST(Glob_NumOfDRMCSteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
allocate(Glob_DRMC(1:Glob_NumOfDRMCSteps))

if (Glob_ProcID==0) then
  do i=1,Glob_NumOfDRMCSteps
    read(1,*) Glob_DRMC(i)%Action(1:9)
  enddo
  do i=1,Glob_NumOfDRMCSteps
    backspace 1
  enddo
  do i=1,Glob_NumOfDRMCSteps
    select case (Glob_DRMC(i)%Action(1:9))
    case('OP_DIPOLE')
      read(1,*) Glob_DRMC(i)%Action(1:9),Glob_DRMC(i)%A,   &
            Glob_DRMC(i)%FileName1(1:Glob_FileNameLength), &
            Glob_DRMC(i)%FileName2(1:Glob_FileNameLength), &
            Glob_DRMC(i)%FileName3(1:Glob_FileNameLength), &
            Glob_DRMC(i)%FileName4(1:Glob_FileNameLength)
        j1=len_trim(Glob_DRMC(i)%FileName1(1:Glob_FileNameLength))
        j2=len_trim(Glob_DRMC(i)%FileName2(1:Glob_FileNameLength))
        j3=len_trim(Glob_DRMC(i)%FileName3(1:Glob_FileNameLength)) 
        j4=len_trim(Glob_DRMC(i)%FileName4(1:Glob_FileNameLength))       
     endselect
  enddo
endif

do i=1,Glob_NumOfDRMCSteps
  do j=1,9
    WorkInt(j)=ichar(Glob_DRMC(i)%Action(j:j))
  enddo
  call MPI_BCAST(WorkInt,9,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,9
     Glob_DRMC(i)%Action(j:j)=char(WorkInt(j))
  enddo
  call MPI_BCAST(Glob_DRMC(i)%A,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_DRMC(i)%B,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_DRMC(i)%C,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_DRMC(i)%D,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_DRMC(i)%E,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_DRMC(i)%F,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)  
  call MPI_BCAST(Glob_DRMC(i)%G,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  call MPI_BCAST(Glob_DRMC(i)%H,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)   
  call MPI_BCAST(Glob_DRMC(i)%Q,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  call MPI_BCAST(Glob_DRMC(i)%R,1,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode) 
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_DRMC(i)%FileName1(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_DRMC(i)%FileName1(j:j)=char(WorkInt(j))
  enddo
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_DRMC(i)%FileName2(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_DRMC(i)%FileName2(j:j)=char(WorkInt(j))
  enddo  
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_DRMC(i)%FileName3(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_DRMC(i)%FileName3(j:j)=char(WorkInt(j))
  enddo  
  do j=1,Glob_FileNameLength
    WorkInt(j)=ichar(Glob_DRMC(i)%FileName4(j:j))
  enddo
  call MPI_BCAST(WorkInt,Glob_FileNameLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do j=1,Glob_FileNameLength
     Glob_DRMC(i)%FileName4(j:j)=char(WorkInt(j))
  enddo  
enddo
close(1)

if (Glob_CurrBasisSize==0) then
   !Stop reading if basis size is zero
   return
endif

allocate(Glob_S0(Glob_CurrBasisSize,Glob_CurrBasisSize))
allocate(Glob_S1(Glob_CurrBasisSize,Glob_CurrBasisSize))
allocate(Glob_c0(Glob_CurrBasisSize))
allocate(Glob_FuncNum0(Glob_CurrBasisSize))
allocate(Glob_NonlinParam0(Glob_npt,Glob_CurrBasisSize))
allocate(Glob_c1(Glob_CurrBasisSize))
allocate(Glob_FuncNum1(Glob_CurrBasisSize))
allocate(Glob_ZIndex(Glob_CurrBasisSize))
allocate(Glob_NonlinParam1(Glob_npt,Glob_CurrBasisSize))
	
if (Glob_ProcID==0) then
  do k=1,Glob_NumOfDRMCSteps
    select case (Glob_DRMC(k)%Action(1:9))
    case('OP_DIPOLE')
        
	! Read S matrix for L=0
        ErrorInDataFile=.false.
       	open(1,file=Glob_DRMC(k)%FileName1,status='old',iostat=OpenFileErr)
  	if (OpenFileErr/=0) then
	   write (*,*) 'Error in DataFileInit: data file not found - ',Glob_DRMC(k)%FileName1
           ErrorInDataFile=.true.
  	endif
        if (ErrorInDataFile) stop
	
        do l=1,Glob_CurrBasisSize*Glob_CurrBasisSize
           read(1,*) i,j,Glob_S0(i,j)
        enddo
	close(1)
	
	! Read S matrix for L=1
        ErrorInDataFile=.false.
       	open(1,file=Glob_DRMC(k)%FileName2,status='old',iostat=OpenFileErr)
  	if (OpenFileErr/=0) then
	   write (*,*) 'Error in DataFileInit: data file not found - ',Glob_DRMC(k)%FileName2
           ErrorInDataFile=.true.
  	endif
        if (ErrorInDataFile) stop
	
        do l=1,Glob_CurrBasisSize*Glob_CurrBasisSize
           read(1,*) i,j,Glob_S1(i,j)
        enddo
	close(1)
	
	! Read Symmetry, Wavefunction, and Eigenvector for L=0
        ErrorInDataFile=.false.
       	open(1,file=Glob_DRMC(k)%FileName3,status='old',iostat=OpenFileErr)
  	if (OpenFileErr/=0) then
	   write (*,*) 'Error in DataFileInit: data file not found - ',Glob_DRMC(k)%FileName3
           ErrorInDataFile=.true.
  	endif
        if (ErrorInDataFile) stop
	
	do l=1,4
           read(1,*) ReadChar
        enddo
	
        read(1,*) ReadChar(1:8),Glob_YOperatorString0
        j=len_trim(Glob_YOperatorString0)
        call writestringadv(6,Glob_YOperatorString0,j)
        do i=1,Glob_YOperatorStringLength
           WorkInt(i)=ichar(Glob_YOperatorString0(i:i))
        enddo
	WorkInt0=WorkInt
	
	do l=1,3
           read(1,*) ReadChar
        enddo
	
        do i=1,Glob_CurrBasisSize
           read(1,*) Glob_FuncNum0(i),Glob_c0(i),ReadChar(1:2),Glob_NonlinParam0(1:Glob_npt,i)
        enddo
	close(1)
	
	! Read Symmetry, Wavefunction, and Eigenvector for L=1
        ErrorInDataFile=.false.
       	open(1,file=Glob_DRMC(k)%FileName4,status='old',iostat=OpenFileErr)
  	if (OpenFileErr/=0) then
	   write (*,*) 'Error in DataFileInit: data file not found - ',Glob_DRMC(k)%FileName4
           ErrorInDataFile=.true.
  	endif
        if (ErrorInDataFile) stop
	
	do l=1,4
           read(1,*) ReadChar
        enddo
        read(1,*) ReadChar(1:8),Glob_YOperatorString1
	
        j=len_trim(Glob_YOperatorString1)
        call writestringadv(6,Glob_YOperatorString1,j)
        do i=1,Glob_YOperatorStringLength
           WorkInt(i)=ichar(Glob_YOperatorString1(i:i))
        enddo
	WorkInt1=WorkInt
	
	do l=1,3
           read(1,*) ReadChar
        enddo
	
        do i=1,Glob_CurrBasisSize
           read(1,*) Glob_FuncNum1(i),Glob_c1(i),ReadChar(1:2),Glob_ZIndex(i),Glob_NonlinParam1(1:Glob_npt,i)
        enddo
	close(1)
	
     endselect			     
  enddo
endif
call MPI_BCAST(Glob_S0,Glob_CurrBasisSize*Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_S0,Glob_CurrBasisSize*Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_S1,Glob_CurrBasisSize*Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c0,Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_c1,Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_NonlinParam0,Glob_npt*Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_NonlinParam1,Glob_npt*Glob_CurrBasisSize,MPI_DPREC,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_FuncNum0,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_FuncNum1,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(Glob_ZIndex,Glob_CurrBasisSize,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(WorkInt0,Glob_YOperatorStringLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
call MPI_BCAST(WorkInt1,Glob_YOperatorStringLength,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
do i=1,Glob_YOperatorStringLength
   Glob_YOperatorString0(i:i)=char(WorkInt0(i))
   Glob_YOperatorString1(i:i)=char(WorkInt1(i))
enddo

end subroutine ReadIOFile

subroutine DataInitForAYoungOp(YOpInput,AreYsIdentical)
!Subroutine DataInitForAYoungOp initializes Young operators and on the outpis it gives
!Glob_NumYTerms0, Glob_YCoeff0, Glob_YMatr0, and Glob_NumYTerms1, Glob_YCoeff1, Glob_YMatr1.
!YOpInput = 0 or 1 in the cases L=0 or L=1, respectively.

integer, intent(in)                   :: YOpInput
logical, intent(out)                  :: AreYsIdentical

integer                               :: n,npart
integer                               :: i,j,k,p,q,t,s,w,ii,jj,kk
character(1)                          :: c1
integer                               :: StrLen,NumFactY
integer                               :: TotNumOfYTerms,CurrNumOfTerms
integer                               :: L,R,FirstLPos,LastRPos
integer,allocatable,dimension(:)      :: TempSymCoeff,TempSymCoeff1,NumTermsInYOpFact
integer,allocatable,dimension(:,:,:)  :: TempSymMatr,TempSymMatr1
integer,allocatable,dimension(:,:)    :: Matr1,Matr2,Matr3,Matr4
integer                               :: Coeff,Cf3
logical                               :: AreTermsIdentical
character(Glob_YOperatorStringLength),allocatable,dimension(:) :: YOpStr

if (Glob_ProcID==0) write(*,*) 'Initializing Young operator for L=',YOpInput

select case (YOpInput)
   case(0)
      Glob_YOperatorString=Glob_YOperatorString0
   case(1)
      Glob_YOperatorString=Glob_YOperatorString1
endselect

n=Glob_n
npart=n+1

!Constructing the Young operator based on the content of
!a string variable Glob_YOperatorString

!First we throw away spaces and multiplication signs
!from Glob_YOperatorString
StrLen=len_trim(Glob_YOperatorString)
do i=1,StrLen
  c1=Glob_YOperatorString(i:i)
  if ((c1==' ').or.(c1=='*')) then
    do j=i,StrLen-1
      Glob_YOperatorString(j:j)=Glob_YOperatorString(j+1:j+1)
	enddo
    Glob_YOperatorString(j:j)=' '
  endif
enddo
StrLen=len_trim(Glob_YOperatorString)

!Checking for wrong symbols in Glob_YOperatorString
do i=1,StrLen
  c1=Glob_YOperatorString(i:i)
  if ((c1/='1').and.(c1/='2').and.(c1/='3').and.(c1/='4').and.(c1/='5').and. &
      (c1/='6').and.(c1/='7').and.(c1/='8').and.(c1/='9').and.(c1/='0').and. &
	  (c1/='P').and.(c1/='+').and.(c1/='-').and.(c1/='*').and.(c1/=')').and. &
	  (c1/='(')) then
    write(*,*) 'Error in ProgramDataInit: the Young operator expression'
    write(*,*) 'contains wrong symbols'
    stop
  endif
enddo

!Checking if the number of left and right brackets is the same
!and counting how many brackets there are
L=0
R=0
do i=1,StrLen
  if (Glob_YOperatorString(i:i)==')') R=R+1
  if (Glob_YOperatorString(i:i)=='(') L=L+1
enddo
if (R/=L) then
  write(*,*) 'Error in ProgramDataInit: the numer of left and right brackets in the'
  write(*,*) 'Young operator is different'
  stop
endif

!NumFactY is the number of factors in the Young operator,
!FirstLPos is the position of the first left bracket,
!LastRPos is the position of the last right bracket,
if (R/=0) then
  FirstLPos=scan(Glob_YOperatorString(1:StrLen),'(')
  LastRPos=scan(Glob_YOperatorString(1:StrLen),')',back=.true.)
  NumFactY=R
  i=0
  do j=1,R
	k=0
	i=i+1
    c1=Glob_YOperatorString(i:i)
    do while (c1/='(')
	  if (c1/='*') k=1
	  i=i+1
	  c1=Glob_YOperatorString(i:i)
	enddo 
	if (k==1) NumFactY=NumFactY+1
	do while (c1/=')')
	  i=i+1
	  c1=Glob_YOperatorString(i:i)
	enddo
  enddo
  if (Glob_YOperatorString(StrLen:StrLen)/=')') NumFactY=NumFactY+1
else
  NumFactY=1 
endif

select case (YOpInput)
   case(0)
      Glob_NumFactY0=NumFactY
   case(1)
      Glob_NumFactY1=NumFactY
endselect

!Splitting Glob_YOperatorString into an array of smaller
!strings, YOpStr. Each column of this array will contain just
!one factor, with no brackets. A '+' or a '-' sign is added in 
!front of the first term in a factor if needed. Multiplication
!signs are dropped.
allocate(YOpStr(NumFactY))
YOpStr=' '
if (R==0) then
  c1=Glob_YOperatorString(1:1)
  if ((c1/='+').or.(c1/='-')) then
    YOpStr(1)(1:1)='+'
	YOpStr(1)(2:StrLen+1)=Glob_YOperatorString(1:StrLen)
  else
	YOpStr(1)(1:StrLen)=Glob_YOperatorString(1:StrLen) 
  endif
else
  i=1
  k=1
  p=i
  q=0
  c1=Glob_YOperatorString(i:i)
  if ((c1/='(').and.(c1/='+').and.(c1/='-')) then
    q=1
    YOpStr(k)(1:1)='+'
  endif
  do while (Glob_YOperatorString(i:i)/='(')
    i=i+1
  enddo
  if (i>1) then
    YOpStr(k)(p+q:i-1+q)=Glob_YOperatorString(p:i-1)
    !if (YOpStr(k)(i-1+q:i-1+q)=='*') YOpStr(k)(i-1+q:i-1+q)=' '
    k=k+1
  endif
  do j=1,R
    i=i+1
	p=i
	q=0
    c1=Glob_YOperatorString(i:i)
    if ((c1/=')').and.(c1/='+').and.(c1/='-')) then
      q=1
      YOpStr(k)(1:1)='+'
	endif
	do while (Glob_YOperatorString(i:i)/=')')
      i=i+1
	enddo
    YOpStr(k)(1+q:i+q-p)=Glob_YOperatorString(p:i-1)
    k=k+1
	i=i+1
    c1=Glob_YOperatorString(i:i)
	if (c1=='*') then
      i=i+1
      c1=Glob_YOperatorString(i:i)
	endif
	if ((c1/='(').and.(c1/=' '))  then
	  p=i
      YOpStr(k)(1:1)='+'	
	  do while ((c1/='(').and.(c1/=' '))
	    i=i+1
        c1=Glob_YOperatorString(i:i)
      enddo
      YOpStr(k)(2:i+1-p)=Glob_YOperatorString(p:i-1)
	  k=k+1
    endif
  enddo
endif

! Checking of similarity of Young operators. If these are the same then AreYsIdentical=.true.
AreYsIdentical=.false.
if ((YOpInput==1).and.(Glob_NumFactY0==NumFactY)) then
    if(all(Glob_YOpStr0==YOpStr)) then
       AreYsIdentical=.true.
       deallocate(YOpStr)
       if (Glob_ProcID==0) write(*,*) 'There are same Young operators!'
       return
    end if
endif

!Counting how many terms there are in each factor of the Young
!operator, as well as the total number of terms in the nonsimplified
!Young operator
allocate(NumTermsInYOpFact(NumFactY))
TotNumOfYTerms=1
do k=1,NumFactY
  j=0
  do i=1,Glob_YOperatorStringLength
    if ((YOpStr(k)(i:i)=='+').or.(YOpStr(k)(i:i)=='-')) j=j+1
  enddo
  NumTermsInYOpFact(k)=j
  TotNumOfYTerms=TotNumOfYTerms*j
enddo
if (Glob_ProcID==0) then
  write(*,*)  'Total number of terms in the nonsimplified Y operator:',TotNumOfYTerms
endif

select case (YOpInput)
   case(0)
      allocate(Glob_YOpStr0(NumFactY))
      Glob_YOpStr0=YOpStr
      allocate(Glob_NumTermsInYOpFact0(NumFactY))
      Glob_NumTermsInYOpFact0=NumTermsInYOpFact
   case(1)
      allocate(Glob_NumTermsInYOpFact1(NumFactY))
      Glob_NumTermsInYOpFact1=NumTermsInYOpFact
endselect

!Multiplying all factors in YOpStr and placing actual matrices
!and coefficients in arrays Glob_YMatr and Glob_YCoeff
!One should remember one important fact here: a product of
!of actual pair permutation operators corresponds to the reversed
!product of matrices that act on the matrix of nonlinear parameters.
!Thus, when doing multiplication we will simultaneously be changing
!the order of permutation matrices.

allocate(Matr1(1:n,1:n))
allocate(Matr2(1:n,1:n))
allocate(Matr3(1:n,1:n))
allocate(Matr4(1:n,1:n))

CurrNumOfTerms=NumTermsInYOpFact(NumFactY)
allocate(TempSymCoeff(CurrNumOfTerms)) 
allocate(TempSymMatr(n,n,CurrNumOfTerms))  
do j=NumFactY,1,-1
  !reading the current factor
  k=0
  i=1
  c1=YOpStr(j)(i:i)
  p=i
  do while (c1/=' ')
    i=i+1
    c1=YOpStr(j)(i:i)
    do while ((c1/='P').and.(c1/='+').and.(c1/='-').and.(i<Glob_YOperatorStringLength))
      i=i+1
      c1=YOpStr(j)(i:i)
    enddo
    if (i-p>1) then
      read(YOpStr(j)(p:i-1),*) Coeff
    else
      if (YOpStr(j)(i-1:i-1)=='+') then
        Coeff=1
	  else
        Coeff=-1
	  endif
    endif 
    Matr1=Glob_Transposit(1:n,1:n,1,1)
    do while (c1=='P')
      read(YOpStr(j)(i+1:i+1),*) p
      read(YOpStr(j)(i+2:i+2),*) q
	  Matr2(1:n,1:n)=Glob_Transposit(1:n,1:n,p,q)
	  Matr4(1:n,1:n)=Matr1(1:n,1:n)
	  do ii=1,n
	    do jj=1,n
	      w=0
	      do kk=1,n
	        w=w+Matr2(ii,kk)*Matr4(kk,jj)
	      enddo
	      Matr1(ii,jj)=w
	    enddo
	  enddo
	  i=i+3
      c1=YOpStr(j)(i:i)
    enddo
    k=k+1
    p=i
    if (j/=NumFactY) then
      if (k==1) then
	    Matr3(1:n,1:n)=Matr1(1:n,1:n)
        Cf3=Coeff
	  else
	    do s=1,t
          Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
          q=t*(k-1)+s
	      do ii=1,n
	        do jj=1,n
	          w=0
	          do kk=1,n
	            w=w+Matr2(ii,kk)*Matr1(kk,jj)
	          enddo
	          TempSymMatr(ii,jj,q)=w
	        enddo
	      enddo
	      TempSymCoeff(q)=Coeff*TempSymCoeff(s)
	    enddo
      endif
	else
      TempSymMatr(1:n,1:n,k)=Matr1(1:n,1:n)
      TempSymCoeff(k)=Coeff
    endif
  enddo
  if (j/=NumFactY) then
    do s=1,t
      Matr2(1:n,1:n)=TempSymMatr(1:n,1:n,s)
	  do ii=1,n
	     do jj=1,n
	        w=0
	        do kk=1,n
	          w=w+Matr2(ii,kk)*Matr3(kk,jj)
	        enddo
	        TempSymMatr(ii,jj,s)=w
	     enddo
	  enddo      
	  TempSymCoeff(s)=Cf3*TempSymCoeff(s)
    enddo
  endif
  !mark the identical terms (adding their coefficients
  !and setting all of them but one to zero)
  t=CurrNumOfTerms
  do i=1,CurrNumOfTerms
    if (TempSymCoeff(i)==0) cycle
    do s=i+1,CurrNumOfTerms
      if (TempSymCoeff(s)==0) cycle    
      if (all(TempSymMatr(1:n,1:n,i)==TempSymMatr(1:n,1:n,s))) then
        TempSymCoeff(i)=TempSymCoeff(i)+TempSymCoeff(s)
        if (TempSymCoeff(i)==0) t=t-1
	    TempSymCoeff(s)=0
	    t=t-1
	  endif
    enddo
  enddo     
  !reallocate arrays containing symmetry terms
  !to allow for multiplication by the next factor  
  if (j/=1) then   
    allocate(TempSymCoeff1(t)) 
    allocate(TempSymMatr1(n,n,t)) 
    s=0
    do i=1,CurrNumOfTerms
      if (TempSymCoeff(i)/=0) then
        s=s+1
        TempSymCoeff1(s)=TempSymCoeff(i)
        TempSymMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
      endif
    enddo
    CurrNumOfTerms=t*NumTermsInYOpFact(j-1)
    deallocate(TempSymCoeff)
    deallocate(TempSymMatr)
    allocate(TempSymCoeff(CurrNumOfTerms)) 
    allocate(TempSymMatr(n,n,CurrNumOfTerms))     
    TempSymCoeff(1:t)=TempSymCoeff1(1:t)
    TempSymMatr(1:n,1:n,1:t)=TempSymMatr1(1:n,1:n,1:t)
    deallocate(TempSymCoeff1)
    deallocate(TempSymMatr1)     
  endif    
enddo

Glob_NumYTerms=t
select case (YOpInput)
   case(0)
      Glob_NumYTerms0=Glob_NumYTerms
      allocate(Glob_YCoeff0(Glob_NumYTerms)) 
      allocate(Glob_YMatr0(n,n,Glob_NumYTerms))
   case(1)
      Glob_NumYTerms1=Glob_NumYTerms
      allocate(Glob_YCoeff1(Glob_NumYTerms)) 
      allocate(Glob_YMatr1(n,n,Glob_NumYTerms))
endselect
s=0
do i=1,CurrNumOfTerms
  if (TempSymCoeff(i)/=0) then
    s=s+1
    select case (YOpInput)
      case(0)
         Glob_YCoeff0(s)=TempSymCoeff(i)
         Glob_YMatr0(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
      case(1)
         Glob_YCoeff1(s)=TempSymCoeff(i)
         Glob_YMatr1(1:n,1:n,s)=TempSymMatr(1:n,1:n,i)
    endselect
  endif
enddo
deallocate(TempSymCoeff)
deallocate(TempSymMatr)
if (Glob_ProcID==0) then
  write(*,*)  'Total number of terms in the simplified Y operator:   ',Glob_NumYTerms
endif

deallocate(Matr4)
deallocate(Matr3)
deallocate(Matr2)
deallocate(Matr1)
deallocate(NumTermsInYOpFact)
deallocate(YOpStr)

end subroutine DataInitForAYoungOp

subroutine ProgramDataInit()
!Subroutine ProgramDataInit initializes some data needed for
!calculations. It should be called at the start of the program,
!right after reading input/output file.

integer                               :: n,npart
integer                               :: i,j,k,p,q,t,s,w
integer                               :: pi,pj,pt,ps
logical                               :: AreYsIdentical
real(dprec)                           :: mk,mi,m0
integer,allocatable,dimension(:)      :: IdentParticleSet
integer,allocatable,dimension(:,:)    :: IdentPseudoPartPairSet

if (Glob_ProcID==0) write(*,*) 'Initializing program data'

n=Glob_n
npart=n+1

!Constructing Glob_MassMatrix
allocate(Glob_MassMatrix(n,n))
Glob_MassMatrix(1:n,1:n)=ONEHALF/Glob_Mass(1)
do i=1,n
  Glob_MassMatrix(i,i)=Glob_MassMatrix(i,i)+ONEHALF/Glob_Mass(i+1)
enddo

!Determine the components of vector Glob_bvc
!that is used in evaluation of particle densities
allocate(Glob_bvc(n,npart))
Glob_MassTotal=sum(Glob_Mass(1:npart))
do i=1,npart
  Glob_bvc(1:n,i)=-Glob_Mass(2:n+1)/Glob_MassTotal
enddo
do i=2,npart
  Glob_bvc(i-1,i)=Glob_bvc(i-1,i)+ONE
enddo

!Determine the mass and the index of the the lightest particle 
!(reference particle excluded). 
!and its index
k=0
mk=2*Glob_MassTotal
do i=1,n
  if (Glob_Mass(i+1)<mk) then
    k=i
    mk=Glob_Mass(i+1)
  endif    
enddo  

m0=Glob_Mass(1)
!alpha = sqrt( 0.5 * (m_0^3 + m_k^3)/(m_0*m_k*(m_0 + m_k)^2) )
Glob_dmva2 = (m0**3 + mk**3)/(TWO*m0*mk*(m0+mk)**2) 
!Glob_dmvB(i,i) = (beta^2 + gamma_i^2)/(alpha^2 * M_ii) - M_ii
Glob_dmvB(1:Glob_MaxAllowedNumOfPseudoParticles,1:Glob_MaxAllowedNumOfPseudoParticles)=ZERO
do i=1,n
  mi=Glob_Mass(i+1)
  Glob_dmvB(i,i)=( (m0**3+mi**3)*mk*(m0+mk)**2 - (m0**3+mk**3)*mi*(m0+mi)**2 ) / ( TWO*(m0+mi)*(m0**3+mk**3)*m0*mi**2 )
enddo  
Glob_dmvM(1:Glob_MaxAllowedNumOfPseudoParticles,1:Glob_MaxAllowedNumOfPseudoParticles)=ZERO
Glob_dmvM(1:n,1:n)=Glob_MassMatrix(1:n,1:n)
Glob_dmvMB=Glob_dmvM+Glob_dmvB

!Constructing all Pij transposition matrices 
!
!  P1i  (i/=1) has the following view:
!
!                    i-1
!
!      |  1   0   0  -1   0   0   0  |
!      |                             |
!      |  0   1   0  -1   0   0   0  |
!      |                             |
!      |  0   0   1  -1   0   0   0  |            
!      |                             |          This particular
!P1i = |  0   0   0  -1   0   0   0  |  i-1     matrix is P15 for the 
!      |                             |          case of 8 particles
!      |  0   0   0  -1   1   0   0  |
!      |                             |
!      |  0   0   0  -1   0   1   0  |
!      |                             |
!      |  0   0   0  -1   0   0   1  |
!
!
!  Pij has the following form:
!
!                i-1         j-1
!
!      |  1   0   0   0   0   0   0  |
!      |                             |
!      |  0   1   0   0   0   0   0  |
!      |                             |
!      |  0   0   0   0   0   1   0  |  i-1
!      |                             |
!Pij = |  0   0   0   1   0   0   0  |          This particular  
!      |                             |          matrix is P47 for the
!      |  0   0   0   0   1   0   0  |          case of 8 particles
!      |                             |
!      |  0   0   1   0   0   0   0  |  j-1
!      |                             |
!      |  0   0   0   0   0   0   1  |
!
allocate(Glob_Transposit(n,n,npart,npart))
!First set all of them to be unit matrices
Glob_Transposit(1:n,1:n,1:npart,1:npart)=0
do i=1,npart
  do j=1,npart
    do k=1,n
      Glob_Transposit(k,k,i,j)=1	  
	enddo
  enddo
enddo
!Now continue depending on type of transposition (P1i or Pij)
do i=2,npart
  Glob_Transposit(1:n,i-1,1,i)=-1
enddo
do i=2,npart
  do j=i+1,npart
    Glob_Transposit(i-1,i-1,i,j)=0
	Glob_Transposit(j-1,j-1,i,j)=0
	Glob_Transposit(j-1,i-1,i,j)=1
	Glob_Transposit(i-1,j-1,i,j)=1
    Glob_Transposit(i-1,i-1,j,i)=0
	Glob_Transposit(j-1,j-1,j,i)=0
	Glob_Transposit(j-1,i-1,j,i)=1
	Glob_Transposit(i-1,j-1,j,i)=1
  enddo
enddo

!Calculate Glob_YCoeff and Glob_YMatr for both case L=0 and L=1.
!The order of calling is important, it must be 0,1.
do i=0,1
call DataInitForAYoungOp(i,AreYsIdentical)
enddo
if(AreYsIdentical) then
   Glob_NumFactY1=Glob_NumFactY0
   allocate(Glob_NumTermsInYOpFact1(Glob_NumFactY1))
   Glob_NumTermsInYOpFact1=Glob_NumTermsInYOpFact0
   Glob_NumYTerms1=Glob_NumYTerms0
   allocate(Glob_YCoeff1(Glob_NumYTerms1)) 
   allocate(Glob_YMatr1(n,n,Glob_NumYTerms1))
   Glob_YCoeff1=Glob_YCoeff0
   Glob_YMatr1=Glob_YMatr0
endif

!Now we determine which particles are identical. This determination
!is based on the input values of masses and charges only. The information
!about the sets of identical particles may be needed for proper
!symmetrization of expectation values of operators that involve
!two-paricle quantities (such as interparticle distances).
!The set of particles to which particle i belongs is labelled by IdentParticleSet(i)
!If IdentParticleSet(i)=IdentParticleSet(j) then it means that
!particles i and j are identical. The largest value in IdentParticleSet
!gives the total number of identical particle sets 
allocate(IdentParticleSet(npart))
IdentParticleSet(1)=1
k=1
do i=2,npart
  s=0
  j=0
  do while ((j<i-1).and.(s==0))
    j=j+1
    if (j>1) then
      if ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge(j-1)==Glob_PseudoCharge(i-1))) then
        IdentParticleSet(i)=IdentParticleSet(j)
		s=1
	  endif
    else 
	  !j=1 case
      if ((Glob_Mass(j)==Glob_Mass(i)).and.(Glob_PseudoCharge0==Glob_PseudoCharge(i-1))) then
        IdentParticleSet(i)=IdentParticleSet(j)
		s=1
	  endif
	endif
  enddo	
  if (s==0) then
    k=k+1
    IdentParticleSet(i)=k
  endif
enddo
Glob_NumOfIdentPartSets=maxval(IdentParticleSet(1:npart))

!Below we determine which pairs of pseudoparticles are identical.
!The information about this is stored in array IdentPseudoPartPairSet(1:n,1:n)
!Diagonal elements do not actually designate pairs of pseudoparticles but
!rather a single pseudoparticle, which corresponds to a certain pair of particles. 
!If IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(k,l) then
!it means these pairs ij and kl should be equivalent.
!The largest value of array IdentPseudoPartPairSet gives the number
!of nonequivalent pairs.
allocate(IdentPseudoPartPairSet(1:n,1:n))
IdentPseudoPartPairSet(1:n,1:n)=0
k=0
do i=1,n
  do j=i,n   
    if (i==j) then
      pi=1; pj=j+1
	else
      pi=i+1; pj=j+1
	endif
    w=0
    do s=1,i
	  if (s==i) then
		q=j-1
      else
        q=n
	  endif
      do t=s,q
	    if (w==1) cycle
        if (s==t) then
          ps=1; pt=t+1
		else
          ps=s+1; pt=t+1
		endif
        if ((IdentParticleSet(ps)==IdentParticleSet(pi)).and. &
		    (IdentParticleSet(pt)==IdentParticleSet(pj))) then
          w=1
          IdentPseudoPartPairSet(i,j)=IdentPseudoPartPairSet(s,t)
        endif
	  enddo
	enddo
	if (w==0) then
      k=k+1
      IdentPseudoPartPairSet(i,j)=k
	endif
	enddo
enddo
Glob_NumOfNoneqvPairSets=maxval(IdentPseudoPartPairSet(1:n,1:n))

!Now we create arrays Glob_NumOfPartInIdentPartSet
!and Glob_IdentPartList
allocate(Glob_NumOfPartInIdentPartSet(Glob_NumOfIdentPartSets))
allocate(Glob_IdentPartList(npart,Glob_NumOfIdentPartSets))
Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets)=0
Glob_IdentPartList(1:npart,1:Glob_NumOfIdentPartSets)=0
do i=1,npart
  k=IdentParticleSet(i)
  Glob_NumOfPartInIdentPartSet(k)=Glob_NumOfPartInIdentPartSet(k)+1
  Glob_IdentPartList(Glob_NumOfPartInIdentPartSet(k),k)=i
enddo

!Create arrays Glob_NumOfPairsInEqvPairSet and
!Glob_EqvPairList
allocate(Glob_NumOfPairsInEqvPairSet(Glob_NumOfNoneqvPairSets))
allocate(Glob_EqvPairList(2,n*(n+1)/2,Glob_NumOfNoneqvPairSets))
Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets)=0
Glob_EqvPairList(1:2,1:n*(n+1)/2,1:Glob_NumOfNoneqvPairSets)=0
do i=1,n
  do j=i,n
    k=IdentPseudoPartPairSet(i,j)
    Glob_NumOfPairsInEqvPairSet(k)=Glob_NumOfPairsInEqvPairSet(k)+1
    Glob_EqvPairList(1,Glob_NumOfPairsInEqvPairSet(k),k)=i
    Glob_EqvPairList(2,Glob_NumOfPairsInEqvPairSet(k),k)=j
  enddo
enddo

deallocate(IdentPseudoPartPairSet)
deallocate(IdentParticleSet)

end subroutine ProgramDataInit

subroutine ComputeExpValL0L1()

!Subroutine ComputeExpValL0L1 computes expectation value
!               < L=0 | Operator | L=1 >
!of the corrections of Hamiltonian for the case Glob_CurrDRMCStep with
!basis functions whose number ranges from 1 to Glob_CurrBasisSize.
!The solution is Glob_ExpVals(Glob_CurrDRMCStep).
!
!Expressions:
!
!       cs = Glob_CurrDRMCStep
!       k,l = 1,Glob_CurrBasisSize
!       N0 = c0 * S0 * c0
!       N1 = c1 * S1 * c1
!       Glob_ExpVals(cs) = (Sum_{k,l} c0_k * H_kl * Glob_c1(l))/sqrt(N0*N1)
!       i = 1,Glob_NumYTerms0 / j = 1,Glob_NumYTerms1
!       H_kl = Sum_{i,j} YCoeff0_i * H_klij * YCoeff1_j
!       H_klij = < f_ki | Operator | f_lj >/sqrt(<f_k|f_k>*<f_l|f_l>)

integer     :: n,np
real(dprec) :: ExpVal,ExpValLoc
real(dprec) :: Hklij,Hkl

integer     :: k,l,i,j,indx
real(dprec) :: temp0,temp0Loc,temp1,temp1Loc

n=Glob_n
np=Glob_np

if(Glob_DRMC(Glob_CurrDRMCStep)%A==0) then
   ExpValLoc=ZERO
   temp0Loc=ZERO
   temp1Loc=ZERO
   indx=0
   do k=1,Glob_CurrBasisSize
      do l=1,Glob_CurrBasisSize
        indx=indx+1
        if(mod(indx,Glob_NumOfProcs)==Glob_ProcID) then
          Hkl=ZERO
          temp0Loc=temp0Loc+Glob_c0(k)*Glob_S0(k,l)*Glob_c0(l)
          temp1Loc=temp1Loc+Glob_c1(k)*Glob_S1(k,l)*Glob_c1(l)
          do i=1,Glob_NumYTerms0
            do j=1,Glob_NumYTerms1
   	       call MatrixElementsL0L1(Glob_ZIndex(l),Glob_NonlinParam0(1:np,k), &
   	            Glob_NonlinParam1(1:np,l),Glob_YMatr0(1:n,1:n,i),Glob_YMatr1(1:n,1:n,j),Hklij)
   	       Hkl=Hkl+Glob_YCoeff0(i)*Hklij*Glob_YCoeff1(j)
            enddo
          enddo
          ExpValLoc=ExpValLoc+Glob_c0(k)*Hkl*Glob_c1(l)
        endif
      enddo
   enddo
   call MPI_ALLREDUCE(ExpValLoc,ExpVal,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
   call MPI_ALLREDUCE(temp0Loc,temp0,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
   call MPI_ALLREDUCE(temp1Loc,temp1,1,MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
   
   temp0=sqrt(temp0*temp1)
   Glob_ExpVals(Glob_CurrDRMCStep)=ExpVal/temp0
else
   if (Glob_ProcID==0) write(*,*) 'Glob_DRMC(Glob_CurrDRMCStep)%A is not 0!'
endif

end subroutine ComputeExpValL0L1

end module workproc
