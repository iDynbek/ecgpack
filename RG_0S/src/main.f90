program main

  use workproc
#ifdef USE_CUDA
  use gpu_backend
#endif
  implicit none

!Local variables
  integer      i,iw,Kstart,Kstop,Kstep,OpenFileErr,OptimizationType
  real(8)      r8

!These variables are used to set random number generators
  integer RNSeedSize
  integer,allocatable :: Seed(:)
  integer       :: seedval, seedios
  character(32) :: seedstr
  character(32) :: gcstr
  integer       :: gcflag, gcios
  logical       :: do_gradcheck

!Initialize MPI
  call MPI_INIT(Glob_MPIErrCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Glob_ProcID,Glob_MPIErrCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_NumOfProcs,Glob_MPIErrCode)
  call system_clock(Glob_WallTick0)   !program-start wall tick for GROWTRACE

  if (Glob_ProcID==0) then
    write (*,*) 'Program ',Glob_BasisType,' has started'
    write (*,*) 'Number of parallel processes running ',Glob_NumOfProcs
    write (*,*)
  endif

#ifdef USE_CUDA
  call gpu_backend_init()   !env-select the GPU backend + bring up the device context (collective)
#endif
  call ReadIOFile()
  if (Glob_IsOptCycleScripted) call ReadBlackList()
  call ProgramDataInit()

!Seed the random number generators
  call random_seed()
  call random_seed(size=RNSeedSize)
  allocate(Seed(RNSeedSize))
  call random_seed(get=Seed(1:RNSeedSize))
!ECG_SEED=<int>: deterministic RNG (all of random_number + drnor_start derive from
!Seed(1)), so generation runs are reproducible and CPU-vs-GPU growth is comparable.
!Unset -> seed from the clock, as before.
  seedval=-1
  if (Glob_ProcID==0) then
    call get_environment_variable('ECG_SEED', seedstr, status=seedios)
    if (seedios==0) read(seedstr,*,iostat=seedios) seedval
    if (seedios/=0) seedval=-1
  endif
  call MPI_BCAST(seedval,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  if (seedval>=0) then
    do i=1,RNSeedSize          !fill ALL seed elements deterministically -- setting
      Seed(i)=seedval+37*(i-1)+1   !only Seed(1) leaves the rest OS-random (random_seed
    enddo                          !(get=) filled them from entropy) -> non-reproducible
    if (Glob_ProcID==0) write(*,'(1x,a,i0)') 'Deterministic RNG seed ECG_SEED=',seedval
  else
    call system_clock(count=Seed(1))
  endif
  Seed(1:RNSeedSize)=Seed(1:RNSeedSize)+Glob_ProcID
  call random_seed(put=Seed(1:RNSeedSize))
  call random_number(r8)
  call random_number(r8)
  call random_number(r8)
  call random_number(r8)
  r8=drnor_start(nint(r8*25000)+Glob_ProcID)

!Empty swap file as it may contain some garbage left
!after last run (if there was a failure)
  if ((Glob_ProcID==0).and.Glob_UseSwapFile) then
    open(1,file=Glob_SwapFileName,form='unformatted',status='replace',iostat=OpenFileErr)
    if (OpenFileErr==0) then
      write(1) 'Swap file is empty'
      close(1)
    endif
  endif

  OptimizationType=1

!ECG_GRADCHECK: if set, run the finite-difference gradient-correctness check on the
!loaded basis (subroutine GradCheck) instead of the scripted BBOP steps, then fall
!through to the normal finalize (profiling summary + MPI_FINALIZE).
  do_gradcheck=.false.
  if (Glob_ProcID==0) then
    call get_environment_variable('ECG_GRADCHECK', gcstr, status=gcios)
    if (gcios==0) do_gradcheck=.true.
  endif
  gcflag=0
  if (do_gradcheck) gcflag=1
  call MPI_BCAST(gcflag,1,MPI_INTEGER,0,MPI_COMM_WORLD,Glob_MPIErrCode)
  do_gradcheck=(gcflag==1)

  if (do_gradcheck) then
    call GradCheck()
  else

  do i=1,Glob_NumOfBBOPSteps
    call linalg_setparam(Glob_CurrBasisSize)

    Glob_CurrBBOPStep=i
    Glob_ApproxEnergy=Glob_CurrEnergy*Glob_InvItParameter

    select case (Glob_BBOP(i)%Action)

    case('BASIS_ENL')
      Kstart=Glob_BBOP(i)%A
      Kstop=Glob_BBOP(i)%B
      Kstep=Glob_BBOP(i)%C
      if (Kstop>Glob_CurrBasisSize) then
        if (Kstart<=Glob_CurrBasisSize+1) then
          Kstart=Glob_CurrBasisSize+1
          select case (Glob_BBOP(i)%GSEPSolutionMethod)
          case('G')
            call BasisEnlG(Kstart,Kstop,Kstep,Glob_BBOP(i)%D,OptimizationType, &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R)
          case('I')
            call BasisEnlI(Kstart,Kstop,Kstep,Glob_BBOP(i)%D,OptimizationType, &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%Q,Glob_BBOP(i)%R)
          endselect
        else
          if (Glob_ProcID==0) then
            write(*,*) 'Error in main: incorrect BBOP step ',i
            write(*,*) 'One or more parameters in BASIS_ENL are incorrect'
          endif
        endif
      endif

    case('OPT_CYCLE')
      if (Glob_ProcID==0) then
        if ((Glob_BBOP(i)%C>Glob_BBOP(i)%A).or.(Glob_BBOP(i)%B>Glob_BBOP(i)%C)) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'One or more parameters in OPT_CYCLE are incorrect'
        endif
      endif
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%B>0).and.  &
          (Glob_BBOP(i)%C<=Glob_CurrBasisSize)) then
        if (Glob_History(Glob_CurrBasisSize)%CyclesDone<Glob_BBOP(i)%F) then
          select case (Glob_BBOP(i)%GSEPSolutionMethod)
          case('G')
            call OptCycleG(Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,  &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
                           Glob_BBOP(i)%H)
          case('I')
            call OptCycleI(Glob_BBOP(i)%A,Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,  &
                           Glob_BBOP(i)%E,Glob_BBOP(i)%F,Glob_BBOP(i)%G,Glob_BBOP(i)%Q,Glob_BBOP(i)%R, &
                           Glob_BBOP(i)%H)
          endselect
        endif
      endif

    case('FULL_OPT1')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%B>0).and.  &
          (Glob_BBOP(i)%C<=Glob_CurrBasisSize)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call FullOpt1G(Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q, &
                         Glob_BBOP(i)%R,real(Glob_BBOP(i)%E,4),real(Glob_BBOP(i)%F,4), &
                         Glob_BBOP(i)%FileName1)
        case('I')
          call FullOpt1I(Glob_BBOP(i)%B,Glob_BBOP(i)%C,Glob_BBOP(i)%D,Glob_BBOP(i)%Q, &
                         Glob_BBOP(i)%R,real(Glob_BBOP(i)%E,4),real(Glob_BBOP(i)%F,4), &
                         Glob_BBOP(i)%FileName1)
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'One or more parameters in FULL_OPT1 are incorrect'
        endif
      endif

    case('ELIM_LCFN')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call EliminateLittleContribFunc(Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1, &
                                          Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in ELIM_LCFN'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'One or more parameters in ELIM_LCFN are incorrect'
        endif
      endif

    case('ELIM_LND1')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call EliminateLinDepFunc(Glob_BBOP(i)%Q,Glob_BBOP(i)%FileName1,Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in ELIM_LND1'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'One or more parameters in ELIM_LND1 are incorrect'
        endif
      endif

    case('SEPR_LND1')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call SeparateLinDepFunc(Glob_BBOP(i)%Q,Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1, &
                                  Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in SEPR_LND1'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'One or more parameters in SEPR_LND1 are incorrect'
        endif
      endif

    case('SEPR_FLCF')
      if ((Glob_BBOP(i)%A==Glob_CurrBasisSize).and.(Glob_BBOP(i)%Q>ZERO)) then
        select case (Glob_BBOP(i)%GSEPSolutionMethod)
        case('G')
          call SeparateFuncLargeCoeff(Glob_BBOP(i)%Q,Glob_BBOP(i)%R,Glob_BBOP(i)%FileName1, &
                                      Glob_ElimRoutPrintSpec)
        case('I')
          if (Glob_ProcID==0) write(*,*) 'Sorry, GSEP soluton method I does not work in SEPR_LND1'
        endselect
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'One or more parameters in SEPR_FLCF are incorrect'
        endif
      endif

    case('EXPC_VALS')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        call ExpectationValues(Glob_BBOP(i)%Action,1,Glob_FileNameNone,Glob_FileNameNone,Glob_FileNameNone, &
                               Glob_FileNameNone,Glob_BBOP(i)%GSEPSolutionMethod)
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'Second parameter in EXPC_VALS is incorrect'
        endif
      endif

    case('DENSITIES')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        call ExpectationValues(Glob_BBOP(i)%Action,1,Glob_BBOP(i)%FileName1,Glob_BBOP(i)%FileName2, &
                               Glob_BBOP(i)%FileName3,Glob_BBOP(i)%FileName4,Glob_BBOP(i)%GSEPSolutionMethod)
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'Second parameter in DENSITIES is incorrect'
        endif
      endif

    case('MOMT_DENS')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        call ExpectationValues(Glob_BBOP(i)%Action,1,Glob_BBOP(i)%FileName1,Glob_BBOP(i)%FileName2, &
                               Glob_BBOP(i)%FileName3,Glob_BBOP(i)%FileName4,Glob_BBOP(i)%GSEPSolutionMethod)
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'Second parameter in MOMT_DENS is incorrect'
        endif
      endif

    case('SAVE_FILE')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        if (Glob_ProcID==0) then
          iw=len_trim(Glob_BBOP(i)%FileName1(1:Glob_FileNameLength))
          write(*,*)
          write(*,*) 'Saving basis in file ',Glob_BBOP(i)%FileName1(1:iw),'...'
          call SaveResults(Filename=Glob_BBOP(i)%FileName1,Sort='no')
          write(*,*) ' done'
          write(*,*)
        endif
      endif

    case('SAVE_HSWF')
      if (Glob_BBOP(i)%A==Glob_CurrBasisSize) then
        call SaveHSWF(Glob_BBOP(i)%FileName1,Glob_BBOP(i)%FileName2, &
                      Glob_BBOP(i)%FileName3,Glob_BBOP(i)%FileName4, &
                      Glob_BBOP(i)%GSEPSolutionMethod)
      else
        if (Glob_ProcID==0) then
          write(*,*) 'Error in main: incorrect BBOP step ',i
          write(*,*) 'Second parameter in SAVE_HSEV is incorrect'
        endif
      endif

    endselect

    !PROFILING: after each BBOP step, report basis size K and the
    !cumulative / last-single wall time in matrix-element builds vs
    !the eigensolve (added for the CUDA-acceleration study).
    if (Glob_ProcID==0) then
      write(*,'(1x,a,i4,1x,a,a,i6,a,f10.3,a,f10.3,a,e11.3,a,e11.3)') &
        'PROF step',i,Glob_BBOP(i)%Action,' K=',Glob_CurrBasisSize, &
        ' | cumME=',Glob_TimeME,'s cumEIG=',Glob_TimeEIG, &
        's | lastME=',Glob_LastME,'s lastEIG=',Glob_LastEIG
      !GROWTRACE now emitted per accepted function from inside
      !BasisEnlG/BasisEnlI (see globvars::GrowTrace) so a single
      !continuous BASIS_ENL step is traced without a warm-start reset.
    endif
  enddo

  endif   !do_gradcheck branch

  if (Glob_ProcID==0) then
    if (Glob_UseSwapFile) then
      !Empty swap file as it may contain the data saved after
      !the last basis building and optimization program step
      open(1,file=Glob_SwapFileName,form='unformatted',status='replace',iostat=OpenFileErr)
      if (OpenFileErr==0) then
        write(1) 'Swap file is empty'
        close(1)
        write(*,*)
        write(*,*) 'Swap file has been cleaned up'
      endif
    endif
    write(*,*) ' '
    write(*,*) 'Basis Building and Optimization Program is completed'
    write(*,*) 'Program has stopped'

    write(*,*)
    write(*,*) '================ PROFILING SUMMARY ================'
    write(*,'(1x,a,f12.3,a,i12,a)') 'Matrix-element build : ',Glob_TimeME, &
        ' s over ',Glob_CntME,' calls'
    write(*,'(1x,a,f12.3,a,i12,a)') 'Eigensolve (G:DSYGVX/I:GSEPIIS): ',Glob_TimeEIG, &
        ' s over ',Glob_CntEIG,' calls'
    if (Glob_TimeME+Glob_TimeEIG>0.0_wp) &
      write(*,'(1x,a,f6.2,a,f6.2,a)') 'Share                : ME ', &
        100*Glob_TimeME/(Glob_TimeME+Glob_TimeEIG),'%   EIG ', &
        100*Glob_TimeEIG/(Glob_TimeME+Glob_TimeEIG),'%'
    write(*,*) '==================================================='
    write(*,*)
  endif

#ifdef USE_CUDA
  call gpu_finalize()   !release cuSOLVER handle + device buffers, reset device
#endif

  call MPI_FINALIZE(Glob_MPIErrCode)

end program main
