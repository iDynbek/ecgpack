module gpu_backend
!High-level orchestration of the CUDA backend. Everything GPU-specific that is
!NOT a raw C binding lives here: run-time selection (env vars), the device
!context lifecycle, the chunked matrix-element builds, and the GPU eigensolver.
!
!This module is compiled ONLY with -DUSE_CUDA, so the shared Fortran sources
!(globvars, matform, workproc, main) carry only a thin, #ifdef-guarded dispatch
!and are otherwise identical to a CPU-only build. The chunk loops call back into
!matform's StoreHS/StoreHSD through procedure arguments, so this module does not
!depend on matform (which would be circular) yet reuses its storage routines.
  use globvars      !Glob_* state, MPI symbols, wp / MPI_WP (via wp_def)
  use matelem_gpu   !ISO_C_BINDING interfaces to matelem_cuda.cu / eigen_cuda.cu
#ifdef USE_CUF
  use matelem_cuf   !SPIKE: native CUDA Fortran energy kernel (replaces the C++ one)
#endif
  implicit none
  private

  !Public surface used by the shared sources (every call site is #ifdef USE_CUDA):
  public :: gpu_backend_init      !collective one-time startup            (main)
  public :: gpu_active            !.true. => use the CUDA ME backend      (matform)
  public :: gpu_eig_active        !.true. => use the cuSOLVER eigensolver (workproc)
  public :: gpu_build_HS          !energy matrix-element build            (matform)
  public :: gpu_build_HS_deriv    !energy+gradient matrix-element build   (matform)
  public :: gpu_dsygvx            !re-exported cuSOLVER eigensolve         (workproc)

  logical, save :: use_me    = .false.   !ECG_GPU=1
  logical, save :: use_eig   = .false.   !ECG_GPU_EIG=1 (implies use_me)
  integer, save :: batch_cap = 16384     !ECG_GPU_BATCH: max pairs per GPU call

  !Signatures of matform's StoreHS / StoreHSD, matched exactly so those routines
  !can be handed in as callbacks. The chunk loops call them one pair at a time.
  abstract interface
    subroutine store_hs_cb(i,j,Hij,Sij)
      import :: wp
      integer,  intent(in) :: i,j
      real(wp), intent(in) :: Hij,Sij
    end subroutine store_hs_cb
    subroutine store_hsd_cb(i,j,Hij,Sij,Di,Dj)
      import :: wp, Glob_npt
      integer,  intent(in) :: i,j
      real(wp), intent(in) :: Hij,Sij
      real(wp), intent(in) :: Di(2*Glob_npt),Dj(2*Glob_npt)
    end subroutine store_hsd_cb
  end interface

contains

  subroutine gpu_backend_init()
  !Read the env knobs on rank 0 and broadcast them so every rank agrees (the
  !on/off flag and the chunk size must match, or the per-chunk collective
  !ALLREDUCEs would deadlock/mismatch). Bring up the device context if enabled.
  !Called once, collectively, at startup.
    character(64) :: val
    integer       :: ios, tmp
    if (Glob_ProcID==0) then
      call get_environment_variable('ECG_GPU', val, status=ios)
      use_me = (ios==0) .and. (val(1:1)=='1')
      call get_environment_variable('ECG_GPU_EIG', val, status=ios)
      use_eig = use_me .and. (ios==0) .and. (val(1:1)=='1')
      call get_environment_variable('ECG_GPU_BATCH', val, status=ios)
      if (ios==0) then
        read(val,*,iostat=ios) tmp
        if ((ios==0).and.(tmp>0)) batch_cap = tmp
      endif
      if (use_me)  write(*,'(1x,a)') 'GPU matrix-element backend ENABLED (ECG_GPU=1)'
      if (use_eig) write(*,'(1x,a)') 'GPU eigensolver (cuSOLVER) ENABLED (ECG_GPU_EIG=1)'
    endif
    call MPI_BCAST(use_me,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(use_eig,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    call MPI_BCAST(batch_cap, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Glob_MPIErrCode)
    if (use_me) call gpu_init(Glob_ProcID)   !device ctx; registers atexit(finalize)
  end subroutine gpu_backend_init

  logical function gpu_active()
    gpu_active = use_me
  end function gpu_active

  logical function gpu_eig_active()
    gpu_eig_active = use_eig
  end function gpu_eig_active

  subroutine gpu_build_HS(Nmin,Nmax,store)
  !Energy-only matrix-element build over the basis pair triangle [Nmin,Nmax],
  !processed in fixed-size chunks so neither host nor device memory scales with
  !the whole K(K+1)/2 pair list. Each pair's summed H/S is handed back via the
  !store callback (matform's StoreHS). Every rank uses the same batch_cap, so the
  !per-chunk ALLREDUCE stays collective.
    integer, intent(in)    :: Nmin,Nmax
    procedure(store_hs_cb) :: store
    integer               :: n,np,k,l,npairs,ipair,nf,ip,batch
    integer, allocatable  :: k_list(:),l_list(:)
    real(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    n=Glob_n; np=Glob_np
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    batch  = min(npairs, batch_cap)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (energy): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    allocate(k_list(batch),l_list(batch),Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    ipair=0; nf=0
    do k=Nmin,Nmax
      do l=k,1,-1
        ipair=ipair+1; nf=nf+1
        k_list(nf)=k; l_list(nf)=l
        if ((nf==batch).or.(ipair==npairs)) then
#ifdef USE_CUF
          call cuf_compute_matelem_batch( &                              !SPIKE: CUDA Fortran energy kernel
              Glob_NonlinParam(1:np,1:Nmax), Nmax, k_list, l_list, nf, &
              Glob_YHYMatr(1,1,1), Glob_YHYCoeff, Glob_NumYHYTerms, &
              n, np, Glob_NumOfProcs, Glob_ProcID, &
              Glob_PiRaised3n2, Glob_MassMatrix(1:n,1:n), &
              Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
              Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
              Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
              Hout, Sout)
#else
          call gpu_compute_matelem_batch( &
              Glob_NonlinParam(1:np,1:Nmax), Nmax, k_list, l_list, nf, &
              Glob_YHYMatr(1,1,1), Glob_YHYCoeff, Glob_NumYHYTerms, &  !whole array by element ref (a (:,:,1) section makes nvfortran copy only term 1)
              n, np, Glob_NumOfProcs, Glob_ProcID, &
              Glob_PiRaised3n2, Glob_MassMatrix(1:n,1:n), &
              Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
              Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
              Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
              Hout, Sout)
#endif
          call MPI_ALLREDUCE(Hout,Hout_r,nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout,Sout_r,nf,MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          do ip=1,nf
            call store(k_list(ip),l_list(ip),Hout_r(ip),Sout_r(ip))
          enddo
          nf=0
        endif
      enddo
    enddo
    deallocate(k_list,l_list,Hout,Sout,Hout_r,Sout_r)
  end subroutine gpu_build_HS

  subroutine gpu_build_HS_deriv(Nmin,Nmax,store)
  !Energy+gradient build. The (2*np x pairs) gradient buffers are the memory
  !wall at large K, so the same chunking applies. Each pair's H/S and gradient
  !slabs are handed back via the store callback (matform's StoreHSD).
    integer, intent(in)     :: Nmin,Nmax
    procedure(store_hsd_cb) :: store
    integer               :: n,np,npt2,k,l,npairs,ipair,nf,ip,batch
    integer, allocatable  :: k_list(:),l_list(:),grad_l(:)
    real(wp),allocatable  :: Hout(:),Sout(:),Hout_r(:),Sout_r(:)
    real(wp),allocatable  :: Dkout(:,:),Dlout(:,:),Dkout_r(:,:),Dlout_r(:,:)
    n=Glob_n; np=Glob_np; npt2=np*2
    npairs = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
    batch  = min(npairs, batch_cap)
    if ((Glob_ProcID==0).and.(npairs>batch)) &
      write(*,'(1x,a,i0,a,i0,a,i0,a)') 'GPU ME (grad): ',npairs, &
        ' pairs in ',(npairs+batch-1)/batch,' chunks of ',batch,' max'
    allocate(k_list(batch),l_list(batch),grad_l(batch))
    allocate(Hout(batch),Sout(batch),Hout_r(batch),Sout_r(batch))
    allocate(Dkout(npt2,batch),Dlout(npt2,batch),Dkout_r(npt2,batch),Dlout_r(npt2,batch))
    ipair=0; nf=0
    do k=Nmin,Nmax
      do l=k,1,-1
        ipair=ipair+1; nf=nf+1
        k_list(nf)=k; l_list(nf)=l
        if ((l>Glob_nfru).and.(l/=k)) then
          grad_l(nf)=1
        else
          grad_l(nf)=0
        endif
        if ((nf==batch).or.(ipair==npairs)) then
          Dkout(1:npt2,1:nf)=ZERO
          Dlout(1:npt2,1:nf)=ZERO
#ifdef USE_CUF
          call cuf_compute_matelem_deriv_batch( &                         !SPIKE: CUDA Fortran gradient kernel
              Glob_NonlinParam(1:np,1:Nmax), Nmax, k_list, l_list, grad_l, nf, &
              Glob_YHYMatr(1,1,1), Glob_YHYCoeff, Glob_NumYHYTerms, &
              n, np, Glob_NumOfProcs, Glob_ProcID, &
              Glob_PiRaised3n2, Glob_MassMatrix(1:n,1:n), &
              Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
              Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
              Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
              Hout, Sout, Dkout, Dlout)
#else
          call gpu_compute_matelem_and_deriv_batch( &
              Glob_NonlinParam(1:np,1:Nmax), Nmax, k_list, l_list, grad_l, nf, &
              Glob_YHYMatr(1,1,1), Glob_YHYCoeff, Glob_NumYHYTerms, &  !whole array by element ref (see note in gpu_build_HS)
              n, np, Glob_NumOfProcs, Glob_ProcID, &
              Glob_PiRaised3n2, Glob_MassMatrix(1:n,1:n), &
              Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
              Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
              Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
              Hout, Sout, Dkout, Dlout)
#endif
          call MPI_ALLREDUCE(Hout, Hout_r, nf,      MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Sout, Sout_r, nf,      MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Dkout,Dkout_r,nf*npt2, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          if (Glob_nfo>1) &
          call MPI_ALLREDUCE(Dlout,Dlout_r,nf*npt2, MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          do ip=1,nf
            call store(k_list(ip),l_list(ip),Hout_r(ip),Sout_r(ip), &
                       Dkout_r(1:npt2,ip),Dlout_r(1:npt2,ip))
          enddo
          nf=0
        endif
      enddo
    enddo
    deallocate(k_list,l_list,grad_l)
    deallocate(Hout,Sout,Hout_r,Sout_r)
    deallocate(Dkout,Dlout,Dkout_r,Dlout_r)
  end subroutine gpu_build_HS_deriv

end module gpu_backend
