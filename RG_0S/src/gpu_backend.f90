module gpu_backend
!The complete GPU backend for RG_0S -- both the high-level orchestration AND the
!native CUDA Fortran device code, in one module. Compiled ONLY with -DUSE_CUDA,
!so the shared Fortran sources (globvars, matform, workproc, main) carry only a
!thin, #ifdef-guarded dispatch and are otherwise identical to a CPU-only build.
!
!Two layers live here, separated by the section banners below:
!  (1) ORCHESTRATION -- run-time selection (env vars), the device-context
!      lifecycle, and the MPI-parallel chunked matrix-element builds. The chunk
!      loops call back into matform's StoreHS/StoreHSD through procedure arguments,
!      so this module does not depend on matform (which would be circular) yet
!      reuses its storage routines.
!  (2) DEVICE CODE -- the attributes(global) matrix-element kernels (which call the
!      SHARED MatrixElements device routine in matelem.f90, so the physics lives
!      once, not duplicated), the device lifecycle, and the cuSOLVER eigensolver.
!      Physics constants are passed as kernel arguments -- device code cannot read
!      host module globals.
  use globvars        !Glob_* state, MPI symbols, wp / MPI_WP (via wp_def)
  use matelem,   only: MatrixElements
  use cudafor
  use cublas
  use cusolverDn
  implicit none
  private

  !Public surface used by the shared sources (every call site is #ifdef USE_CUDA):
  public :: gpu_backend_init      !collective one-time startup            (main)
  public :: gpu_active            !.true. => use the CUDA ME backend      (matform)
  public :: gpu_eig_active        !.true. => use the cuSOLVER eigensolver (workproc)
  public :: gpu_build_HS          !energy matrix-element build            (matform)
  public :: gpu_build_HS_deriv    !energy+gradient matrix-element build   (matform)
  public :: gpu_dsygvx            !cuSOLVER eigensolve                    (workproc)

  logical, save :: use_me    = .false.   !ECG_GPU=1
  logical, save :: use_eig   = .false.   !ECG_GPU_EIG=1 (implies use_me)
  integer, save :: batch_cap = 16384     !ECG_GPU_BATCH: max pairs per GPU call

  integer, parameter :: NN  = Glob_AllowedNumOfPseudoParticles
  integer, parameter :: NNP = NN*(NN+1)/2   !max packed vech length; Dk/Dl are 2*np long (np=n(n+1)/2)
  integer, parameter :: CUF_BLK = 128       !threads/block. The term loop STRIDES (see the kernels), so
                                            !block size is independent of NumYHYTerms. The shared
                                            !MatrixElements costs ~254 regs/thread on device; 254*128 <
                                            !65536 (V100 per-block reg limit), so heavy atoms (large
                                            !NumYHYTerms, e.g. Carbon=720) launch instead of failing.
                                            !One-thread-per-term busted this at ~258 terms.

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

  ! ==========================================================================
  !  (1) ORCHESTRATION -- env selection, lifecycle, chunked MPI builds
  ! ==========================================================================

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
    if (use_me) call gpu_init(Glob_ProcID)   !device ctx
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
          call cuf_compute_matelem_batch( &
              Glob_NonlinParam(1:np,1:Nmax), Nmax, k_list, l_list, nf, &
              Glob_YHYMatr(1,1,1), Glob_YHYCoeff, Glob_NumYHYTerms, &
              n, np, Glob_NumOfProcs, Glob_ProcID, &
              Glob_PiRaised3n2, Glob_MassMatrix(1:n,1:n), &
              Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
              Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
              Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
              Hout, Sout)
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
          call cuf_compute_matelem_deriv_batch( &
              Glob_NonlinParam(1:np,1:Nmax), Nmax, k_list, l_list, grad_l, nf, &
              Glob_YHYMatr(1,1,1), Glob_YHYCoeff, Glob_NumYHYTerms, &
              n, np, Glob_NumOfProcs, Glob_ProcID, &
              Glob_PiRaised3n2, Glob_MassMatrix(1:n,1:n), &
              Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
              Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
              Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
              Hout, Sout, Dkout, Dlout)
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

  ! ==========================================================================
  !  (2) DEVICE CODE -- kernels, host launchers, lifecycle, eigensolver
  !  (all private except gpu_dsygvx)
  ! ==========================================================================

  !Energy kernel: one block per (k,l) pair; threads STRIDE over the permutation
  !terms. Each thread calls the shared MatrixElements (grad=false) for its terms.
  attributes(global) subroutine me_energy_kernel(NonlinParam, np, Kmax, &
      k_list, l_list, YHYMatr, YHYCoeff, nterms, n, nprocs, procid, &
      mass, charge, charge0, sqrtpi, pir3n2, attr, rep, repp, repm, Hout, Sout)
    integer, value   :: np, Kmax, nterms, n, nprocs, procid
    real(wp)         :: NonlinParam(np,Kmax)
    integer          :: k_list(*), l_list(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    real(wp)         :: mass(n,n), charge(n)
    real(wp), value  :: charge0, sqrtpi, pir3n2, attr, rep, repp, repm
    real(wp)         :: Hout(*), Sout(*)
    real(wp), shared :: sh_H, sh_S
    integer  :: pair_idx, j, qbase, k0, l0, istat
    real(wp) :: Hkl, Skl, coeff
    real(wp) :: dDk(2*NNP), dDl(2*NNP)   !unused energy-path gradient outputs

    pair_idx = blockIdx%x
    if (threadIdx%x==1) then
      sh_H=0.0_wp; sh_S=0.0_wp
    endif
    call syncthreads()

    qbase = (pair_idx-1)*nterms - 1          !matches Fortran q=(i-1)*nterms-1
    k0    = k_list(pair_idx)                 !pair-only: hoist out of the term loop
    l0    = l_list(pair_idx)
    do j = threadIdx%x, nterms, blockDim%x   !STRIDED: block size independent of nterms
      if (mod(qbase+j, nprocs) == procid) then
        call MatrixElements(n, np, NonlinParam(1,k0), NonlinParam(1,l0), &
                            YHYMatr((j-1)*n*n+1), mass, charge, charge0, &
                            sqrtpi, pir3n2, attr, rep, repp, repm, &
                            Hkl, Skl, dDk, dDl, .false., .false.)
        coeff = YHYCoeff(j)
        istat = atomicadd(sh_H, coeff*Hkl)
        istat = atomicadd(sh_S, coeff*Skl)
      endif
    enddo
    call syncthreads()

    if (threadIdx%x==1) then
      Hout(pair_idx)=sh_H; Sout(pair_idx)=sh_S
    endif
  end subroutine me_energy_kernel

  !Host launcher for the energy build. Stages inputs to the device, launches, copies back.
  subroutine cuf_compute_matelem_batch(NonlinParam, Kmax, k_list, l_list, npairs, &
      YHYMatr, YHYCoeff, nterms, n, np, nprocs, procid, &
      pir3n2, mass, charge, charge0, attr, rep, repp, repm, Hout, Sout)
    integer,  intent(in) :: Kmax, npairs, nterms, n, np, nprocs, procid
    real(wp), intent(in) :: NonlinParam(np,Kmax)
    integer,  intent(in) :: k_list(npairs), l_list(npairs)
    real(wp), intent(in) :: YHYMatr(n*n*nterms), YHYCoeff(nterms)
    real(wp), intent(in) :: pir3n2, charge0, attr, rep, repp, repm
    real(wp), intent(in) :: mass(n,n), charge(n)
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    real(wp), device, allocatable :: d_Nonlin(:,:), d_YHY(:), d_coeff(:), d_H(:), d_S(:)
    real(wp), device, allocatable :: d_mass(:,:), d_charge(:)
    integer,  device, allocatable :: d_k(:), d_l(:)
    real(wp) :: sqrtpi
    integer  :: blk, istat

    sqrtpi = sqrt(4.0_wp*atan(1.0_wp))
    allocate(d_Nonlin(np,Kmax), d_YHY(n*n*nterms), d_coeff(nterms))
    allocate(d_mass(n,n), d_charge(n))
    allocate(d_k(npairs), d_l(npairs), d_H(npairs), d_S(npairs))
    d_Nonlin = NonlinParam; d_YHY = YHYMatr; d_coeff = YHYCoeff
    d_mass = mass; d_charge = charge
    d_k = k_list; d_l = l_list

    blk = min(nterms, CUF_BLK)
    call me_energy_kernel<<<npairs, blk>>>(d_Nonlin, np, Kmax, d_k, d_l, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, &
        d_mass, d_charge, charge0, sqrtpi, pir3n2, attr, rep, repp, repm, d_H, d_S)
    istat = cudaGetLastError()
    if (istat /= cudaSuccess) write(*,'(1x,a,a)') &
        'FATAL cuf_compute_matelem_batch: energy kernel launch failed: ',trim(cudaGetErrorString(istat))

    Hout = d_H; Sout = d_S
    deallocate(d_Nonlin, d_YHY, d_coeff, d_mass, d_charge, d_k, d_l, d_H, d_S)
  end subroutine cuf_compute_matelem_batch

  !Energy+gradient kernel: one block per (k,l) pair; threads STRIDE over the terms.
  !Each thread calls the shared MatrixElements (grad_k=true, grad_l=flag) and
  !atomic-accumulates H,S and the gradient slabs Dk,Dl into shared memory.
  attributes(global) subroutine me_grad_kernel(NonlinParam, np, Kmax, &
      k_list, l_list, grad_l_flag, YHYMatr, YHYCoeff, nterms, n, nprocs, procid, &
      mass, charge, charge0, sqrtpi, pir3n2, attr, rep, repp, repm, &
      Hout, Sout, Dkout, Dlout)
    integer, value   :: np, Kmax, nterms, n, nprocs, procid
    real(wp)         :: NonlinParam(np,Kmax)
    integer          :: k_list(*), l_list(*), grad_l_flag(*)
    real(wp)         :: YHYMatr(*), YHYCoeff(*)
    real(wp)         :: mass(n,n), charge(n)
    real(wp), value  :: charge0, sqrtpi, pir3n2, attr, rep, repp, repm
    real(wp)         :: Hout(*), Sout(*), Dkout(*), Dlout(*)   !Dk/Dl slabs (2*np,npairs) flattened
    real(wp), shared :: sh_H, sh_S
    real(wp), shared :: sh_Dk(2*NNP), sh_Dl(2*NNP)
    integer  :: pair_idx, j, qbase, k0, l0, istat, q, npt2, nthr
    real(wp) :: Hkl, Skl, coeff, Dk(2*NNP), Dl(2*NNP)
    logical  :: gl

    pair_idx = blockIdx%x
    nthr     = blockDim%x
    npt2     = 2*np
    if (threadIdx%x==1) then
      sh_H=0.0_wp; sh_S=0.0_wp
    endif
    do q=threadIdx%x,npt2,nthr    !stride-init the shared gradient slabs
      sh_Dk(q)=0.0_wp; sh_Dl(q)=0.0_wp
    enddo
    call syncthreads()

    qbase = (pair_idx-1)*nterms - 1
    k0    = k_list(pair_idx)                !pair-only: hoist out of the term loop
    l0    = l_list(pair_idx)
    gl    = (grad_l_flag(pair_idx)==1)
    do j = threadIdx%x, nterms, nthr        !STRIDED: block size independent of nterms
      if (mod(qbase+j, nprocs) == procid) then
        call MatrixElements(n, np, NonlinParam(1,k0), NonlinParam(1,l0), &
                            YHYMatr((j-1)*n*n+1), mass, charge, charge0, &
                            sqrtpi, pir3n2, attr, rep, repp, repm, &
                            Hkl, Skl, Dk, Dl, .true., gl)
        coeff = YHYCoeff(j)
        istat = atomicadd(sh_H, coeff*Hkl)
        istat = atomicadd(sh_S, coeff*Skl)
        do q=1,npt2
          istat = atomicadd(sh_Dk(q), coeff*Dk(q))
          if (gl) istat = atomicadd(sh_Dl(q), coeff*Dl(q))
        enddo
      endif
    enddo
    call syncthreads()

    if (threadIdx%x==1) then
      Hout(pair_idx)=sh_H; Sout(pair_idx)=sh_S
    endif
    do q=threadIdx%x,npt2,nthr    !stride-write the slabs back to global
      Dkout((pair_idx-1)*npt2+q)=sh_Dk(q)
      Dlout((pair_idx-1)*npt2+q)=sh_Dl(q)
    enddo
  end subroutine me_grad_kernel

  !Host launcher for the gradient build.
  subroutine cuf_compute_matelem_deriv_batch(NonlinParam, Kmax, k_list, l_list, grad_l_flag, &
      npairs, YHYMatr, YHYCoeff, nterms, n, np, nprocs, procid, &
      pir3n2, mass, charge, charge0, attr, rep, repp, repm, Hout, Sout, Dkout, Dlout)
    integer,  intent(in) :: Kmax, npairs, nterms, n, np, nprocs, procid
    real(wp), intent(in) :: NonlinParam(np,Kmax)
    integer,  intent(in) :: k_list(npairs), l_list(npairs), grad_l_flag(npairs)
    real(wp), intent(in) :: YHYMatr(n*n*nterms), YHYCoeff(nterms)
    real(wp), intent(in) :: pir3n2, charge0, attr, rep, repp, repm
    real(wp), intent(in) :: mass(n,n), charge(n)
    real(wp), intent(out):: Hout(npairs), Sout(npairs)
    real(wp), intent(out):: Dkout(2*np*npairs), Dlout(2*np*npairs)
    real(wp), device, allocatable :: d_Nonlin(:,:), d_YHY(:), d_coeff(:), d_H(:), d_S(:)
    real(wp), device, allocatable :: d_mass(:,:), d_charge(:), d_Dk(:), d_Dl(:)
    integer,  device, allocatable :: d_k(:), d_l(:), d_gl(:)
    real(wp) :: sqrtpi
    integer  :: npt2, blk, istat

    npt2   = 2*np
    sqrtpi = sqrt(4.0_wp*atan(1.0_wp))
    allocate(d_Nonlin(np,Kmax), d_YHY(n*n*nterms), d_coeff(nterms))
    allocate(d_mass(n,n), d_charge(n))
    allocate(d_k(npairs), d_l(npairs), d_gl(npairs), d_H(npairs), d_S(npairs))
    allocate(d_Dk(npt2*npairs), d_Dl(npt2*npairs))
    d_Nonlin = NonlinParam; d_YHY = YHYMatr; d_coeff = YHYCoeff
    d_mass = mass; d_charge = charge
    d_k = k_list; d_l = l_list; d_gl = grad_l_flag

    blk = min(nterms, CUF_BLK)
    call me_grad_kernel<<<npairs, blk>>>(d_Nonlin, np, Kmax, d_k, d_l, d_gl, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, &
        d_mass, d_charge, charge0, sqrtpi, pir3n2, attr, rep, repp, repm, &
        d_H, d_S, d_Dk, d_Dl)
    istat = cudaGetLastError()
    if (istat /= cudaSuccess) write(*,'(1x,a,a)') &
        'FATAL cuf_compute_matelem_deriv_batch: grad kernel launch failed: ',trim(cudaGetErrorString(istat))

    Hout = d_H; Sout = d_S; Dkout = d_Dk; Dlout = d_Dl
    deallocate(d_Nonlin, d_YHY, d_coeff, d_mass, d_charge, d_k, d_l, d_gl, d_H, d_S, d_Dk, d_Dl)
  end subroutine cuf_compute_matelem_deriv_batch

  !Device lifecycle. gpu_init picks device = rank % nGPUs; called from gpu_backend_init.
  subroutine gpu_init(local_rank)
    integer :: local_rank
    integer :: ndev, dev, istat
    type(cudaDeviceProp) :: prop
    istat = cudaGetDeviceCount(ndev)
    if (ndev <= 0) then
      write(*,*) 'gpu_init: no CUDA devices'; return
    endif
    dev   = mod(local_rank, ndev)
    istat = cudaSetDevice(dev)
    istat = cudaGetDeviceProperties(prop, dev)
    write(*,'(1x,a,i0,a,i0,a,a)') 'gpu_init: rank ',local_rank,' -> device ',dev,' ',trim(prop%name)
  end subroutine gpu_init

  subroutine gpu_finalize()
    integer :: istat
    istat = cudaDeviceReset()
  end subroutine gpu_finalize

  !cuSOLVER generalized symmetric eigensolver. Solves H x = lambda S x (itype 1),
  !returns the iwhich-th eigenvalue and, if jobz_vec/=0, its eigenvector. Handle +
  !device buffers cached across calls.
  subroutine gpu_dsygvx(jobz_vec, n, H, ldH, S, ldS, iwhich, eval, Z, info_out)
    integer  :: jobz_vec, n, ldH, ldS, iwhich, info_out
    real(wp) :: H(ldH,n), S(ldS,n), eval, Z(*)
    real(wp), device, allocatable, save :: dA(:,:), dB(:,:), dW(:), dwork(:)
    integer,  device, allocatable, save :: dinfo(:)
    type(cusolverDnHandle), save :: hdl
    logical, save :: created = .false.
    integer, save :: cap_n = 0, cap_work = 0
    integer :: istat, lwork, meig, jobz, uplo
    real(wp) :: hA(n,n), hB(n,n), hZ(n), hW1(1)   !host staging (device<->host via cudaMemcpy)
    integer  :: hinfo(1)

    if (.not. created) then
      istat = cusolverDnCreate(hdl)
      allocate(dinfo(1))
      created = .true.
    endif
    if (n > cap_n) then
      if (allocated(dA)) deallocate(dA, dB, dW)
      allocate(dA(n,n), dB(n,n), dW(n)); cap_n = n
    endif
    hA(1:n,1:n) = H(1:n,1:n)      !host de-stride (ld -> n), then explicit H2D
    hB(1:n,1:n) = S(1:n,1:n)
    istat = cudaMemcpy(dA, hA, n*n)   !array-section device assignment ICEs nvfortran -> use cudaMemcpy
    istat = cudaMemcpy(dB, hB, n*n)

    jobz = CUSOLVER_EIG_MODE_NOVECTOR
    if (jobz_vec /= 0) jobz = CUSOLVER_EIG_MODE_VECTOR
    uplo = CUBLAS_FILL_MODE_UPPER

    istat = cusolverDnDsygvdx_bufferSize(hdl, CUSOLVER_EIG_TYPE_1, jobz, &
        CUSOLVER_EIG_RANGE_I, uplo, n, dA, n, dB, n, 0.0_wp, 0.0_wp, &
        iwhich, iwhich, meig, dW, lwork)
    if (lwork > cap_work) then
      if (allocated(dwork)) deallocate(dwork)
      allocate(dwork(lwork)); cap_work = lwork
    endif
    istat = cusolverDnDsygvdx(hdl, CUSOLVER_EIG_TYPE_1, jobz, &
        CUSOLVER_EIG_RANGE_I, uplo, n, dA, n, dB, n, 0.0_wp, 0.0_wp, &
        iwhich, iwhich, meig, dW, dwork, lwork, dinfo(1))
    istat = cudaDeviceSynchronize()

    istat = cudaMemcpy(hinfo, dinfo, 1); info_out = hinfo(1)   !D2H
    istat = cudaMemcpy(hW1, dW, 1);      eval     = hW1(1)     !selected eigenvalue (first of meig)
    if (jobz_vec /= 0) then
      istat = cudaMemcpy(hZ, dA, n)      !D2H first column of A (the eigenvector)
      Z(1:n) = hZ(1:n)
    endif
  end subroutine gpu_dsygvx

end module gpu_backend
