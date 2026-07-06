module matelem_cuf
!SPIKE (single-source CUDA Fortran): GPU matrix-element kernels that call the
!SHARED MatrixElements device routine in matelem.f90 -- NO duplicated physics.
!The per-element math lives once, in matelem.f90, compiled as host,device.
!Compiled only with -DUSE_CUF (nvfortran -cuda).
!
!Constants are passed as kernel arguments (mass/charge as contiguous device
!arrays, scalars by value) rather than via constant memory, which avoids any
!leading-dimension mismatch when handing them to MatrixElements' mass(n,n).
  use cudafor
  use cublas
  use cusolverDn
  use matelem,  only: MatrixElements
  use globvars, only: wp, Glob_AllowedNumOfPseudoParticles
  implicit none
  private
  public :: cuf_compute_matelem_batch, cuf_compute_matelem_deriv_batch
  public :: gpu_init, gpu_finalize, gpu_dsygvx

  integer, parameter :: NN  = Glob_AllowedNumOfPseudoParticles
  integer, parameter :: NNP = NN*(NN+1)/2       !max packed vech length; Dk/Dl are 2*np long (np=n(n+1)/2)

contains

  !===========================================================================
  ! Energy kernel: one block per (k,l) pair, one thread per permutation term.
  ! Each thread calls the shared MatrixElements (grad=false) for its term.
  !===========================================================================
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
    j        = threadIdx%x
    if (j==1) then
      sh_H=0.0_wp; sh_S=0.0_wp
    endif
    call syncthreads()

    if (j <= nterms) then
      qbase = (pair_idx-1)*nterms - 1        !matches Fortran q=(i-1)*nterms-1
      if (mod(qbase+j, nprocs) == procid) then
        k0 = k_list(pair_idx)
        l0 = l_list(pair_idx)
        call MatrixElements(n, np, NonlinParam(1,k0), NonlinParam(1,l0), &
                            YHYMatr((j-1)*n*n+1), mass, charge, charge0, &
                            sqrtpi, pir3n2, attr, rep, repp, repm, &
                            Hkl, Skl, dDk, dDl, .false., .false.)
        coeff = YHYCoeff(j)
        istat = atomicadd(sh_H, coeff*Hkl)
        istat = atomicadd(sh_S, coeff*Skl)
      endif
    endif
    call syncthreads()

    if (j==1) then
      Hout(pair_idx)=sh_H; Sout(pair_idx)=sh_S
    endif
  end subroutine me_energy_kernel

  !===========================================================================
  ! Host launcher. Same signature as the C gpu_compute_matelem_batch so it
  ! drops into gpu_build_HS. Stages inputs to the device, launches, copies back.
  !===========================================================================
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

    sqrtpi = sqrt(4.0_wp*atan(1.0_wp))
    allocate(d_Nonlin(np,Kmax), d_YHY(n*n*nterms), d_coeff(nterms))
    allocate(d_mass(n,n), d_charge(n))
    allocate(d_k(npairs), d_l(npairs), d_H(npairs), d_S(npairs))
    d_Nonlin = NonlinParam; d_YHY = YHYMatr; d_coeff = YHYCoeff
    d_mass = mass; d_charge = charge
    d_k = k_list; d_l = l_list

    call me_energy_kernel<<<npairs, nterms>>>(d_Nonlin, np, Kmax, d_k, d_l, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, &
        d_mass, d_charge, charge0, sqrtpi, pir3n2, attr, rep, repp, repm, d_H, d_S)

    Hout = d_H; Sout = d_S
    deallocate(d_Nonlin, d_YHY, d_coeff, d_mass, d_charge, d_k, d_l, d_H, d_S)
  end subroutine cuf_compute_matelem_batch

  !===========================================================================
  ! Energy+gradient kernel: one block per (k,l) pair, one thread per term.
  ! Each thread calls the shared MatrixElements (grad_k=true, grad_l=flag) and
  ! atomic-accumulates H,S and the gradient slabs Dk,Dl into shared memory.
  !===========================================================================
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
    j        = threadIdx%x
    nthr     = blockDim%x
    npt2     = 2*np
    if (j==1) then
      sh_H=0.0_wp; sh_S=0.0_wp
    endif
    do q=j,npt2,nthr          !stride-init the shared gradient slabs
      sh_Dk(q)=0.0_wp; sh_Dl(q)=0.0_wp
    enddo
    call syncthreads()

    if (j <= nterms) then
      qbase = (pair_idx-1)*nterms - 1
      if (mod(qbase+j, nprocs) == procid) then
        k0 = k_list(pair_idx)
        l0 = l_list(pair_idx)
        gl = (grad_l_flag(pair_idx)==1)
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
    endif
    call syncthreads()

    if (j==1) then
      Hout(pair_idx)=sh_H; Sout(pair_idx)=sh_S
    endif
    do q=j,npt2,nthr          !stride-write the slabs back to global
      Dkout((pair_idx-1)*npt2+q)=sh_Dk(q)
      Dlout((pair_idx-1)*npt2+q)=sh_Dl(q)
    enddo
  end subroutine me_grad_kernel

  !===========================================================================
  ! Host launcher for the gradient path. Same signature as the C
  ! gpu_compute_matelem_and_deriv_batch so it drops into gpu_build_HS_deriv.
  !===========================================================================
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
    integer  :: npt2

    npt2   = 2*np
    sqrtpi = sqrt(4.0_wp*atan(1.0_wp))
    allocate(d_Nonlin(np,Kmax), d_YHY(n*n*nterms), d_coeff(nterms))
    allocate(d_mass(n,n), d_charge(n))
    allocate(d_k(npairs), d_l(npairs), d_gl(npairs), d_H(npairs), d_S(npairs))
    allocate(d_Dk(npt2*npairs), d_Dl(npt2*npairs))
    d_Nonlin = NonlinParam; d_YHY = YHYMatr; d_coeff = YHYCoeff
    d_mass = mass; d_charge = charge
    d_k = k_list; d_l = l_list; d_gl = grad_l_flag

    call me_grad_kernel<<<npairs, nterms>>>(d_Nonlin, np, Kmax, d_k, d_l, d_gl, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, &
        d_mass, d_charge, charge0, sqrtpi, pir3n2, attr, rep, repp, repm, &
        d_H, d_S, d_Dk, d_Dl)

    Hout = d_H; Sout = d_S; Dkout = d_Dk; Dlout = d_Dl
    deallocate(d_Nonlin, d_YHY, d_coeff, d_mass, d_charge, d_k, d_l, d_gl, d_H, d_S, d_Dk, d_Dl)
  end subroutine cuf_compute_matelem_deriv_batch

  !===========================================================================
  ! Device lifecycle (replaces gpu_init_/gpu_finalize_ in matelem_cuda.cu).
  !===========================================================================
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

  !===========================================================================
  ! cuSOLVER generalized symmetric eigensolver (replaces eigen_cuda.cu).
  ! Solves H x = lambda S x (itype 1), returns the iwhich-th eigenvalue and,
  ! if jobz_vec/=0, its eigenvector. Handle + device buffers cached across calls.
  !===========================================================================
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

end module matelem_cuf
