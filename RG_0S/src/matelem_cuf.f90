module matelem_cuf
!SPIKE: native CUDA Fortran energy-path matrix-element backend.
!
!An alternative to the CUDA C++ energy kernel in matelem_cuda.cu, written in
!nvfortran CUDA Fortran: same block-per-pair / thread-per-permutation-term layout
!and atomicAdd reduction, but the per-element math is a direct device port of the
!CPU Fortran (matelem.f90 MatrixElements, energy path, lines 68-236) rather than a
!hand C++ translation. Physics constants live in device `constant` memory, copied
!once from the Glob_* values. Compiled only with -DUSE_CUF (nvfortran -cuda).
!
!Scope: ENERGY path only (Hkl/Skl, no gradients). The gradient build and the
!eigensolver still use the C++/cuSOLVER path, so this drops in behind gpu_build_HS
!for A/B against matelem_cuda.cu on energy-only workloads (e.g. EXPC_VALS).
  use cudafor
  use globvars, only: wp, Glob_AllowedNumOfPseudoParticles
  implicit none
  private
  public :: cuf_compute_matelem_batch

  integer, parameter :: NN = Glob_AllowedNumOfPseudoParticles

  !--- Physics constants in device constant memory (copied once at first use) ---
  real(wp), constant :: cu_pir3n2, cu_sqrtpi, cu_charge0
  real(wp), constant :: cu_attr, cu_rep, cu_repp, cu_repm
  real(wp), constant :: cu_mass(NN,NN)
  real(wp), constant :: cu_charge(NN)
  logical, save :: const_loaded = .false.

contains

  !===========================================================================
  ! Device: scaled charge product (port of ScaledChargeProd, matelem.f90:1615)
  !===========================================================================
  attributes(device) real(wp) function scaled_charge_prod(q1,q2)
    real(wp), value :: q1,q2
    real(wp) :: x
    x = q1*q2
    if (x < 0.0_wp) then
      scaled_charge_prod = x*cu_attr
    else if ((q1 > 0.0_wp) .and. (q2 > 0.0_wp)) then
      scaled_charge_prod = x*cu_rep*cu_repp
    else
      scaled_charge_prod = x*cu_rep*cu_repm
    endif
  end function scaled_charge_prod

  !===========================================================================
  ! Device: energy-only matrix element for one (k,l) pair and one permutation
  ! term P (n x n). Direct port of matelem.f90 MatrixElements energy path.
  !===========================================================================
  attributes(device) subroutine me_energy_dev(n, vechLk, vechLl, P, Hkl, Skl)
    integer, value        :: n
    real(wp), intent(in)  :: vechLk(*), vechLl(*)
    real(wp), intent(in)  :: P(n,n)
    real(wp), intent(out) :: Hkl, Skl
    real(wp) :: Lk(NN,NN), Ll(NN,NN), Ak(NN,NN), tAl(NN,NN), tAkl(NN,NN)
    real(wp) :: W1(NN,NN), W2(NN,NN), inv_tAkl(NN,NN), inv_tAkltAlM(NN,NN)
    real(wp) :: det_tAkl, Tkl, Vkl, temp1, temp3, temp4, temp5
    integer  :: i,j,k,indx

    !Build Lk, Ll (lower-triangular) from the packed vech parameters
    indx = 0
    do i=1,n
      do j=i,n
        indx = indx+1
        Lk(i,j)=0.0_wp; Lk(j,i)=vechLk(indx)
        Ll(i,j)=0.0_wp; Ll(j,i)=vechLl(indx)
      enddo
    enddo

    !Ak = Lk Lk^T,  tAl = Ll Ll^T
    do i=1,n
      do j=i,n
        temp1=0.0_wp
        do k=1,i
          temp1=temp1+Lk(i,k)*Lk(j,k)
        enddo
        Ak(i,j)=temp1; Ak(j,i)=temp1
        temp1=0.0_wp
        do k=1,i
          temp1=temp1+Ll(i,k)*Ll(j,k)
        enddo
        tAl(i,j)=temp1; tAl(j,i)=temp1
      enddo
    enddo

    !tAl = P^T Al P ;  tAkl = Ak + tAl
    do i=1,n
      do j=1,n
        temp1=0.0_wp
        do k=1,n
          temp1=temp1+P(k,j)*tAl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    do i=1,n
      do j=i,n
        temp1=0.0_wp
        do k=1,n
          temp1=temp1+W1(i,k)*P(k,j)
        enddo
        tAl(i,j)=temp1; tAl(j,i)=temp1
        tAkl(i,j)=Ak(i,j)+temp1; tAkl(j,i)=tAkl(i,j)
      enddo
    enddo

    !Cholesky of tAkl into the lower triangle of W1
    det_tAkl=1.0_wp
    do i=1,n
      do j=i,n
        temp1=tAkl(i,j)
        do k=i-1,1,-1
          temp1=temp1-W1(i,k)*W1(j,k)
        enddo
        if (i==j) then
          W1(i,i)=sqrt(temp1); det_tAkl=det_tAkl*temp1
        else
          W1(j,i)=temp1/W1(i,i); W1(i,j)=0.0_wp
        endif
      enddo
    enddo

    !Invert the Cholesky factor, then inv_tAkl = W1^T W1
    do i=1,n
      W1(i,i)=1.0_wp/W1(i,i)
      do j=i+1,n
        temp1=0.0_wp
        do k=i,j-1
          temp1=temp1-W1(j,k)*W1(k,i)
        enddo
        W1(j,i)=temp1/W1(j,j)
      enddo
    enddo
    do i=1,n
      do j=i,n
        temp1=0.0_wp
        do k=j,n
          temp1=temp1+W1(k,i)*W1(k,j)
        enddo
        inv_tAkl(i,j)=temp1; inv_tAkl(j,i)=temp1
      enddo
    enddo

    !Overlap
    Skl = cu_pir3n2/(det_tAkl*sqrt(det_tAkl))

    !W2 = inv_tAkl tAl ;  inv_tAkltAlM = W2 M
    do i=1,n
      do j=1,n
        temp1=0.0_wp
        do k=1,n
          temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
        enddo
        W2(j,i)=temp1
      enddo
    enddo
    do i=1,n
      do j=1,n
        temp1=0.0_wp
        do k=1,n
          temp1=temp1+W2(j,k)*cu_mass(k,i)
        enddo
        inv_tAkltAlM(j,i)=temp1
      enddo
    enddo

    !Kinetic energy Tkl = 6 Skl tr[inv_tAkltAlM Ak]
    Tkl=0.0_wp
    do i=1,n
      temp1=0.0_wp
      do k=1,n
        temp1=temp1+inv_tAkltAlM(i,k)*Ak(k,i)
      enddo
      Tkl=Tkl+temp1
    enddo
    Tkl=6.0_wp*Skl*Tkl

    !Potential energy Vkl
    temp1=(2.0_wp/cu_sqrtpi)*Skl
    Vkl=0.0_wp
    do i=1,n
      temp3=inv_tAkl(i,i)
      temp4=sqrt(temp3)
      temp5=temp1/temp4
      Vkl=Vkl+scaled_charge_prod(cu_charge(i),cu_charge0)*temp5
    enddo
    do i=1,n
      do j=i+1,n
        temp3=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
        temp4=sqrt(temp3)
        temp5=temp1/temp4
        Vkl=Vkl+scaled_charge_prod(cu_charge(i),cu_charge(j))*temp5
      enddo
    enddo

    Hkl=Tkl+Vkl
  end subroutine me_energy_dev

  !===========================================================================
  ! Global kernel: one block per (k,l) pair, one thread per permutation term.
  ! Matches the C++ matelem_batch_energy_kernel (incl. the MPI term split).
  !===========================================================================
  attributes(global) subroutine me_energy_kernel(NonlinParam, np, Kmax, &
      k_list, l_list, YHYMatr, YHYCoeff, nterms, n, nprocs, procid, Hout, Sout)
    integer, value       :: np, Kmax, nterms, n, nprocs, procid
    real(wp)             :: NonlinParam(np,Kmax)
    integer              :: k_list(*), l_list(*)
    real(wp)             :: YHYMatr(*), YHYCoeff(*)
    real(wp)             :: Hout(*), Sout(*)
    real(wp), shared     :: sh_H, sh_S
    integer  :: pair_idx, j, qbase, k0, l0, istat
    real(wp) :: Hkl, Skl, coeff

    pair_idx = blockIdx%x
    j        = threadIdx%x            !1-based term index
    if (j==1) then
      sh_H=0.0_wp; sh_S=0.0_wp
    endif
    call syncthreads()

    if (j <= nterms) then
      qbase = (pair_idx-1)*nterms - 1     !matches Fortran q=(i-1)*nterms-1
      if (mod(qbase+j, nprocs) == procid) then
        k0 = k_list(pair_idx)
        l0 = l_list(pair_idx)
        !P for term j starts at YHYMatr((j-1)*n*n + 1); pass as n x n
        call me_energy_dev(n, NonlinParam(1,k0), NonlinParam(1,l0), &
                           YHYMatr((j-1)*n*n+1), Hkl, Skl)
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
  ! Host launcher. Same argument list as the C gpu_compute_matelem_batch, so it
  ! drops into gpu_build_HS. Copies constants once, stages inputs to the device,
  ! launches, and copies the per-pair partial sums back (ALLREDUCE happens on the
  ! host side, as in the C++ path).
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
    !device mirrors
    real(wp), device, allocatable :: d_Nonlin(:,:), d_YHY(:), d_coeff(:), d_H(:), d_S(:)
    integer,  device, allocatable :: d_k(:), d_l(:)
    integer :: threads

    if (.not. const_loaded) then
      cu_pir3n2 = pir3n2; cu_sqrtpi = sqrt(4.0_wp*atan(1.0_wp))
      cu_charge0 = charge0
      cu_attr = attr; cu_rep = rep; cu_repp = repp; cu_repm = repm
      cu_mass(1:n,1:n) = mass(1:n,1:n)
      cu_charge(1:n)   = charge(1:n)
      const_loaded = .true.
    endif

    allocate(d_Nonlin(np,Kmax), d_YHY(n*n*nterms), d_coeff(nterms))
    allocate(d_k(npairs), d_l(npairs), d_H(npairs), d_S(npairs))
    d_Nonlin = NonlinParam
    d_YHY    = YHYMatr
    d_coeff  = YHYCoeff
    d_k      = k_list
    d_l      = l_list

    threads = nterms
    call me_energy_kernel<<<npairs, threads>>>(d_Nonlin, np, Kmax, d_k, d_l, &
        d_YHY, d_coeff, nterms, n, nprocs, procid, d_H, d_S)

    Hout = d_H
    Sout = d_S
    deallocate(d_Nonlin, d_YHY, d_coeff, d_k, d_l, d_H, d_S)
  end subroutine cuf_compute_matelem_batch

end module matelem_cuf
