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
  use matelem,  only: MatrixElements
  use globvars, only: wp, Glob_AllowedNumOfPseudoParticles
  implicit none
  private
  public :: cuf_compute_matelem_batch

  integer, parameter :: NN = Glob_AllowedNumOfPseudoParticles

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
    real(wp) :: dDk(2*NN), dDl(2*NN)   !unused energy-path gradient outputs

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

end module matelem_cuf
