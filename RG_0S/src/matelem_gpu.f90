module matelem_gpu
!ISO_C_BINDING interfaces to the CUDA backend (matelem_cuda.cu).
!Active only when compiled with -DUSE_CUDA.
use iso_c_binding
implicit none

interface

  subroutine gpu_init(local_rank) bind(C, name="gpu_init_")
    use iso_c_binding
    integer(c_int), intent(inout) :: local_rank
  end subroutine gpu_init

  subroutine gpu_finalize() bind(C, name="gpu_finalize_")
  end subroutine gpu_finalize

  !Batched energy-only kernel (ComputeMatElem GPU path).
  !NonlinParam(np,Nmax): all basis function parameters, Fortran col-major.
  !k_list(npairs), l_list(npairs): 1-indexed pair indices.
  !Hout(npairs), Sout(npairs): this rank's partial H/S sums (before ALLREDUCE).
  subroutine gpu_compute_matelem_batch( &
      NonlinParam, Nmax, k_list, l_list, npairs, &
      YHYMatr, YHYCoeff, NumYHYTerms, n, np, &
      NumOfProcs, ProcID, &
      piraised3n2, massmat, charge, charge0, &
      attr_scale, rep_scale, rep_plus, rep_minus, &
      Hout, Sout) &
      bind(C, name="gpu_compute_matelem_batch_")
    use iso_c_binding
    real(c_double), intent(in)  :: NonlinParam(*)
    integer(c_int), intent(in)  :: Nmax
    integer(c_int), intent(in)  :: k_list(*), l_list(*)
    integer(c_int), intent(in)  :: npairs
    real(c_double), intent(in)  :: YHYMatr(*), YHYCoeff(*)
    integer(c_int), intent(in)  :: NumYHYTerms, n, np
    integer(c_int), intent(in)  :: NumOfProcs, ProcID
    real(c_double), intent(in)  :: piraised3n2
    real(c_double), intent(in)  :: massmat(*), charge(*)
    real(c_double), intent(in)  :: charge0
    real(c_double), intent(in)  :: attr_scale, rep_scale, rep_plus, rep_minus
    real(c_double), intent(out) :: Hout(*), Sout(*)
  end subroutine gpu_compute_matelem_batch

  !Batched energy+gradient kernel (ComputeMatElemAndDeriv GPU path).
  !grad_l_list(npairs): 1 if grad_l for that pair, 0 otherwise.
  !Dkout(2*np,npairs), Dlout(2*np,npairs): gradient partial sums.
  subroutine gpu_compute_matelem_and_deriv_batch( &
      NonlinParam, Nmax, k_list, l_list, grad_l_list, npairs, &
      YHYMatr, YHYCoeff, NumYHYTerms, n, np, &
      NumOfProcs, ProcID, &
      piraised3n2, massmat, charge, charge0, &
      attr_scale, rep_scale, rep_plus, rep_minus, &
      Hout, Sout, Dkout, Dlout) &
      bind(C, name="gpu_compute_matelem_and_deriv_batch_")
    use iso_c_binding
    real(c_double), intent(in)  :: NonlinParam(*)
    integer(c_int), intent(in)  :: Nmax
    integer(c_int), intent(in)  :: k_list(*), l_list(*), grad_l_list(*)
    integer(c_int), intent(in)  :: npairs
    real(c_double), intent(in)  :: YHYMatr(*), YHYCoeff(*)
    integer(c_int), intent(in)  :: NumYHYTerms, n, np
    integer(c_int), intent(in)  :: NumOfProcs, ProcID
    real(c_double), intent(in)  :: piraised3n2
    real(c_double), intent(in)  :: massmat(*), charge(*)
    real(c_double), intent(in)  :: charge0
    real(c_double), intent(in)  :: attr_scale, rep_scale, rep_plus, rep_minus
    real(c_double), intent(out) :: Hout(*), Sout(*)
    real(c_double), intent(out) :: Dkout(*), Dlout(*)
  end subroutine gpu_compute_matelem_and_deriv_batch

end interface

end module matelem_gpu
