module matelem_acc
!OpenACC energy-path matrix element, callable from a device parallel loop.
!This is the single-source-of-truth alternative to the hand-written CUDA
!kernel: the physics stays in ordinary Fortran, marked !$acc routine seq so
!nvfortran can compile it for the GPU. Used only when built with -DUSE_OPENACC.
!
!Ports the energy path of MatrixElements (matelem.f90 lines ~68-242). All
!inputs are passed as arguments (no module globals, no system_clock) so the
!routine is self-contained on the device. Gradient path not included in this
!prototype.
use globvars
implicit none

!Fixed maximum pseudoparticle count (matches Glob_MaxAllowedNumOfPseudoParticles)
integer, parameter :: ACC_NN = 7

contains

subroutine MatrixElements_energy_acc(vechLk, vechLl, P, n, np, &
    piraised3n2, MassMatrix, charge, charge0, &
    attr_scale, rep_scale, rep_plus, rep_minus, Hkl, Skl)
!$acc routine seq
real(dprec), intent(in)  :: vechLk(*), vechLl(*)
real(dprec), intent(in)  :: P(n,n)
integer,     intent(in)  :: n, np
real(dprec), intent(in)  :: piraised3n2
real(dprec), intent(in)  :: MassMatrix(n,n)
real(dprec), intent(in)  :: charge(n), charge0
real(dprec), intent(in)  :: attr_scale, rep_scale, rep_plus, rep_minus
real(dprec), intent(out) :: Hkl, Skl

!Local fixed-size scratch (so the routine needs no heap on the device)
real(dprec) :: Lk(ACC_NN,ACC_NN), Ll(ACC_NN,ACC_NN)
real(dprec) :: Ak(ACC_NN,ACC_NN), tAl(ACC_NN,ACC_NN), tAkl(ACC_NN,ACC_NN)
real(dprec) :: W1(ACC_NN,ACC_NN), W2(ACC_NN,ACC_NN)
real(dprec) :: inv_tAkl(ACC_NN,ACC_NN), inv_tAkltAlM(ACC_NN,ACC_NN)
real(dprec) :: temp1, temp3, temp4, temp5
real(dprec) :: det_tAkl, Tkl, Vkl, cprod
integer     :: i, j, k, indx

!--- Build Lk, Ll from vech (lower-triangular packing) ---
indx = 0
do i = 1, n
  do j = i, n
    indx = indx + 1
    Lk(i,j) = ZERO; Lk(j,i) = vechLk(indx)
    Ll(i,j) = ZERO; Ll(j,i) = vechLl(indx)
  enddo
enddo

!--- Ak = Lk*Lk^T, tAl = Ll*Ll^T (before permutation) ---
do i = 1, n
  do j = i, n
    temp1 = ZERO
    do k = 1, i
      temp1 = temp1 + Lk(i,k)*Lk(j,k)
    enddo
    Ak(i,j) = temp1; Ak(j,i) = temp1
    temp1 = ZERO
    do k = 1, i
      temp1 = temp1 + Ll(i,k)*Ll(j,k)
    enddo
    tAl(i,j) = temp1; tAl(j,i) = temp1
  enddo
enddo

!--- tAl = P^T*Al*P, tAkl = Ak + tAl ---
do i = 1, n
  do j = 1, n
    temp1 = ZERO
    do k = 1, n
      temp1 = temp1 + P(k,j)*tAl(k,i)
    enddo
    W1(j,i) = temp1
  enddo
enddo
do i = 1, n
  do j = i, n
    temp1 = ZERO
    do k = 1, n
      temp1 = temp1 + W1(i,k)*P(k,j)
    enddo
    tAl(i,j) = temp1; tAl(j,i) = temp1
    tAkl(i,j) = Ak(i,j) + temp1; tAkl(j,i) = tAkl(i,j)
  enddo
enddo

!--- Cholesky of tAkl, factor stored in lower triangle of W1 ---
det_tAkl = ONE
do i = 1, n
  do j = i, n
    temp1 = tAkl(i,j)
    do k = i-1, 1, -1
      temp1 = temp1 - W1(i,k)*W1(j,k)
    enddo
    if (i == j) then
      W1(i,i) = sqrt(temp1)
      det_tAkl = det_tAkl*temp1
    else
      W1(j,i) = temp1/W1(i,i)
      W1(i,j) = ZERO
    endif
  enddo
enddo

!--- Invert Cholesky factor in place ---
do i = 1, n
  W1(i,i) = ONE/W1(i,i)
  do j = i+1, n
    temp1 = ZERO
    do k = i, j-1
      temp1 = temp1 - W1(j,k)*W1(k,i)
    enddo
    W1(j,i) = temp1/W1(j,j)
  enddo
enddo

!--- inv_tAkl = W1^T * W1 ---
do i = 1, n
  do j = i, n
    temp1 = ZERO
    do k = j, n
      temp1 = temp1 + W1(k,i)*W1(k,j)
    enddo
    inv_tAkl(i,j) = temp1; inv_tAkl(j,i) = temp1
  enddo
enddo

!--- Overlap ---
Skl = piraised3n2/(det_tAkl*sqrt(det_tAkl))

!--- W2 = inv_tAkl * tAl ---
do i = 1, n
  do j = 1, n
    temp1 = ZERO
    do k = 1, n
      temp1 = temp1 + inv_tAkl(j,k)*tAl(k,i)
    enddo
    W2(j,i) = temp1
  enddo
enddo

!--- inv_tAkltAlM = W2 * M ---
do i = 1, n
  do j = 1, n
    temp1 = ZERO
    do k = 1, n
      temp1 = temp1 + W2(j,k)*MassMatrix(k,i)
    enddo
    inv_tAkltAlM(j,i) = temp1
  enddo
enddo

!--- Kinetic energy ---
Tkl = ZERO
do i = 1, n
  temp1 = ZERO
  do k = 1, n
    temp1 = temp1 + inv_tAkltAlM(i,k)*Ak(k,i)
  enddo
  Tkl = Tkl + temp1
enddo
Tkl = SIX*Skl*Tkl

!--- Potential energy ---
temp1 = (TWO/SQRTPI)*Skl
Vkl = ZERO
do i = 1, n
  temp3 = inv_tAkl(i,i)
  temp4 = sqrt(temp3)
  temp5 = temp1/temp4
  ! ScaledChargeProd(charge(i), charge0)
  cprod = charge(i)*charge0
  if (cprod < ZERO) then
    cprod = cprod*attr_scale
  else if (charge(i) > ZERO .and. charge0 > ZERO) then
    cprod = cprod*rep_scale*rep_plus
  else
    cprod = cprod*rep_scale*rep_minus
  endif
  Vkl = Vkl + cprod*temp5
enddo
do i = 1, n
  do j = i+1, n
    temp3 = inv_tAkl(i,i) + inv_tAkl(j,j) - inv_tAkl(j,i) - inv_tAkl(j,i)
    temp4 = sqrt(temp3)
    temp5 = temp1/temp4
    cprod = charge(i)*charge(j)
    if (cprod < ZERO) then
      cprod = cprod*attr_scale
    else if (charge(i) > ZERO .and. charge(j) > ZERO) then
      cprod = cprod*rep_scale*rep_plus
    else
      cprod = cprod*rep_scale*rep_minus
    endif
    Vkl = Vkl + cprod*temp5
  enddo
enddo

Hkl = Tkl + Vkl

end subroutine MatrixElements_energy_acc

end module matelem_acc
