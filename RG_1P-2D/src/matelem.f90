module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=1 and L=2 Gaussians.

use globvars
implicit none

contains



SUBROUTINE OverLapElementS1(m_k, vechLk, P, Skk)
!==============================================================================================
! This subroutine calculates the diagonal overlap for L=1 correlated Gaussians:
!
! fk = z{m_k} * exp[-r'(Lk*Lk')r]
!
!                  <P*fk|P*fk>               ||L_kk||^3       v_k' * inv_tAkk * tv_k
! <psi_k|psi_k> = -------------- = 2^(3n/2) ------------- * --------------------------
!                    <fk|fk>                 |tAkk|^{3/2}     v_k' * inv_Akk * v_k
!
! where tAkk = Ak + P'*Ak*P,  tv_k = P'*v_k
!
!----------------------------------------------------------------------------------------------
! Input:
!   m_k      :: integer that determines which z-component is in the premultiplier
!   vechLk   :: Array of length n(n+1)/2 of exponential parameters
!   P        :: The symmetry permutation matrix of size n x n
!
! Output:
!   Skk      :: <psi_k|psi_k>
!==============================================================================================

 IMPLICIT NONE

 !---------------------------------------------------------------------------
 ! ARGUMENTS
 !---------------------------------------------------------------------------
 INTEGER, INTENT(IN)      :: m_k
 REAL(dprec), INTENT(IN)  :: vechLk(Glob_np)
 REAL(dprec), INTENT(IN)  :: P(Glob_n, Glob_n)
 REAL(dprec), INTENT(OUT) :: Skk

 !---------------------------------------------------------------------------
 ! LOCAL VARIABLES
 !---------------------------------------------------------------------------
 INTEGER                  :: n
 INTEGER                  :: i, j, k, indx

 REAL(dprec)              :: Lk(Glob_n, Glob_n)
 REAL(dprec)              :: tAk(Glob_n, Glob_n), tAkk(Glob_n, Glob_n)
 REAL(dprec)              :: inv_tAkk(Glob_n, Glob_n), inv_Akk(Glob_n, Glob_n)
 REAL(dprec)              :: W1(Glob_n, Glob_n)

 REAL(dprec)              :: tvk(Glob_n)

 REAL(dprec)              :: temp1
 REAL(dprec)              :: det_Lk, det_tAkk, tau3

 !===========================================================================
 n = Glob_n

 !---------------------------------------------------------------------------
 ! 1. Build lower-triangular Lk from vechLk
 !---------------------------------------------------------------------------
 indx = 0
 DO i = 1, n
     DO j = i, n
         indx = indx + 1
         Lk(i, j) = ZERO
         Lk(j, i) = vechLk(indx)
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 2. Compute inv_Akk using Cholesky factor Lk
 !    inv_Akk = (1/2) * Lk^{-T} * Lk^{-1} = (2*Ak)^{-1} = Akk^{-1}
 !---------------------------------------------------------------------------
 DO i = 1, n
     W1(i, i) = ONE / Lk(i, i)
     DO j = i + 1, n
         temp1 = ZERO
         DO k = i, j - 1
             temp1 = temp1 - Lk(j, k) * W1(k, i)
         END DO
         W1(j, i) = temp1 / Lk(j, j)
     END DO
 END DO

 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = j, n
             temp1 = temp1 + W1(k, i) * W1(k, j)
         END DO
         inv_Akk(i, j) = ONEHALF * temp1
         inv_Akk(j, i) = ONEHALF * temp1
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 3. Compute Ak = Lk * Lk'
 !---------------------------------------------------------------------------
 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = 1, i
             temp1 = temp1 + Lk(i, k) * Lk(j, k)
         END DO
         tAk(i, j) = temp1
         tAk(j, i) = temp1
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 4. Apply symmetry permutation and form tAkk = Ak + P' * Ak * P
 !---------------------------------------------------------------------------
 ! W1 = P' * Ak
 DO i = 1, n
     DO j = 1, n
         temp1 = ZERO
         DO k = 1, n
             temp1 = temp1 + P(k, j) * tAk(k, i)
         END DO
         W1(j, i) = temp1
     END DO
 END DO

 ! tAkk = Ak + P' * Ak * P
 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = 1, n
             temp1 = temp1 + W1(i, k) * P(k, j)
         END DO
         tAkk(i, j) = tAk(i, j) + temp1
         tAkk(j, i) = tAkk(i, j)
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 5. Determinant of Lk (product of diagonal elements)
 !---------------------------------------------------------------------------
 det_Lk = ONE
 DO i = 1, n
    det_Lk = det_Lk * Lk(i, i)
 END DO

 !---------------------------------------------------------------------------
 ! 6. Cholesky factorization of tAkk -> W1, and determinant
 !---------------------------------------------------------------------------
 det_tAkk = ONE
 DO i = 1, n
     DO j = i, n
         temp1 = tAkk(i, j)
         DO k = i - 1, 1, -1
             temp1 = temp1 - W1(i, k) * W1(j, k)
         END DO
         IF (i == j) THEN
             W1(i, i) = SQRT(temp1)

              IF (temp1 <= 0.0_dprec) THEN
                 WRITE(*,*) 'ERROR: tAkk not SPD at i=',i
                 ERROR STOP
              END IF      

             det_tAkk = det_tAkk * temp1
         ELSE
             W1(j, i) = temp1 / W1(i, i)
             W1(i, j) = ZERO
         END IF
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 7. Invert tAkk using its Cholesky factor (stored in W1)
 !---------------------------------------------------------------------------
 DO i = 1, n
     W1(i, i) = ONE / W1(i, i)
     DO j = i + 1, n
         temp1 = ZERO
         DO k = i, j - 1
             temp1 = temp1 - W1(j, k) * W1(k, i)
         END DO
         W1(j, i) = temp1 / W1(j, j)
     END DO
 END DO

 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = j, n
             temp1 = temp1 + W1(k, i) * W1(k, j)
         END DO
         inv_tAkk(i, j) = temp1
         inv_tAkk(j, i) = temp1
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 8. Compute tau3 = v_k' * inv_tAkk * tv_k
 !    where v_k = e_{m_k} (unit vector) and tv_k = P(m_k, :)
 !---------------------------------------------------------------------------
 DO i = 1, n
     tvk(i) = P(m_k, i)
 END DO

 tau3 = ZERO
 DO i = 1, n
     tau3 = tau3 + inv_tAkk(m_k, i) * tvk(i)
 END DO

 !---------------------------------------------------------------------------
 ! 9. Final overlap calculation
 !    Skk = 2^{3n/2} * ||L_kk||^3 / |tAkk|^{3/2} * tau3 / inv_Akk(m_k, m_k)
 !---------------------------------------------------------------------------
 temp1 = (det_Lk * det_Lk) / det_tAkk
 temp1 = ABS(temp1) * SQRT(ABS(temp1))

 Skk = Glob_2raised3n2 * temp1 * tau3 / inv_Akk(m_k, m_k)


END SUBROUTINE OverLapElementS1




SUBROUTINE OverLapElementS2(ml_1, ml_2, vechLl, P, Sll)
!
! The subroutine was copied from "RG_2D" code.
! This subroutine computes symmetry adapted matrix element with
! two real L=2 correlated Gaussians:
!
!       —                    —
!      | (v_k1' r)*(w_k1' r)+ | 
! fk = | (v_k2' r)*(w_k2' r)+ | exp[-r'(Lk*Lk')r]
!      | (v_k3' r)*(w_k3' r)  |
!       —                    —
!
!       —                                                     —
!      |                                                       |
!     =| x(ml_1)*x(ml_2) + y(ml_1)*y(ml_2) - 2 z(ml_1)*z(ml_2) |  exp[-r'(Lk*Lk')r]
!      |                                                       |
!       —                                                     —
!
!        m_k and mm_k are integers between 1 and n 
!         (n is the number ofpseudoparticles).
!        Symmetry adaption is applied to the ket using permutation matrix P
!
!
!                  <P*fk|P*fl>                 || L_{ll} ||^3      
! <psi_k|psi_l> = -------------- =  2^(3*n/2) ---------------- 
!                    <fk|fl>                    det_tAkl^3/2     
!
!                                     (tgamma1 * tgamma2 + tgamma5 * tgamma6)
!                                * ---------------------------------------------
!                                       (gamma1 * gamma2 + gamma5 * gamma6)
! 
!                gamma1 = v_k'*inv_tAkl*v_l
!                gamma2 = w_k'*inv_tAkl*w_l
!                gamma5 = v_k'*inv_tAkl*w_l
!                gamma6 = w_k'*inv_tAkl*v_l
!----------------------------------------------------------------------------------------------
! Input:
!  ml_1, ml_2      :: integers that determine which psuedoparticles carry l=1 momentum
!  vechLl         :: Array of length (n(n+1)/2) of exponential parameters.
!  P              :: The symmetry permutation matrix of size n x n
!
! Output:
!  Sll	          :: Overlap matrix element (normalized)
!
!===============================================================================================

 IMPLICIT NONE
 !---------------------------------------------------------------------------
 ! ARGUMENTS
 !---------------------------------------------------------------------------
 INTEGER, INTENT(IN)      :: ml_1, ml_2
 REAL(dprec), INTENT(IN)  :: vechLl(Glob_np)
 REAL(dprec), INTENT(IN)  :: P(Glob_n, Glob_n)
 REAL(dprec), INTENT(OUT) :: Sll

 !---------------------------------------------------------------------------
 ! LOCAL VARIABLES
 !---------------------------------------------------------------------------
 INTEGER                  :: i, j, k, indx
 INTEGER                  :: n
 ! Vectors for Bra/Ket components
 REAL(dprec)              :: tvl(Glob_n), twl(Glob_n)
 REAL(dprec)              :: tvl_inv_tAll(Glob_n), twl_inv_tAll(Glob_n)
 ! Matrices
 REAL(dprec)              :: Ll(Glob_n, Glob_n)
 REAL(dprec)              :: tAl(Glob_n, Glob_n), tAll(Glob_n, Glob_n)
 REAL(dprec)              :: inv_tAll(Glob_n, Glob_n), inv_All(Glob_n, Glob_n)
 REAL(dprec)              :: W1(Glob_n, Glob_n), W2(Glob_n, Glob_n)
 ! Scalars
 REAL(dprec)              :: temp1, temp2
 REAL(dprec)              :: det_tAll, det_Ll
 REAL(dprec)              :: tgamma1, tgamma2, tgamma5, tgamma6
 REAL(dprec)              :: gamma1, gamma2, gamma5, gamma6
 REAL(dprec)              :: tm, m

 !===========================================================================
 ! INITIALIZATION
 !===========================================================================
 n = Glob_n

 !---------------------------------------------------------------------------
 ! 1. Construct Lower-Triangular Matrix Ll from Packed Vector
 !---------------------------------------------------------------------------
 indx = 0
 DO i = 1, n
     DO j = i, n
         indx = indx + 1
         Ll(i, j) = ZERO
         Ll(j, i) = vechLl(indx)
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 2. Compute Inverse of All (Unpermuted)
 !    We first compute Ll^-1, then construct inv_All = 0.5 * (L L^T)^-1
 !---------------------------------------------------------------------------
 
 ! Invert Lower Triangular Matrix Ll -> Store in W1
 DO i = 1, n
     W1(i, i) = ONE / Ll(i, i)
     DO j = i + 1, n
         temp1 = ZERO
         DO k = i, j - 1
             temp1 = temp1 - Ll(j, k) * W1(k, i)
         END DO
         W1(j, i) = temp1 / Ll(j, j)
     END DO
 END DO 

 ! Construct inv_All = 0.5 * W1^T * W1
 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = j, n
             temp1 = temp1 + W1(k, i) * W1(k, j)
         END DO
         inv_All(i, j) = ONEHALF * temp1
         inv_All(j, i) = ONEHALF * temp1
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 3. Compute Matrix Al = Ll * Ll^T
 !---------------------------------------------------------------------------
 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = 1, i
             temp1 = temp1 + Ll(i, k) * Ll(j, k)
         END DO 
         tAl(i, j) = temp1
         tAl(j, i) = temp1
     END DO
 END DO

 !---------------------------------------------------------------------------
 ! 4. Apply symmetry permutation and form tAll = Al + P' * Al * P
 !---------------------------------------------------------------------------
 ! W1 = P' * Al
 DO i = 1, n
     DO j = 1, n
         temp1 = ZERO
         DO k = 1, n
             temp1 = temp1 + P(k, j) * tAl(k, i)
         END DO
         W1(j, i) = temp1
     END DO
 END DO

 ! tAll = Al + P' * Al * P
 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = 1, n
             temp1 = temp1 + W1(i, k) * P(k, j)
         END DO
         tAll(i, j) = tAl(i, j) + temp1
         tAll(j, i) = tAll(i, j)
     END DO
 END DO


 !---------------------------------------------------------------------------
 ! 5. Determinants and Inversion of tAll
 !---------------------------------------------------------------------------
 
 ! Determinant of Ll is the product of diagonal elements
 det_Ll = ONE
 DO i = 1, n
    det_Ll = det_Ll * Ll(i, i)
 END DO

 ! Cholesky Decomposition of tAll -> Store in W1
 ! Also accumulate the determinant of tAll
 det_tAll = ONE
 DO i = 1, n
     DO j = i, n
         temp1 = tAll(i, j)
         DO k = i - 1, 1, -1
             temp1 = temp1 - W1(i, k) * W1(j, k)
         END DO
         
         IF (i == j) THEN
             W1(i, i) = SQRT(temp1)

              IF (temp1 <= 0.0_dprec) THEN
                 WRITE(*,*) 'ERROR: tAll not SPD at i=',i
                 ERROR STOP
              END IF       

             det_tAll = det_tAll * temp1
         ELSE
             W1(j, i) = temp1 / W1(i, i)
             W1(i, j) = ZERO
         END IF
     END DO
 END DO

 ! Invert tAll using its Cholesky factors (W1)
 DO i = 1, n
     W1(i, i) = ONE / W1(i, i)
     DO j = i + 1, n
         temp1 = ZERO
         DO k = i, j - 1
             temp1 = temp1 - W1(j, k) * W1(k, i)
         END DO
         W1(j, i) = temp1 / W1(j, j)
     END DO
 END DO 

 ! Construct the full inverse matrix inv_tAll
 DO i = 1, n
     DO j = i, n
         temp1 = ZERO
         DO k = j, n
             temp1 = temp1 + W1(k, i) * W1(k, j)
         END DO
         inv_tAll(i, j) = temp1
         inv_tAll(j, i) = temp1
     END DO
 END DO  

 !---------------------------------------------------------------------------
 ! 6. Compute Gamma Terms (Vector-Matrix-Vector products)
 !---------------------------------------------------------------------------
 
 ! Extract transformed vectors from Permutation Matrix
 DO i = 1, n
     tvl(i) = P(ml_1, i)  ! "v" vector
     twl(i) = P(ml_2, i)  ! "w" vector
 END DO

!  ! Calculate Matrix-Vector Products: vector' * inv_tAll
!  DO i = 1, n
!      temp1 = ZERO
!      temp2 = ZERO
!      DO j = 1, n
!          temp1 = temp1 + tvl(j) * inv_tAll(j, i)
!          temp2 = temp2 + twl(j) * inv_tAll(j, i)
!      END DO
!      tvl_inv_tAll(i) = temp1
!      twl_inv_tAll(i) = temp2
!  END DO

 ! Calculate Numerator Gammas (tgamma)
 ! Note: Since k=l, Bra = Ket. 
 ! tgamma1 = v' * A^-1 * v
 ! tgamma2 = w' * A^-1 * w
 tgamma1 = ZERO
 tgamma2 = ZERO
 tgamma5 = ZERO
 tgamma6 = ZERO
 DO i = 1, n
   ! tgamma1 = e_{ml1}' * inv_tAll * P'*e_{ml1}
   !         = (unpermuted bra)  *  (permuted ket)
   tgamma1 = tgamma1 + inv_tAll(ml_1, i) * tvl(i)
   
   ! tgamma2 = e_{ml2}' * inv_tAll * P'*e_{ml2}
   tgamma2 = tgamma2 + inv_tAll(ml_2, i) * twl(i)
   
   ! tgamma5 = e_{ml1}' * inv_tAll * P'*e_{ml2}
   tgamma5 = tgamma5 + inv_tAll(ml_1, i) * twl(i)
   
   ! tgamma6 = e_{ml2}' * inv_tAll * P'*e_{ml1}
   tgamma6 = tgamma6 + inv_tAll(ml_2, i) * tvl(i)
 END DO

 ! Calculate Denominator Gammas (from inv_All direct access)
 gamma1 = inv_All(ml_1, ml_1)
 gamma2 = inv_All(ml_2, ml_2)
 gamma5 = inv_All(ml_1, ml_2)
 gamma6 = inv_All(ml_2, ml_1)
 ! Combine Gammas
 tm = tgamma1 * tgamma2 + tgamma5 * tgamma6
 m  = gamma1  * gamma2  + gamma5  * gamma6

 !---------------------------------------------------------------------------
 ! 7. Final Overlap Calculation
 !---------------------------------------------------------------------------
 
 ! Determinant Ratio: (|Ll| * |Ll|) / |tAll|
 temp1 = (det_Ll * det_Ll) / det_tAll
 
 ! Apply power 3/2 (mathematically equivalent to |L|^3 / |A|^1.5)
 temp1 = ABS(temp1) * SQRT(ABS(temp1))
 
 ! Final Result
 Sll = Glob_2raised3n2 * temp1 * tm / m



END SUBROUTINE OverLapElementS2




SUBROUTINE MatrixElemenTranDipoleMoment(mk, vechLk, Pk, ml_1, ml_2, vechLl, Pl, &
                                        TranDipolLength_kl, TranDipolVelocity_kl)

                                       !This subroutine computes symmetry adapted matrix element with
 !a real L=1 and a real L=2 correlated Gaussians:
 !
 !fk = z_{mk} exp[-r'(Lk*Lk')r] 
 !
 !       —                    —
 !      | (v_k1' r)*(w_k1' r)+ | 
 ! fk = | (v_k2' r)*(w_k2' r)+ | exp[-r'(Lk*Lk')r]
 !      | (v_k3' r)*(w_k3' r)  |
 !       —                    —
 !
 !       —                                                     —
 !      |                                                       |
 !    = | x(ml_1)*x(ml_2) + y(ml_1)*y(ml_2) - 2 z(ml_1)*z(ml_2) |  exp[-r'(Lk*Lk')r]
 !      |                                                       |
 !       —                                                     —
 !  where m_l,m_k and mm_k are some integer between 1 and n (n is the number of 
 !  pseuDOparticles). Symmetry adaption is applied to the bra using 
 !  permutation matrix Pk and to the ket using permutation matrix Pl.
 !
 ! 
 !  Pl fl = (Pl z_{m_l}) exp[-r' {Pl'*(Ll*Ll')*Pl} r] 
 !
 !
 !       —                          —
 !      | (Pk v_k1' r)*(Pk w_k1' r)+ | 
 ! fk = | (Pk v_k2' r)*(Pk w_k2' r)+ | exp[-r'(Pk' Lk*Lk' Pk')r]
 !      | (Pk v_k3' r)*(Pk w_k3' r)  |
 !       —                          —
 ! this subroutine computes the following expressions
 ! which is a part of transition dipole momentum calculation.
 !
 !
 !
 !====================================================================
 !			Transition dipole integral in Lenght form
 !====================================================================
 !                               
 !                      ---  
 !                      \                  m_i         < Pk*fk | z_i | Pl*fl >
 ! TranDipolLength_kl = /  (q_i - Q_tot * -----)*-------------------------------------------
 !                      --                  m0      Sqrt( <fk|fk> ) * Sqrt( <fl|fl> )
 !                      i             
 !                                         
 !                                 
 !      < Pk*fk | z_i | Pl*fl >               
 ! ----------------------------------------- =  
 !   Sqrt( <fk|fk> ) * Sqrt( <fl|fl> )        
 !           
 !
 !               2^(3*n/2)     (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5       
 !            - ----------- * ----------------------------------------  
 !                sqrt(3)                (det_tAkl)^1.5                       
 !
 !                tgamma1 * gamma2 + gamma5 * gamma6   
 !             * --------------------------------------
 !                        Sqrt(vl'*inv_Akk*vl)         
 !
 !______________________________________________________________________________________________
 
 !Input:     
 !	m_l				    :: integer that determine which z-component is in the premultiplier of the Gaussian
 !	vechLk, vechLl     :: Arrays of length (n(n+1)/2) of exponential parameters. 
 !	Pk, Pl             :: The symmetry permutation matrices of size n x n
 
 
 !Output:
 !	TDkl               :: Matrix element (normalized)
 
 
 !==============================================================================================

  IMPLICIT NONE
  !---------------------------------------------------------------------------
  ! ARGUMENTS
  !---------------------------------------------------------------------------
  INTEGER, INTENT(IN)      :: mk, ml_1, ml_2
  REAL(dprec), INTENT(IN)  :: vechLk(Glob_np), vechLl(Glob_np)
  REAL(dprec), INTENT(IN)  :: Pk(Glob_n, Glob_n), Pl(Glob_n, Glob_n)
  REAL(dprec), INTENT(OUT) :: TranDipolLength_kl, TranDipolVelocity_kl
  !---------------------------------------------------------------------------
  ! PARAMETERS
  !---------------------------------------------------------------------------
  INTEGER, PARAMETER       :: nn  = Glob_MaxAllowedNumOfPseudoParticles
  INTEGER, PARAMETER       :: nnp = nn * (nn + 1) / 2
  !---------------------------------------------------------------------------
  ! LOCAL VARIABLES
  !---------------------------------------------------------------------------
  INTEGER                  :: n, np, Qtotal
  INTEGER                  :: i, j, k, indx
  REAL(dprec)              :: Lk(nn, nn), Ll(nn, nn)
  REAL(dprec)              :: tAk(nn, nn), tAl(nn, nn), tAkl(nn, nn)
  REAL(dprec)              :: W1(nn, nn), W2(nn, nn)
  
  REAL(dprec)              :: inv_Akk(nn, nn), inv_All(nn, nn), inv_tAkl(nn, nn)
  REAL(dprec)              :: tvk(Glob_np)               ! Transformed Bra vector v
  REAL(dprec)              :: tvl(Glob_np), twl(Glob_np)      ! Transformed Ket vectors v, w
  
  REAL(dprec)              :: tvkinv_tAkl(Glob_np)       ! v_k' * A^-1
  REAL(dprec)              :: inv_tAkl_twl(Glob_np)      ! A^-1 * w_l
  REAL(dprec)              :: inv_tAkl_tvl(Glob_np)      ! A^-1 * v_l
  REAL(dprec)              :: det_Lk, det_Ll, det_tAkl
  REAL(dprec)              :: temp01, temp02, temp1, temp2
  
  REAL(dprec)              :: tgamma1, tgamma5, tgamma2, tgamma6     ! Invariant Scalars
  REAL(dprec)              :: gamma1, gamma2, gamma5, gamma6         ! For normalization
  REAL(dprec)              :: tm1, tm3, tm_num, m_den, m1, m3
  !===========================================================================
  ! INITIALIZATION
  !===========================================================================
  n  = Glob_n
  np = Glob_np
  !---------------------------------------------------------------------------
  ! 1. Build Matrices Lk and Ll from vector storage
  !---------------------------------------------------------------------------
  indx = 0
  DO i = 1, n
      DO j = i, n
          indx = indx + 1
          Lk(i, j) = ZERO
          Lk(j, i) = vechLk(indx)
          Ll(i, j) = ZERO
          Ll(j, i) = vechLl(indx)
      END DO
  END DO
  !---------------------------------------------------------------------------
  ! 2. Compute A matrices: A = L * L^T
  !---------------------------------------------------------------------------
  DO i = 1, n
      DO j = i, n
          ! Compute Ak
          temp1 = ZERO
          DO k = 1, i
              temp1 = temp1 + Lk(i, k) * Lk(j, k)
          END DO 
          tAk(i, j) = temp1
          tAk(j, i) = temp1
          
          ! Compute Al
          temp1 = ZERO
          DO k = 1, i
              temp1 = temp1 + Ll(i, k) * Ll(j, k)
          END DO 
          tAl(i, j) = temp1
          tAl(j, i) = temp1
      END DO
  END DO
  !---------------------------------------------------------------------------
  ! 3. Apply Symmetry Permutations
  !    tAl_trans = P * Al * P^T  (using intermediate W1, W2)
  !---------------------------------------------------------------------------
  DO i = 1, n
      DO j = 1, n
          temp1 = ZERO
          temp2 = ZERO
          DO k = 1, n
              temp1 = temp1 + Pl(k, j) * tAl(k, i)
              temp2 = temp2 + tAk(j, k) * Pk(k, i)
          END DO
          W1(j, i) = temp1
          W2(j, i) = temp2
      END DO
  END DO

  ! Sum to get Total Exponent Matrix A_kl = Ak' + Al'
  DO i = 1, n  
      DO j = i, n
          temp1 = ZERO
          temp2 = ZERO
          DO k = 1, n
              temp1 = temp1 + W1(j, k) * Pl(k, i)
              temp2 = temp2 + Pk(k, j) * W2(k, i)
          END DO
          tAl(j, i)  = temp1
          tAl(i, j)  = temp1
          tAk(j, i)  = temp2
          tAk(i, j)  = temp2  
          tAkl(j, i) = temp1 + temp2
          tAkl(i, j) = temp1 + temp2
      END DO
  END DO
  !---------------------------------------------------------------------------
  ! 4. Determinants and Inversion of tAkl
  !---------------------------------------------------------------------------
  ! Determinants of L (product of diagonals)
  det_Lk = ONE
  det_Ll = ONE

  DO i = 1, n
     det_Lk = det_Lk * Lk(i, i)
     det_Ll = det_Ll * Ll(i, i)
  END DO

  ! Cholesky Decomposition of tAkl -> Store in W1
  det_tAkl = ONE
  DO i = 1, n
      DO j = i, n
          temp1 = tAkl(i, j)
          DO k = i - 1, 1, -1
              temp1 = temp1 - W1(i, k) * W1(j, k)
          END DO
          IF (i == j) THEN
              W1(i, i) = SQRT(temp1)

              IF (temp1 <= 0.0_dprec) THEN
                 WRITE(*,*) 'ERROR: tAkl not SPD at i=',i,' temp1=',temp1,' mk=',mk,' ml1=',ml_1,' ml2=',ml_2
                 ERROR STOP
              END IF

              det_tAkl = det_tAkl * temp1
          ELSE
              W1(j, i) = temp1 / W1(i, i)
              W1(i, j) = ZERO
          END IF
      END DO
  END DO

  ! Invert tAkl using Cholesky factor
  DO i = 1, n
      W1(i, i) = ONE / W1(i, i)
      DO j = i + 1, n
          temp1 = ZERO
          DO k = i, j - 1
              temp1 = temp1 - W1(j, k) * W1(k, i)
          END DO
          W1(j, i) = temp1 / W1(j, j)
      END DO
  END DO 

  ! Construct inv_tAkl from inverse Cholesky factors
  DO i = 1, n
      DO j = i, n
          temp1 = ZERO
          DO k = j, n
              temp1 = temp1 + W1(k, i) * W1(k, j)
          END DO
          inv_tAkl(i, j) = temp1
          inv_tAkl(j, i) = temp1
      END DO
  END DO  
  !---------------------------------------------------------------------------
  ! 5. Compute Transformed Vectors (Bra and Ket)
  !---------------------------------------------------------------------------
  DO i = 1, n
      tvk(i) = Pk(mk, i)     ! Bra vector v_k
      tvl(i) = Pl(ml_1, i)   ! Ket vector v_l
      twl(i) = Pl(ml_2, i)   ! Ket vector w_l
  END DO
  !---------------------------------------------------------------------------
  ! 6. Compute Inverse Normalization Matrices (inv_Akk, inv_All)
  !    Uses 0.5 * (L L^T)^-1
  !---------------------------------------------------------------------------
  DO i = 1, n
      W1(i, i) = ONE / Lk(i, i)
      W2(i, i) = ONE / Ll(i, i)
      DO j = i + 1, n
          temp1 = ZERO
          temp2 = ZERO
          DO k = i, j - 1
              temp1 = temp1 - Lk(j, k) * W1(k, i)
              temp2 = temp2 - Ll(j, k) * W2(k, i)
          END DO
          W1(j, i) = temp1 / Lk(j, j)
          W2(j, i) = temp2 / Ll(j, j)
      END DO
  END DO 

  DO i = 1, n
      DO j = i, n
          temp1 = ZERO
          temp2 = ZERO
          DO k = j, n
              temp1 = temp1 + W1(k, i) * W1(k, j)
              temp2 = temp2 + W2(k, i) * W2(k, j)       
          END DO
          inv_Akk(i, j) = ONEHALF * temp1
          inv_Akk(j, i) = ONEHALF * temp1
          inv_All(i, j) = ONEHALF * temp2
          inv_All(j, i) = ONEHALF * temp2  
      END DO
  END DO

  !---------------------------------------------------------------------------
  ! 7. Calculate Normalization Denominator Terms
  !---------------------------------------------------------------------------
  gamma1 = inv_All(ml_1, ml_1)    ! v_l * inv_All * v_l
  gamma2 = inv_All(ml_2, ml_2)    ! w_l * inv_All * w_l
  gamma5 = inv_All(ml_1, ml_2)    ! v_l * inv_All * w_l
  gamma6 = inv_All(ml_2, ml_1)    ! w_l * inv_All * v_l

  m1 = gamma1 * gamma2
  m3 = gamma5 * gamma6
  m_den = SQRT(m1 + m3) 

  IF(Verbose==2 .OR. Verbose==3) THEN

     IF (m_den /= m_den) THEN
        WRITE(*,*)
        WRITE(*,*) 'WARNING: m_den is NaN'
        ERROR STOP
     END IF

     IF(Verbose==3) THEN
        WRITE(*,*)
        WRITE(*,'(4X,A31,ES22.10)') 'Denominator (m_den): ', m_den
     END IF
     
   END IF

  !---------------------------------------------------------------------------
  ! 8. Calculate Constant Prefactor
  !---------------------------------------------------------------------------
  ! Determinant part: |Lk|^1.5 * |Ll|^1.5 / |Akl|^1.5
  temp1 = ABS(det_Ll * det_Lk) / det_tAkl

  IF(Verbose==2 .OR. Verbose==3) THEN

     IF (temp1 /= temp1) THEN
        WRITE(*,*) 'WARNING: |Lk|^1.5*|Ll|^1.5 / |Akl|^1.5 is NaN'
        ERROR STOP
     END IF

     IF(Verbose==3) &
        WRITE(*,'(4X,A31,ES22.10)') '|Lk|^1.5*|Ll|^1.5 / |Akl|^1.5: ', temp1

  END IF

  temp1 = temp1 * SQRT(temp1)

  ! Constants: 2^(3n/2) / sqrt(3)
  temp1 = temp1 * Glob_2raised3n2
  temp1 = temp1 / SQRT(Three)
  
  !  1 / sqrt(v_k * inv_Akk * v_k)
  ! Note: inv_Akk(mk, mk) is the diagonal element corresponding to v_k
  temp1 = temp1 / SQRT(inv_Akk(mk, mk))
  TranDipolLength_kl = temp1 / m_den

      IF(Verbose==2 .OR. Verbose==3) THEN

         IF (SQRT(inv_Akk(mk, mk)) /= SQRT(inv_Akk(mk, mk))) THEN
            WRITE(*,*) 'WARNING: SQRT(inv_Akk(mk, mk)) is NaN'
            ERROR STOP
         END IF

         IF(Verbose==3) &
            WRITE(*,'(4X,A31,ES22.10)') 'SQRT(inv_Akk(mk, mk)): ', SQRT(inv_Akk(mk, mk))

       END IF
  !---------------------------------------------------------------------------
  ! 9. PRE-CALCULATION FOR DIPOLE SUMMATION
  !---------------------------------------------------------------------------
  
  ! Calculate tv_k * A^-1 (Used for Scalar Gammas 1 & 5)
  DO i = 1, n
      temp1 = ZERO
      DO j = 1, n
          temp1 = temp1 + tvk(j) * inv_tAkl(j, i)
      END DO
      tvkinv_tAkl(i) = temp1
  END DO

  ! Calculate Gammas (Invariant in summation loop)
  ! gamma1 = tv_k * A^-1 * tv_l
  ! gamma5 = tv_k * A^-1 * tw_l
  tgamma1 = ZERO
  tgamma5 = ZERO
  DO j = 1, n
      tgamma1 = tgamma1 + tvkinv_tAkl(j) * tvl(j)
      tgamma5 = tgamma5 + tvkinv_tAkl(j) * twl(j)
  END DO

  ! Calculate Matrix-Vector Products for Loop (Vector Gammas 2 & 6)
  ! inv_tAkl_twl = A^-1 * w_l
  ! inv_tAkl_tvl = A^-1 * v_l
  DO i = 1, n
      temp1 = ZERO
      temp2 = ZERO
      DO j = 1, n
          temp1 = temp1 + inv_tAkl(i, j) * twl(j)
          temp2 = temp2 + inv_tAkl(i, j) * tvl(j)
      END DO
      inv_tAkl_twl(i) = temp1 
      inv_tAkl_tvl(i) = temp2
  END DO

  !---------------------------------------------------------------------------
  ! 10. DIPOLE MOMENT SUMMATION LOOP
  !     Sum over coordinate i (pseudo-particles)
  !---------------------------------------------------------------------------
  Qtotal = Glob_PseudoCharge0
  DO i = 1, n
     Qtotal = Qtotal + Glob_PseudoCharge(i)
  END DO
  
  temp01 = ZERO
  DO i = 1, n
      ! gamma2 corresponds to the i-th component of (A^-1 * tw_l)
      ! This effectively represents: w_k * A^-1 * tw_l where w_k = e_i
      tgamma2 = inv_tAkl_twl(i) 
      
      ! gamma6 corresponds to the i-th component of (A^-1 * tv_l)
      tgamma6 = inv_tAkl_tvl(i) 
      ! Formula: gamma1*gamma2 + gamma5*gamma6
      tm1 = tgamma1 * tgamma2
      tm3 = tgamma5 * tgamma6
      tm_num = (tm1 + tm3) 
      
         IF(Verbose==2 .OR. Verbose==3) THEN
  
            IF (tm_num /= tm_num) THEN
               WRITE(*,*) 'WARNING: tm_num is NaN'
               ERROR STOP
            END IF
  
            IF(Verbose==3) &
               WRITE(*,'(4X,A31,ES22.10)') 'tm_num: ', tm_num
  
          END IF
      
      ! Charge and Mass term: (q_i - Q_tot * m_i / m_0)
      temp02 = Glob_PseudoCharge(i) - Qtotal * Glob_Mass(i + 1) / Glob_MassTotal
      ! Accumulate
      temp01 = temp01 + temp02 * tm_num     
  END DO

  !---------------------------------------------------------------------------
  ! 11. Final Result
  !---------------------------------------------------------------------------
  TranDipolLength_kl   = - TranDipolLength_kl * temp01
  TranDipolVelocity_kl = TranDipolLength_kl


    IF(Verbose==2 .OR. Verbose==3) THEN

       IF (TranDipolLength_kl /= TranDipolLength_kl) THEN
          WRITE(*,*) 'WARNING: TranDipolLength_kl is NaN'
          ERROR STOP
       END IF

       IF(Verbose==3) THEN
          WRITE(*,'(4X,A31,ES22.10)') 'TranDipolLength_kl: ', TranDipolLength_kl
          WRITE(*,*)
       END IF

     END IF


     
END SUBROUTINE MatrixElemenTranDipoleMoment


end module matelem




