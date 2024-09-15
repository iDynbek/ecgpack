module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=0 and L=1 Gaussians.
use globvars
implicit none

contains

subroutine overlapMatrixElementsLP(m_k, vechLk, P, Skk)
	!This subroutine computes symmetry adapted matrix element with
	!two real L=1 correlated Gaussians:
	!
	!fk = z_{m_k} exp[-r'(Lk*Lk')r]
	!
	!where m_k is some integer between 1 and n (n is the number of
	!pseudoparticles). Symmetry adaption is applied to the ket using
	!permutation matrix P
	!
	!Input:
	!   m_k :: integer that determine which z-component is in the
	!                premultiplier of the Gaussian
	!   vechLk :: Array of length (n(n+1)/2) of
	!     exponential parameters.
	!   P  :: The symmetry permutation matrix of size n x n
	!Output:
	!   Skk	 ::	Overlap matrix element (normalized)
	
	!Arguments
	integer,intent(in)          :: m_k
	real(dprec),intent(in)      :: vechLk(Glob_np)
	real(dprec),intent(in)      :: P(Glob_n,Glob_n)
	real(dprec),intent(out)     :: Skk
	
	!Parameters (These are needed to declare static arrays. Using static
	!arrays makes the function call a little faster in comparison with
	!the case when arrays are dynamically allocated in stack)
	integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
	integer,parameter :: nnp=nn*(nn+1)/2
	
	!Local variables
	integer           n, np
	integer           tvl(nn)
	real(dprec)       Lk(nn,nn),Ll(nn,nn)
	real(dprec)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
	real(dprec)       inv_tAkl(nn,nn)
	real(dprec)       W1(nn,nn)
	real(dprec)       inv_tAkltvl(nn),vkinv_tAkl(nn)
	real(dprec)       temp1
	real(dprec)       det_tAkl
	real(dprec)       tau3
	integer           i,j,k,q,t,indx
	
	n=Glob_n
	np=Glob_np
	!First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
	indx=0
	do i=1,n
	  do j=i,n
		indx=indx+1
		Lk(i,j)=ZERO
		Lk(j,i)=vechLk(indx)
		Ll(i,j)=Lk(i,j)
		Ll(j,i)=Lk(j,i)
	  enddo
	enddo
	
	do i=1,n
	  do j=i,n
		temp1=ZERO
		do k=1,i
		  temp1=temp1+Lk(i,k)*Lk(j,k)
		enddo
		Ak(i,j)=temp1
		Ak(j,i)=temp1
		temp1=ZERO
		do k=1,i
		  temp1=temp1+Ll(i,k)*Ll(j,k)
		enddo
		tAl(i,j)=temp1
		tAl(j,i)=temp1
	  enddo
	enddo
	
	!Then we permute elements of Al to account for
	!the action of the permutation matrix
	!tAl=P'*Al*P
	!We also form matrix tAkl=Ak+tAl
	do i=1,n
	  do j=1,n
		temp1=ZERO
		do k=1,n
		   temp1=temp1+P(k,j)*tAl(k,i)
		enddo
		W1(j,i)=temp1
	  enddo
	enddo
	do i=1,n
	  do j=i,n
		temp1=ZERO
		do k=1,n
		   temp1=temp1+W1(i,k)*P(k,j)
		enddo
		  tAl(i,j)=temp1
		tAl(j,i)=temp1
		tAkl(i,j)=Ak(i,j)+temp1
		tAkl(j,i)=tAkl(i,j)
	  enddo
	enddo
	
	
	!After this we can do Cholesky factorization of tAkl.
	!The Cholesky factor will be temporarily stored in the
	!lower triangle of W1
	det_tAkl=ONE
	do i=1,n
	  do j=i,n
		temp1=tAkl(i,j)
		do k=i-1,1,-1
		  temp1=temp1-W1(i,k)*W1(j,k)
		enddo
		if (i==j) then
		  W1(i,i)=sqrt(temp1)
		  det_tAkl=det_tAkl*temp1
		else
		  W1(j,i)=temp1/W1(i,i)
		  W1(i,j)=ZERO
		endif
	  enddo
	enddo
	
	!Inverting tAkl using its Cholesky factor (stored in W1)
	!and placing the result into inv_tAkl
	do i=1,n
	  W1(i,i)=ONE/W1(i,i)
	  do j=i+1,n
		temp1=ZERO
		do k=i,j-1
		  temp1=temp1-W1(j,k)*W1(k,i)
		enddo
		W1(j,i)=temp1/W1(j,j)
	  enddo
	enddo
	
	do i=1,n
	  do j=i,n
		 temp1=ZERO
		 do k=j,n
		   temp1=temp1+W1(k,i)*W1(k,j)
		 enddo
		 inv_tAkl(i,j)=temp1
		 inv_tAkl(j,i)=temp1
	   enddo
	enddo
	
	!Computing tvl=P'*vl
	do i=1,n
	  tvl(i)=P(m_k,i)
	enddo
	
	!Compute inv_tAkltvl = inv_tAkl * tvl
	do i=1,n
	  temp1=ZERO
	  do j=1,n
		temp1=temp1+inv_tAkl(j,i)*tvl(j)
	  enddo
	  inv_tAkltvl(i)=temp1
	enddo
	
	!Compute vkinv_tAkl=vk'*inv_tAkl
	do i=1,n
	  vkinv_tAkl(i)=inv_tAkl(m_k,i)
	enddo
	
	!Compute tau3=vkinv_tAkl*tvl
	tau3=ZERO
	do i=1,n
	  tau3=tau3+vkinv_tAkl(i)*tvl(i)
	enddo
	
	!Evaluating overlap
	!temp1=abs(det_Ll*det_Lk)/det_tAkl
	temp1=det_tAkl*sqrt(det_tAkl)
	!Skl=Glob_2raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
	Skk=Glob_Piraised3n2*tau3/(TWO*temp1)
	
	
	end subroutine overlapMatrixElementsLP
	


subroutine spinPreCalc(n, nFactorial, parityFactor, SSNCmassChargeCoefficient, &
	SOmassChargeCoefficient, AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, &
	ketMatrix, spatialYoung0, spatialYoung1, &
	SSNCspinME, SiMinusME, SiPlusME, SziME, spinFreeME)
	use spinStuff
	implicit none
  
	!input vars:
	character(len = maxLen), intent(in) :: spatialYoung0, spatialYoung1
	integer, intent(in) :: n, nFactorial
  
	!output vars:
	real(dprec), dimension(n, n), intent(out) :: SSNCmassChargeCoefficient

	real(dprec), dimension(nFactorial), intent(out) :: parityFactor
	real(dprec), dimension(n, n, 4), intent(out) :: SOmassChargeCoefficient, AMMmassChargeCoefficient, &
	AMMFinmassChargeCoefficient
	real(dprec), dimension(n, n, nFactorial), intent(out) :: ketMatrix, SSNCspinME
	real(dprec), dimension(nFactorial, 2), intent(out) :: spinFreeME
	real(kind = dprec), dimension(n, nFactorial), intent(out) :: SiMinusME, SiPlusME, SziME 
	
	! local variables
	integer :: i, j, k, l, m
	character(len = maxLen) :: mySpatialYoung0, mySpatialYoung1
	integer, dimension(nFactorial) :: parities
	integer, dimension(n, n, nFactorial) :: allPermutations
	real(kind = dprec), dimension(:), allocatable :: finalSpinFunction0, finalSpinFunction1
	integer, dimension(:, :), allocatable :: primitives0, primitives1
	integer :: numberOfPrimitives0, numberOfPrimitives1

  
	SOmassChargeCoefficient = ZERO
	do i = 1, n
	  SOmassChargeCoefficient(i, i, 1) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / Glob_Mass(i + 1) * &
	  (ONE / Glob_Mass(i + 1) + TWO / Glob_Mass(1))
	enddo
  
	do i = 1, n
	  SOmassChargeCoefficient(i, i, 2) = -Glob_PseudoCharge0 * Glob_PseudoCharge(i) / &
	  (Glob_Mass(i + 1) * Glob_Mass(1))
	enddo
  
	do i = 1, n
	  do j = 1, n
		SOmassChargeCoefficient(i, j, 3) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
		(Glob_Mass(i + 1) * Glob_Mass(j + 1))
		SSNCmassChargeCoefficient(i, j) = Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
           (Glob_Mass(i + 1) * Glob_Mass(j + 1))
	  enddo
	enddo
  
	do i = 1, n
	  do j = 1, n
		SOmassChargeCoefficient(i, j, 4) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
		(Glob_Mass(i + 1)**TWO)
	  enddo
	enddo
  
	AMMmassChargeCoefficient = ZERO
	AMMfinmassChargeCoefficient = ZERO
	do i = 1, n
	  AMMmassChargeCoefficient(i, i, 1) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / Glob_Mass(i + 1)**TWO
	  AMMFinmassChargeCoefficient(i, i, 1) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / Glob_Mass(i + 1) * &
	  (ONE / Glob_Mass(1) + ONE / Glob_Mass(i + 1))
	enddo
  
	do i = 1, n
		AMMFinmassChargeCoefficient(i, i, 2) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / &
		(Glob_Mass(i + 1) * Glob_Mass(1))
	enddo

	do i = 1, n
	  do j = 1, n
		AMMmassChargeCoefficient(i, j, 3) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
		(Glob_Mass(i + 1) * Glob_Mass(j + 1))
		AMMFinmassChargeCoefficient(i, j, 3) = AMMmassChargeCoefficient(i, j, 3)
	  enddo
	enddo
  
	do i = 1, n
	  do j = 1, n
		AMMmassChargeCoefficient(i, j, 4) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
		(Glob_Mass(i + 1)**TWO)
		AMMFinmassChargeCoefficient(i, j, 4) = AMMmassChargeCoefficient(i, j, 4)
	  enddo
	enddo
  
	! now we deal with the spin stuff
	! rename the particles
	mySpatialYoung0 = spatialYoung0
	do i = 1, maxLen
	  if (mySpatialYoung0(i:i) == 'P') then
		read(mySpatialYoung0(i + 1:i + 1), *) k
		read(mySpatialYoung0(i + 2:i + 2), *) j
		write(mySpatialYoung0(i + 1:i + 1), '(i1)') k - 1
		write(mySpatialYoung0(i + 2:i + 2), '(i1)') j - 1
	  endif
	enddo

	mySpatialYoung1 = spatialYoung1
	do i = 1, maxLen
	  if (mySpatialYoung1(i:i) == 'P') then
		read(mySpatialYoung1(i + 1:i + 1), *) k
		read(mySpatialYoung1(i + 2:i + 2), *) j
		write(mySpatialYoung1(i + 1:i + 1), '(i1)') k - 1
		write(mySpatialYoung1(i + 2:i + 2), '(i1)') j - 1
	  endif
	enddo

	call generatePermutationMatrices(allPermutations, n, nFactorial, parities)
	
	!singlet WF
	call getSpinFunction(n, nFactorial, mySpatialYoung1, allPermutations, parities, &
	finalSpinFunction0, primitives0, numberOfPrimitives0)

	!triplet WF
	call getSpinFunction(n, nFactorial, mySpatialYoung0, allPermutations, parities, &
	finalSpinFunction1, primitives1, numberOfPrimitives1)


	call getSpinOpMeanValues(n, nFactorial, allPermutations, finalSpinFunction0, finalSpinFunction1, &
		primitives0, primitives1, numberOfPrimitives0, numberOfPrimitives1, &
		spinFreeME, SSNCspinME, SiMinusME, SiPlusME, SziME)

	ketMatrix = ZERO
	do i = 1, nFactorial
	  do k = 1, n
		do l = 1, n
		! note the transposition here
		ketMatrix(k, l, i) = real(allPermutations(l, k, i), kind=dprec)
		enddo
	  enddo
  
	enddo
  
	do i = 1, nFactorial
	  parityFactor(i) = real(parities(i), kind=dprec)
	enddo
  
  end subroutine spinPreCalc


  subroutine spinDependentMatrixElements(selectTransition, m_k, m_l, vechLk, vechLl, Pket, &
	SOspinME, SSNCspinME, SSNCmassChargeCoefficient, SOmassChargeCoefficient, &
	AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, SSNCkl, SO1kl, SO2kl, &
	AMM1kl, AMM2kl, AMM1finkl, AMM2finkl)
 !This subroutine computes symmetry adapted matrix element
 !with two real L=1 correlated Gaussians. These matrix element
 !is used in calculations of expectation values.

 !Input:
 !   m_k, m_l :: integers that determine which z-component is in the
 !       premultiplier of the Gaussian
 !   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.

 !Output (all matrix elements are computed with normilized functions):

 !   SSNCkl :: Non-contact spin-spin term (without the factor of alpha**2)
 !   SO1kl, SO2kl  :: Spin-Orbit corrections (without the factor of alpha**2)
 !         1 and 2 stay for spin-same orbit and spin-another orbit contributions
 !   AMM1kl, AMM2kl  :: AMM corrections (without the factor of alpha**2)
 !         1 and 2 stay for spin-same orbit and spin-another orbit contributions

 !Arguments
 integer,intent(in)	      :: selectTransition
 integer,intent(in)       :: m_k, m_l
 real(dprec),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
 real(dprec),intent(in)   :: Pket(Glob_n,Glob_n)

 real(dprec), intent(out)  :: SO1kl, SO2kl, AMM1kl, AMM2kl, AMM1finkl, AMM2finkl, SSNCkl
 !Parameters (These are needed to declare static arrays. Using static
 !arrays makes the function call a little faster in comparison with
 !the case when arrays are dynamically allocated in stack)
 integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
 integer,parameter :: nnp=nn*(nn+1)/2
 real(dprec),intent(in)   :: SSNCspinME(Glob_n, Glob_n), &
 							 SOspinME(Glob_n), &
							 SOmassChargeCoefficient(Glob_n, Glob_n, 4), &
							 AMMmassChargeCoefficient(Glob_n, Glob_n, 4), &
							 AMMFinmassChargeCoefficient(Glob_n, Glob_n, 4), &
							 SSNCmassChargeCoefficient(Glob_n, Glob_n)

 !Local variables
 integer           n, np
 integer           tvk(nn),tvl(nn)
 real(dprec)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
 real(dprec)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
 real(dprec)       inv_tAkl(nn,nn)


 real(dprec)       W1(nn,nn)
 real(dprec)       temp1, temp2, det_tAkl
 integer :: i, j, k, indx


 integer :: pm_k, pm_l ! new non-zero components of v_k and v_l
 real(dprec) :: commonFactorSO, commonFactorSS, gamma, gamma_diag, &
 				jiVl, jiAlAklinvVk, jiAlAklinvVl, jiAklinvVk, jiAklinvVl, &
				jjAlAklinvVl, jjAlAklinvVk, jjAklinvVk, jjAklinvVl, jjVl, localEps

 integer :: indexI, indexJ ! indeces enumerating particles from H_SO and AMM operators

 localEps = 1.d-14 ! if the corresponding spin mean value is less then localEps, we don't calculate the spatial part

 ! basically copy-paste from the old ExpecVals subroutine
 n=Glob_n
 np=Glob_np

 !First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
 indx=0
 do i=1,n
   do j=i,n
	 indx=indx+1
	 Lk(i,j)=ZERO
	 Lk(j,i)=vechLk(indx)
	 Ll(i,j)=ZERO
	 Ll(j,i)=vechLl(indx)
   enddo
 enddo

 do i=1,n
   do j=i,n
	 temp1=ZERO
	   do k=1,i
		 temp1=temp1+Lk(i,k)*Lk(j,k)
	   enddo
	   tAk(i,j)=temp1
	   tAk(j,i)=temp1
	 temp1=ZERO
	   do k=1,i
		 temp1=temp1+Ll(i,k)*Ll(j,k)
	   enddo
	   tAl(i,j)=temp1
	   tAl(j,i)=temp1
   enddo
 enddo

 !Then we permute elements of Al to account for
 !the action of the permutation matrix
 !  tAl=Pket'*Al*Pket

 !We also form matrix tAkl=tAk+tAl
 do i=1,n
   do j=1,n
	 temp1=ZERO
	 temp2=ZERO
	 do k=1,n
		temp1=temp1+Pket(k,j)*tAl(k,i)
	 enddo
	 W1(j,i)=temp1
   enddo
 enddo
 !tAl=W1*Pket

 do i=1,n
   do j=i,n
	 temp1=ZERO
	 temp2=ZERO
	 do k=1,n
		temp1=temp1+W1(j,k)*Pket(k,i)
	   enddo
	   tAl(j,i)=temp1
	   tAl(i,j)=temp1
	   tAkl(j,i)=temp1+tAk(j,i)
	   tAkl(i,j)=temp1+tAk(i,j)
   enddo
 enddo

 !After this we can do Cholesky factorization of tAkl.
 !The Cholesky factor will be temporarily stored in the
 !lower triangle of W1
 det_tAkl=ONE
 do i=1,n
   do j=i,n
	 temp1=tAkl(i,j)
	 do k=i-1,1,-1
	   temp1=temp1-W1(i,k)*W1(j,k)
	 enddo
	 if (i==j) then
	   W1(i,i)=sqrt(temp1)
	   det_tAkl=det_tAkl*temp1
	 else
	   W1(j,i)=temp1/W1(i,i)
	   W1(i,j)=ZERO
	 endif
   enddo
 enddo

 !Inverting tAkl using its Cholesky factor (stored in W1)
 !and placing the result into inv_tAkl
 do i=1,n
   W1(i,i)=ONE/W1(i,i)
   do j=i+1,n
	 temp1=ZERO
	 do k=i,j-1
	   temp1=temp1-W1(j,k)*W1(k,i)
	 enddo
	 W1(j,i)=temp1/W1(j,j)
   enddo
 enddo

 do i=1,n
   do j=i,n
	  temp1=ZERO
	  do k=j,n
		temp1=temp1+W1(k,i)*W1(k,j)
	  enddo
	  inv_tAkl(i,j)=temp1
	  inv_tAkl(j,i)=temp1
	enddo
 enddo

 ! new code from here

 ! new m_k and m_l
 ! new v_l = (P TRANSPOSED) * v_l
 pm_k = m_k

 pm_l = m_l
 do i = 1, n
   if (abs(Pket(m_l, i) - 1.d0) < 1.d-13) then ! for integers it would be == 1
	 pm_l = i
	 exit
   endif
 enddo

 !common factor
 commonFactorSO = ZERO
 commonFactorSS = ZERO
 if (selectTransition == 1) then
 	commonFactorSO = TWO * Glob_Piraised3n2 / (SQRTPI * det_tAkl * sqrt(det_tAkl))
	commonFactorSS = TWO * Glob_Piraised3n2 / (SQRTPI * det_tAkl * sqrt(det_tAkl))
 else if (selectTransition == 2) then
	commonFactorSO = -TWO * Glob_Piraised3n2 / (SQRTPI * det_tAkl * sqrt(det_tAkl))
 endif

 SO1kl = ZERO
 SO2kl = ZERO

 AMM1kl = ZERO
 AMM2kl = ZERO
 AMM1finkl = ZERO
 AMM2finkl = ZERO

 SSNCkl = ZERO


 do indexI = 1, n

   if (abs(SOspinME(indexI)) < localEps) cycle
 

   ! gamma diagonal coefficient
   gamma_diag = ONE / sqrt(inv_tAkl(indexI, indexI))

   gamma = gamma_diag ! for spin-same-orbit

   ! kronecker deltas
   jiVl = ZERO
   if (pm_l == indexI) then
	 jiVl = ONE
   endif

   ! calculating all the traces we need
   ! tr(Axy') is computed as (y, Ax) everywhere
   ! variable names: jiAlAklinvVk = (j^i, A_l A_{kl}^(-1) v_k) (names doesnt account for permutations)


   jiAlAklinvVl = ZERO
   do i = 1, n
	 jiAlAklinvVl = jiAlAklinvVl + tAl(indexI, i) * inv_tAkl(i, pm_l)
   enddo

   jiAlAklinvVk = ZERO
   do i = 1, n
	 jiAlAklinvVk = jiAlAklinvVk + tAl(indexI, i) * inv_tAkl(i, pm_k)
   enddo

   jiAklinvVl = inv_tAkl(indexI, pm_l)

   jiAklinvVk = inv_tAkl(indexI, pm_k)

   ! diagonal (spin-same orbit) matrix element
   temp1 =  &
   gamma**3 / THREE * (jiVl * jiAklinvVk + &
   jiAlAklinvVk * jiAklinvVl - &
   jiAklinvVk * jiAlAklinvVl)


	SO1kl = SO1kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
	AMM1kl = AMM1kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
	AMM1finkl = AMM1finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 1) * temp1


   ! these traces are needed for spin-other orbit contribution

   do indexJ = 1, n
	 if (indexI == indexJ) cycle

	 gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
	 inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))

	 ! we need more traces

	 jjAlAklinvVl = ZERO
	 do i = 1, n
	   jjAlAklinvVl = jjAlAklinvVl + tAl(indexJ, i) * inv_tAkl(i, pm_l)
	 enddo

	 jjAlAklinvVk = ZERO
	 do i = 1, n
	   jjAlAklinvVk = jjAlAklinvVk + tAl(indexJ, i) * inv_tAkl(i, pm_k)
	 enddo

	 jjAklinvVk = inv_tAkl(indexJ, pm_k)
	 
	 jjAklinvVl = inv_tAkl(indexJ, pm_l)

	 ! kronecker deltas
	 jjVl = ZERO
	 if (pm_l == indexJ) then
	   jjVl = ONE
	 endif


	 ! NOTE for this term we need gamma_diag (tr[Cklinv J_ii]), not gamma (tr[Cklinv J_ij])
	 temp1 = &
	 gamma_diag**3 / THREE *  (jiAklinvVk * jjVl + &
	 jjAlAklinvVk * jiAklinvVl - &
	 jiAklinvVk * jjAlAklinvVl)
	 

	 SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
	 AMM2finkl = AMM2finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 2) * temp1
	 
	 temp1 = &
	 gamma**3 / THREE * (jjVl * (jjAklinvVk - jiAklinvVk) + &
	 jjAlAklinvVk * (jjAklinvVl - jiAklinvVl) + &
	 jjAlAklinvVl * (jiAklinvVk - jjAklinvVk)) 


	SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
	AMM2kl = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
	AMM2finkl = AMM2finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexJ, 3) * temp1

	 
	 temp1 = &
	 gamma**3 / THREE * (jiVl * (jiAklinvVk - jjAklinvVk) + &
	 jiAlAklinvVk * (jiAklinvVl - jjAklinvVl) + &
	 jiAlAklinvVl * (jjAklinvVk - jiAklinvVk)) 


	SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
	AMM2kl = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
	AMM2finkl = AMM2finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexJ, 4) * temp1
	

   enddo ! indexJ cycle

 enddo ! indexI cycle


 ! SS term (separate loop for not to interfere with the case Siz == 0 at previous loop)
 
 do indexI = 1,n
	do indexJ = indexI + 1,n !(we need only indexJ > indexI)
	   if (abs(SSNCspinME(indexI,indexJ)) < localEps)  cycle
	   gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
	   inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))
 
	   jjAklinvVk = inv_tAkl(indexJ, pm_k)
	   jjAklinvVl = inv_tAkl(indexJ, pm_l)
	   jiAklinvVl = inv_tAkl(indexI, pm_l)
	   jiAklinvVk = inv_tAkl(indexI, pm_k)

	   temp1 = &
		  (gamma**5 / 15._dprec) * ( jiAklinvVk * (jiAklinvVl  - jjAklinvVl) + &
		  jjAklinvVk * (jjAklinvVl - jiAklinvVl) ) !additional factor of 1/sqrt(6) is taken from spin part 

		SSNCkl = SSNCkl + SSNCspinME(indexI, indexJ) * SSNCmassChargeCoefficient(indexI, indexJ) * temp1
	   
	enddo !indexJ loop
 enddo !indexI loop
   
 SSNCkl = SSNCkl * commonFactorSS
 SO1kl = SO1kl * commonFactorSO
 SO2kl = SO2kl * commonFactorSO
 AMM1kl = AMM1kl * commonFactorSO
 AMM2kl = AMM2kl * commonFactorSO
 AMM1finkl = AMM1finkl * commonFactorSO
 AMM2finkl = AMM2finkl * commonFactorSO


end subroutine spinDependentMatrixElements

end module matelem