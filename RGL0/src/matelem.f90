module matelem
!Module matelem contains subroutines for computing 
!matrix elements with real L=0 Gaussians without normalization.
use globvars
implicit none

contains


subroutine MatrixElements(vechLk, vechLl, P, &
               Hkl, Skl, Dk, Dl, grad_k, grad_l)
!This subroutine computes symmetry adapted matrix element with 
!two real L=0 correlated Gaussians:
! 
!fk =  exp[-r'(Lk*Lk')r] 
!
!Symmetry adaption is applied to the ket using 
!permutation matrices Glob_YHYMatr(:,:,1:Glob_NumYHYTerms)
!
!Input:     
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters. 
!   P   :: The symmetry permutation matrix of size n x n
!   grad_k, grad_l :: Gradient flags
!   grad_k=.true.  means that dHkldvechLk, dSkldvechLk need to be computed. 
!   grad_l=.true.  means that dHkldvechLl, dSkldvechLl need to be computed.
!Output:
!   Hkl	 ::	Hamiltonian term (normalized)
!   Skl	 ::	Overlap matrix element (normalized) 
!   Dk,Dl:: derivatives of normalized Hkl and Skl wrt vechLk
!           and vechLl respectively. They are ordered in the 
!           following manner:
!           Dk=(dHkldvechLk,dSkldvechLk)
!           Dl=(dHkldvechLl,dSkldvechLl)


!Arguments
real(dprec),intent(in)      :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec),intent(in)      :: P(Glob_n,Glob_n)
real(dprec),intent(out)     :: Skl,Hkl
real(dprec),intent(out)     :: Dk(2*Glob_np),Dl(2*Glob_np)
logical,intent(in)          :: grad_k, grad_l

!Parameters (These are needed to declare static arrays. Using static 
!arrays makes the function call a little faster in comparison with 
!the case when arrays are dynamically allocated in stack)
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
integer           n, np
real(dprec)       dHkldvechLk(nnp), dHkldvechLl(nnp)
real(dprec)       dSkldvechLk(nnp), dSkldvechLl(nnp)
real(dprec)       Lk(nn,nn), Ll(nn,nn), inv_Lk(nn,nn), inv_Ll(nn,nn)
real(dprec)       Ak(nn,nn), tAl(nn,nn), tAkl(nn,nn)
real(dprec)       inv_tAkl(nn,nn), inv_ttAkl(nn,nn)
real(dprec)       inv_tAkltAlM(nn,nn)
real(dprec)       tr_inv_tAklJij32(nn,nn)
real(dprec)       F(nn,nn),G(nn,nn)
real(dprec)       tQ(nn,nn),ttQ(nn,nn)
real(dprec)       WVcLk(nnp), WVcLl(nnp)
real(dprec)       W1(nn,nn), W2(nn,nn), W3(nn,nn)
real(dprec)       temp1, temp2, temp3, temp4, temp5
real(dprec)       det_Lk, det_Ll, det_tAkl
real(dprec)       Tkl, Vkl
integer           i, j, k, kk, kkk, q, t, indx

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

!The determinants of Lk and Ll are just
!the products of their diagonal elements
!det_Lk=ONE
!det_Ll=ONE
!do i=1,n
!  det_Lk=det_Lk*Lk(i,i)
!  det_Ll=det_Ll*Ll(i,i)
!enddo

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

!Evaluating overlap

!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2raised3n2*temp1*sqrt(temp1)
Skl=Glob_Piraised3n2/(det_tAkl*sqrt(det_tAkl))  !new line

!Doing multiplication W2=inv_tAkl*tAl
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
    enddo
    W2(j,i)=temp1
  enddo
enddo

!Doing multiplication inv_tAkltAlM=inv_tAkl*tAl*M=W2*M
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+W2(j,k)*Glob_MassMatrix(k,i)
    enddo
    inv_tAkltAlM(j,i)=temp1
  enddo
enddo

!Computing kinetic energy, Tkl=tr[inv_tAkltAlM*Ak]
Tkl=ZERO
do i=1,n
  temp1=ZERO
  do k=1,n
    temp1=temp1+inv_tAkltAlM(i,k)*Ak(k,i)
  enddo
  Tkl=Tkl+temp1
enddo
Tkl=SIX*Skl*Tkl

!Evaluating potential energy, Vkl, and tr[invCkl*Jij]^(-3/2)
!The lower triangle of array trinvCklJij32
!will contain the corresponding quantities.
temp1=(TWO/SQRTPI)*Skl
Vkl=ZERO
do i=1,n
  temp3=inv_tAkl(i,i)
  temp4=sqrt(temp3)
  tr_inv_tAklJij32(i,i)=1/(temp4*temp3)
  temp5=temp1/temp4
  Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*temp5
enddo
do i=1,n
  do j=i+1,n
    temp3=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    temp4=sqrt(temp3)
    tr_inv_tAklJij32(j,i)=1/(temp4*temp3)
	temp5=temp1/temp4
    Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*temp5
  enddo
enddo

Hkl=Tkl+Vkl

!Now we start computing the gradient of Skl

if (grad_k) then
  !Multiplying inv_tAkl and Lk and storing it in W1
  do i=1,n
    do j=1,n
	  temp1=ZERO
	  do k=j,n
        temp1=temp1+inv_tAkl(i,k)*Lk(k,j)
	  enddo
      W1(i,j)=temp1
    enddo
  enddo
  !I delete not-necessary parts of the code
  
  !Inverting Lk. The inverse of a lower triangular matrix is a lower
  !triangular matrix. The inverse of Lk will be stored in inv_Lk. The upper
  !triangle of inv_Lk is set to be zero.
 ! do i=1,n
 !   do j=i,n
 !     inv_Lk(j,i)=Lk(j,i) 
 !   enddo
 ! enddo
 ! do i=1,n
 !   inv_Lk(i,i)=1/inv_Lk(i,i)
 !   do j=i+1,n
!	  inv_Lk(i,j)=ZERO
!      temp1=ZERO
!      do k=i,j-1
!        temp1=temp1-inv_Lk(j,k)*inv_Lk(k,i)
!      enddo
!      inv_Lk(j,i)=temp1/inv_Lk(j,j)
!    enddo
!  enddo 
  !storing dSkldvechLk
  !temp1=Skl*THREE/TWO
  temp1=-Skl*THREE  ! new line
  indx=0
  do i=1,n
    do j=i,n
	  indx=indx+1
      !dSkldvechLk(indx)=(inv_Lk(i,j)-TWO*W1(j,i))*temp1
      dSkldvechLk(indx)=W1(j,i)*temp1  !new line
    enddo
  enddo 
endif !end if (grad_k)

if (grad_l) then
  !calculating inv_ttAkl: inv_ttAkl=P*inv_tAkl*P'
  do i=1,n
    do j=1,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+P(i,k)*inv_tAkl(k,j) 
	  enddo
      W1(i,j)=temp1
	enddo
  enddo
  do i=1,n
    do j=i,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+W1(i,k)*P(j,k)
	  enddo
	  inv_ttAkl(i,j)=temp1
      inv_ttAkl(j,i)=temp1
	enddo
  enddo
  !Multiplying inv_ttAkl and Ll and storing the result in W1
  do i=1,n
    do j=1,n
	  temp1=ZERO
	  do k=1,n
        temp1=temp1+inv_ttAkl(i,k)*Ll(k,j)
	  enddo
      W1(i,j)=temp1
	enddo
  enddo
  !I delete not-necessary parts of the code
  
  !Inverting Ll. The inverse of a lower triangular matrix 
  !is a lower triangular matrix. The inverse of Ll will be 
  !stored in inv_Ll. The upper triangle of inv_Ll is not set to 
  !be zero.
  !do i=1,n
  !  do j=i,n
  !    inv_Ll(j,i)=Ll(j,i) 
  !  enddo
  !enddo
  !do i=1,n
  !  inv_Ll(i,i)=1/inv_Ll(i,i)
  !  do j=i+1,n
!	  inv_Ll(i,j)=ZERO
!      temp1=ZERO
!      do k=i,j-1
!        temp1=temp1-inv_Ll(j,k)*inv_Ll(k,i)
!      enddo
 !     inv_Ll(j,i)=temp1/inv_Ll(j,j)
 !   enddo
 ! enddo  
  !storing dSkldvechLl
!  temp1=Skl*THREE/TWO
  temp1=-Skl*THREE  !new line
  indx=0
  do i=1,n
	do j=i,n
	  indx=indx+1
    !  dSkldvechLl(indx)=(inv_Ll(i,j)-TWO*W1(j,i))*temp1
      dSkldvechLl(indx)=W1(j,i)*temp1  !new line
	enddo
  enddo 
endif !end if (grad_l)
  
!Gradient of Tkl

if (grad_l) then
  !Compute matrix G=P*inv_tAkl*Ak*M*Ak*inv_tAkl*P' 
  do i=1,n
    do j=1,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+Ak(i,k)*Glob_MassMatrix(k,j)
	  enddo
	  W1(i,j)=temp1
    enddo
  enddo
  do i=1,n
    do j=i,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+W1(i,k)*Ak(k,j)
	  enddo
	  W2(i,j)=temp1
	  W2(j,i)=temp1
    enddo
  enddo
  do i=1,n
    do j=1,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+inv_tAkl(i,k)*W2(k,j)
	  enddo
	  W1(i,j)=temp1
    enddo
  enddo
  do i=1,n
    do j=i,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+W1(i,k)*inv_tAkl(k,j)
	  enddo
	  W2(i,j)=temp1
	  W2(j,i)=temp1
    enddo
  enddo
  do i=1,n
    do j=1,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+P(i,k)*W2(k,j) 
	  enddo
      W1(i,j)=temp1 
	enddo
  enddo
  do i=1,n
    do j=i,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+W1(i,k)*P(j,k)
	  enddo
	  G(i,j)=temp1
      G(j,i)=temp1
	enddo
  enddo
  !Computing the product G*Ll and 
  !putting it in array W1.
  do i=1,n
    do j=1,n
      temp1=ZERO
	  do k=j,n
        temp1=temp1+G(i,k)*Ll(k,j)
	  enddo
	  W1(i,j)=temp1
    enddo
  enddo
  !storing dTkldvechLl
  temp1=Tkl/Skl
  temp2=12*Skl
  indx=0
  do i=1,n
    do j=i,n
	  indx=indx+1
      dHkldvechLl(indx)=temp1*dSkldvechLl(indx)+temp2*W1(j,i)
    enddo
  enddo
endif !end if (grad_l)

if (grad_k) then
  !Compute matrix F=inv_tAkl*Al*M*Al*inv_tAkl=inv_tAklAlM*Al*inv_tAkl 
  do i=1,n
    do j=1,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+inv_tAkltAlM(i,k)*tAl(k,j)
	  enddo
	  W1(i,j)=temp1
    enddo
  enddo
  do i=1,n
    do j=i,n
      temp1=ZERO
	  do k=1,n
        temp1=temp1+W1(i,k)*inv_tAkl (k,j)
	  enddo
	  F(i,j)=temp1
	  F(j,i)=temp1
    enddo
  enddo
  !Computing the product F*Lk and 
  !putting it in array W1.
  do i=1,n
    do j=1,n
      temp1=ZERO
	  do k=j,n
        temp1=temp1+F(i,k)*Lk(k,j)
	  enddo
	  W1(i,j)=temp1
    enddo
  enddo
  !storing dTkldvechLk
  temp1=Tkl/Skl
  temp2=12*Skl
  indx=0
  do i=1,n
    do j=i,n
	  indx=indx+1
      dHkldvechLk(indx)=temp1*dSkldvechLk(indx)+temp2*W1(j,i)
    enddo
  enddo
endif !end if (grad_k)


!Gradient of Vkl

if (grad_l.OR.grad_k) then
  !First we set to zero some work arrays
  if (grad_k) WVcLk(1:np)=ZERO
  if (grad_l) WVcLl(1:np)=ZERO
  !Now we can proceed
  do i=1,n
    do j=1,i
      !Evaluating tQ=inv_tAkl*Jij*inv_tAkl
	  if (i==j) then
        do k=1,n
          tQ(k,k)=inv_tAkl(k,i)*inv_tAkl(i,k)
          do kk=k+1,n
	         temp1=inv_tAkl(k,i)*inv_tAkl(i,kk)
             tQ(k,kk)=temp1
	         tQ(kk,k)=temp1
	      enddo
        enddo
      else
        do k=1,n
          tQ(k,k)=inv_tAkl(k,i)*inv_tAkl(i,k)+inv_tAkl(k,j)*inv_tAkl(j,k) &
		          -inv_tAkl(k,i)*inv_tAkl(j,k)-inv_tAkl(k,j)*inv_tAkl(i,k)
          do kk=k+1,n
	         temp1=inv_tAkl(k,i)*inv_tAkl(i,kk)+inv_tAkl(k,j)*inv_tAkl(j,kk) &
		          -inv_tAkl(k,i)*inv_tAkl(j,kk)-inv_tAkl(k,j)*inv_tAkl(i,kk)
             tQ(k,kk)=temp1
	         tQ(kk,k)=temp1
	      enddo
        enddo
	  endif
	  if (grad_k) then
        !Multiplying tQ by Lk
	    !and storing the results in the lower
	    !triangle of W1. 	 
	    do k=1,n
          do kk=1,k
	         temp2=ZERO
	         do kkk=kk,n
               temp2=temp2+tQ(k,kkk)*Lk(kkk,kk)
		     enddo
		     W1(k,kk)=temp2
	      enddo
	    enddo
      endif 
	  if (grad_l) then
	    !Evaluating ttQ=P*tQ*P'
        do k=1,n
          do kk=1,n
            temp2=ZERO
	        do kkk=1,n
              temp2=temp2+P(k,kkk)*tQ(kkk,kk) 
	        enddo
            W2(k,kk)=temp2 
	      enddo
        enddo
        do k=1,n
          do kk=k,n
            temp2=ZERO
	        do kkk=1,n
              temp2=temp2+W2(k,kkk)*P(kk,kkk)
	        enddo
	        ttQ(k,kk)=temp2
	        ttQ(kk,k)=temp2
	      enddo
        enddo
        !Multiplying ttQ by Ll
	    !and storing the results in the lower
	    !triangle of W2. 	 
	    do k=1,n
          do kk=1,k
	         temp2=ZERO
	         do kkk=kk,n
               temp2=temp2+ttQ(k,kkk)*Ll(kkk,kk)
		     enddo
		     W2(k,kk)=temp2
	      enddo
	    enddo
      endif 

	  !Calculating ij-terms of the sums in the expressions for
	  !the dVkldvechLk and dVkldvechLl 
      if (i==j) then
	temp1=ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(i))*tr_inv_tAklJij32(i,i)
      else
	temp1=ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*tr_inv_tAklJij32(i,j)           
	  endif
	  if (grad_k) then
	    indx=0
	    do k=1,n
		  indx=indx+1
		  WVcLk(indx)=WVcLk(indx)+temp1*W1(k,k)
	      do kk=k+1,n
		     indx=indx+1
			 WVcLk(indx)=WVcLk(indx)+temp1*W1(kk,k)
		  enddo
        enddo 
	  endif
	  if (grad_l) then
	    indx=0
	    do k=1,n
		  indx=indx+1
		  WVcLl(indx)=WVcLl(indx)+temp1*W2(k,k)
	      do kk=k+1,n
		     indx=indx+1
			 WVcLl(indx)=WVcLl(indx)+temp1*W2(kk,k)
		  enddo
        enddo 
	  endif      
	enddo
  enddo
  !Multiplying by common factors and getting the final 
  !result for the gradient of Vkl
  temp1=(TWO/SQRTPI)*Skl
  if (grad_k) dHkldvechLk(1:np)=dHkldvechLk(1:np)+(Vkl/Skl)*dSkldvechLk(1:np)+temp1*WVcLk(1:np)
  if (grad_l) dHkldvechLl(1:np)=dHkldvechLl(1:np)+(Vkl/Skl)*dSkldvechLl(1:np)+temp1*WVcLl(1:np)
endif

!Packing derivatives into the output array
if (grad_k) then
  Dk(1:np)=dHkldvechLk(1:np)
  Dk(np+1:2*np)=dSkldvechLk(1:np)
endif
if (grad_l) then
  Dl(1:np)=dHkldvechLl(1:np)
  Dl(np+1:2*np)=dSkldvechLl(1:np)
endif


end subroutine MatrixElements



subroutine MatrixElementsForExpcVals(vechLk, vechLl, Pbra, Pket, &
           Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
           MVkl, drach_MVkl, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, del2kl, prvalkl, &
           wf2originkl, NumCFGridPoints, CFGrid, CFkl, NumDensGridPoints, DensGrid, Denskl, &
           AreCorrFuncNeeded, ArePartDensNeeded, AreMCorrFuncNeeded, AreMPartDensNeeded)
!This subroutine computes symmetry adapted matrix elements 
!with two real L=0 correlated Gaussians. These matrix elements
!are used in calculations of expectation values.
!Symmetry adaption is applied to the bra and ket using permutation matrices Pbra and Pket
!
!Input:     
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters. 
!   Pbra    :: The symmetry permutation matrix of size n x n that is applied to bra
!   Pket    :: The symmetry permutation matrix of size n x n that is applied to ket
!Output (all matrix elements are computed with normilized functions):
!   Hkl	     ::	Hamiltonian
!   Skl	     ::	Overlap
!   Tkl      :: Kinetic energy
!   Vkl      :: Potential energy
!   rm2kl    :: 1/r_i^2, 1/r_{ij}^2
!   rmkl     :: 1/r_i, 1/r_{ij}
!   rkl      :: r_i, r_{ij}
!   r2kl     :: r_i^2, r_{ij}^2
! deltarkl   :: delta(r_i), delta(r_{ij})
! drach_deltarkl:: Drachmanized delta(r_i), delta(r_{ij})
!   MVkl     :: Mass-velocity correction (without the factor of alpha**2)
! drach_MVkl :: Drachminized mass-velocity correction (without the factor of alpha**2)
! Darwinkl  :: Darwin correction (without the factor of alpha**2)
! drach_Darwinkl:: Drachmanized Darwin correction (without the factor of alpha**2)
!   OOkl    :: Orbit-Orbit correction (without the factor of alpha**2)
! rmrmkl    :: 1/(r_{ij}*r_{pq})
! del2kl    :: delta(r_{ij})delta(r_{pq}) when r_{ij}/=r_{pq}
! prvalkl   :: P(1/r^3_ij) - principal values of matric element 1/r^3_ij  (appears in the Araki-Sucker term for QED correction)                
! wf2originkl:: n-particle density at all-particle coalescence point (absolute square of the wave function at the origin)
!NumCFGridPoints   :: Number of grid points for correlation function calculations
!CFGrid            :: Array containing grid points where matrix elements of 
!                     correlation functions should be computed   
!CFkl              :: Matrix elements of correlation functions
!NumDensGridPoints :: Number of grid points for particle density calculations
!DensGrid          :: Array containing grid points where matrix elements of 
!                     particle densities should be computed
!Denskl            :: Matrix elements of particle densities
!AreCorrFuncNeeded :: flag indicating whether matrix elements of correlation 
!                     functions need be computed
!ArePartDensNeeded :: flag indicating whether matrix elements of particle
!                     densities need be computed


!Arguments
real(dprec),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec),intent(in)   :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
real(dprec),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl
real(dprec),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
real(dprec),intent(out)  :: deltarkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: drach_deltarkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: prvalkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: rmrmkl(Glob_n,Glob_n,Glob_n,Glob_n)
real(dprec),intent(out)  :: del2kl(Glob_n,Glob_n,Glob_n,Glob_n)
real(dprec),intent(out)  :: wf2originkl
integer,intent(in)       :: NumCFGridPoints,NumDensGridPoints
real(dprec),intent(in)   :: CFGrid(NumCFGridPoints),DensGrid(NumDensGridPoints)
real(dprec),intent(out)  :: CFkl(Glob_n*(Glob_n+1)/2,NumCFGridPoints)
real(dprec),intent(out)  :: Denskl(Glob_n+1,NumDensGridPoints)
logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded,AreMCorrFuncNeeded, AreMPartDensNeeded

!Parameters (These are needed to declare static arrays. Using static 
!arrays makes the function call a little faster in comparison with 
!the case when arrays are dynamically allocated in stack)
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles

!Local variables
integer           n, np
real(dprec)       Lk(nn,nn), Ll(nn,nn), inv_Lk(nn,nn), inv_Ll(nn,nn)
real(dprec)       tAk(nn,nn), tAl(nn,nn), tAkl(nn,nn), tAkl_copy(nn,nn)
real(dprec)       inv_tAk(nn,nn), inv_tAl(nn,nn), inv_tAkl(nn,nn), inv_ttAkl(nn,nn)
real(dprec)       inv_tAkltAlM(nn,nn), inv_invtAkinvtAl(nn,nn)
real(dprec)       tr_inv_tAklJij32(nn,nn)
real(dprec)       W1(nn,nn), W2(nn,nn), W3(nn,nn), W4(nn,nn), W5(nn,nn)
real(dprec)       temp1, temp2, temp3, temp4, temp5, temp6
real(dprec)       tr1, tr2, tr3, tr4, tr5
real(dprec)       det_Lk, det_Ll, det_tAkl, det_tAk, det_tAl, det_invtAkinvtAl
integer           i,j,k,kk,kkk,indx,p,q
real(dprec)       TrAJ(nn,nn),sqrtTrAJ(nn,nn),TrAJAJ(nn,nn,nn,nn),MTrAJ(nn,nn),sqrtMTrAJ(nn,nn)
real(dprec)       Mass_For_Darwin(0:nn)
real(dprec)       V2kl, MSkl

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

!Then we permute elements of Ak and Al to account for 
!the action of the permutation matrix
!  tAl=Pket'*Al*Pket
!  tAk=Pbra'*Ak*Pbra
!We also form matrix tAkl=tAk+tAl
do i=1,n
  do j=1,n
	temp1=ZERO
	temp2=ZERO
    do k=1,n
       temp1=temp1+Pket(k,j)*tAl(k,i)
       temp2=temp2+tAk(j,k)*Pbra(k,i)
	enddo
	W1(j,i)=temp1
	W2(j,i)=temp2
  enddo
enddo
!tAl=W1*Pket
!tAk=Pbra'*W2
do i=1,n  
  do j=i,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
	   temp1=temp1+W1(j,k)*Pket(k,i)
	   temp2=temp2+Pbra(k,j)*W2(k,i)
	enddo
  	tAl(j,i)=temp1
	tAl(i,j)=temp1
  	tAk(j,i)=temp2
	tAk(i,j)=temp2	
	tAkl(j,i)=temp1+temp2
	tAkl(i,j)=temp1+temp2
  enddo
enddo

!I will delete not-nesessary parts of the code
!The determinants of Lk and Ll are just
!the products of their diagonal elements
!det_Lk=ONE
!det_Ll=ONE
!do i=1,n
!  det_Lk=det_Lk*Lk(i,i)
!  det_Ll=det_Ll*Ll(i,i)
!enddo

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

!Do Cholesky factorization of tAk and tAk and then invert
if (AreMCorrFuncNeeded.or.AreMPartDensNeeded) then
det_tAk=ONE
det_tAl=ONE
do i=1,n
  do j=i,n
    temp1=tAk(i,j)
    temp2=tAl(i,j)
    do k=i-1,1,-1
      temp1=temp1-W1(i,k)*W1(j,k)
      temp2=temp2-W2(i,k)*W2(j,k)
    enddo
    if (i==j) then
      W1(i,i)=sqrt(temp1)
      det_tAk=det_tAk*temp1
      W2(i,i)=sqrt(temp2)
      det_tAl=det_tAl*temp2
    else
      W1(j,i)=temp1/W1(i,i)
      W1(i,j)=ZERO
      W2(j,i)=temp2/W2(i,i)
      W2(i,j)=ZERO
    endif
  enddo
enddo

!Inverting tAk and tAl using its Cholesky factors (stored in W1, W2)
!and placing the result into inv_tAk, inv_tAl
do i=1,n
  W1(i,i)=ONE/W1(i,i)
  W2(i,i)=ONE/W2(i,i)
  do j=i+1,n
    temp1=ZERO
    temp2=ZERO
    do k=i,j-1
      temp1=temp1-W1(j,k)*W1(k,i)
      temp2=temp2-W2(j,k)*W2(k,i)
    enddo
    W1(j,i)=temp1/W1(j,j)
    W2(j,i)=temp2/W2(j,j)
  enddo
enddo 

do i=1,n
  do j=i,n
     temp1=ZERO
     temp2=ZERO
     do k=j,n
       temp1=temp1+W1(k,i)*W1(k,j)
       temp2=temp2+W2(k,i)*W2(k,j)
     enddo
     inv_tAk(i,j)=temp1
     inv_tAl(i,j)=temp2
	 inv_tAk(j,i)=temp1
   inv_tAl(j,i)=temp2
   enddo
enddo  

!Now calculate inv_invtAkinvtAl
det_invtAkinvtAl=ONE
do i=1,n
  do j=i,n
    temp1=inv_tAk(i,j)+inv_tAl(i,j)
    do k=i-1,1,-1
      temp1=temp1-W1(i,k)*W1(j,k)
    enddo
    if (i==j) then
      W1(i,i)=sqrt(temp1)
      det_invtAkinvtAl=det_invtAkinvtAl*temp1
    else
      W1(j,i)=temp1/W1(i,i)
      W1(i,j)=ZERO
    endif
  enddo
enddo

!Inverting invtAk+invtAl
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
     inv_invtAkinvtAl(i,j)=temp1
	 inv_invtAkinvtAl(j,i)=temp1
   enddo
enddo  
endif

!Evaluating overlap
!temp2=abs(det_Ll*det_Lk)
!temp1=temp2/det_tAkl
!Skl=Glob_2raised3n2*temp1*sqrt(temp1)
!wf2originkl=Glob_2raised3n2*(temp2*sqrt(temp2))/(PI**(THREE*n/TWO))
wf2originkl=ONE
Skl=Glob_Piraised3n2/(det_tAkl*sqrt(det_tAkl))  !new line  

if(AreMCorrFuncNeeded.or.AreMPartDensNeeded) then
  temp1=1/det_tAk/det_tAl/det_invtAkinvtAl
  MSkl=Glob_Piraised3n2*temp1*sqrt(temp1)
endif

!Doing multiplication W2=inv_tAkl*tAl
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
    enddo
    W2(j,i)=temp1
  enddo
enddo

!Doing multiplication inv_tAkltAlM=inv_tAkl*tAl*M=W2*M
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+W2(j,k)*Glob_MassMatrix(k,i)
    enddo
    inv_tAkltAlM(j,i)=temp1
  enddo
enddo

!Computing kinetic energy, Tkl=tr[inv_tAkltAlM*Ak]
Tkl=ZERO
do i=1,n
  temp1=ZERO
  do k=1,n
    temp1=temp1+inv_tAkltAlM(i,k)*tAk(k,i)
  enddo
  Tkl=Tkl+temp1
enddo
Tkl=SIX*Skl*Tkl

!Evaluating potential energy, Vkl, 
!(1/r_{ij})_kl, (r_{ij})_kl, (r_{ij}^2)_kl
!and delta(r_{ij})_kl
temp1=(TWO/SQRTPI)*Skl
temp2=THREEHALF*Skl
temp3=Skl/(PI*SQRTPI)
Vkl=ZERO
do i=1,n
  TrAJ(i,i)=inv_tAkl(i,i)
  sqrtTrAJ(i,i)=sqrt(TrAJ(i,i))
  temp4=sqrtTrAJ(i,i)
  temp5=TrAJ(i,i)
  rmkl(i,i)=temp1/temp4
  rkl(i,i)=temp1*temp4
  r2kl(i,i)=temp2*temp5
  deltarkl(i,i)=temp3/(temp4*temp5)
  prvalkl(i,i)=(temp1/(temp4*temp5))*(Glob_EulerConst+log(temp5))
  Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*rmkl(i,i)
enddo
do i=1,n
  do j=i+1,n
    TrAJ(i,j)=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    TrAJ(j,i)=TrAJ(i,j)
    sqrtTrAJ(j,i)=sqrt(TrAJ(j,i))
    sqrtTrAJ(i,j)=sqrtTrAJ(j,i) 
    temp4=sqrtTrAJ(j,i)
    temp5=TrAJ(j,i)    
    rmkl(i,j)=temp1/temp4
    rmkl(j,i)=rmkl(i,j)
    rkl(i,j)=temp1*temp4
    rkl(j,i)=rkl(i,j)
    r2kl(i,j)=temp2*temp5
    r2kl(j,i)=r2kl(i,j)
    deltarkl(i,j)=temp3/(temp4*temp5)
    deltarkl(j,i)=deltarkl(i,j)
    prvalkl(i,j)=(temp1/(temp4*temp5))*(Glob_EulerConst+log(temp5))
    prvalkl(j,i)=prvalkl(i,j)
    Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*rmkl(i,j)
  enddo
enddo
Hkl=Tkl+Vkl

if (AreMCorrFuncNeeded) then
  do i=1,n
    MTrAJ(i,i)=inv_invtAkinvtAl(i,i)*4
    sqrtMTrAJ(i,i)=sqrt(MTrAJ(i,i))
  enddo
  do i=1,n
    do j=i+1,n
      MTrAJ(i,j)=(inv_invtAkinvtAl(i,i)+inv_invtAkinvtAl(j,j)-inv_invtAkinvtAl(j,i)-inv_invtAkinvtAl(j,i))*4
      MTrAJ(j,i)=MTrAJ(i,j)
      sqrtMTrAJ(j,i)=sqrt(MTrAJ(j,i))
      sqrtMTrAJ(i,j)=sqrtMTrAJ(j,i) 
    enddo
  enddo
end if

!Evaluating tr[inv_tAkl Jij inv_tAkl Jpq]
do i=1,n
  temp1=inv_tAkl(i,i)*inv_tAkl(i,i)
  TrAJAJ(i,i,i,i)=temp1
  do p=i+1,n
    do q=p+1,n
      temp2=inv_tAkl(p,i)-inv_tAkl(q,i)
      temp1=temp2*temp2
      TrAJAJ(i,i,p,q)=temp1
      TrAJAJ(i,i,q,p)=temp1
      TrAJAJ(p,q,i,i)=temp1
      TrAJAJ(q,p,i,i)=temp1
    enddo
  enddo
  do j=i+1,n
    temp1=inv_tAkl(j,i)*inv_tAkl(j,i)
    TrAJAJ(j,j,i,i)=temp1
    TrAJAJ(i,i,j,j)=temp1
    do p=i,n
      temp2=inv_tAkl(p,i)-inv_tAkl(p,j)
      temp1=temp2*temp2
      TrAJAJ(i,j,p,p)=temp1
      TrAJAJ(j,i,p,p)=temp1
      TrAJAJ(p,p,i,j)=temp1
      TrAJAJ(p,p,j,i)=temp1
      do q=p+1,n
        temp2=inv_tAkl(q,i)-inv_tAkl(p,i)-inv_tAkl(q,j)+inv_tAkl(p,j)
        temp1=temp2*temp2
        TrAJAJ(i,j,p,q)=temp1
        TrAJAJ(j,i,p,q)=temp1
        TrAJAJ(i,j,q,p)=temp1
        TrAJAJ(j,i,q,p)=temp1
        TrAJAJ(p,q,i,j)=temp1
        TrAJAJ(p,q,j,i)=temp1
        TrAJAJ(q,p,i,j)=temp1
        TrAJAJ(q,p,j,i)=temp1
      enddo
    enddo
  enddo
enddo

!This is a slow old version of the previous loop that computes TrAJAJ
!do i=1,n
!  do j=1,n
!    do p=1,n
!      do q=1,n
!        W1(1:n,1:n)=ZERO; W1(i,j)=-ONE; W1(j,i)=-ONE; W1(i,i)=ONE; W1(j,j)=ONE
!        W2(1:n,1:n)=ZERO; W2(p,q)=-ONE; W2(q,p)=-ONE; W2(p,p)=ONE; W2(q,q)=ONE;
!        W3(1:n,1:n)=matmul(inv_tAkl(1:n,1:n),W1(1:n,1:n))
!        W4(1:n,1:n)=matmul(inv_tAkl(1:n,1:n),W2(1:n,1:n))
!        W5(1:n,1:n)=matmul(W3(1:n,1:n),W4(1:n,1:n))
!        TrAJAJ(i,j,p,q)=trace(n,W5)
!      enddo
!    enddo
!  enddo
!enddo

!Evaluating (1/r_{ij}*1/r_{pq}))_kl
temp2=4*Skl/PI
temp3=2*Skl
do i=1,n
  do j=i,n
    do p=i,n
      do q=p,n
        if (((p==i).and.(q==j)).or.((p==j).and.(q==i))) then
          temp1=temp3/TrAJ(i,j)
          rmrmkl(i,j,p,q)=temp1
          rmrmkl(j,i,p,q)=temp1
          rmrmkl(i,j,q,p)=temp1
          rmrmkl(j,i,q,p)=temp1
          rmrmkl(p,q,i,j)=temp1
          rmrmkl(p,q,j,i)=temp1
          rmrmkl(q,p,i,j)=temp1
          rmrmkl(q,p,j,i)=temp1
        else
          if (TrAJAJ(i,j,p,q)>ZERO) then
            temp4=sqrt(TrAJAJ(i,j,p,q))
            temp5=temp4/(sqrtTrAJ(i,j)*sqrtTrAJ(p,q))
            temp1=temp2*asin(temp5)/temp4
          else
            temp1=temp2/(sqrtTrAJ(i,j)*sqrtTrAJ(p,q))
          endif  
          rmrmkl(i,j,p,q)=temp1
          rmrmkl(j,i,p,q)=temp1
          rmrmkl(i,j,q,p)=temp1
          rmrmkl(j,i,q,p)=temp1
          rmrmkl(p,q,i,j)=temp1
          rmrmkl(p,q,j,i)=temp1
          rmrmkl(q,p,i,j)=temp1
          rmrmkl(q,p,j,i)=temp1
        endif
      enddo  
    enddo
  enddo
enddo

!Evaluating <delta(r_{ij})delta(r_{pq})>_kl
del2kl(1:n,1:n,1:n,1:n)=ZERO
temp1=Skl/(PI**3)
do i=1,n
  do j=i,n
    do p=i,n
      do q=p,n
        if (.not.((p==i).and.(q==j))) then
          temp2=trAJ(i,j)*trAJ(p,q)-trAJAJ(i,j,p,q)
          temp3=temp1/(temp2*sqrt(temp2))
          del2kl(i,j,p,q)=temp3
          del2kl(i,j,q,p)=temp3
          del2kl(j,i,p,q)=temp3
          del2kl(j,i,q,p)=temp3
          del2kl(p,q,i,j)=temp3
          del2kl(p,q,j,i)=temp3
          del2kl(q,p,i,j)=temp3
          del2kl(q,p,j,i)=temp3
        endif
      enddo 
    enddo
  enddo
enddo


!Extracting rm2kl from rmrmkl
do i=1,n
  do j=1,n
    rm2kl(j,i)=rmrmkl(j,i,j,i)
  enddo
enddo

!Evaluating drachmanized delta(r_{ij})_kl
!Computing tAk*M*tAl and placing it in W2
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+Glob_MassMatrix(i,k)*tAl(k,j)
	enddo
    W1(i,j)=temp1
  enddo
enddo
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+tAk(i,k)*W1(k,j)
	enddo
	W2(i,j)=temp1
  enddo
enddo
!W3(i,j) will contain all matrix elements r'(tAk*M*tAl)r / rij
call ME_rXr_over_rij_all(W2,inv_tAkl,rmkl,TrAJ,W3)
!Loop that computes all drachmanized delta(r_{ij})_kl  as well as V^2_kl
V2kl=ZERO
do p=1,n
  do q=p,n
    temp1=ZERO
    do i=1,n
      temp1=temp1+ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(i))*rmrmkl(p,q,i,i)
      do j=i+1,n
        temp1=temp1+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*rmrmkl(p,q,i,j)
      enddo
    enddo
    if (p==q) then
      temp4=2*PI*Glob_MassMatrix(p,p)
      temp5=ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(p))
    else
      temp4=2*PI*(Glob_MassMatrix(p,p)+Glob_MassMatrix(q,q) &
        -Glob_MassMatrix(q,p)-Glob_MassMatrix(q,p))
      temp5=ScaledChargeProd(Glob_PseudoCharge(p),Glob_PseudoCharge(q))  
    endif
    drach_deltarkl(q,p)=(Glob_CurrEnergy*rmkl(q,p)-temp1-4*W3(q,p))/temp4
    drach_deltarkl(p,q)=drach_deltarkl(q,p)
    V2kl=V2kl+temp5*temp1
  enddo
enddo

!Evaluating Orbit-Orbit (OO) matrix element (without the factor of alpha**2)
OOkl=ZERO
!First double loop for OO
do i=1,n
  do j=1,n    
    tr1=tAl(j,i)
    tr2=tr1 
    tr3=3*tAl(j,j)
    !W1 = Al Eij Al + Al Ejj Eji Al 
    do p=1,n
      do q=1,n
        W1(q,p)=tAl(q,i)*tAl(p,j)+tAl(q,j)*tAl(p,i)
      enddo
    enddo
    !W1 = W1 + Akl Ejj Al Eij 
    do p=1,n
      W1(p,j)=W1(p,j)+tAkl(p,j)*tAl(j,i)      
    enddo
    !W1 = W1 + Eji Al Ejj Al + tr3 Eji Al 
    do p=1,n
      W1(j,p)=W1(j,p)+tAl(i,j)*tAl(p,j)+tr3*tAl(p,i)      
    enddo    
    !W2 = Akl Ejj Al
    do p=1,n
      do q=1,n
        W2(p,q)=tAkl(p,j)*tAl(q,j)
      enddo
    enddo
    !W3 = Eji Al
    W3(1:n,1:n)=ZERO
    do p=1,n
      W3(j,p)=tAl(p,i)
    enddo
    !compute integrals
    temp1=ME_rXr_over_rij(W1,j,j,inv_tAkl,rmkl(j,j),TrAJ(j,j))
    temp2=ME_rXr_rYr_over_rij(W2,W3,j,j,inv_tAkl,rmkl(j,j),TrAJ(j,j)) 
    temp3=-6*(tr1+tr2)*rmkl(j,j)+4*temp1-8*temp2
    OOkl=OOkl-temp3*ScaledChargeProd(Glob_PseudoCharge(j),Glob_PseudoCharge0)/Glob_Mass(j+1)
  enddo
enddo
OOkl=OOkl/Glob_Mass(1)

!Second double loop for OO
do i=1,n
  do j=i+1,n
    tr1=tAl(j,i)
    tr2=tr1 
    tr3=3*tAl(j,j)
    !W1 = Al Eij Al + Al Ejj Eji Al 
    do p=1,n
      do q=1,n
        W1(q,p)=tAl(q,i)*tAl(p,j)+tAl(q,j)*tAl(p,i)
      enddo
    enddo
    !W1 = W1 + Akl Ejj Al (Eij - Eii) 
    do p=1,n
      W1(p,j)=W1(p,j)+tAkl(p,j)*tAl(j,i)
      W1(p,i)=W1(p,i)-tAkl(p,j)*tAl(j,i)    
    enddo
    !W1 = W1 + (Eji - Eii) Al Ejj Al + tr3 (Eji - Eii) Al 
    do p=1,n
      temp1=tAl(i,j)*tAl(p,j)+tr3*tAl(p,i)
      W1(j,p)=W1(j,p)+temp1
      W1(i,p)=W1(i,p)-temp1    
    enddo    
    !W2 = Akl Ejj Al
    do p=1,n
      do q=1,n
        W2(p,q)=tAkl(p,j)*tAl(q,j)
      enddo
    enddo
    !W3 = (Eji - Eii) Al
    W3(1:n,1:n)=ZERO
    do p=1,n
      W3(j,p)=tAl(p,i)
      W3(i,p)=-tAl(p,i)
    enddo
    !compute integrals
    temp1=ME_rXr_over_rij(W1,i,j,inv_tAkl,rmkl(i,j),TrAJ(i,j))
    temp2=ME_rXr_rYr_over_rij(W2,W3,i,j,inv_tAkl,rmkl(i,j),TrAJ(i,j)) 
    temp3=-6*(tr1+tr2)*rmkl(i,j)+4*temp1-8*temp2
    OOkl=OOkl+ &
      temp3*ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))/(Glob_Mass(i+1)*Glob_Mass(j+1))
  enddo
enddo
OOkl=OOkl/2

!Evaluating mass-velocity matrix element
!W3=J*tAk
do p=1,n
  temp1=ZERO
  do q=1,n
    temp1=temp1+tAk(q,p)
  enddo
  do q=1,n
    W3(q,p)=temp1
  enddo
enddo
!tr1=tr[J*tAk]
tr1=ZERO
do p=1,n
  tr1=tr1+W3(p,p)
enddo
!W1=tAk*J*tAk
do p=1,n
  temp1=ZERO
  do i=1,n
    temp1=temp1+tAk(p,i)*W3(i,p)
  enddo
  W1(p,p)=temp1
  do q=p+1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+tAk(q,j)*W3(j,p)
    enddo
    W1(q,p)=temp1
    W1(p,q)=temp1
  enddo
enddo
!W3=J*tAl
do p=1,n
  temp1=ZERO
  do q=1,n
    temp1=temp1+tAl(q,p)
  enddo
  do q=1,n
    W3(q,p)=temp1
  enddo
enddo
!tr2=tr[J*tAk]
tr2=ZERO
do p=1,n
  tr2=tr2+W3(p,p)
enddo
!W2=tAl*J*tAl
do p=1,n
  temp1=ZERO
  do i=1,n
    temp1=temp1+tAl(p,i)*W3(i,p)
  enddo
  W2(p,p)=temp1
  do q=p+1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+tAl(q,j)*W3(j,p)
    enddo
    W2(q,p)=temp1
    W2(p,q)=temp1
  enddo
enddo
!W3 = tr1*W2+tr2*W1 = tr[J*tAk]*tAl*J*tAl + tr[J*tAl]*tAk*J*tAk 
do p=1,n
  W3(p,p)=tr1*W2(p,p)+tr2*W1(p,p)
  do q=p+1,n
    W3(q,p)=tr1*W2(q,p)+tr2*W1(q,p)
    W3(p,q)=W3(q,p)
  enddo
enddo
!W4=inv_tAkl*W1
do p=1,n
  do q=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_tAkl(q,j)*W1(j,p)
    enddo  
    W4(q,p)=temp1
  enddo
enddo
!tr4=tr[W4]
tr4=ZERO
do p=1,n
  tr4=tr4+W4(p,p)
enddo
!W5=inv_tAkl*W2
do p=1,n
  do q=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_tAkl(q,j)*W2(j,p)
    enddo  
    W5(q,p)=temp1
  enddo
enddo
!tr5=tr[W5]
tr5=ZERO
do p=1,n
  tr5=tr5+W5(p,p)
enddo
!temp3=tr[W4*W5]
temp3=ZERO
do p=1,n
  do q=1,n
    temp3=temp3+W4(p,q)*W5(q,p)
  enddo
enddo
!temp2=tr[inv_tAkl*W3]
temp2=ZERO
do p=1,n
  do q=1,n
    temp2=temp2+inv_tAkl(p,q)*W3(q,p)
  enddo
enddo

MVkl=36*Skl*(tr4*tr5 + (TWO/THREE)*temp3 - temp2 + tr1*tr2) &
      /(Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(1))

!sum for mass-velocity from 1 to n 
do i=1,n
  !W1=tAk*Jii*tAk
  do p=1,n
    W1(p,p)=tAk(p,i)*tAk(i,p)
    do q=p+1,n
      W1(q,p)=tAk(q,i)*tAk(i,p)
      W1(p,q)=W1(q,p)
    enddo
  enddo
  !tr1=tr[tAk*Jii]
  tr1=tAk(i,i)
  !W2=tAl*Jii*tAl
  do p=1,n
    W2(p,p)=tAl(p,i)*tAl(i,p)
    do q=p+1,n
      W2(q,p)=tAl(q,i)*tAl(i,p)
      W2(p,q)=W2(q,p)
    enddo
  enddo
  !tr2=tr[tAl*Jii]
  tr2=tAl(i,i)  
  !temp2=tr[inv_tAkl*(tr1*W2+tr2*W1)]
  temp2=ZERO
  do p=1,n
    do q=1,n
      temp2=temp2+inv_tAkl(p,q)*(tr1*W2(q,p)+tr2*W1(q,p))
    enddo
  enddo
  !W4=inv_tAkl*W1
  do p=1,n
    do q=1,n
      temp1=ZERO
      do j=1,n
        temp1=temp1+inv_tAkl(q,j)*W1(j,p)
      enddo  
      W4(q,p)=temp1
    enddo
  enddo
  !tr4=tr[W4]
  tr4=ZERO
  do p=1,n
    tr4=tr4+W4(p,p)
  enddo
  !W5=inv_tAkl*W2
  do p=1,n
    do q=1,n
      temp1=ZERO
      do j=1,n
        temp1=temp1+inv_tAkl(q,j)*W2(j,p)
      enddo  
      W5(q,p)=temp1
    enddo
  enddo
  !tr5=tr[W5]
  tr5=ZERO
  do p=1,n
    tr5=tr5+W5(p,p)
  enddo
  !temp3=tr[W4*W5]
  temp3=ZERO
  do p=1,n
    do q=1,n
      temp3=temp3+W4(p,q)*W5(q,p)
    enddo
  enddo  
  MVkl=MVkl+36*Skl*(tr4*tr5 + (TWO/THREE)*temp3 - temp2 + tr1*tr2) &
      /(Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1))
enddo

MVkl=-MVkl/8

!Evaluation of the drachmanized mass-velocity
drach_MVkl= ME_dXd_dYd(Glob_dmvM,Glob_dmvMB,inv_tAkl,tAk,tAl,Skl) &
 - V2kl - Glob_CurrEnergy*Glob_CurrEnergy*Skl + 2*Glob_CurrEnergy*Vkl &
 + Glob_CurrEnergy*ME_dXd(Glob_dmvB,inv_tAkl,tAl,Skl)
!W1(i,j) will contain matrix elements <phi_k| (1/r_{ij})(nabla_r'*(M-B)*nabla_r) |phi_l>
call ME_1_over_rij_dXd_all(Glob_dmvB,inv_tAkl,tAl,rmkl,TrAJ,W1) 
do i=1,n 
  drach_MVkl=drach_MVkl-ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(i))*W1(i,i)                       
  do j=i+1,n
    drach_MVkl=drach_MVkl-ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*W1(j,i)                        
  enddo
enddo 
drach_MVkl = drach_MVkl*Glob_dmva2 + MVkl

!Evaluation of the Darwin correction
Mass_For_Darwin(0)=Glob_Mass(1)
Mass_For_Darwin(1:n)=Glob_Mass(2:n+1)
Darwinkl=ZERO
do i=1,n
  Darwinkl=Darwinkl+(   &
     ONE/(Mass_For_Darwin(0)*Mass_For_Darwin(0)) &
    +ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
    )*ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(i))*deltarkl(i,i)
enddo
do i=1,n
  do j=1,n
    if(j/=i) then
      Darwinkl=Darwinkl+   &
        ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
       *ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*deltarkl(i,j)
    endif   
  enddo  
enddo
Darwinkl=-Darwinkl*PI/2
!Evaluation of the drachmanized Darwin correction
drach_Darwinkl=ZERO
do i=1,n
  drach_Darwinkl=drach_Darwinkl+(   &
     ONE/(Mass_For_Darwin(0)*Mass_For_Darwin(0)) &
    +ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
    )*ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(i))*drach_deltarkl(i,i)
enddo
do i=1,n
  do j=1,n
    if(j/=i) then
      drach_Darwinkl=drach_Darwinkl+   &
        ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
       *ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*drach_deltarkl(i,j)
    endif   
  enddo  
enddo
drach_Darwinkl=-drach_Darwinkl*PI/2
 

!Evaluation of correlation functions
if (AreCorrFuncNeeded) then
  temp1=Skl/(PI*SQRTPI)
  p=0
  do i=1,n
    do j=i,n
      p=p+1
      temp3=temp1/(sqrtTrAJ(j,i)*TrAJ(j,i))
      do k=1,NumCFGridPoints
        temp2=CFGrid(k)*CFGrid(k)
        CFkl(p,k)=temp3*exp(-temp2/TrAJ(j,i))
        !CFkl(p,k)=temp2*temp3*exp(-temp2/TrAJ(j,i))  !Multiplied by \xi^2
      enddo
    enddo
  enddo
endif

if (ArePartDensNeeded) then
  temp1=Skl/(PI*SQRTPI)
  do i=1,n+1
    temp3=ZERO
    do p=1,n
      temp3=temp3+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tAkl(p,p)
      do q=p+1,n
        temp3=temp3+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tAkl(q,p)
      enddo
    enddo
    temp4=temp1/(sqrt(temp3)*temp3)
    do k=1,NumDensGridPoints
      temp2=DensGrid(k)*DensGrid(k)
      Denskl(i,k)=temp4*exp(-temp2/temp3)
      !Denskl(i,k)=temp2*temp4*exp(-temp2/temp3) !Multiplied by \xi^2
    enddo
  enddo
endif

if (AreMCorrFuncNeeded) then
  temp1=MSkl/(PI*SQRTPI)
  p=0
  do i=1,n
    do j=i,n
      p=p+1
      temp3=temp1/(sqrtMTrAJ(j,i)*MTrAJ(j,i))
      do k=1,NumCFGridPoints
        temp2=CFGrid(k)*CFGrid(k)
        CFkl(p,k)=temp3*exp(-temp2/MTrAJ(j,i))
        !CFkl(p,k)=temp2*temp3*exp(-temp2/TrAJ(j,i))  !Multiplied by \xi^2
      enddo
    enddo
  enddo
end if

if (AreMPartDensNeeded) then
  temp1=MSkl/(PI*SQRTPI)
  do i=1,n+1
    temp3=ZERO
    do p=1,n
      temp3=temp3+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(p,p)
      do q=p+1,n
        temp3=temp3+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(q,p)
      enddo
    enddo
    temp3=temp3*4
    temp4=temp1/(sqrt(temp3)*temp3)
    do k=1,NumDensGridPoints
      temp2=DensGrid(k)*DensGrid(k)
      Denskl(i,k)=temp4*exp(-temp2/temp3)
      !Denskl(i,k)=temp2*temp4*exp(-temp2/temp3) !Multiplied by \xi^2
    enddo
  enddo
end if
end subroutine MatrixElementsForExpcVals


function trace(k,M)
real(dprec) trace
integer k
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
real(dprec) M(nn,nn)
integer i
trace=ZERO
do i=1,k
  trace=trace+M(i,i)
enddo
end function trace


function ScaledChargeProd(q1,q2)
real(dprec) ScaledChargeProd,q1,q2,x
x=q1*q2
if (x<0.0_dprec) then
  ScaledChargeProd=x*Glob_AttractionScalingParam
else
  if ((q1>0.0_dprec).and.(q2>0.0_dprec)) then
    ScaledChargeProd=x*Glob_RepulsionScalingParam*Glob_RepulsionScalingParamPlus
  else 
    ScaledChargeProd=x*Glob_RepulsionScalingParam*Glob_RepulsionScalingParamMinus
  endif 
endif
end function ScaledChargeProd


function ME_rXr_over_rij(X,i,j,inv_tAkl,ME_1_over_rij,TrAJ)
!Function ME_rXr_over_rij computes the following
!matrix element with real L=0 Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>
!Here X is an arbitrary (i.e. nonsymmetric) real matrix.
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n real matrix
!           i,j :: indices denoting i and j.
!        inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element 
!          TrAJ :: the value of Tr[inv_tAkl Jij]
!Note that n=Glob_n and nn=Glob_MaxAllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced. 
real(dprec) ME_rXr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles

!Arguments:
real(dprec)  X(nn,nn),inv_tAkl(nn,nn),ME_1_over_rij,TrAJ
integer      i,j
!Local variables
real(dprec)  AXi(nn),AXj(nn)
real(dprec)  TrAX,TrAXAJ
integer      n,m,p,q

n=Glob_n

TrAX=ZERO
do m=1,n
  do p=1,n
    TrAX=TrAX+inv_tAkl(m,p)*X(p,m)
  enddo
enddo

do m=1,n
  AXi(m)=ZERO
  do p=1,n
    AXi(m)=AXi(m)+inv_tAkl(i,p)*X(p,m)
  enddo
enddo
if (j/=i) then
  do m=1,n
    AXj(m)=ZERO
    do p=1,n
      AXj(m)=AXj(m)+inv_tAkl(j,p)*X(p,m)
    enddo
  enddo
endif

if (i==j) then
  TrAXAJ=ZERO
  do m=1,n
    TrAXAJ=TrAXAJ+AXi(m)*inv_tAkl(m,i)
  enddo
else
  TrAXAJ=ZERO
  do m=1,n
    TrAXAJ=TrAXAJ+(AXi(m)-AXj(m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
  enddo  
endif

ME_rXr_over_rij=ME_1_over_rij*(3*TrAX-TrAXAJ/TrAJ)/2

end function ME_rXr_over_rij


subroutine ME_rXr_over_rij_all(X,inv_tAkl,rmkl,TrAJ,ME)
!Subroutine ME_rXr_over_rij_all computes the following
!matrix elements with real L=0 Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>    for all combinations i,j=1..n at the same time
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
!          rmkl :: the values of <phi_k| 1/r_{ij} |phi_l> matrix element 
!          TrAJ :: the values of Tr[inv_tAkl Jij]
!Output:
!            ME :: n x n real matrix where all computed matric elements are returned    
!Note that n=Glob_n and nn=Glob_MaxAllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced. 
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles

!Arguments:
real(dprec)  X(nn,nn),inv_tAkl(nn,nn),rmkl(Glob_n,Glob_n),TrAJ(nn,nn),ME(nn,nn)
!Local variables
integer      i,j
real(dprec)  AXi(nn),AXj(nn)
real(dprec)  TrAX,TrAXAJ
integer      n,m,p,q

n=Glob_n

TrAX=ZERO
do m=1,n
  do p=1,n
    TrAX=TrAX+inv_tAkl(m,p)*X(p,m)
  enddo
enddo

do i=1,n
  do j=i,n    
    do m=1,n
      AXi(m)=ZERO
      do p=1,n
        AXi(m)=AXi(m)+inv_tAkl(i,p)*X(p,m)
      enddo
    enddo
    if (j/=i) then
      do m=1,n
        AXj(m)=ZERO
        do p=1,n
          AXj(m)=AXj(m)+inv_tAkl(j,p)*X(p,m)
        enddo
      enddo
    endif

    if (i==j) then
      TrAXAJ=ZERO
      do m=1,n
        TrAXAJ=TrAXAJ+AXi(m)*inv_tAkl(m,i)
      enddo
    else
      TrAXAJ=ZERO
      do m=1,n
        TrAXAJ=TrAXAJ+(AXi(m)-AXj(m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
      enddo  
    endif

    ME(i,j)=rmkl(i,j)*(3*TrAX-TrAXAJ/TrAJ(i,j))/2
    ME(j,i)=ME(i,j)
    
  enddo
enddo

end subroutine ME_rXr_over_rij_all


function ME_rXr_rYr_over_rij(X,Y,i,j,inv_tAkl,ME_1_over_rij,TrAJ)
!Function ME_rXr_rYr_over_rij computes the following
!matrix element with real L=0 Gaussians phi_k and phi_l:
!<phi_k| (r'Xr)(r'Yr)/r_{ij} |phi_l>
!Here X and Y are arbitrary (i.e. nonsymmetric) real matrices.
!Index i can be equal to j. In the latter case
!<phi_k| (r'Xr)(r'Yr)/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!           i,j :: indices denoting i and j.
!        inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element 
!          TrAJ :: the value of Tr[inv_tAkl Jij]
!Note that n=Glob_n and nn=Glob_MaxAllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced. 
real(dprec) ME_rXr_rYr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles

!Arguments:
real(dprec)  X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),ME_1_over_rij,TrAJ
integer      i,j
!Local variables
real(dprec)  Ys(nn,nn),Xs(nn,nn),AX(nn,nn),AY(nn,nn),AXAYi(nn),AXAYj(nn)
real(dprec)  TrAX,TrAY,TrAXAY,TrAXAJ,TrAYAJ,TrAXAYAJ
integer      n,m,p,q

n=Glob_n

do p=1,n
  Ys(p,p)=Y(p,p)
  do q=p+1,n
    Ys(p,q)=(Y(p,q)+Y(q,p))/2
    Ys(q,p)=Ys(p,q)
  enddo
enddo
do p=1,n
  Xs(p,p)=X(p,p)
  do q=p+1,n
    Xs(p,q)=(X(p,q)+X(q,p))/2
    Xs(q,p)=Xs(p,q)
  enddo
enddo

do p=1,n
  do q=1,n
    AY(p,q)=ZERO
    do m=1,n
      AY(p,q)=AY(p,q)+inv_tAkl(p,m)*Ys(m,q)
    enddo
  enddo
enddo
do p=1,n
  do q=1,n
    AX(p,q)=ZERO
    do m=1,n
      AX(p,q)=AX(p,q)+inv_tAkl(p,m)*Xs(m,q)
    enddo
  enddo
enddo

do m=1,n
  AXAYi(m)=ZERO
  do p=1,n
    AXAYi(m)=AXAYi(m)+AX(i,p)*AY(p,m)
  enddo
enddo
if (j/=i) then
  do m=1,n
    AXAYj(m)=ZERO
    do p=1,n
      AXAYj(m)=AXAYj(m)+AX(j,p)*AY(p,m)
    enddo
  enddo
endif

TrAY=ZERO
TrAX=ZERO
do m=1,n
  TrAY=TrAY+AY(m,m)
  TrAX=TrAX+AX(m,m)  
enddo
TrAXAY=ZERO
do m=1,n
  do p=1,n
    TrAXAY=TrAXAY+AX(m,p)*AY(p,m)
  enddo
enddo

if (i==j) then
  TrAXAJ=ZERO
  do m=1,n
    TrAXAJ=TrAXAJ+AX(i,m)*inv_tAkl(m,i)
  enddo
  TrAYAJ=ZERO
  do m=1,n
    TrAYAJ=TrAYAJ+AY(i,m)*inv_tAkl(m,i)
  enddo
  TrAXAYAJ=ZERO
  do m=1,n
    TrAXAYAJ=TrAXAYAJ+AXAYi(m)*inv_tAkl(m,i)
  enddo  
else
  TrAXAJ=ZERO
  do m=1,n
    TrAXAJ=TrAXAJ+(AX(i,m)-AX(j,m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
  enddo  
  TrAYAJ=ZERO
  do m=1,n
    TrAYAJ=TrAYAJ+(AY(i,m)-AY(j,m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
  enddo   
  TrAXAYAJ=ZERO
  do m=1,n
    TrAXAYAJ=TrAXAYAJ+(AXAYi(m)-AXAYj(m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
  enddo  
endif

ME_rXr_rYr_over_rij=ME_1_over_rij*(   &
     9*TrAY*TrAX + 6*TrAXAY           & 
     - ( 3*(TrAY*TrAXAJ + TrAX*TrAYAJ) + 4*TrAXAYAJ )/TrAJ   &
     + 3*TrAXAJ*TrAYAJ/(TrAJ*TrAJ)    &
)/4 

end function ME_rXr_rYr_over_rij


function ME_dXd_dYd(X,Y,inv_tAkl,tAk,tAl,Skl)
real(dprec) ME_dXd_dYd
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
real(dprec) X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),tAk(nn,nn),tAl(nn,nn)
real(dprec) Skl,M1(nn,nn),M2(nn,nn),Z1(nn,nn),Z2(nn,nn)
integer p,q,k,n
real(dprec) tr1,tr2,tr3,tr4,tr5,t1,t2
n=Glob_n
!M1(1:n,1:n)=matmul(X(1:n,1:n),tAk(1:n,1:n))
!tr1=trace(M1)
!M2(1:n,1:n)=matmul(Y(1:n,1:n),tAl(1:n,1:n))
!tr2=trace(M2)
tr1=ZERO
tr2=ZERO
do p=1,n
  do q=1,n
    t1=ZERO
    t2=ZERO
    do k=1,n
      t1=t1+X(q,k)*tAk(k,p)   
      t2=t2+Y(q,k)*tAl(k,p)
    enddo    
    M1(q,p)=t1
    M2(q,p)=t2
  enddo 
  tr1=tr1+M1(p,p)
  tr2=tr2+M2(p,p)
enddo 
!Z1=tAk*M1
!Z2=tAl*M2
do p=1,n
  do q=1,n
    t1=ZERO
    t2=ZERO
    do k=1,n
      t1=t1+tAk(q,k)*M1(k,p)  
      t2=t2+tAl(q,k)*M2(k,p)
    enddo    
    Z1(q,p)=t1
    Z2(q,p)=t2
  enddo 
enddo 
!M1(1:n,1:n)=matmul(inv_tAkl(1:n,1:n),Z1(1:n,1:n))
!tr3=trace(M1)
!M2(1:n,1:n)=matmul(inv_tAkl(1:n,1:n),Z2(1:n,1:n))
!tr4=trace(M2)
tr3=ZERO
tr4=ZERO
do p=1,n
  do q=1,n
    t1=ZERO
    t2=ZERO
    do k=1,n
      t1=t1+inv_tAkl(q,k)*Z1(k,p)   
      t2=t2+inv_tAkl(q,k)*Z2(k,p)
    enddo    
    M1(q,p)=t1
    M2(q,p)=t2
  enddo 
  tr3=tr3+M1(p,p)
  tr4=tr4+M2(p,p)
enddo
!tr5=trace(M1*M2)
tr5=ZERO
do p=1,n
  do k=1,n
    tr5=tr5+M1(p,k)*M2(k,p)   
  enddo    
enddo
ME_dXd_dYd=(24*tr5+36*(tr3*tr4-tr1*tr4-tr2*tr3+tr1*tr2))*Skl
end function ME_dXd_dYd


function ME_dXd(X,inv_tAkl,tAl,Skl)
real(dprec) ME_dXd
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
real(dprec) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn),Skl,M(nn,nn),Z(nn,nn)
integer p,q,k,n
real(dprec) tr1,tr2,t
n=Glob_n
!M(1:n,1:n)=matmul(tAl(1:n,1:n),X(1:n,1:n))
!tr1=trace(M)
tr1=ZERO
do p=1,n
  do q=1,n
    t=ZERO  
    do k=1,n
      t=t+tAl(q,k)*X(k,p)    
    enddo    
    M(q,p)=t
  enddo 
  tr1=tr1+M(p,p)
enddo   
!Z(1:n,1:n)=matmul(M(1:n,1:n),tAl(1:n,1:n))
do p=1,n
  do q=1,n
    t=ZERO  
    do k=1,n
      t=t+M(q,k)*tAl(k,p)    
    enddo    
    Z(q,p)=t
  enddo    
enddo
!M=matmul(inv_tAkl(1:n,1:n),Z(1:n,1:n))
!tr2=trace(M)
tr2=ZERO
do p=1,n
  do q=1,n
    t=ZERO  
    do k=1,n
      t=t+inv_tAkl(q,k)*Z(k,p)    
    enddo    
    M(q,p)=t
  enddo 
  tr2=tr1+M(p,p)
enddo   
!Evaluating matrix elements
ME_dXd=6*(tr2-tr1)*Skl
end function ME_dXd


function ME_1_over_rij_dXd(X,i,j,inv_tAkl,tAl,ME_1_over_rij,TrAJ)
!Function ME_1_over_rij_dXd computes the following
!matrix element with real L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!Here X and Y are arbitrary (i.e. nonsymmetric) real matrices.
!Index i can be equal to j. In the latter case
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>is computed
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!           i,j :: indices denoting i and j.
!      inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element 
!          TrAJ :: the value of Tr[inv_tAkl Jij]
!Note that n=Glob_n and nn=Glob_MaxAllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced. 
!Input parameters:
real(dprec) ME_1_over_rij_dXd
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
real(dprec) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn),ME_1_over_rij,TrAJ
!Local variables:
integer i,j,p,q,k,n
real(dprec) M(nn,nn),Z(nn,nn)
real(dprec) tr1,t

n=Glob_n
!Compute Z=tAl*X, tr1=trace[Z], and M=tAl*X*tAl
tr1=ZERO
do p=1,n
  do q=1,n
    t=ZERO  
    do k=1,n
      t=t+tAl(q,k)*X(k,p)    
    enddo    
    Z(q,p)=t
  enddo 
  tr1=tr1+Z(p,p)
enddo    
do p=1,n
  do q=1,n
    t=ZERO  
    do k=1,n
      t=t+Z(q,k)*tAl(k,p)    
    enddo    
    M(q,p)=t
  enddo    
enddo
!Compute the matrix element
ME_1_over_rij_dXd=4*ME_rXr_over_rij(M,i,j,inv_tAkl,ME_1_over_rij,TrAJ)-6*tr1*ME_1_over_rij
end function ME_1_over_rij_dXd


subroutine ME_1_over_rij_dXd_all(X,inv_tAkl,tAl,rmkl,TrAJ,ME)
!Subroutine ME_1_over_rij_dXd computes 
!matrix elements with real L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!for all combinations of i and j indexes.    
!Here X and Y are arbitrary (i.e. nonsymmetric) real matrices.
!When i is equal to j it stands for the matrix elements
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!           i,j :: indices denoting i and j.
!      inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
!          rmkl :: the values of <phi_k| 1/r_{ij} |phi_l> matrix element 
!          TrAJ :: the values of Tr[inv_tAkl Jij]
!Output:
!            ME :: n x n real matrix where all computed matric elements are returned 
!Note that n=Glob_n and nn=Glob_MaxAllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced. 
!Input parameters:
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
real(dprec) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn),rmkl(Glob_n,Glob_n),TrAJ(nn,nn)
real(dprec) ME(nn,nn)
!Local variables:
integer i,j,p,q,k,n
real(dprec) M(nn,nn),Z(nn,nn)
real(dprec) tr1,t

n=Glob_n
!Compute Z=tAl*X, tr1=trace[Z], and M=tAl*X*tAl
tr1=ZERO
do p=1,n
  do q=1,n
    t=ZERO  
    do k=1,n
      t=t+tAl(q,k)*X(k,p)    
    enddo    
    Z(q,p)=t
  enddo 
  tr1=tr1+Z(p,p)
enddo    
do p=1,n
  do q=1,n
    t=ZERO  
    do k=1,n
      t=t+Z(q,k)*tAl(k,p)    
    enddo    
    M(q,p)=t
  enddo    
enddo
!Compute all matrix elements <phi_k| r'Mr/r_{ij} |phi_l> 
call ME_rXr_over_rij_all(M,inv_tAkl,rmkl,TrAJ,Z)
t=6*tr1
do i=1,n
  do j=i,n
    ME(j,i)=4*Z(j,i)-t*rmkl(j,i)
    ME(i,j)=ME(j,i)
  enddo    
enddo    

end subroutine ME_1_over_rij_dXd_all

!!Old, slow, and simple version
!function ME_1_over_rij_dXd(X,i,j,inv_tAkl,tAl,ME_1_over_rij,TrCJ)
!real(dprec) ME_1_over_rij_dXd
!integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!real(dprec) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn),ME_1_over_rij,TrCJ
!integer i,j,n
!real(dprec) M(nn,nn),Z(nn,nn)
!real(dprec) tr1
!n=Glob_n
!M(1:n,1:n)=matmul(tAl(1:n,1:n),matmul(X(1:n,1:n),tAl(1:n,1:n)))
!Z(1:n,1:n)=matmul(tAl(1:n,1:n),X(1:n,1:n))
!tr1=trace(n,Z)
!ME_1_over_rij_dXd=4*ME_rXr_over_rij(M,i,j,inv_tAkl,ME_1_over_rij,TrCJ)-6*tr1*ME_1_over_rij
!end function ME_1_over_rij_dXd


subroutine spinPreCalc(n, nFactorial, parityFactor, SSFmassChargeCoefficient, &
  AnihMassChargeCoefficient, ketMatrix, spatialYoung, &
  positronPosition, numberOfSpinFunctions, spinFreeME, SiSjME)
  use spinStuff
  implicit none

  character(len = maxLen), intent(in) :: spatialYoung
  integer, intent(in) :: n, nFactorial

  real(dprec), dimension(nFactorial), intent(out) :: parityFactor
  real(dprec), dimension(n, n, nFactorial), intent(out) :: ketMatrix
  real(dprec), dimension(n, n), intent(out) :: SSFmassChargeCoefficient, AnihMassChargeCoefficient
  integer, intent(out) :: positronPosition, numberOfSpinFunctions

  real(dprec), dimension(nFactorial), intent(out) :: spinFreeME
  real(kind = dprec), dimension(n, n, 2, nFactorial), intent(out) :: SiSjME

  ! local variables
  integer :: i, j, k, l, m
  character(len = maxLen) :: mySpatialYoung
  integer, dimension(nFactorial) :: parities
  integer, dimension(n, n, nFactorial) :: allPermutations

  SSFmassChargeCoefficient = ZERO
  do i = 1, n
    do j = 1, n
      SSFmassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
      (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * EIGHT * PI / THREE
    enddo
  enddo

  AnihMassChargeCoefficient = ZERO
  do i = 1, n
    do j = 1, n
      AnihMassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
      (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * TWO * PI
    enddo
  enddo


  ! now we deal with the spin stuff

  ! rename the particles
  mySpatialYoung = spatialYoung
  do i = 1, maxLen
    if (mySpatialYoung(i:i) == 'P') then
      read(mySpatialYoung(i + 1:i + 1), *) k
      read(mySpatialYoung(i + 2:i + 2), *) j
      write(mySpatialYoung(i + 1:i + 1), '(i1)') k - 1
      write(mySpatialYoung(i + 2:i + 2), '(i1)') j - 1
    endif
  enddo

  call getSpinOperatorsMeanValues(n, nFactorial, mySpatialYoung, positronPosition, numberOfSpinFunctions, &
  allPermutations, parities, spinFreeME, SiSjMe)

  ketMatrix = ZERO
  do i = 1, nFactorial

    do k = 1, n
      do l = 1, n

        ketMatrix(k, l, i) = real(allPermutations(l, k, i))
        ! note the transposition here

      enddo
    enddo

  enddo

  do i = 1, nFactorial
    parityFactor(i) = real(parities(i))
  enddo

end subroutine


subroutine overlapMatrixElements(vechLk, P, Skk)
!This subroutine computes symmetry adapted matrix element with 
!two real L=0 correlated Gaussians:
! 
!fk =  exp[-r'(Lk*Lk')r] 
!
!Symmetry adaption is applied to the ket using 
!permutation matrices Glob_YHYMatr(:,:,1:Glob_NumYHYTerms)
!
!Input:     
!   vechLk :: Array of length (n(n+1)/2) of exponential parameters. 
!   P   :: The symmetry permutation matrix of size n x n
!Output:
!   Skk	 ::	Overlap matrix element (normalized)


!Arguments
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
real(dprec)       Lk(nn,nn), Ll(nn,nn)
real(dprec)       Ak(nn,nn), tAl(nn,nn), tAkl(nn,nn)
real(dprec)       inv_tAkl(nn,nn)
real(dprec)       W1(nn,nn)
real(dprec)       temp1
real(dprec)       det_tAkl
integer           i, j, k, indx

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

!Evaluating overlap

!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2raised3n2*temp1*sqrt(temp1)
Skk=Glob_Piraised3n2/(det_tAkl*sqrt(det_tAkl))  !new line


end subroutine overlapMatrixElements


end module matelem
