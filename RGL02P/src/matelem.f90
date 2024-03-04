module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=0 and L=1 Gaussians.
use globvars
implicit none

contains

subroutine OverLapElementS0(vechLk, P, Skk)

!
! This subroutine computes symmetry adapted matrix element <P*fk|P*fk>
! with two real L=0 correlated Gaussians

! fk =  exp[-r'(Lk*Lk')r]
!
! Symmetry adaption is applied to the bra using 
! permutation matrix Pk and to the ket using permutation matrix Pl.
!
! Input :     
!		vechLl		::	Arrays of length (n(n+1)/2) of exponential parameters.
!		P	  	::	The symmetry permutation matrix of size n x n
!
! Output:
!		Skk			::	< P*fk | P*fk >
!
!                      <P*fk|Y*fk>               abs(det_Lk)^3
!   < P*fk | P*fk > = --------------= 2^(3*n/2) ---------------- 
!                        <fk|fk>                  det_tAkk^3/2


!Arguments
real(dprec),intent(in)   :: vechLk(Glob_np)
real(dprec),intent(in)   :: P(Glob_n,Glob_n)
real(dprec),intent(out)  :: Skk


!Parameters (These are needed to declare static arrays. Using static 
!arrays makes the function call a little faster in comparison with 
!the case when arrays are dynamically allocated in stack)
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles

!Local variables
real(dprec),allocatable,dimension(:)         :: vechLl

integer           n, np
real(dprec)       Lk(nn,nn), Ll(nn,nn), inv_Lk(nn,nn), inv_Ll(nn,nn)
real(dprec)       tAk(nn,nn), tAl(nn,nn), tAkl(nn,nn), Ak(nn,nn)
real(dprec)       inv_tAkl(nn,nn)
real(dprec)       W1(nn,nn), W2(nn,nn)
real(dprec)       temp1, temp2
real(dprec)       det_Lk, det_Ll, det_tAkl
integer           i,j,k,indx,q


n=Glob_n
np=Glob_np

allocate(vechLl(n))
vechLl=vechLk


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
det_Lk=ONE
det_Ll=ONE
do i=1,n
  det_Lk=det_Lk*Lk(i,i)
  det_Ll=det_Ll*Ll(i,i)
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

!Evaluating overlap
temp2=abs(det_Ll*det_Lk)
temp1=temp2/det_tAkl
Skk=Glob_2raised3n2*temp1*sqrt(temp1)
!Skk=Glob_Piraised3n2/(det_tAkl*sqrt(det_tAkl))  !new line

end subroutine OverLapElementS0



subroutine OverLapElementS1(m_l,vechLl, P, Sll)

! this subroutine calculates <P*fk|P*fk>
!
! two real L=1 correlated Gaussians
! fl = z{m_l} * exp[-r'(Lk*Lk')r]
!
!                  <P*fl|P*fl>               abs(det_Ll)^3      v_l'*inv_tAll*tv_l 
! <psi_l|psi_l> = --------------= 2^(3*n/2) ---------------- * -------------------- 
!                    <fl|fl>                  det_tAll^3/2      v_l'*inv_All*v_l
!
!where m_l is some integer between 1 and n (n is the number of pseudoparticles). 
!
!
!Input :     
!	m_l			::	integers that determine which z-component is in the
!					premultiplier of the Gaussian
!	vechLl		::	Arrays of length (n(n+1)/2) of exponential parameters.
!   P           ::  The symmetry permutation matrix of size n x n 
!Output:
!	Sll			::	< psi_l | psi_l >


!Arguments
integer,intent(in)       :: m_l
real(dprec),intent(in)   :: vechLl(Glob_np)
real(dprec),intent(in)   :: P(Glob_n,Glob_n)
real(dprec),intent(out)  :: Sll


!Parameters (These are needed to declare static arrays. Using static 
!arrays makes the function call a little faster in comparison with 
!the case when arrays are dynamically allocated in stack)
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
real(dprec),allocatable,dimension(:)         :: vechLk

integer           n,np,m_k
integer           tvk(nn),tvl(nn)
real(dprec)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
real(dprec)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn),Ak(nn,nn)
real(dprec)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn), inv_tAkltAl(nn,nn)
real(dprec)       W1(nn,nn),W2(nn,nn)
real(dprec)       inv_tAkltvl(nn),vkinv_tAkl(nn)
real(dprec)       temp1,temp2
real(dprec)       det_Lk, det_Ll, det_tAkl,tau3
integer           i,j,k,indx


n=Glob_n
np=Glob_np

allocate(vechLk(n))
vechLk=vechLl
m_k=m_l

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
det_Lk=ONE
det_Ll=ONE
do i=1,n
  det_Lk=det_Lk*Lk(i,i)
  det_Ll=det_Ll*Ll(i,i)
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

!Finding the inverse of Akk and All using their Cholesky factors
!The result is placed in inv_Akk and inv_All
do i=1,n
  W1(i,i)=ONE/Lk(i,i)
  W2(i,i)=ONE/Ll(i,i)  
  do j=i+1,n
    temp1=ZERO
    temp2=ZERO
    do k=i,j-1
      temp1=temp1-Lk(j,k)*W1(k,i)
      temp2=temp2-Ll(j,k)*W2(k,i)
    enddo
    W1(j,i)=temp1/Lk(j,j)
    W2(j,i)=temp2/Ll(j,j)
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
     inv_Akk(i,j)=ONEHALF*temp1
	 inv_Akk(j,i)=ONEHALF*temp1
     inv_All(i,j)=ONEHALF*temp2
	 inv_All(j,i)=ONEHALF*temp2	 
   enddo
enddo  

!Computing tvl=P'*vl 
do i=1,n
  tvl(i)=P(m_l,i)
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
temp1=abs(det_Ll*det_Lk)/det_tAkl
Sll=Glob_2raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))



end subroutine OverLapElementS1

subroutine OverlapMatrixElementsLS(vechLk, P, Skk)

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
!Output:
!   Skk	 ::	Overlap matrix element 



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
	real(dprec)       Lk(nn,nn), Ll(nn,nn), inv_Lk(nn,nn), inv_Ll(nn,nn)
	real(dprec)       Ak(nn,nn), tAl(nn,nn), tAkl(nn,nn)
	real(dprec)       inv_tAkl(nn,nn)
	real(dprec)       W1(nn,nn)
	real(dprec)       temp1
	real(dprec)       det_Lk, det_Ll, det_tAkl
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
		Ll(j,i)=vechLk(indx)
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
	
	!Evaluating overlap
	Skk=Glob_Piraised3n2/(det_tAkl*sqrt(det_tAkl))  !new line

end subroutine OverlapMatrixElementsLS



subroutine overlapMatrixElementsLP(m_k, mm_k, vechLk, P, Skk)
	!This subroutine computes symmetry adapted matrix element with
	!two real L=2 correlated Gaussians:
	!
	!fk = (v'_k * r) (w'_k * r) exp[-r'(Lk*Lk')r]
	!
	!m_k and mm_k are integers between 1 and n (n is the number of
	!pseudoparticles). Symmetry adaption is applied to the ket using
	!permutation matrix P
	!
	!Input:
	!   m_k, mm_k :: integers that determine which psuedoparticles carry l=1 momentum
	!   vechLk :: Array of length (n(n+1)/2) of
	!     exponential parameters.
	!   P  :: The symmetry permutation matrix of size n x n
	!Output:
	!   Skk	 ::	Overlap matrix element (normalized)
	
	!Arguments
	integer,intent(in)          :: m_k, mm_k
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
	real(dprec)       vl(nn), bl(nn), vkinv_tAkl(nn), bkinv_tAkl(nn)
	real(dprec)       Lk(nn,nn),Ll(nn,nn)
	real(dprec)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
	real(dprec)       inv_tAkl(nn,nn)
	real(dprec)       W1(nn,nn)
	real(dprec)       temp1, temp2
	real(dprec)       det_tAkl
	real(dprec)       tau3, tau33, tau333, tau334, m, m1, m3
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
	  Ll(i,j)=ZERO
	  Ll(j,i)=vechLk(indx)
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
	!temp1=ZERO
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
	
	!Computing vl=P'*vl, bl=P'*bl 
	do i=1,n
	  vl(i)=P(m_k,i)
	  bl(i)=P(mm_k,i)
	enddo
	
	!Compute vkinv_tAkl=vk'*inv_tAkl, bkinv_tAkl=bk'*inv_tAkl
	do i=1,n
	  vkinv_tAkl(i)=inv_tAkl(m_k,i)
	  bkinv_tAkl(i)=inv_tAkl(mm_k,i)
	enddo
	
	!Compute tau3=vkinv_tAkl*vl, tau33=bkinv_tAkl*bl
	tau3=ZERO
	tau33=ZERO
	tau333=ZERO
	tau334=ZERO
	do i=1,n
	  tau3=tau3+vkinv_tAkl(i)*vl(i)
	  tau33=tau33+bkinv_tAkl(i)*bl(i)
	  tau333=tau333+vkinv_tAkl(i)*bl(i)
	  tau334=tau334+bkinv_tAkl(i)*vl(i)
	enddo
	m1=tau3*tau33
	m3=tau333*tau334
	m=m1-m3  !look formula 40 in document
	
	!Evaluating overlap
	!temp1=ZERO
	temp1=FOUR*det_tAkl*sqrt(det_tAkl)
	Skk=Glob_Piraised3n2*m/temp1
	
	end subroutine overlapMatrixElementsLP

subroutine MatrixElemenTranDipoleMoment(ml, vechLk, vechLl, Pk, Pl, TranDipolLength_kl,TranDipolVelocity_kl)

!This subroutine computes symmetry adapted matrix element with
!a real L=0 and a real L=1 correlated Gaussians:
!
!fk = exp[-r'(Lk*Lk')r]  
!fl = z_{ml} exp[-r'(Ll*Ll')r] 
!
!where m_l is some integer between 1 and n (n is the number of 
!pseudoparticles). Symmetry adaption is applied to the bra using 
!permutation matrix Pk and to the ket using permutation matrix Pl.

!Pk fk = exp[-r' {Pk'*(Lk*Lk')*Pk} r]  
!Pl fl = (Pl z_{m_l}) exp[-r' {Pl'*(Ll*Ll')*Pl} r] 

! this subroutine computes the following expressions
! which is a part of transition dipole momentum calculation.
!
!
!
!====================================================================
!				Transition dipole integral in Lenght form
!====================================================================
!
!
!                                            m_i         < Pk*fk0 | z_i | Pl*fl1 >
! TranDipolLength_kl = SUM{i}(q_i - Q_tot * -----)*-------------------------------------------
!                                            m0      Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )
!
!
!      < Pk*fk0 | z_i | Pl*fl1 >               2^(3*n/2)     (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5        vi'*inv_tAkl*vl
! ----------------------------------------- = ----------- * ---------------------------------------- * ----------------------- 
!   Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )       sqrt(2)                (det_tAkl)^1.5                   Sqrt(vl'*inv_All*vl)   
!           
!
!
!
!
!====================================================================
!				Transition dipole integral in Lenght form
!====================================================================
!
!
!                                    < Pk*fk0 | P(z_i) | Pl*fl1 >
! TranDipolLength_kl = SUM{i}*-------------------------------------------
!                                   Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )
!
!
!     < Pk*fk0 | P(z_i) | Pl*fl1 >             2^(3*n/2)     (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5       vi'*tAk*inv_tAkl*vl
! ----------------------------------------- = ----------- * ---------------------------------------- * ----------------------- 
!   Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )       sqrt(2)                (det_tAkl)^1.5                    Sqrt(vl'*inv_All*vl)   
!           
!




!Input:     
!	m_l				::	integer that determine which z-component is in the
!						premultiplier of the Gaussian
!	vechLk, vechLl	::	Arrays of length (n(n+1)/2) of exponential parameters. 
!	Pk, Pl			::	The symmetry permutation matrices of size n x n
!
!Output:
!	TDkl  			::	Matrix element (normalized)


!Arguments
integer, intent(in)          :: ml
real(dprec), intent(in)      :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec), intent(in)      :: Pk(Glob_n,Glob_n), Pl(Glob_n,Glob_n)
real(dprec), intent(out)     :: TranDipolLength_kl, TranDipolVelocity_kl


!Parameters (These are needed to declare static arrays. Using static 
!arrays makes the function call a little faster in comparison with 
!the case when arrays are dynamically allocated in stack)
integer,parameter            :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer,parameter            :: nnp=nn*(nn+1)/2

!Local variables
integer           :: n, np, Qtotal
real(dprec)       :: Lk(nn,nn),Ll(nn,nn)
real(dprec)       :: inv_Lk(nn,nn),inv_Ll(nn,nn)
real(dprec)       :: det_Lk,det_Ll
real(dprec)       :: tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
real(dprec)       :: inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn),tAk_inv_tAkl(nn,nn)
real(dprec)       :: det_tAkl
real(dprec)       :: tvl(nn),vi(nn),vi_tAk(nn),vi_tAk_inv_tAkl(nn)

integer           :: i,j,k,indx,ii
real(dprec)       :: temp0, temp1, temp2, temp3, charge_mass0
real(dprec)       :: W1(nn,nn),W2(nn,nn)

n=Glob_n
np=Glob_np

!First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
indx=0
Do i=1,n
	Do j=i,n
		indx=indx+1
		Lk(i,j)=ZERO
		Lk(j,i)=vechLk(indx)
		Ll(i,j)=ZERO
		Ll(j,i)=vechLl(indx)
	EndDo
EndDo

Do i=1,n
	Do j=i,n
		temp1=ZERO
		Do k=1,i
			temp1=temp1+Lk(i,k)*Lk(j,k)
		EndDo 
		tAk(i,j)=temp1
		tAk(j,i)=temp1
		temp1=ZERO
		Do k=1,i
			temp1=temp1+Ll(i,k)*Ll(j,k)
		EndDo 
		tAl(i,j)=temp1
		tAl(j,i)=temp1
	EndDo
EndDo

!Then we permute elements of Ak and Al to account for 
!the action of the permutation matrix
!  tAl=Pl'*Al*Pl
!  tAk=Pk'*Ak*Pk
!We also form matrix tAkl=tAk+tAl
Do i=1,n
	Do j=1,n
		temp1=ZERO
		temp2=ZERO
		Do k=1,n
			temp1=temp1+Pl(k,j)*tAl(k,i)
			temp2=temp2+tAk(j,k)*Pk(k,i)
		EndDo
		W1(j,i)=temp1
		W2(j,i)=temp2
	EndDo
EndDo
!tAl=W1*Pl
!tAk=Pk'*W2
Do i=1,n  
  Do j=i,n
    temp1=ZERO
    temp2=ZERO
    Do k=1,n
	   temp1=temp1+W1(j,k)*Pl(k,i)
	   temp2=temp2+Pk(k,j)*W2(k,i)
	EndDo
  	tAl(j,i)=temp1
	tAl(i,j)=temp1
  	tAk(j,i)=temp2
	tAk(i,j)=temp2	
	tAkl(j,i)=temp1+temp2
	tAkl(i,j)=temp1+temp2
  EndDo
EndDo

!The determinants of Lk and Ll are just
!the products of their diagonal elements
det_Lk=ONE
det_Ll=ONE
Do i=1,n
	det_Lk=det_Lk*Lk(i,i)
	det_Ll=det_Ll*Ll(i,i)
EndDo

!After this we can do Cholesky factorization of tAkl.
!The Cholesky factor will be temporarily stored in the 
!lower triangle of W1
det_tAkl=ONE
Do i=1,n
	Do j=i,n
		temp1=tAkl(i,j)
		Do k=i-1,1,-1
			temp1=temp1-W1(i,k)*W1(j,k)
		EndDo
		IF (i==j) then
			W1(i,i)=sqrt(temp1)
			det_tAkl=det_tAkl*temp1
		Else
			W1(j,i)=temp1/W1(i,i)
			W1(i,j)=ZERO
		EndIF
	EndDo
EndDo

!Inverting tAkl using its Cholesky factor (stored in W1)
!and placing the result into inv_tAkl
Do i=1,n
	W1(i,i)=ONE/W1(i,i)
	Do j=i+1,n
		temp1=ZERO
		Do k=i,j-1
			temp1=temp1-W1(j,k)*W1(k,i)
		EndDo
		W1(j,i)=temp1/W1(j,j)
	EndDo
EndDo 

Do i=1,n
	Do j=i,n
		temp1=ZERO
		Do k=j,n
			temp1=temp1+W1(k,i)*W1(k,j)
		EndDo
		inv_tAkl(i,j)=temp1
		inv_tAkl(j,i)=temp1
	EndDo
EndDo  

!Finding the inverse of Akk and All using their Cholesky factors
!The result is placed in inv_Akk and inv_All
Do i=1,n
	! W1(i,i)=ONE/Lk(i,i)
	W2(i,i)=ONE/Ll(i,i)  
	Do j=i+1,n
		! temp1=ZERO
		temp2=ZERO
		Do k=i,j-1
			! temp1=temp1-Lk(j,k)*W1(k,i)
			temp2=temp2-Ll(j,k)*W2(k,i)
		EndDo
		! W1(j,i)=temp1/Lk(j,j)
		W2(j,i)=temp2/Ll(j,j)
	EndDo
EndDo 

Do i=1,n
	Do j=i,n
		! temp1=ZERO
		temp2=ZERO
		Do k=j,n
			! temp1=temp1+W1(k,i)*W1(k,j)
			temp2=temp2+W2(k,i)*W2(k,j)       
		EndDo
		! inv_Akk(i,j)=ONEHALF*temp1
		! inv_Akk(j,i)=ONEHALF*temp1
		inv_All(i,j)=ONEHALF*temp2
		inv_All(j,i)=ONEHALF*temp2	 
	EndDo
EndDo

!tvl=Pl'*vl
Do i=1,n
	tvl(i)=Pl(ml,i)
EndDo


!====================================================
!			Evaluating Matrix Elements				
!====================================================


!         (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5       
! temp1= ----------------------------------------
!                    (det_tAkl)^1.5
 
 
temp1 = abs(det_Ll*det_Lk) / det_tAkl
temp1 = temp1 * sqrt(abs(temp1))



!                       2^(3*n/2)            temp1
! TranDipolLength_kl = ----------- * -----------------------
!                        sqrt(2)      Sqrt(vl'*inv_All*vl) 


TranDipolLength_kl = Glob_2raised3n2 * temp1 / sqrt(TWO*inv_All(ml,ml))
TranDipolVelocity_kl = TranDipolLength_kl * TWO



!
!                                                                m_i
! TranDipolLength_kl = TranDipolLength_kl *SUM{i}(q_i - Q_tot * -----)* vi'*inv_tAkl*vl
!                                                                m0
!
Qtotal = Glob_PseudoCharge0

Do i=1,n
	Qtotal = Qtotal + Glob_PseudoCharge(i)
EndDo



temp1 = ZERO

Do i=1,n 							!pseudo-particles

	temp2 = ZERO
	Do j=1,n 						!trace elements
		temp2 = temp2 + inv_tAkl(i,j)* tvl(j)
	EndDo


	temp0 = Glob_PseudoCharge(i) - Qtotal*Glob_Mass(i+1)/Glob_MassTotal
	temp1 = temp1 + temp0 * temp2

EndDo

TranDipolLength_kl = TranDipolLength_kl * temp1

!
!                                                      
! TranDipolLength_kl = TranDipolLength_kl * SUM{i} vi' * tAk inv_tAkl * vl
!                                                      
!

temp0 = ZERO
temp1 = ZERO
temp2 = ZERO
temp3 = ZERO

charge_mass0 = Glob_PseudoCharge0 / Glob_Mass(1)

vi_tAk = ZERO
vi_tAk_inv_tAkl = ZERO

Do i=1,n

	vi = ZERO
	vi(i) = ONE
	
	temp0 = charge_mass0 - Glob_PseudoCharge(i)/Glob_Mass(i+1)
	
	Do k=1,n

		temp1 = ZERO
		Do j=1,n
			temp1 = temp1 + vi(j) * tAk(j,k) 
		EndDo
		vi_tAk(k) = temp1
		
	EndDo


	Do k=1,n

		temp2 = ZERO
		Do j=1,n
			temp2 = temp2 + vi_tAk(j) * inv_tAkl(j,k) 
		EndDo
		vi_tAk_inv_tAkl(k) = temp2
		
	EndDo
	

	Do k=1,n

		temp3 = temp3 + vi_tAk_inv_tAkl(k) * tvl(k) * temp0
		
	EndDo
EndDo

TranDipolVelocity_kl = TranDipolVelocity_kl * temp3


end subroutine MatrixElemenTranDipoleMoment

subroutine spinPreCalc(n, nFactorial, parityFactor, SOmassChargeCoefficient, AMMmassChargeCoefficient, &
	ketMatrix, spatialYoung0, spatialYoung1, SiMinusME, SiPlusME, SziME, spinFreeME)
	use spinStuff
	implicit none
  
	!input vars:
	character(len = maxLen), intent(in) :: spatialYoung0, spatialYoung1
	integer, intent(in) :: n, nFactorial
  
	!output vars:
	real(dprec), dimension(nFactorial), intent(out) :: parityFactor
	real(dprec), dimension(n, n, 4), intent(out) :: SOmassChargeCoefficient, AMMmassChargeCoefficient
	real(dprec), dimension(n, n, nFactorial), intent(out) :: ketMatrix
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
	  enddo
	enddo
  
	do i = 1, n
	  do j = 1, n
		SOmassChargeCoefficient(i, j, 4) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
		(Glob_Mass(i + 1)**TWO)
	  enddo
	enddo
  
	AMMmassChargeCoefficient = ZERO
	do i = 1, n
	  AMMmassChargeCoefficient(i, i, 1) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / Glob_Mass(i + 1)**TWO
	enddo
  
	do i = 1, n
	  do j = 1, n
		AMMmassChargeCoefficient(i, j, 3) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
		(Glob_Mass(i + 1) * Glob_Mass(j + 1))
	  enddo
	enddo
  
	do i = 1, n
	  do j = 1, n
		AMMmassChargeCoefficient(i, j, 4) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
		(Glob_Mass(i + 1)**TWO)
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
	
	call getSpinFunction(n, nFactorial, mySpatialYoung0, allPermutations, parities, &
	finalSpinFunction0, primitives0, numberOfPrimitives0)

	call getSpinFunction(n, nFactorial, mySpatialYoung1, allPermutations, parities, &
	finalSpinFunction1, primitives1, numberOfPrimitives1)

	call getSpinOpMeanValues(n, nFactorial, allPermutations, finalSpinFunction0, finalSpinFunction1, primitives0, primitives1, &
	numberOfPrimitives0, numberOfPrimitives1, spinFreeME, SiMinusME, SiPlusME, SziME)

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


  subroutine spinDependentMatrixElements(m_k, mm_k, vechLk, vechLl, Pket, &
	SOspinCoeff, SOmassChargeCoefficient, AMMmassChargeCoefficient, SO1kl, SO2kl, &
	AMM1kl, AMM2kl)
  !This subroutine computes symmetry adapted off-diagonal SO matrix element
  !between 2p(P) and S Gaussians. This matrix element
  !is used in calculations of expectation values.
  
  !Input:
  !   m_k,mm_k :: integers that determine which x or y-components is in the
  !                premultiplier of the Gaussian
  !   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.
  
  !Output (all matrix elements are computed with normilized functions):
  
  !   SO1kl, SO2kl  :: Spin-Orbit corrections (without the factor of alpha**2)
  !         1 and 2 stay for spin-same orbit and spin-another orbit contributions
  !   AMM1kl, AMM2kl  :: AMM corrections (without the factor of alpha**2)
  !         1 and 2 stay for spin-same orbit and spin-another orbit contributions
  
  !Input vars:
  integer,intent(in)       :: m_k, mm_k
  real(dprec),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
  real(dprec),intent(in)   :: Pket(Glob_n,Glob_n)
  real(dprec),intent(in)   :: SOspinCoeff(Glob_n), &
  SOmassChargeCoefficient(Glob_n, Glob_n, 4), AMMmassChargeCoefficient(Glob_n, Glob_n, 4)

  !Output vars:
  real(dprec), intent(out)  :: SO1kl, SO2kl, AMM1kl, AMM2kl

  !Parameters (These are needed to declare static arrays. Using static
  !arrays makes the function call a little faster in comparison with
  !the case when arrays are dynamically allocated in stack)
  integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
  integer,parameter :: nnp=nn*(nn+1)/2

  !Local variables
  integer           n, np
  integer           tvk(nn),tvl(nn)
  real(dprec)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
  real(dprec)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
  real(dprec)       inv_tAkl(nn,nn)
  
  
  real(dprec)       W1(nn,nn)
  real(dprec)       temp1, temp2, det_tAkl
  integer :: i, j, k, indx
  
  
  integer :: pm_k, pmm_k ! new non-zero components of v_k and v_l
  real(dprec) :: commonFactor, gamma, gamma_diag, localEps
  
  ! V-quantities 
  real(dprec) :: jiAlAklinvVk, jiAklinvVk, jjAlAklinvVk, jjAklinvVk
  
  ! W-quantities
  real(dprec) :: jiAklinvWk, jiAlAklinvWk, jjAlAklinvWk, jjAklinvWk
  
  integer :: indexI, indexJ ! indices enumerating particles from H_SO and AMM operators
  
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
  pm_k = m_k
  pmm_k = mm_k
  !common factor (sqrt(ONEHALF) - for consistent normalization with Skl)
  commonFactor = TWO / (THREE * sqrt(THREE)) * Glob_Piraised3n2 / (SQRTPI * det_tAkl * sqrt(det_tAkl))
  
  SO1kl = ZERO
  SO2kl = ZERO
  
  AMM1kl = ZERO
  AMM2kl = ZERO
  
  do indexI = 1, n
	if (abs(SOspinCoeff(indexI)) < localEps) cycle

   	! gamma diagonal coefficient
   	gamma_diag = ONE / sqrt(inv_tAkl(indexI, indexI))
  
  	! calculating all the traces we need
   	! tr(Axy') is computed as (y, Ax) everywhere
   	! variable names: jiAlAklinvVk = (j^i, A_l A_{kl}^(-1) v_k) (names doesnt account for permutations)
  
  
   	jiAlAklinvVk = ZERO
   	do i = 1, n
		jiAlAklinvVk = jiAlAklinvVk + tAl(indexI, i) * inv_tAkl(i, pm_k)
   	enddo
   	jiAlAklinvWk = ZERO
   	do i = 1, n
		jiAlAklinvWk = jiAlAklinvWk + tAl(indexI, i) * inv_tAkl(i, pmm_k)
   	enddo
   
   	jiAklinvVk = inv_tAkl(indexI, pm_k)
   	jiAklinvWk = inv_tAkl(indexI, pmm_k)
  
   	! I term -> diagonal (spin-same orbit) matrix element (f[ii, ii])
   	temp1 = gamma_diag ** 3 *  (jiAklinvVk * jiAlAklinvWk - jiAklinvWk * jiAlAklinvVk) 
   	SO1kl = SO1kl + SOspinCoeff(indexI) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
   	AMM1kl = AMM1kl + SOspinCoeff(indexI) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
   

   ! these traces are needed for spin-other-orbit contribution and SSNC
   	do indexJ = 1, n
		if (indexI == indexJ) cycle
  
	 	gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
	 	inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))

  
	 	jjAlAklinvVk = ZERO
	 	do i = 1, n
	   		jjAlAklinvVk = jjAlAklinvVk + tAl(indexJ, i) * inv_tAkl(i, pm_k)
	 	enddo

	 	jjAlAklinvWk = ZERO
	 	do i = 1, n
	   		jjAlAklinvWk = jjAlAklinvWk + tAl(indexJ, i) * inv_tAkl(i, pmm_k)
	 	enddo
  
	 	jjAklinvVk = inv_tAkl(indexJ, pm_k)
	 	jjAklinvWk = inv_tAkl(indexJ, pmm_k)
  
	 	!! II term -> f[ii, ij]
	 	temp1 = gamma ** 3 * (jiAklinvVk * jjAlAklinvWk - jiAklinvWk * jjAlAklinvVk)
	   	SO2kl = SO2kl + SOspinCoeff(indexI) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
	 
  
	 	!! III term -> f[ij, jj]
		temp1 = gamma ** 3 * (jjAlAklinvWk * (jjAklinvVk - jiAklinvVk) + &
		jjAlAklinvVk * (jiAklinvWk - jjAklinvWk))
	   	SO2kl = SO2kl + SOspinCoeff(indexI) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
	   	AMM2kl = AMM2kl + SOspinCoeff(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
  
  
	 	!! IV term -> f[ij, ii] (i <->j of the III term)
	 	temp1 = gamma ** 3 * (jiAlAklinvWk * (jiAklinvVk - jjAklinvVk) + &
		jiAlAklinvVk * (jjAklinvWk - jiAklinvWk))
	    SO2kl = SO2kl + SOspinCoeff(indexI) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
	    AMM2kl = AMM2kl + SOspinCoeff(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
  
   enddo ! indexJ cycle
  enddo ! indexI cycle
  
  SO1kl = SO1kl * commonFactor 
  SO2kl = SO2kl * commonFactor
  AMM1kl = AMM1kl * commonFactor
  AMM2kl = AMM2kl * commonFactor
  
  
  end subroutine spinDependentMatrixElements


end module matelem




