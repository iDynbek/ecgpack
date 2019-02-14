module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=0 and L=1 Gaussians.
use globvars
implicit none

contains


subroutine MatrixElementsL0L1(ml, vechLk, vechLl, Pk, Pl, Hklij)

!This subroutine computes symmetry adapted matrix element with
!a real L=0 and a real L=1 correlated Gaussians:
!
!fk = exp[-r'(Lk*Lk')r]  
!fl = z_{ml} exp[-r'(Ll*Ll')r] 
!
!where m_l is some integer between 1 and n (n is the number of 
!pseudoparticles). Symmetry adaption is applied to the bra using 
!permutation matrix Pk and to the ket using permutation matrix Pl.
!
!Pk fk = exp[-r' {Pk'*(Lk*Lk')*Pk} r]  
!Pl fl = (Pl z_{m_l}) exp[-r' {Pl'*(Ll*Ll')*Pl} r] 
!
!Input:     
!   m_l            :: integer that determine which z-component is in the
!                     premultiplier of the Gaussian
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters. 
!   Pk, Pl         :: The symmetry permutation matrices of size n x n
!
!Output:
!   Hklij	   :: Matrix element (normalized)

!Arguments
integer,intent(in)          :: ml
real(dprec),intent(in)      :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec),intent(in)      :: Pk(Glob_n,Glob_n), Pl(Glob_n,Glob_n)
real(dprec),intent(out)     :: Hklij

!Parameters (These are needed to declare static arrays. Using static 
!arrays makes the function call a little faster in comparison with 
!the case when arrays are dynamically allocated in stack)
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
integer           :: n, np
real(dprec)       :: Lk(nn,nn),Ll(nn,nn)
real(dprec)       :: inv_Lk(nn,nn),inv_Ll(nn,nn)
real(dprec)       :: det_Lk,det_Ll
real(dprec)       :: tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
real(dprec)       :: inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn)
real(dprec)       :: det_tAkl
integer           :: tml

integer           :: i,j,k,indx
real(dprec)       :: temp1, temp2
real(dprec)       :: W1(nn,nn),W2(nn,nn)

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
!  tAl=Pl'*Al*Pl
!  tAk=Pk'*Ak*Pk
!We also form matrix tAkl=tAk+tAl
do i=1,n
  do j=1,n
	temp1=ZERO
	temp2=ZERO
    do k=1,n
       temp1=temp1+Pl(k,j)*tAl(k,i)
       temp2=temp2+tAk(j,k)*Pk(k,i)
	enddo
	W1(j,i)=temp1
	W2(j,i)=temp2
  enddo
enddo
!tAl=W1*Pl
!tAk=Pk'*W2
do i=1,n  
  do j=i,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
	   temp1=temp1+W1(j,k)*Pl(k,i)
	   temp2=temp2+Pk(k,j)*W2(k,i)
	enddo
  	tAl(j,i)=temp1
	tAl(i,j)=temp1
  	tAk(j,i)=temp2
	tAk(i,j)=temp2	
	tAkl(j,i)=temp1+temp2
	tAkl(i,j)=temp1+temp2
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

!Computing tml, the index of nonzero (=1) component of vector tvl=Pl'*vl,
!where vl(i) = 1 if i=ml and 0 otherwise
do i=1,n
   if(abs(Pl(ml,i))==1) tml=i
enddo

!Evaluating Matrix Elements
select case (Glob_DRMC(Glob_CurrDRMCStep)%Action(1:9))
   case('OP_DIPOLE')
      temp1=(abs(det_Ll*det_Lk)/det_tAkl)**THREEHALF
      Hklij=SIX*temp1/sqrt(inv_All(ml,ml))
      temp1=ZERO
      do i=1,n
         temp1=temp1+Glob_PseudoCharge(i)*inv_tAkl(tml,i)
      enddo
      Hklij=Hklij*temp1
endselect

end subroutine MatrixElementsL0L1

end module matelem
