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



subroutine MatrixElemenTranDipoleMoment(ml, vechLk, vechLl, Pk, Pl, TDkl)

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

! this subroutine computes the dollowing expression
! which is a part of transition dipole momentum calculation.
!                                 m_i         < Pk*fk0 | z_i | Pl*fl1 >
! TDkl{i} = SUM{i}(q_i - Q_tot * -----)*-------------------------------------------
!                                 m0    Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )
!
!
!        < Pk*fk0 | z_i | Pl*fl1 >                 2^(3*n/2)  
! --------------------------------------------- = ----------- * 
!     Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )         sqrt(2)
!
!  (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5         vi'*inv_tAkl*vl
! ---------------------------------------- * -----------------------
!            (det_tAkl)^1.5                   Sqrt(vl'*inv_All*vl)              

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
real(dprec), intent(out)     :: TDkl

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
real(dprec)       :: inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn)
real(dprec)       :: det_tAkl
real(dprec)       :: tvl(nn)

integer           :: i,j,k,indx
real(dprec)       :: temp0, temp1, temp2
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

!Evaluating Matrix Elements

!         (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5       
! temp1= ----------------------------------------
!                    (det_tAkl)^1.5
                             
temp1=abs(det_Ll*det_Lk)/det_tAkl
temp1=temp1*sqrt(abs(temp1))

!        2^(3*n/2)            temp1
! TDkl = ----------- * -----------------------
!         sqrt(2)      Sqrt(vl'*inv_All*vl) 
 
TDkl=Glob_2raised3n2*temp1/sqrt(TWO*inv_All(ml,ml))


!                                        m_i
! TDkl{ij} = TDkl *SUM{i}(q_i - Q_tot * -----)* vi'*inv_tAkl*vl
!                                         m0

Qtotal=Glob_PseudoCharge0
Do i=1,n
	Qtotal=Qtotal+Glob_PseudoCharge(i)
EndDo

temp1=ZERO
Do i=1,n 						!pseudo-particles
	temp2=ZERO
	Do j=1,n 					!trace elements
		temp2=temp2+inv_tAkl(j,i)*tvl(j)
	EndDo
	! temp1=temp1+Glob_PseudoCharge(i)*temp2
	temp0=Glob_PseudoCharge(i)-Qtotal*Glob_Mass(i+1)/Glob_Mass(1)
	temp1=temp1+temp0*temp2
EndDo

TDkl=TDkl*temp1

end subroutine MatrixElemenTranDipoleMoment


end module matelem




