module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=1 Gaussians without normalization
use globvars
implicit none

contains

subroutine MatrixElementsL1(m_k, m_l, vechLk, vechLl, P, &
               Hkl, Skl, Dk, Dl, grad_k, grad_l)
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
!   m_k, m_l :: integers that determine which z-component is in the
!                premultiplier of the Gaussian
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of
!     exponential parameters.
!   P  :: The symmetry permutation matrix of size n x n
!   grad_k, grad_l :: Gradient flags
!   grad_k=.true.  means that dHkldvechLk, dSkldvechLk need to be computed.
!   grad_l=.true.  means that dHkldvechLl, dSkldvechLl need to be computed.
!Output:
!   Hkl	 ::	Hamiltonian term (normalized)
!   Skl	 ::	Overlap matrix element (normalized)
!   Dk,Dl:: derivatives of normalized Hkl and Skl wrt Paramk
!           and Paraml respectively. They are ordered in the
!           following manner:
!           Dk=(dHkldvechLk,dSkldvechLk)
!           Dl=(dHkldvechLl,dSkldvechLl)


!Arguments
integer,intent(in)          :: m_k,m_l
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
integer           tvl(nn)
real(dprec)       dHkldvechLk(nnp), dHkldvechLl(nnp)
real(dprec)       dSkldvechLk(nnp), dSkldvechLl(nnp)
real(dprec)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
real(dprec)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
real(dprec)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn)
real(dprec)       inv_tAkltAl(nn,nn),inv_tAkltAlM(nn,nn)
real(dprec)       inv_tAklAk(nn,nn),inv_tAklAkM(nn,nn)
real(dprec)       eta1(nn,nn),sqrt_eta1(nn,nn),eta2(nn,nn),Rkl(nn,nn)
real(dprec)       W1(nn,nn),W2(nn,nn),W3(nn,nn),W4(nn,nn)
real(dprec)       twosym_tFkl(nn,nn),two_Fkk(nn,nn),two_Fll(nn,nn),twosym_tGkl(nn,nn)
real(dprec)       tKkl(nn,nn),tUkl(nn,nn),tWkl(nn,nn)
real(dprec)       twosym_tQkl(nn,nn),twosym_tDkl(nn,nn)
real(dprec)       inv_tAkltvl(nn),vkinv_tAkl(nn),vkinv_tAkltAlM(nn)
real(dprec)       u1(nn),u2(nn),u3(nn)
real(dprec)       temp1, temp2, temp3, temp4, temp5, temp6
real(dprec)       det_Lk, det_Ll, det_tAkl
real(dprec)       tau1,tau2,tau3,tau11
real(dprec)       Tkl, Vkl
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

!I delete non-necessary parts of the code

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

!I delete non-necessary parts of the code

!Finding the inverse of Akk and All using their Cholesky factors
!The result is placed in inv_Akk and inv_All
!do i=1,n
!  W1(i,i)=ONE/Lk(i,i)
!  W2(i,i)=ONE/Ll(i,i)
!  do j=i+1,n
!    temp1=ZERO
!    temp2=ZERO
!    do k=i,j-1
!      temp1=temp1-Lk(j,k)*W1(k,i)
!      temp2=temp2-Ll(j,k)*W2(k,i)
!    enddo
!    W1(j,i)=temp1/Lk(j,j)
!    W2(j,i)=temp2/Ll(j,j)
!  enddo
!enddo

!do i=1,n
!  do j=i,n
!     temp1=ZERO
!     temp2=ZERO
!    do k=j,n
!       temp1=temp1+W1(k,i)*W1(k,j)
!       temp2=temp2+W2(k,i)*W2(k,j)
!     enddo
!     inv_Akk(i,j)=ONEHALF*temp1
!	 inv_Akk(j,i)=ONEHALF*temp1
!     inv_All(i,j)=ONEHALF*temp2
!	 inv_All(j,i)=ONEHALF*temp2
!   enddo
!enddo

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
!temp1=abs(det_Ll*det_Lk)/det_tAkl
temp1=det_tAkl*sqrt(det_tAkl)
!Skl=Glob_2raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
Skl=Glob_Piraised3n2*tau3/(TWO*temp1)

!Doing multiplication inv_tAkltAl=inv_tAkl*tAl
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
    enddo
    inv_tAkltAl(j,i)=temp1
  enddo
enddo

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkltAl(j,k)*Glob_MassMatrix(k,i)
    enddo
    inv_tAkltAlM(j,i)=temp1
  enddo
enddo

!Computing tau1=tr[inv_tAkltAlM*Ak]
tau1=ZERO
do i=1,n
  temp1=ZERO
  do k=1,n
    temp1=temp1+inv_tAkltAlM(i,k)*Ak(k,i)
  enddo
  tau1=tau1+temp1
enddo
!Computing tau2 = vk'*inv_tAkltAlM*Ak*inv_tAkltvl
!We do it by multiplying twice the row-vector on the left
!by a matrix on the right and computing a dot product in the end.

!vkinv_tAkltAlM'=vk'*inv_tAkltAlM
do i=1,n
  vkinv_tAkltAlM(i)=inv_tAkltAlM(m_k,i)
enddo
!u1=vkinv_tAkltAlM'*Ak
!tau2=u1'*inv_tAkltvl (storage for u1 as such is not needed, we use temp1=u1(i))
tau2=ZERO
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+vkinv_tAkltAlM(j)*Ak(j,i)
  enddo
  tau2=tau2+temp1*inv_tAkltvl(i)
enddo

!Evaluating the kinetic energy
Tkl=Skl*(SIX*tau1+FOUR*tau2/tau3)

!Evaluating eta1(i,j), sqrt_eta1(i,j), eta2(i,j), Rkl(i,j),
!and the potential energy. Notice that only the lower triangles
!of eta1, sqrt_eta1, eta2, and Rkl are filled.
Vkl=ZERO
temp1=Skl*(TWO/SQRTPI)
do i=1,n
  temp2=inv_tAkl(i,i)
  temp3=sqrt(temp2)
  eta1(i,i)=temp2
  sqrt_eta1(i,i)=temp3
  !Getting row m_k of matrix inv_tAkl*Jii*inv_tAkl
  !as only this row is needed to compute eta2(i,i)
  do k=1,n
    u1(k)=inv_tAkl(i,m_k)*inv_tAkl(k,i)
  enddo
  temp4=ZERO
  do k=1,n
    temp4=temp4+u1(k)*tvl(k)
  enddo
  eta2(i,i)=temp4
  Rkl(i,i)=temp1*(ONE-temp4/(THREE*temp2*tau3))/temp3
  Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*Rkl(i,i)
enddo
do i=1,n
  do j=i+1,n
    temp2=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    temp3=sqrt(temp2)
    eta1(j,i)=temp2
    sqrt_eta1(j,i)=temp3
    !Getting row m_k of matrix inv_tAkl*Jij*inv_tAkl
    !as only this row is needed to compute eta2(i,i)
    do k=1,n
      u1(k)=(inv_tAkl(i,m_k)-inv_tAkl(j,m_k))*(inv_tAkl(k,i)-inv_tAkl(k,j))
    enddo
    temp4=ZERO
    do k=1,n
      temp4=temp4+u1(k)*tvl(k)
    enddo
    eta2(j,i)=temp4
    Rkl(j,i)=temp1*(ONE-temp4/(THREE*temp2*tau3))/temp3
    Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*Rkl(j,i)
  enddo
enddo

Hkl=Tkl+Vkl

!Now we start computing the gradient of Skl

if (grad_k.or.grad_l) then
  !Evaluating matrix tKkl = inv_tAkltvl * vkinv_tAkl'
  !which will be used a lot below.
  do i=1,n
    do j=1,n
      tKkl(i,j)=inv_tAkltvl(i)*vkinv_tAkl(j)
    enddo
  enddo
  !Evaluating twosym_tFkl = tFkl + tFkl' , where tFkl = (3/2)*inv_tAkl + tKkl/tau3
  do i=1,n
    do j=1,i
      twosym_tFkl(i,j)=THREE*inv_tAkl(j,i)+(tKkl(i,j)+tKkl(j,i))/tau3
      twosym_tFkl(j,i)=twosym_tFkl(i,j)
    enddo
  enddo
endif

if (grad_k) then

  !I delete non-necessary parts of the code
  !Evaluating two_Fkk = 2*Fkk, where
  !Fkk = (3/2) * inv_Akk + (inv_Akk * vk * vk' * inv_Akk)/(vk' * inv_Akk * vk)
  !do i=1,n
  !  u1(i)=inv_Akk(m_k,i)
  !enddo
  !temp1=TWO/u1(m_k)
  !do i=1,n
  !  do j=1,i
  !     two_Fkk(i,j)=THREE*inv_Akk(i,j)+temp1*u1(i)*u1(j)
  !     two_Fkk(j,i)=two_Fkk(i,j)
  !  enddo
  !enddo
  !Evaluating Skl*vech((two_Fkk-twosym_tFkl)*Lk)'
  !Evaluating -Skl*vech((twosym_tFkl)*Lk)'  !new line
  indx=0
  do i=1,n
    do j=i,n
      temp1=ZERO
      do k=i,n
        temp1=temp1-twosym_tFkl(k,j)*Lk(k,i)
      enddo
      indx=indx+1
      Dk(Glob_np+indx)=Skl*temp1
    enddo
  enddo
endif

if (grad_l) then
  !Evaluating twosym_tGkl = P * twosym_tFkl *P'
  do i=1,n
    do j=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+P(i,k)*twosym_tFkl(k,j)
      enddo
      W1(i,j)=temp1
    enddo
  enddo
  do i=1,n
    do j=1,i
      temp1=ZERO
      do k=1,n
        temp1=temp1+W1(i,k)*P(j,k)
      enddo
      twosym_tGkl(i,j)=temp1
      twosym_tGkl(j,i)=temp1
    enddo
  enddo

  !I delete non-necessary parts of the code
  !Evaluating two_Fll = 2*Fll, where
  !Fll = (3/2) * inv_All + (inv_All * vl * vl' * inv_All)/(vl' * inv_All * vl)
  !do i=1,n
  !  u1(i)=inv_All(m_l,i)
  !enddo
  !temp1=TWO/u1(m_l)
  !do i=1,n
  !  do j=1,i
  !     two_Fll(i,j)=THREE*inv_All(i,j)+temp1*u1(i)*u1(j)
  !     two_Fll(j,i)=two_Fll(i,j)
  !  enddo
  !enddo
  !Evaluating Skl*vech((two_Fll-twosym_tGkl)*Ll)'
  !Evaluating -Skl*vech((twosym_tGkl)*Ll)'  !new line
  indx=0
  do i=1,n
    do j=i,n
      temp1=ZERO
      do k=i,n
        temp1=temp1-twosym_tGkl(k,j)*Ll(k,i)
      enddo
      indx=indx+1
      Dl(Glob_np+indx)=Skl*temp1
    enddo
  enddo
endif

!Gradient of Tkl

if (grad_k) then
  !Computing W1=6*inv_tAkltAlM*inv_tAkltAl'
  do i=1,n
    do j=i,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAkltAlM(j,k)*inv_tAkltAl(i,k)
      enddo
      W1(j,i)=SIX*temp1
      W1(i,j)=W1(j,i)
    enddo
  enddo
  !Computing u1'=vkinv_tAkltAlM'*inv_tAkltAl'
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+vkinv_tAkltAlM(j)*inv_tAkltAl(i,j)
    enddo
    u1(i)=temp1
  enddo
  !Computing u2=inv_tAkltAlM*Ak*inv_tAkltvl
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+Ak(i,j)*inv_tAkltvl(j)
    enddo
    u3(i)=temp1
  enddo
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_tAkltAlM(i,j)*u3(j)
    enddo
    u2(i)=temp1
  enddo
  !Computing matrix tUkl=W1+(4/tau3)*(inv_tAkltvl*u1'-u2*vkinv_tAkl')+(4*tau2/tau3^2)*tKkl
  temp1=FOUR/tau3
  temp2=temp1*tau2/tau3
  do i=1,n
    do j=1,n
      tUkl(j,i)=W1(j,i)+temp1*(inv_tAkltvl(j)*u1(i)-u2(j)*vkinv_tAkl(i))+temp2*tKkl(j,i)
    enddo
  enddo
  !Evaluating (Tkl/Skl)*dSkldvechLk' + Skl*vech((tUkl+tUkl')*Lk)'
  temp4=Tkl/Skl
  indx=0
  do i=1,n
    do j=i,n
      temp1=ZERO
      do k=i,n
        temp1=temp1+(tUkl(k,j)+tUkl(j,k))*Lk(k,i)
      enddo
      indx=indx+1
      Dk(indx)=Skl*temp1+temp4*Dk(Glob_np+indx)
    enddo
  enddo
endif

if (grad_l) then
  !Computing inv_tAklAk=inv_tAkl*Ak
  do i=1,n
    do j=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAkl(j,k)*Ak(i,k)
      enddo
      inv_tAklAk(j,i)=temp1
    enddo
  enddo
  !Computing inv_tAklAkM=inv_tAklAk*M
  do i=1,n
    do j=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAklAk(j,k)*Glob_MassMatrix(i,k)
      enddo
      inv_tAklAkM(j,i)=temp1
    enddo
  enddo
  !Computing W1=6*inv_tAklAkM*inv_tAklAk'
  do i=1,n
    do j=i,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAklAkM(j,k)*inv_tAklAk(i,k)
      enddo
      W1(j,i)=SIX*temp1
      W1(i,j)=W1(j,i)
    enddo
  enddo
  !Computing u1=inv_tAklAkM*inv_tAklAk'*tvl
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_tAklAk(j,i)*tvl(j)
    enddo
    u3(i)=temp1
  enddo
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_tAklAkM(i,j)*u3(j)
    enddo
    u1(i)=temp1
  enddo
  !Computing u2'=vkinv_tAkltAlM'*inv_tAklAk'
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+vkinv_tAkltAlM(j)*inv_tAklAk(i,j)
    enddo
    u2(i)=temp1
  enddo
  !Computing tWkl = P*(
  !     W1 + (4/tau3)*(u1*vkinv_tAkl'-inv_tAkltvl*u2') + (4*tau2/tau3^2)*tKkl
  !                   )*P'
  temp1=FOUR/tau3
  temp2=temp1*tau2/tau3
  do i=1,n
    do j=1,n
      W3(j,i)=W1(j,i)+temp1*(u1(j)*vkinv_tAkl(i)-inv_tAkltvl(j)*u2(i))+temp2*tKkl(j,i)
    enddo
  enddo
  do i=1,n
    do j=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+P(i,k)*W3(k,j)
      enddo
      W2(i,j)=temp1
    enddo
  enddo
  do i=1,n
    do j=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+W2(i,k)*P(j,k)
      enddo
      tWkl(i,j)=temp1
    enddo
  enddo
  !Evaluating (Tkl/Skl)*dSkldvechLl' + Skl*vech((tWkl+tWkl')*Ll)'
  temp4=Tkl/Skl
  indx=0
  do i=1,n
    do j=i,n
      temp1=ZERO
      do k=i,n
        temp1=temp1+(tWkl(k,j)+tWkl(j,k))*Ll(k,i)
      enddo
      indx=indx+1
      Dl(indx)=Skl*temp1+temp4*Dl(Glob_np+indx)
    enddo
  enddo
endif

!Dk(1:Glob_np)=ZERO
!Dl(1:Glob_np)=ZERO

!Gradient of Vkl
if (grad_k.or.grad_l) then
  do i=1,n
    !Computing twosym_tQkl(i,i)=tQkl(i,i)+tQkl(i,i)'
    temp1=(TWO/SQRTPI)/(eta1(i,i)*sqrt_eta1(i,i))
    temp2=ONE-eta2(i,i)/(eta1(i,i)*tau3)
    temp3=ONETHIRD/tau3
    temp4=ONETHIRD*eta2(i,i)/(tau3*tau3)
    do t=1,n
      do q=t,n
        temp5=inv_tAkl(q,i)*inv_tAkl(i,t)
        temp6=inv_tAkl(q,i)*(tKkl(i,t)+tKkl(t,i))+inv_tAkl(t,i)*(tKkl(i,q)+tKkl(q,i))
        twosym_tQkl(q,t)=temp1*(temp2*temp5+temp3*temp6-temp4*(tKkl(q,t)+tKkl(t,q)))
        twosym_tQkl(t,q)=twosym_tQkl(q,t)
      enddo
    enddo
    temp5=ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)
    if (grad_k) then
      !Evaluating (Rkl(i,i)/Skl)*dSkldvechLk' + Skl*vech(twosym_tQkl*Lk)'
      !and updating Dk
      temp2=Rkl(i,i)/Skl
      indx=0
      do t=1,n
        do q=t,n
          temp1=ZERO
          do k=t,n
            temp1=temp1+twosym_tQkl(k,q)*Lk(k,t)
          enddo
          indx=indx+1
          Dk(indx)=Dk(indx)+temp5*(temp2*Dk(Glob_np+indx)+Skl*temp1)
        enddo
      enddo
    endif
    if (grad_l) then
      !Computing twosym_tDkl = P * twosym_tQkl * P'
      do t=1,n
        do q=1,n
          temp1=ZERO
          do k=1,n
            temp1=temp1+P(q,k)*twosym_tQkl(k,t)
          enddo
          W1(q,t)=temp1
        enddo
      enddo
      do t=1,n
        do q=t,n
          temp1=ZERO
          do k=1,n
            temp1=temp1+W1(q,k)*P(t,k)
          enddo
          twosym_tDkl(q,t)=temp1
          twosym_tDkl(t,q)=temp1
        enddo
      enddo
      !Evaluating (Rkl(i,i)/Skl)*dSkldvechLl' + Skl*vech(twosym_tDkl*Ll)'
      !and updating Dl
      temp2=Rkl(i,i)/Skl
      indx=0
      do t=1,n
        do q=t,n
          temp1=ZERO
          do k=t,n
            temp1=temp1+twosym_tDkl(k,q)*Ll(k,t)
          enddo
          indx=indx+1
          Dl(indx)=Dl(indx)+temp5*(temp2*Dl(Glob_np+indx)+Skl*temp1)
        enddo
      enddo
    endif
  enddo
  do i=1,n
    do j=i+1,n
      !Computing twosym_tQkl(j,i)=tQkl(j,i)+tQkl(j,i)'
      temp1=(TWO/SQRTPI)/(eta1(j,i)*sqrt_eta1(j,i))
      temp2=ONE-eta2(j,i)/(eta1(j,i)*tau3)
      temp3=ONETHIRD/tau3
      temp4=ONETHIRD*eta2(j,i)/(tau3*tau3)
      do t=1,n
        do q=t,n
          temp5=(inv_tAkl(q,i)-inv_tAkl(q,j))*(inv_tAkl(i,t)-inv_tAkl(j,t))
          temp6=(inv_tAkl(q,i)-inv_tAkl(q,j))*(tKkl(i,t)-tKkl(j,t)+tKkl(t,i)-tKkl(t,j))+ &
                (tKkl(q,i)-tKkl(q,j)+tKkl(i,q)-tKkl(j,q))*(inv_tAkl(t,i)-inv_tAkl(t,j))
          twosym_tQkl(q,t)=temp1*(temp2*temp5+temp3*temp6-temp4*(tKkl(q,t)+tKkl(t,q)))
          twosym_tQkl(t,q)=twosym_tQkl(q,t)
        enddo
      enddo
      temp5=ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))
      if (grad_k) then
        !Evaluating (Rkl(j,i)/Skl)*dSkldvechLk' + Skl*vech(twosym_tQkl*Lk)'
        !and updating Dk
        temp2=Rkl(j,i)/Skl
        indx=0
        do t=1,n
          do q=t,n
            temp1=ZERO
            do k=t,n
              temp1=temp1+twosym_tQkl(k,q)*Lk(k,t)
            enddo
            indx=indx+1
            Dk(indx)=Dk(indx)+temp5*(temp2*Dk(Glob_np+indx)+Skl*temp1)
          enddo
        enddo
      endif
      if (grad_l) then
        !Computing twosym_tDkl = P * twosym_tQkl * P'
        do t=1,n
          do q=1,n
            temp1=ZERO
            do k=1,n
              temp1=temp1+P(q,k)*twosym_tQkl(k,t)
            enddo
            W1(q,t)=temp1
          enddo
        enddo
        do t=1,n
          do q=t,n
            temp1=ZERO
            do k=1,n
              temp1=temp1+W1(q,k)*P(t,k)
            enddo
            twosym_tDkl(q,t)=temp1
            twosym_tDkl(t,q)=temp1
          enddo
        enddo
        !Evaluating (Rkl(j,i)/Skl)*dSkldvechLl' + Skl*vech(twosym_tQkl*Ll)'
        !and updating Dl
        temp2=Rkl(j,i)/Skl
        indx=0
        do t=1,n
          do q=t,n
            temp1=ZERO
            do k=t,n
              temp1=temp1+twosym_tDkl(k,q)*Ll(k,t)
            enddo
            indx=indx+1
            Dl(indx)=Dl(indx)+temp5*(temp2*Dl(Glob_np+indx)+Skl*temp1)
          enddo
        enddo
      endif
    enddo
  enddo
endif

end subroutine MatrixElementsL1


subroutine MatrixElementsL1ForExpcVals(m_k, m_l, vechLk, vechLl, Pbra, Pket, &
           Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
           MVkl, drach_MVkl, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, prvalkl, &
           NumCFGridPoints, CFGrid, CFkl, NumDensGridPoints, DensGrid, Denskl, &
           AreCorrFuncNeeded, ArePartDensNeeded, AreMCorrFuncNeeded, AreMPartDensNeeded)
!This subroutine computes symmetry adapted matrix elements
!with two real L=1 correlated Gaussians. These matrix elements
!are used in calculations of expectation values.
!Symmetry adaptation is applied to the bra and ket using permutation matrices Pbra and Pket
!
!Input:
!   m_k, m_l :: integers that determine which z-component is in the
!       premultiplier of the Gaussian
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.
!   Pbra :: The symmetry permutation matrix of size n x n that is applied to bra
!   Pket :: The symmetry permutation matrix of size n x n that is applied to ket
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
! prvalkl   :: P(1/r^3_ij) - principal values of matric element 1/r^3_ij  (appears in the Araki-Sucker term for QED correction)
!NumCFGridPoints    :: Number of grid points for correlation function calculations
!CFGrid             :: Array containing grid points where matrix elements of
!                      correlation functions should be computed
!CFkl               :: Matrix elements of correlation functions
!NumDensGridPoints  :: Number of grid points for particle density calculations
!DensGrid           :: Array containing grid points where matrix elements of
!                      particle densities should be computed
!Denskl             :: Matrix elements of particle densities
!AreCorrFuncNeeded  :: flag indicating whether matrix elements of correlation
!                      functions need be computed
!ArePartDensNeeded  :: flag indicating whether matrix elements of particle
!                      densities need be computed
!AreMCorrFuncNeeded :: flag indicating whether matrix elements of momentum correlation
!                      functions need be computed
!AreMPartDensNeeded :: flag indicating whether matrix elements of particle
!                      momentum densities need be computed

!Arguments
integer,intent(in)       :: m_k,m_l
real(dprec),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec),intent(in)   :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
real(dprec),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl
real(dprec),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
real(dprec),intent(out)  :: deltarkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: drach_deltarkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: prvalkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: rmrmkl(Glob_n,Glob_n,Glob_n,Glob_n)
integer,intent(in)       :: NumCFGridPoints,NumDensGridPoints
real(dprec),intent(in)   :: CFGrid(2,NumCFGridPoints),DensGrid(2,NumDensGridPoints)
real(dprec),intent(out)  :: CFkl(Glob_n*(Glob_n+1)/2,NumCFGridPoints)
real(dprec),intent(out)  :: Denskl(Glob_n+1,NumDensGridPoints)
logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded,AreMCorrFuncNeeded,AreMPartDensNeeded

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
integer           n,np
integer           tvk(nn),tvl(nn)
real(dprec)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
real(dprec)       inv_tAk(nn,nn),inv_tAl(nn,nn),tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
real(dprec)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn), inv_tAkltAl(nn,nn)
real(dprec)       inv_invtAkinvtAl(nn,nn),tvkinv_tAk(nn),inv_tAltvl(nn)
real(dprec)       eta2(nn,nn),meta2(nn,nn),inv_tAkltAlM(nn,nn)
real(dprec)       W1(nn,nn),W2(nn,nn),W3(nn,nn),W4(nn,nn),W5(nn,nn),W6(nn,nn),W7(nn,nn)
real(dprec)       inv_tAkltvl(nn),tvkinv_tAkl(nn),tvkinv_tAkltAlM(nn),u1(nn)
real(dprec)       tvkinv_tAkinv_invtAkinvtAl(nn),inv_invtAkinvtAlinv_tAltvl(nn)
real(dprec)       temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9
real(dprec)       temp10,temp11,temp12,temp13,temp14,threshold,tr1, tr2, tr3, tr4
real(dprec)       det_Lk, det_Ll, det_tAkl, det_tAk, det_tAl, det_invtAkinvtAl
real(dprec)       tau1,tau2,tau3,inv_tau3 ,V2kl, tau4, MSkl
integer           i,j,k,t,indx,p,q
real(dprec)       TrAJ(nn,nn),sqrtTrAJ(nn,nn),TrAJAJ(nn,nn,nn,nn),MTrAJ(nn,nn),sqrtMTrAJ(nn,nn)
real(dprec)       jAj(nn,nn,nn,nn),jAtvl(nn,nn),tvkAj(nn,nn),Mass_For_Darwin(0:nn)

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

!I delete non-necessary parts of the code
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

!I delete non-necessary parts of the code
!Finding the inverse of Akk and All using their Cholesky factors
!The result is placed in inv_Akk and inv_All
!do i=1,n
!  W1(i,i)=ONE/Lk(i,i)
!  W2(i,i)=ONE/Ll(i,i)
!  do j=i+1,n
!    temp1=ZERO
!    temp2=ZERO
!    do k=i,j-1
!      temp1=temp1-Lk(j,k)*W1(k,i)
!      temp2=temp2-Ll(j,k)*W2(k,i)
!    enddo
!    W1(j,i)=temp1/Lk(j,j)
!    W2(j,i)=temp2/Ll(j,j)
!  enddo
!enddo

!do i=1,n
!  do j=i,n
!     temp1=ZERO
!     temp2=ZERO
!     do k=j,n
!       temp1=temp1+W1(k,i)*W1(k,j)
!       temp2=temp2+W2(k,i)*W2(k,j)
!     enddo
!     inv_Akk(i,j)=ONEHALF*temp1
!	 inv_Akk(j,i)=ONEHALF*temp1
!     inv_All(i,j)=ONEHALF*temp2
!	 inv_All(j,i)=ONEHALF*temp2
!   enddo
!enddo

!Computing tvk=Pbra'*vk and tvl=Pket'*vl
do i=1,n
  tvk(i)=Pbra(m_k,i)
  tvl(i)=Pket(m_l,i)
enddo



!Compute inv_tAkltvl = inv_tAkl * tvl
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+inv_tAkl(j,i)*tvl(j)
  enddo
  inv_tAkltvl(i)=temp1
enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAkl(j,i)
  enddo
  tvkinv_tAkl(i)=temp1
enddo


!Compute tau3=tvkinv_tAkl*tvl
tau3=ZERO
do i=1,n
  tau3=tau3+tvkinv_tAkl(i)*tvl(i)
enddo

if (AreMCorrFuncNeeded.or.AreMPartDensNeeded)  then
  !Compute tvkinv_tAk
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+tvk(j)*inv_tAk(j,i)
    enddo
    tvkinv_tAk(i)=temp1
  enddo

  !Compute inv_tAltvl
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_tAl(j,i)*tvl(j)
    enddo
    inv_tAltvl(i)=temp1
  enddo

  !Compute inv_invtAkinvtAlinv_tAltvl
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_invtAkinvtAl(j,i)*inv_tAltvl(j)
    enddo
    inv_invtAkinvtAlinv_tAltvl(i)=temp1
  enddo

  !Compute tvkinv_tAkinv_invtAkinvtAl
  do i=1,n
    temp1=ZERO
    do j=1,n
      temp1=temp1+inv_invtAkinvtAl(j,i)*tvkinv_tAk(j)
    enddo
    tvkinv_tAkinv_invtAkinvtAl(i)=temp1
  enddo

  !Compute tau4 = tvkinv_tAk*inv_invtAkinvtAl*inv_tAltvl
  tau4=ZERO
  do i=1,n
    do j=1,n
      tau4=tau4+tvkinv_tAk(i)*inv_invtAkinvtAl(i,j)*inv_tAltvl(j)
    enddo
  enddo
end if

!Evaluating overlap
!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
temp1=det_tAkl*sqrt(det_tAkl)
Skl=Glob_Piraised3n2*tau3/(TWO*temp1)


if(AreMCorrFuncNeeded.or.AreMPartDensNeeded) then
  temp1=det_tAk*det_tAl*det_invtAkinvtAl
  temp2=temp1*sqrt(temp1)
  MSkl=Glob_Piraised3n2*tau4/(TWO*temp2)
endif

!Doing multiplication inv_tAkltAl=inv_tAkl*tAl, inv_tAkltAk=inv_tAkl*tAk
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
    enddo
    inv_tAkltAl(j,i)=temp1
  enddo
enddo

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkltAl(j,k)*Glob_MassMatrix(k,i)
    enddo
    inv_tAkltAlM(j,i)=temp1
  enddo
enddo

!Computing tau1=tr[inv_tAkltAlM*tAk]
tau1=ZERO
do i=1,n
  temp1=ZERO
  do k=1,n
    temp1=temp1+inv_tAkltAlM(i,k)*tAk(k,i)
  enddo
  tau1=tau1+temp1
enddo
!Computing tau2 = tvk'*inv_tAkltAlM*tAk*inv_tAkltvl
!We do it by multiplying twice the row-vector on the left
!by a matrix on the right and computing a dot product in the end.

!tvkinv_tAkltAlM'=tvk'*inv_tAkltAlM
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAkltAlM(j,i)
  enddo
  tvkinv_tAkltAlM(i)=temp1
enddo

!u1=tvkinv_tAkltAlM'*tAk
!tau2=u1'*inv_tAkltvl (storage for u1 as such is not needed, we use temp1=u1(i))
tau2=ZERO
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkltAlM(j)*tAk(j,i)
  enddo
  tau2=tau2+temp1*inv_tAkltvl(i)
enddo

!Evaluating the kinetic energy
Tkl=Skl*(SIX*tau1+FOUR*tau2/tau3)

!Evaluating TrAJ(i,j)=Tr[inv_tAkl*Jij], sqrtTrAJ(i,j)=sqrt(TrAJ(i,j)),
!           eta2(i,j)=tvk'*inv_tAkl*Jij*inv_tAkl*tvl,
!Vkl, (1/r_{ij}^2)_kl, (1/r_{ij})_kl, (r_{ij})_kl, (r_{ij}^2)_kl
!and delta(r_{ij})_kl
Vkl=ZERO
temp5=Skl*TWO
temp1=temp5/SQRTPI
temp8=Skl/(PI*SQRTPI)
do i=1,n
  temp2=inv_tAkl(i,i)
  TrAJ(i,i)=temp2
  temp3=sqrt(temp2)
  sqrtTrAJ(i,i)=temp3
  !u1'=tvk'*inv_tAkl*Jii*inv_tAkl
  do q=1,n
    temp4=ZERO
    do k=1,n
      temp4=temp4+tvk(k)*inv_tAkl(k,i)*inv_tAkl(q,i)
    enddo
    u1(q)=temp4
  enddo
  !eta2=u1'*tvl
  temp4=ZERO
  do k=1,n
    temp4=temp4+u1(k)*tvl(k)
  enddo
  eta2(i,i)=temp4
  temp7=temp4/(temp2*tau3)
  temp6=temp7/THREE
  temp9=ONE-temp7
  rm2kl(i,i)=temp5*(ONE-2*temp6)/temp2
  rmkl(i,i)=temp1*(ONE-temp6)/temp3
  Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*rmkl(i,i)
  rkl(i,i)=temp1*(ONE+temp6)*temp3
  r2kl(i,i)=Skl*THREEHALF*(ONE+2*temp6)*temp2
  temp10=temp8/(temp2*temp3)
  deltarkl(i,i)=temp10*temp9
  prvalkl(i,i)=PI*temp10*( TWO*temp9*(Glob_EulerConst+log(temp2)) + temp7*FOUR/THREE )
  !prvalkl(i,i)=(temp1/(temp2*temp3))*temp9*(Glob_EulerConst+log(temp2))+(FOUR*Skl/3)*temp7/(temp2*temp3)
  !prvalkl(i,i)=(temp9*(Glob_EulerConst+log(temp2))/SQRTPI+temp7*TWO/THREE)*temp5/(temp2*temp3)
  !prvalkl(i,i)=(temp1/(temp2*temp3))*(1-temp7)*(Glob_EulerConst+log(temp2))
  !prvalkl(i,i)=TWO*PI*(Glob_EulerConst+log(temp2))*deltarkl(i,i)
enddo
do i=1,n
  do j=i+1,n
    temp2=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    TrAJ(i,j)=temp2
    TrAJ(j,i)=temp2
    temp3=sqrt(temp2)
    sqrtTrAJ(j,i)=temp3
    sqrtTrAJ(i,j)=temp3
    !u1'=tvk'*inv_tAkl*Jij*inv_tAkl
    do q=1,n
      temp4=ZERO
	  temp10=ZERO
      do k=1,n
        temp4=temp4+tvk(k)*(inv_tAkl(k,i)-inv_tAkl(k,j))*(inv_tAkl(q,i)-inv_tAkl(q,j))
      enddo
      u1(q)=temp4
    enddo
    temp4=ZERO
    do k=1,n
      temp4=temp4+u1(k)*tvl(k)
    enddo
    eta2(j,i)=temp4
    eta2(i,j)=temp4
    temp7=temp4/(temp2*tau3)
    temp6=temp7/THREE
    temp9=ONE-temp7
    rm2kl(j,i)=temp5*(ONE-2*temp6)/temp2
    rm2kl(i,j)=rm2kl(j,i)
    rmkl(j,i)=temp1*(ONE-temp6)/temp3
    rmkl(i,j)=rmkl(j,i)
    Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*rmkl(j,i)
    rkl(j,i)=temp1*(ONE+temp6)*temp3
    rkl(i,j)=rkl(j,i)
    r2kl(j,i)=Skl*THREEHALF*(ONE+2*temp6)*temp2
    r2kl(i,j)=r2kl(j,i)
    temp10=temp8/(temp2*temp3)
    deltarkl(j,i)=temp10*temp9
    deltarkl(i,j)=deltarkl(j,i)
    prvalkl(j,i)=PI*temp10*( TWO*temp9*(Glob_EulerConst+log(temp2)) + temp7*FOUR/THREE )
    !prvalkl(j,i)=(temp1/(temp2*temp3))*temp9*(Glob_EulerConst+log(temp2))+(FOUR*Skl/3)*temp7/(temp2*temp3)
    !prvalkl(j,i)=(temp9*(Glob_EulerConst+log(temp2))/SQRTPI+temp7*TWO/THREE)*temp5/(temp2*temp3)
    !prvalkl(j,i)=(temp1/(temp2*temp3))*(1-temp7)*(Glob_EulerConst+log(temp2))
    !prvalkl(j,i)=TWO*PI*(Glob_EulerConst+log(temp2))*deltarkl(j,i)
    prvalkl(i,j)=prvalkl(j,i)
  enddo
enddo

if (AreMCorrFuncNeeded) then
  do i=1,n
    MTrAJ(i,i)=inv_invtAkinvtAl(i,i)
    sqrtMTrAJ(i,i)=sqrt(MTrAJ(i,i))
    do q=1,n
      temp4=ZERO
      do k=1,n
        temp4=temp4+tvkinv_tAk(k)*inv_invtAkinvtAl(k,i)*inv_invtAkinvtAl(q,i)
      enddo
      u1(q)=temp4
    enddo
    !eta2=u1'*tvl
    temp4=ZERO
    do k=1,n
      temp4=temp4+u1(k)*inv_tAltvl(k)
    enddo
    meta2(i,i)=temp4
  enddo
  do i=1,n
    do j=i+1,n
      MTrAJ(i,j)=inv_invtAkinvtAl(i,i)+inv_invtAkinvtAl(j,j)-inv_invtAkinvtAl(j,i)-inv_invtAkinvtAl(j,i)
      MTrAJ(j,i)=MTrAJ(i,j)
      sqrtMTrAJ(j,i)=sqrt(MTrAJ(j,i))
      sqrtMTrAJ(i,j)=sqrtMTrAJ(j,i)
    enddo
  enddo
end if

Hkl=Tkl+Vkl

!Evaluating tr[inv_tAkl Jij inv_tAkl Jpq] and
!j^{ij}' inv_tAkl j^{p,q}  where j^{ij}=e^i-e^j and
!the only nonzero element of e^i is the i-th element
do i=1,n
  temp2=inv_tAkl(i,i)
  temp1=temp2*temp2
  TrAJAJ(i,i,i,i)=temp1
  jAj(i,i,i,i)=temp2
  do p=i+1,n
    do q=p+1,n
      temp2=inv_tAkl(p,i)-inv_tAkl(q,i)
      temp1=temp2*temp2
      TrAJAJ(i,i,p,q)=temp1
      TrAJAJ(i,i,q,p)=temp1
      TrAJAJ(p,q,i,i)=temp1
      TrAJAJ(q,p,i,i)=temp1
      jAj(i,i,p,q)=temp2
      jAj(i,i,q,p)=-temp2
      jAj(p,q,i,i)=temp2
      jAj(q,p,i,i)=-temp2
    enddo
  enddo
  do j=i+1,n
    temp2=inv_tAkl(j,i)
    temp1=temp2*temp2
    TrAJAJ(j,j,i,i)=temp1
    TrAJAJ(i,i,j,j)=temp1
    jAj(j,j,i,i)=temp2
    jAj(i,i,j,j)=temp2
    do p=i,n
      temp2=inv_tAkl(p,i)-inv_tAkl(p,j)
      temp1=temp2*temp2
      TrAJAJ(i,j,p,p)=temp1
      TrAJAJ(j,i,p,p)=temp1
      TrAJAJ(p,p,i,j)=temp1
      TrAJAJ(p,p,j,i)=temp1
      jAj(i,j,p,p)=temp2
      jAj(j,i,p,p)=-temp2
      jAj(p,p,i,j)=temp2
      jAj(p,p,j,i)=-temp2
      do q=p+1,n
        temp2=inv_tAkl(p,i)-inv_tAkl(q,i)-inv_tAkl(p,j)+inv_tAkl(q,j)
        temp1=temp2*temp2
        TrAJAJ(i,j,p,q)=temp1
        TrAJAJ(j,i,p,q)=temp1
        TrAJAJ(i,j,q,p)=temp1
        TrAJAJ(j,i,q,p)=temp1
        TrAJAJ(p,q,i,j)=temp1
        TrAJAJ(p,q,j,i)=temp1
        TrAJAJ(q,p,i,j)=temp1
        TrAJAJ(q,p,j,i)=temp1
        jAj(i,j,p,q)=temp2
        jAj(j,i,p,q)=-temp2
        jAj(i,j,q,p)=-temp2
        jAj(j,i,q,p)=temp2
        jAj(p,q,i,j)=temp2
        jAj(p,q,j,i)=-temp2
        jAj(q,p,i,j)=-temp2
        jAj(q,p,j,i)=temp2
      enddo
    enddo
  enddo
enddo

!Evaluating vector-matrix-vector products
!j^{ij}' inv_tAkl tvl
!tvk' inv_tAkl j^{ij}
do j=1,n
  do i=1,n
    if (i==j) then
      jAtvl(i,i)=inv_tAkltvl(i)
      tvkAj(i,i)=tvkinv_tAkl(i)
    else
      jAtvl(i,j)=inv_tAkltvl(i)-inv_tAkltvl(j)
      tvkAj(i,j)=tvkinv_tAkl(i)-tvkinv_tAkl(j)
    endif
  enddo
enddo

!!Evaluating tr [inv_tAkl J_{ij} inv_tAkl J_{pq} inv_tAkl V] where V = tvk tvl'
!!and tr [inv_tAkl J_{ij} inv_tAkl J_{pq} inv_tAkl J_{ij} inv_tAkl V]
!do i=1,n
! do j=i,n
!   do p=1,n
!     do q=p,n
!       temp1=tvkAj(i,j)*jAj(i,j,p,q)*jAtvl(p,q)
!       temp2=tvkAj(i,j)*jAj(i,j,p,q)*jAj(p,q,i,j)*jAtvl(i,j)
!       TrAJAJAV(i,j,p,q)=temp1
!       TrAJAJAV(j,i,p,q)=temp1
!       TrAJAJAV(i,j,q,p)=temp1
!       TrAJAJAV(j,i,q,p)=temp1
!       omega(i,j,p,q)=temp2
!       omega(j,i,p,q)=temp2
!       omega(i,j,q,p)=temp2
!       omega(j,i,q,p)=temp2
!     enddo
!   enddo
! enddo
!enddo

!Evaluation of (1/r_{ij}*1/r_{pq}))_kl is not implemented yet
temp1=4*Skl/(3*PI)
do i=1,n
  do j=i,n
    do p=i,n
      do q=p,n   !try q=max(p,j),n - it may speed things up a little
        if (((p==i).and.(q==j)).or.((p==j).and.(q==i))) then
          temp2=rm2kl(i,j)
          rmrmkl(i,j,p,q)=temp2
          rmrmkl(j,i,p,q)=temp2
          rmrmkl(i,j,q,p)=temp2
          rmrmkl(j,i,q,p)=temp2
          rmrmkl(p,q,i,j)=temp2
          rmrmkl(p,q,j,i)=temp2
          rmrmkl(q,p,i,j)=temp2
          rmrmkl(q,p,j,i)=temp2
        else
          !jAj(nn,nn,nn,nn),jAtvl(nn,nn),tvkAj(nn,nn)
          temp2=TrAJ(i,j)          !sigma    a
          temp3=TrAJ(p,q)          !chi      b
          temp4=jAj(p,q,i,j)       !         c
          !tau3                    !zeta     w
          temp5=tvkAj(i,j)         !         f
          temp6=tvkAj(p,q)         !         g
          temp7=jAtvl(i,j)         !         u
          temp8=jAtvl(p,q)         !         v
          temp9 = sqrt(temp2*temp3)
          temp10 = temp4/temp9
          threshold = ONE - temp10*temp10
          if ( threshold  < TEN * EPSILON(threshold) )then
            ! employed expression
            ! expansion_term_2 = 2*Skl* (Sqrt(a*b)*g*u + Sqrt(a*b)*f*v - 3*a*b*tau3) / (3*a*b*Sqrt(a*b)*tau3)
            temp11 = temp9 * temp6 * temp7                                                              ! Sqrt(a*b)*g*u
            temp12 = temp9 * temp5 * temp8                                                              ! Sqrt(a*b)*f*v
            temp13 = THREE * temp2 * temp3 * tau3                                                       ! 3*a*b*tau3
            temp14 = (TWO * Skl) * (temp11 + temp12 - temp13) / (temp13 * temp9)
          else
            temp11=sqrt(ONE-temp10*temp10)                                                                              ! sqrt(1-x^2)
            temp12=(temp6*temp7+temp5*temp8)/(temp2*temp3*tau3)                                                         ! h
            temp13=(temp2*temp6*temp8+temp3*temp5*temp7)/(temp2*temp3*temp9*tau3)                                       ! t
            temp14=(temp1/temp11)*(THREE/temp9-temp13+(temp12-temp10*THREE/temp9)*ftransaux(temp10))
          endif
          rmrmkl(i,j,p,q)=temp14
          rmrmkl(j,i,p,q)=temp14
          rmrmkl(i,j,q,p)=temp14
          rmrmkl(j,i,q,p)=temp14
          rmrmkl(p,q,i,j)=temp14
          rmrmkl(p,q,j,i)=temp14
          rmrmkl(q,p,i,j)=temp14
          rmrmkl(q,p,j,i)=temp14
        endif
      enddo
    enddo
  enddo
enddo

!write(*,*) '===================================='
!write(*,'(1x,a,e23.16)') 'Skl=',Skl
!write(*,*) 'tvk='
!do i=1,n
!    write(*,'(2x,i3)') tvk(i)
!enddo
!write(*,*) 'tvl='
!do i=1,n
!    write(*,'(2x,i3)') tvl(i)
!enddo
!write(*,*) 'tAk='
!do i=1,n
!  do j=1,n
!    write(*,'(2x,f23.20)',advance='no') tAk(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) 'tAl='
!do i=1,n
!  do j=1,n
!    write(*,'(2x,f23.20)',advance='no') tAl(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) 'inv_tAkl='
!do i=1,n
!  do j=1,n
!    write(*,'(2x,f23.20)',advance='no') inv_tAkl(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) 'M='
!do i=1,n
!  do j=1,n
!    write(*,'(2x,f23.20)',advance='no') Glob_MassMatrix(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) '1/(rijrpq)='
!do i=1,n
!  do j=1,n
!    do p=1,n
!      do q=1,n
!        write(*,'(2x,i3,i3,i3,i3,2x,e23.16)') i,j,p,q,rmrmkl(i,j,p,q)
!      enddo
!    enddo
!  enddo
!enddo
!write(*,*) 'ME_d_X_over_rij_d='

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
	temp4=ZERO
	temp5=ZERO
    if (p==q) then
      temp4=2*PI*Glob_MassMatrix(p,p)
      temp5=ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(p))
    else
      temp4=2*PI*(Glob_MassMatrix(p,p)+Glob_MassMatrix(q,q) &
        -Glob_MassMatrix(p,q)-Glob_MassMatrix(p,q))
      temp5=ScaledChargeProd(Glob_PseudoCharge(p),Glob_PseudoCharge(q))
    endif

    !temp2=ME_rXr_over_rij(W2,p,q,inv_tAkl,rmkl(p,q),TrAJ(p,q))
    temp2=ZERO
    temp2=ME_d_X_over_rij_d(p,q,Glob_dmvM,tAk,tAl,inv_tAkl,tvk,tvl, &
    inv_tAkltvl,tvkinv_tAkl,trAJ(p,q),tau3,eta2(p,q),Skl)
    drach_deltarkl(p,q)=(Glob_CurrEnergy*rmkl(p,q)-temp1-temp2)/temp4
    drach_deltarkl(q,p)=drach_deltarkl(p,q)
    V2kl=V2kl+temp5*temp1
  enddo
enddo

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

!Mass-velocity correction
inv_tau3=1/tau3
W1(1:n,1:n)=ONE
temp1=Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(1)
MVkl=ME_dWd2(W1,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)/temp1
W1(1:n,1:n)=ZERO
do i=1,n
  W1(i,i)=ONE
  temp1=Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1)
  MVkl=MVkl+ME_dWd2(W1,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)/temp1
  W1(i,i)=ZERO
enddo
MVkl=-MVkl/8
drach_MVkl=ZERO
drach_MVkl=ME_dWd21(Glob_dmvM,Glob_dmvMB,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl) &
  -V2kl-Glob_CurrEnergy*Glob_CurrEnergy*Skl+2*Glob_CurrEnergy*Vkl+ &
  Glob_CurrEnergy*ME_dXd(Glob_dmvB,tvk,tvl,inv_tAkltvl,inv_tAkl,tAk,tAl,inv_tAkltAl,Skl,tau3)!+&
do i=1,n
  drach_MVkl=drach_MVkl-ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(i)) &
                        *ME_d_X_over_rij_d(i,i,Glob_dmvB,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl, &
                                   tvkinv_tAkl,trAJ(i,i),tau3,eta2(i,i),Skl)
  do j=i+1,n
  drach_MVkl=drach_MVkl-ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j)) &
                        *ME_d_X_over_rij_d(i,j,Glob_dmvB,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl, &
                                   tvkinv_tAkl,trAJ(i,j),tau3,eta2(i,j),Skl)
  enddo
enddo
drach_MVkl = drach_MVkl*Glob_dmva2 + MVkl




!Evaluating Orbit-Orbit (OO) matrix element (without the factor of alpha**2)
OOkl=ZERO
!First double loop for OO
do i=1,n
  do j=1,n
    tr1=tAl(j,i)
    tr2=tAl(i,j)
    tr3=3*tAl(j,j)
    tr4=tvl(j)*tvk(j)
    !W1 = Al Eij Al + Al Ejj Eji Al
    do p=1,n
      do q=1,n
        W1(p,q)=tAl(p,i)*tAl(j,q)+tAl(p,j)*tAl(i,q)
      enddo
    enddo
    !W1 = W1 + Akl Ejj Al Eij
    do p=1,n
      W1(p,j)=W1(p,j)+tAkl(p,j)*tAl(j,i)
    enddo
    !W1 = W1 + Eji Al Ejj Al + tr3 Eji Al
    do q=1,n
      W1(j,q)=W1(j,q)+tAl(i,j)*tAl(j,q)+tr3*tAl(i,q)
    enddo
    !W2 = Akl Ejj Al
    do p=1,n
      do q=1,n
        W2(p,q)=tAkl(p,j)*tAl(j,q)
      enddo
    enddo
    !W3 = Eji Al
    W3(1:n,1:n)=ZERO
    do q=1,n
      W3(j,q)=tAl(i,q)
    enddo

    !W4 = Al(Ejj Eji + Eij Ejj)vl vk' + Al(Eji + Eij)vl vk'
    do p=1,n
      do q=1,n
        W4(p,q) = 2*tAl(p,j)*tvl(i)*tvk(q) + 2*tAl(p,i)*tvl(j)*tvk(q)
      enddo
    enddo
    !W4 = W4 + vl vk' Ejj Al Eij
    do p=1,n
      W4(p,j) = W4(p,j) + tvl(p)*tvk(j)*tAl(j,i)
    enddo
    !W4 = W4 + 2 Eji Al Ejj vl vk' + tr3 Eji vl vk' + Eji vl vk' Ejj Al
    do q=1,n
      W4(j,q) = W4(j,q) + 2*tAl(i,j)*tvl(j)*tvk(q) + tr3*tvl(i)*tvk(q) + tvl(i)*tvk(j)*tAl(j,q)
    enddo

    !W5 = Akl Ejj vl vk' + Al Ejj vl vk' + vl vk' Ejj Al
    do p=1,n
      do q=1,n
        W5(p,q) = tAkl(p,j)*tvl(j)*tvk(q) + tAl(p,j)*tvl(j)*tvk(q) + tvl(p)*tvk(j)*tAl(j,q)
      enddo
    enddo

    !W6 = Akl Ejj Al
    do p=1,n
      do q=1,n
        W6(p,q) = tAkl(p,j)*tAl(j,q)
      enddo
    enddo

    !W7 = Eji vl vk'
    W7(1:n,1:n)=ZERO
    do q=1,n
      W7(j,q) = tvl(i)*tvk(q)
    enddo

    call symmetrize_matrix(W1)
    call symmetrize_matrix(W2)
    call symmetrize_matrix(W3)
    call symmetrize_matrix(W4)
    call symmetrize_matrix(W5)
    call symmetrize_matrix(W6)
    call symmetrize_matrix(W7)
    !compute integrals
    temp1=ME_rXr_over_rij(j,j,W1,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,tau3,Skl)
    temp2=ME_rXr_rYr_over_rij(j,j,W2,W3,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,tau3,Skl)
    temp3=SG_ME_rXr_over_rij(j,j,W3,inv_tAkl,tau3,Skl)
    temp4=ONETHIRD*SG_ME_rXr_over_rij(j,j,W4,inv_tAkl,tau3,Skl)
    temp5=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W3,W5,inv_tAkl,tau3,Skl)
    temp6=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W6,W7,inv_tAkl,tau3,Skl)
    temp7=-6*(tr1+tr2)*rmkl(j,j)+4*temp1-8*temp2-2*tr4*temp3-2*temp4+4*temp5+4*temp6
    OOkl=OOkl-temp7*ScaledChargeProd(Glob_PseudoCharge(j),Glob_PseudoCharge0)/Glob_Mass(j+1)
  enddo
enddo
OOkl=OOkl/Glob_Mass(1)


!Second double loop for OO
do i=1,n
  do j=i+1,n
    tr1=tAl(j,i)
    tr2=tAl(i,j)
    tr3=3*tAl(j,j)
    tr4=tvl(j)*tvk(j)
    !W1 = Al Eji Al + Al Ejj (Eji - Eii) Al
    do p=1,n
      do q=1,n
        W1(p,q)=tAl(p,i)*tAl(j,q)+tAl(p,j)*tAl(i,q)
      enddo
    enddo
    !W1 = W1 + Akl Ejj Al (Eij - Eii)
    do p=1,n
      W1(p,j)=W1(p,j)+tAkl(p,j)*tAl(j,i)
      W1(p,i)=W1(p,i)-tAkl(p,j)*tAl(j,i)
    enddo
    !W1 = W1 + (Eji - Eii) Al Ejj Al + tr3 (Eji - Eii) Al
    do q=1,n
      temp1=tAl(i,j)*tAl(j,q)+tr3*tAl(i,q)
      W1(j,q)=W1(j,q)+temp1
      W1(i,q)=W1(i,q)-temp1
    enddo
    !W2 = Akl Ejj Al
    do p=1,n
      do q=1,n
        W2(p,q)=tAkl(p,j)*tAl(j,q)
      enddo
    enddo
    !W3 = (Eji - Eii) Al
    W3(1:n,1:n)=ZERO
    do q=1,n
      W3(j,q)=tAl(i,q)
      W3(i,q)=-tAl(i,q)
    enddo

    !W4 = Al(Ejj Eji + Eij Ejj - Ejj Eii - Eii Ejj)vl vk' + Al(Eji + Eij)vl vk'
    do p=1,n
      do q=1,n
        W4(p,q) = 2*tAl(p,j)*tvl(i)*tvk(q) + 2*tAl(p,i)*tvl(j)*tvk(q)
      enddo
    enddo
    !W4 = W4 + vl vk' Ejj Al (Eij-Eii)
    do p=1,n
      W4(p,j) = W4(p,j) + tvl(p)*tvk(j)*tAl(j,i)
      W4(p,i) = W4(p,i) - tvl(p)*tvk(j)*tAl(j,i)
    enddo
    !W4 = W4 + 2 (Eji-Eii) Al Ejj vl vk' + tr3 (Eji-Eii) vl vk' + (Eji-Eii) vl vk' Ejj Al
    do q=1,n
      temp1 = 2*tAl(i,j)*tvl(j)*tvk(q) + tr3*tvl(i)*tvk(q) + tvl(i)*tvk(j)*tAl(j,q)
      W4(j,q) = W4(j,q) + temp1
      W4(i,q) = W4(i,q) - temp1
    enddo

    !W5 = Akl Ejj vl vk' + Al Ejj vl vk' + vl vk' Ejj Al
    do p=1,n
      do q=1,n
        W5(p,q) = tAkl(p,j)*tvl(j)*tvk(q) + tAl(p,j)*tvl(j)*tvk(q) + tvl(p)*tvk(j)*tAl(j,q)
      enddo
    enddo

    !W6 = Akl Ejj Al
    do p=1,n
      do q=1,n
        W6(p,q) = tAkl(p,j)*tAl(j,q)
      enddo
    enddo

    !W7 = (Eji - Eii) vl vk'
    W7(1:n,1:n)=ZERO
    do q=1,n
      W7(j,q) = tvl(i)*tvk(q)
      W7(i,q) = -tvl(i)*tvk(q)
    enddo

    call symmetrize_matrix(W1)
    call symmetrize_matrix(W2)
    call symmetrize_matrix(W3)
    call symmetrize_matrix(W4)
    call symmetrize_matrix(W5)
    call symmetrize_matrix(W6)
    call symmetrize_matrix(W7)

    !compute integrals
    temp1=ME_rXr_over_rij(i,j,W1,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,tau3,Skl)
    temp2=ME_rXr_rYr_over_rij(i,j,W2,W3,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,tau3,Skl)
    temp3=SG_ME_rXr_over_rij(i,j,W3,inv_tAkl,tau3,Skl)
    temp4=ONETHIRD*SG_ME_rXr_over_rij(i,j,W4,inv_tAkl,tau3,Skl)
    temp5=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W3,W5,inv_tAkl,tau3,Skl)
    temp6=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W6,W7,inv_tAkl,tau3,Skl)
    temp7=-6*(tr1+tr2)*rmkl(i,j)+4*temp1-8*temp2-2*tr4*temp3-2*temp4+4*temp5+4*temp6
    OOkl=OOkl+&
      temp7*ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))/(Glob_Mass(i+1)*Glob_Mass(j+1))
  enddo
enddo
OOkl=OOkl/2


!Evaluation of correlation functions
if (AreCorrFuncNeeded) then
  temp1=Skl/(PI*SQRTPI)
  p=0
  do i=1,n
    do j=i,n
      p=p+1
      temp2=temp1/(sqrtTrAJ(j,i)*TrAJ(j,i))
      temp3=-1/TrAJ(j,i)
      temp4=2/TrAJ(j,i)
      temp5=eta2(j,i)/(TrAJ(j,i)*tau3)
      do k=1,NumCFGridPoints
        temp6=CFGrid(2,k)*CFGrid(2,k)         !this is \xi_z^2
        temp7=temp6+CFGrid(1,k)*CFGrid(1,k)   !this is  \xi^2
        temp8=ONE+(temp4*temp6-ONE)*temp5
        CFkl(p,k)=temp2*temp8*exp(temp7*temp3)
      enddo
    enddo
  enddo
endif

if (ArePartDensNeeded) then
  temp1=Skl/(PI*SQRTPI)
  do i=1,n+1
    temp2=ZERO
    do p=1,n
      temp2=temp2+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tAkl(p,p)
      do q=p+1,n
        temp2=temp2+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tAkl(q,p)
      enddo
    enddo
    temp3=ZERO
    temp4=ZERO
    do p=1,n
      temp3=temp3+tvkinv_tAkl(p)*Glob_bvc(p,i)
      temp4=temp4+Glob_bvc(p,i)*inv_tAkltvl(p)
    enddo
    temp5=temp3*temp4/(temp2*tau3)
    temp6=-1/temp2
    temp7=2/temp2
    temp8=temp1/(sqrt(temp2)*temp2)
    do k=1,NumDensGridPoints
      temp9=DensGrid(2,k)*DensGrid(2,k)          !this is  \xi_z^2
      temp10=temp9+DensGrid(1,k)*DensGrid(1,k)   !this is -\xi^2
      temp11=ONE+(temp7*temp9-ONE)*temp5
      Denskl(i,k)=temp8*temp11*exp(temp10*temp6)
    enddo
  enddo
endif

if (AreMCorrFuncNeeded) then
  temp1=MSkl/(PI*SQRTPI)
  p=0
  do i=1,n
    do j=i,n
      p=p+1
      temp2=temp1/(sqrtMTrAJ(j,i)*MTrAJ(j,i))/8
      temp3=-1/MTrAJ(j,i)/4
      temp4=1/MTrAJ(j,i)/2
      temp5=meta2(j,i)/(MTrAJ(j,i)*tau4)
      do k=1,NumCFGridPoints
        temp6=CFGrid(2,k)*CFGrid(2,k)         !this is \xi_z^2
        temp7=temp6+CFGrid(1,k)*CFGrid(1,k)   !this is  \xi^2
        temp8=ONE+(temp4*temp6-ONE)*temp5
        CFkl(p,k)=temp2*temp8*exp(temp7*temp3)
      enddo
    enddo
  enddo
endif

if (AreMPartDensNeeded) then
  temp1=MSkl/(PI*SQRTPI)
  do i=1,n+1
    temp2=ZERO
    do p=1,n
      temp2=temp2+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(p,p)
      do q=p+1,n
        temp2=temp2+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(q,p)
      enddo
    enddo
    temp3=ZERO
    temp4=ZERO
    do p=1,n
      temp3=temp3+tvkinv_tAkinv_invtAkinvtAl(p)*Glob_bvc(p,i)
      temp4=temp4+Glob_bvc(p,i)*inv_invtAkinvtAlinv_tAltvl(p)
    enddo
    temp5=temp3*temp4/(temp2*tau4)
    temp6=-1/temp2
    temp7=2/temp2
    temp8=temp1/(sqrt(temp2)*temp2)
    do k=1,NumDensGridPoints
      temp9=DensGrid(2,k)*DensGrid(2,k)          !this is  \xi_z^2
      temp10=temp9+DensGrid(1,k)*DensGrid(1,k)   !this is -\xi^2
      temp11=ONE+(temp7*temp9/4-ONE)*temp5
      Denskl(i,k)=temp8*temp11*exp(temp10*temp6/4)/8
    enddo
  enddo
endif

end subroutine MatrixElementsL1ForExpcVals

subroutine symmetrize_matrix(W)
!subroutine symmetrize_matrix makes an arbitrary square matrix W
!symmetric by the following procedure:
!W = (1/2)*(W + W')
!Input:
!   W :: n x n real matrix

integer, parameter :: nn = Glob_MaxAllowedNumOfPseudoParticles
real(dprec)           W(nn, nn), t
integer               i,j,n

n = Glob_n

do i = 1,n
    do j = i+1,n
	    t=ONEHALF*(W(j,i)+W(i,j))
        W(j,i) = t
        W(i,j) = t
    end do
end do

end subroutine symmetrize_matrix



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


function ME_rXr(X,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_rXr computes the following matrix element:
!<\tilde phi_k| r'Xr |\tilde phi_l>
!Here X is a symmetric real matrix. If matrix X is not symmetric
!then user need to symmetrize it before calling this function.
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   ME_rXr
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),inv_tAkl(nn,nn),inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
integer       i,j,n
real(dprec)   trAX,trAXAtvltvk,temp1,workvec1(nn)

n=Glob_n
trAX=ZERO
do i=1,n
  do j=1,n
    trAX=trAX+inv_tAkl(i,j)*X(j,i)
  enddo
enddo
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkl(j)*X(j,i)
  enddo
  workvec1(i)=temp1
enddo
trAXAtvltvk=ZERO
do i=1,n
  trAXAtvltvk=trAXAtvltvk+workvec1(i)*inv_tAkltvl(i)
enddo

ME_rXr=Skl*(THREEHALF*trAX+trAXAtvltvk*inv_tau3)

end function ME_rXr


function ME_rXr_rYr(X,Y,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_rXr_rYr computes the following matrix element:
!<\tilde phi_k| r'Xr r'Yr |\tilde phi_l>
!Here X and Y are symmetric real matrices. If matrices are not symmetric
!then user needs to symmetrize them before calling this function.
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   ME_rXr_rYr
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
integer       i,j,k,n
real(dprec)   trAX,trAY,trAXAY,temp1,temp2
real(dprec)   trAXAtvltvk,trAYAtvltvk,trAXAYAtvltvk,trAYAXAtvltvk
real(dprec)   workvec1(nn),workvec2(nn),workvec3(nn),workvec4(nn)
real(dprec)   AX(nn,nn),AY(nn,nn)

n=Glob_n
trAX=ZERO
trAY=ZERO
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*X(k,i)
      temp2=temp2+inv_tAkl(j,k)*Y(k,i)
    enddo
    AX(j,i)=temp1
    AY(j,i)=temp2
  enddo
  trAX=trAX+AX(i,i)
  trAY=trAY+AY(i,i)
enddo
trAXAY=ZERO
do i=1,n
  do k=1,n
    trAXAY=trAXAY+AX(i,k)*AY(k,i)
  enddo
enddo

do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkl(j)*X(j,i)
    temp2=temp2+tvkinv_tAkl(j)*Y(j,i)
  enddo
  workvec1(i)=temp1
  workvec2(i)=temp2
enddo
trAXAtvltvk=ZERO
trAYAtvltvk=ZERO
do i=1,n
  trAXAtvltvk=trAXAtvltvk+workvec1(i)*inv_tAkltvl(i)
  trAYAtvltvk=trAYAtvltvk+workvec2(i)*inv_tAkltvl(i)
enddo

do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+workvec1(j)*AY(j,i)
    temp2=temp2+workvec2(j)*AX(j,i)
  enddo
  workvec3(i)=temp1
  workvec4(i)=temp2
enddo
trAXAYAtvltvk=ZERO
trAYAXAtvltvk=ZERO
do i=1,n
  trAXAYAtvltvk=trAXAYAtvltvk+workvec3(i)*inv_tAkltvl(i)
  trAYAXAtvltvk=trAYAXAtvltvk+workvec4(i)*inv_tAkltvl(i)
enddo

ME_rXr_rYr=Skl*( THREEHALF*THREEHALF*trAX*trAY + THREEHALF*trAXAY+    &
  (THREEHALF*(trAX*trAYAtvltvk+trAY*trAXAtvltvk)+trAXAYAtvltvk+trAYAXAtvltvk)*inv_tau3 )

end function ME_rXr_rYr


function ME_dWd2(W,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_dWd2 computes the following matrix element:
!<\tilde phi_k| (\nabla_r' W \nabla_r)^2 |\tilde phi_l>
!Here W is a real symmetric matrix. If matrix W is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!            W  :: n x n real matrix
!           tAk :: n x n real matrix where \tilde Ak is stored
!           tAl :: n x n real matrix where \tilde Al is stored
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   ME_dWd2
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   W(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn)
integer       tvk(nn),tvl(nn)
real(dprec)   inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
integer       i,j,k,n
real(dprec)   temp1,temp2,tuk(nn),tul(nn),tAkWtAk(nn,nn),tAlWtAl(nn,nn)
real(dprec)   WorkMat1(nn,nn),WorkMat2(nn,nn)
real(dprec)   workvec1(nn),workvec2(nn),workvec3(nn),workvec4(nn)
real(dprec)   trWorkMat1,trWorkMat2,trAtAkWtAk,trAtAlWtAl,trAtvltuk,trAtultvk
real(dprec)   trAtAkWtAkAtultvk,trAtAlWtAlAtvltuk,trAtultuk,temp3

n=Glob_n

!Compute  WorkMat1=W*tAk  WorkMat2=W*tAl  and their traces
trWorkMat1=ZERO
trWorkMat2=ZERO
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+W(j,k)*tAk(k,i)            !  tAk(k,i)*W(j,k)
      temp2=temp2+W(j,k)*tAl(k,i)  	  !  tAl(k,i)*W(j,k)
    enddo
    WorkMat1(j,i)=temp1
    WorkMat2(j,i)=temp2
  enddo
  trWorkMat1=trWorkMat1+WorkMat1(i,i)
  trWorkMat2=trWorkMat2+WorkMat2(i,i)
enddo

!Compute tuk' = 4*tvk'*W*tAk + 6*tr[W*tAk]*tvk'
!        tul' = 4*tvl'*W*tAl + 6*tr[W*tAl]*tvl'
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*WorkMat1(j,i)
    temp2=temp2+tvl(j)*WorkMat2(j,i)
  enddo
  tuk(i)=FOUR*temp1+trWorkMat1*tvk(i)*6
  tul(i)=FOUR*temp2+trWorkMat2*tvl(i)*6
enddo

!Compute tAkWtAk=tAk*W*tAk  tAlWtAl=tAl*W*tAl
do i=1,n
  do j=1,i
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+tAk(i,k)*WorkMat1(k,j)
      temp2=temp2+tAl(i,k)*WorkMat2(k,j)
    enddo
    tAkWtAk(i,j)=temp1
    tAkWtAk(j,i)=temp1
    tAlWtAl(i,j)=temp2
    tAlWtAl(j,i)=temp2
  enddo
enddo

!Compute trAtAkWtAk=tr[inv_tAkl*tAk*W*tAk]
!        trAtAlWtAl=tr[inv_tAkl*tAl*W*tAl]
trAtAkWtAk=ZERO
trAtAlWtAl=ZERO
do i=1,n
  do j=1,n
    trAtAkWtAk=trAtAkWtAk+inv_tAkl(i,j)*tAkWtAk(j,i)
    trAtAlWtAl=trAtAlWtAl+inv_tAkl(i,j)*tAlWtAl(j,i)
  enddo
enddo

!Compute trAtvltuk=tr[inv_tAkl*tvl*tuk], trAtultvk=tr[inv_tAkl*tul*tvk]
trAtvltuk=ZERO
trAtultvk=ZERO
do i=1,n
  trAtvltuk=trAtvltuk+tuk(i)*inv_tAkltvl(i)
  trAtultvk=trAtultvk+tvkinv_tAkl(i)*tul(i)
enddo

!Compute
!trAtAkWtAkAtultvk = tr[inv_tAkl*tAkWtAk*inv_tAkl*tul*tvk'] = tvkinv_tAkl'*tAkWtAk*inv_tAkl*tul
!trAtAlWtAlAtvltuk = tr[inv_tAkl*tAlWtAl*inv_tAkl*tvl*tuk'] = tuk'*inv_tAkl*tAlWtAl*inv_tAkltvl
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkl(j)*tAkWtAk(j,i)
    temp2=temp2+tuk(j)*inv_tAkl(j,i)
  enddo
  workvec1(i)=temp1
  workvec2(i)=temp2
enddo
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+workvec1(j)*inv_tAkl(j,i)
    temp2=temp2+workvec2(j)*tAlWtAl(j,i)
  enddo
  workvec3(i)=temp1
  workvec4(i)=temp2
enddo
trAtAkWtAkAtultvk=ZERO
trAtAlWtAlAtvltuk=ZERO
do i=1,n
  trAtAkWtAkAtultvk=trAtAkWtAkAtultvk+workvec3(i)*tul(i)
  trAtAlWtAlAtvltuk=trAtAlWtAlAtvltuk+workvec4(i)*inv_tAkltvl(i)
enddo

!Compute trAtultuk=tr[inv_tAkl*tul*tuk']
trAtultuk=ZERO
do i=1,n
  trAtultuk=trAtultuk+workvec2(i)*tul(i)
enddo

ME_dWd2= 16*ME_rXr_rYr(tAkWtAk,tAlWtAl,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)  &
        +Skl*inv_tau3*(trAtultuk-SIX*(trAtAkWtAk*trAtultvk +trAtAlWtAl*trAtvltuk) &
           -FOUR*(trAtAkWtAkAtultvk+trAtAlWtAlAtvltuk))

end function ME_dWd2



function ME_dWd21(X,Glob_B,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_dWd2 computes the following matrix element:
!<\tilde phi_k| (\nabla_r' X \nabla_r)^2 |\tilde phi_l>
!Here X is a real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!             X :: n x n real matrix (Glob_MassMatrix)
!        Glob_B :: n x n real matrix
!           tAk :: n x n real matrix where \tilde Ak is stored
!           tAl :: n x n real matrix where \tilde Al is stored
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   ME_dWd21
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),Glob_B(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn)
integer       tvk(nn),tvl(nn)
real(dprec)   inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
integer       i,j,k,n
real(dprec)   temp1,temp2,tuk(nn),tul(nn),tAkXtAk(nn,nn),tAlXtAl(nn,nn)
real(dprec)   WorkMat1(nn,nn),WorkMat2(nn,nn)
real(dprec)   workvec1(nn),workvec2(nn),workvec3(nn),workvec4(nn)
real(dprec)   trWorkMat1,trWorkMat2,trAtAkXtAk,trAtAlXtAl,trAtvltuk,trAtultvk
real(dprec)   trAtAkXtAkAtultvk,trAtAlXtAlAtvltuk,trAtultuk,temp3

n=Glob_n

!Compute  WorkMat1=X*tAk  WorkMat2=X*tAl  and their traces
trWorkMat1=ZERO
trWorkMat2=ZERO
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+X(j,k)*tAk(k,i)
      temp2=temp2+Glob_B(j,k)*tAl(k,i)
    enddo
    WorkMat1(j,i)=temp1
    WorkMat2(j,i)=temp2
  enddo
  trWorkMat1=trWorkMat1+WorkMat1(i,i)
  trWorkMat2=trWorkMat2+WorkMat2(i,i)
enddo

!Compute tuk' = 4*tvk'*X*tAk + 6*tr[X*tAk]*tvk'
!        tul' = 4*tvl'*X*tAl + 6*tr[X*tAl]*tvl'

do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*WorkMat1(j,i)
	temp2=temp2+tvl(j)*WorkMat2(j,i)
  enddo
  tuk(i)=FOUR*temp1+trWorkMat1*tvk(i)*6
  tul(i)=FOUR*temp2+trWorkMat2*tvl(i)*6
enddo

!Compute tAkXtAk=tAk*X*tAk  tAlWtAl=tAl*X*tAl
do i=1,n
  do j=1,i
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+tAk(i,k)*WorkMat1(k,j)
      temp2=temp2+tAl(i,k)*WorkMat2(k,j)
    enddo
    tAkXtAk(i,j)=temp1
    tAkXtAk(j,i)=temp1
    tAlXtAl(i,j)=temp2
    tAlXtAl(j,i)=temp2
  enddo
enddo

!Compute trAtAkXtAk=tr[inv_tAkl*tAk*X*tAk]
!        trAtAlXtAl=tr[inv_tAkl*tAl*X*tAl]
trAtAkXtAk=ZERO
trAtAlXtAl=ZERO
do i=1,n
  do j=1,n
    trAtAkXtAk=trAtAkXtAk+inv_tAkl(i,j)*tAkXtAk(j,i)
    trAtAlXtAl=trAtAlXtAl+inv_tAkl(i,j)*tAlXtAl(j,i)
  enddo
enddo

!Compute trAtvltuk=tr[inv_tAkl*tvl*tuk], trAtultvk=tr[inv_tAkl*tul*tvk]
trAtvltuk=ZERO
trAtultvk=ZERO
do i=1,n
  trAtvltuk=trAtvltuk+tuk(i)*inv_tAkltvl(i)
  trAtultvk=trAtultvk+tvkinv_tAkl(i)*tul(i)
enddo

!Compute
!trAtAkXtAkAtultvk = tr[inv_tAkl*tAkXtAk*inv_tAkl*tul*tvk'] = tvkinv_tAkl'*tAkXtAk*inv_tAkl*tul
!trAtAlXtAlAtvltuk = tr[inv_tAkl*tAlXtAl*inv_tAkl*tvl*tuk'] = tuk'*inv_tAkl*tAlXtAl*inv_tAkltvl
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkl(j)*tAkXtAk(j,i)
    temp2=temp2+tuk(j)*inv_tAkl(j,i)
  enddo
  workvec1(i)=temp1
  workvec2(i)=temp2
enddo
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+workvec1(j)*inv_tAkl(j,i)
    temp2=temp2+workvec2(j)*tAlXtAl(j,i)
  enddo
  workvec3(i)=temp1
  workvec4(i)=temp2
enddo
trAtAkXtAkAtultvk=ZERO
trAtAlXtAlAtvltuk=ZERO
do i=1,n
  trAtAkXtAkAtultvk=trAtAkXtAkAtultvk+workvec3(i)*tul(i)
  trAtAlXtAlAtvltuk=trAtAlXtAlAtvltuk+workvec4(i)*inv_tAkltvl(i)
enddo

!Compute trAtultuk=tr[inv_tAkl*tul*tuk']
trAtultuk=ZERO
do i=1,n
  trAtultuk=trAtultuk+workvec2(i)*tul(i)
enddo

ME_dWd21= 16*ME_rXr_rYr(tAkXtAk,tAlXtAl,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)  &
        +Skl*inv_tau3*(trAtultuk-SIX*(trAtAkXtAk*trAtultvk +trAtAlXtAl*trAtvltuk) &
           -FOUR*(trAtAkXtAkAtultvk+trAtAlXtAlAtvltuk))

end function ME_dWd21





function ME_dXd(X,tvk,tvl,inv_tAkltvl,inv_tAkl,tAk,tAl,inv_tAkltAl,Skl,tau3)
real(dprec)   ME_dXd
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
real(dprec)   X(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn),inv_tAkltAl(nn,nn)
integer       i,j,n,k,tvk(nn),tvl(nn)
real(dprec)   inv_tAkltAlX(nn,nn),inv_tAkltAlXtAk(nn,nn),tvkinv_tAkltAlX(nn),inv_tAkltvl(nn)
real(dprec)   temp1, Skl,tau,tau1,tau2,tau3
n=Glob_n
do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkltAl(j,k)*X(k,i)
    enddo
    inv_tAkltAlX(j,i)=temp1
  enddo
enddo
tau1=ZERO
do i=1,n
  temp1=ZERO
  do k=1,n
    temp1=temp1+inv_tAkltAlX(i,k)*tAk(k,i)
  enddo
  tau1=tau1+temp1
enddo
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAkltAlX(j,i)
  enddo
  tvkinv_tAkltAlX(i)=temp1
enddo
tau2=ZERO
do i=1,n
  temp1=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkltAlX(j)*tAk(j,i)
  enddo
  tau2=tau2+temp1*inv_tAkltvl(i)
enddo
ME_dXd=Skl*(SIX*tau1+FOUR*tau2/tau3)
end function ME_dXd





function SG_ME_rXr_over_rij(i,j,X,inv_tAkl,t_V,Skl)
!function SG_ME_rXr_over_rij computes the following matrix element:
!<\tilde psi_k| (r' X r)/r_ij |\tilde psi_l>
!where psi_k = exp[-(r' Ak r)] is a Simple Gaussian wavefunction
!Here X is a some real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   SG_ME_rXr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),inv_tAkl(nn,nn)
integer       i,j
real(dprec)   t_V,Skl
!Local variables:
integer       p,q,n
real(dprec)   temp1,temp2,temp3
real(dprec)   Aj(nn),AjX(nn)
real(dprec)   t_J,t_X,t_XJ

n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
if (i==j) then
 do p=1,n
   Aj(p)=inv_tAkl(p,i)
 enddo
else
 do p=1,n
   Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
 enddo
endif

!Compute AjX'=Aj'*X
do p=1,n
  temp1=ZERO
  do q=1,n
    temp1=temp1+Aj(q)*X(q,p)
  enddo
  AjX(p)=temp1
enddo

!Compute t_J=tr[inv_tAkl*Jij]
if (i==j) then
  t_J=inv_tAkl(i,i)
else
  t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
endif

!Compute t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
t_XJ=ZERO
do p=1,n
  t_XJ=t_XJ+AjX(p)*Aj(p)
enddo

!Compute t_X=tr[inv_tAkl*X]
t_X=ZERO
do p=1,n
  do q=1,n
    t_X=t_X+inv_tAkl(q,p)*X(q,p)
  enddo
enddo

temp1=(TWO/SQRTPI)*Skl
temp2=1/t_V
temp3=1/t_J
SG_ME_rXr_over_rij=temp1*temp2*temp3*sqrt(temp3)*(THREE*t_J*t_X - t_XJ)

end function SG_ME_rXr_over_rij


function SG_ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,t_V,Skl)
!function ME_rXr_rYr_over_rij computes the following matrix element:
!<\tilde phi_k| (r' X r)(r' Y r)/r_ij |\tilde phi_l>
!Here X and Y are some real symmetric matrices. If matrices X or Y are not symmetric
!then user needs to symmetrize them before calling this function.
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   SG_ME_rXr_rYr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn)
integer       i,j
real(dprec)   t_V,Skl
!Local variables:
integer       p,q,s,n
real(dprec)   temp1,temp2,temp3,temp4,temp5,temp6
real(dprec)   AX(nn,nn),AY(nn,nn)
real(dprec)   Aj(nn),AjX(nn),AjY(nn),AXAj(nn),AYAj(nn)
real(dprec)   t_J,t_X,t_Y
real(dprec)   t_XJ,t_YJ,t_XY
real(dprec)   t_XYJ,t_YXJ

n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
if (i==j) then
 do p=1,n
   Aj(p)=inv_tAkl(p,i)
 enddo
else
 do p=1,n
   Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
 enddo
endif

!Compute AjX'=Aj'*X
!    and AjY'=Aj'*Y
do p=1,n
  temp1=ZERO
  temp2=ZERO
  do q=1,n
    temp1=temp1+Aj(q)*X(q,p)
    temp2=temp2+Aj(q)*Y(q,p)
  enddo
  AjX(p)=temp1
  AjY(p)=temp2
enddo

!Compute AX=inv_tAkl*X  t_X=tr[inv_tAkl*X]
!        AY=inv_tAkl*Y  t_Y=tr[inv_tAkl*Y]
t_X=ZERO
t_Y=ZERO
do p=1,n
  do q=1,n
    temp1=ZERO
    temp2=ZERO
    do s=1,n
      temp1=temp1+inv_tAkl(s,q)*X(p,s)
      temp2=temp2+inv_tAkl(s,q)*Y(p,s)
    enddo
    AX(q,p)=temp1
    AY(q,p)=temp2
  enddo
  t_X=t_X+AX(p,p)
  t_Y=t_Y+AY(p,p)
enddo

!Compute t_J=tr[inv_tAkl*Jij]
if (i==j) then
  t_J=inv_tAkl(i,i)
else
  t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
endif

!Compute t_XY=tr[inv_tAkl*X*inv_tAkl*Y]=tr[AX*AY]
!        AXAj=AX*Aj
!        AYAj=AY*Aj
t_XY=ZERO
do p=1,n
  temp3=ZERO
  temp4=ZERO
  do q=1,n
    t_XY=t_XY+AX(p,q)*AY(q,p)
    temp3=temp3+AX(p,q)*Aj(q)
    temp4=temp4+AY(p,q)*Aj(q)
  enddo
  AXAj(p)=temp3
  AYAj(p)=temp4
enddo

!Compute
!t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_YJ=tr[inv_tAkl*Y*inv_tAkl*Jij]=AjY'*Aj
!t_XYJ=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij]=AjX'*AYAj

t_XJ=ZERO
t_YJ=ZERO
t_XYJ=ZERO


do p=1,n
  t_XJ=t_XJ+AjX(p)*Aj(p)
  t_YJ=t_YJ+AjY(p)*Aj(p)
  t_XYJ=t_XYJ+AjX(p)*AYAj(p)
enddo

!Compute t_YXJ=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij]
t_YXJ=t_XYJ

temp1=(TWO/SQRTPI)*Skl
temp2=1/t_V
temp3=1/t_J
SG_ME_rXr_rYr_over_rij=THREE*temp1*temp2*temp3*sqrt(temp3)*(  &
THREEHALF*t_J*t_X*t_Y - ONEHALF*(t_Y*t_XJ + t_X*t_YJ) + &
t_J*t_XY  - ONETHIRD*(t_XYJ + t_YXJ) + ONEHALF*temp3*t_XJ*t_YJ&
)

end function SG_ME_rXr_rYr_over_rij


function ME_rXr_over_rij(i,j,X,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,t_V,Skl)
!function ME_rXr_over_rij computes the following matrix element:
!<\tilde phi_k| (r' X r)/r_ij |\tilde phi_l>
!Here X is a some real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   ME_rXr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),inv_tAkl(nn,nn)
integer       i,j,tvk(nn),tvl(nn)
real(dprec)   inv_tAkltvl(nn),tvkinv_tAkl(nn),t_V,Skl
!Local variables:
integer       p,q,n
real(dprec)   Ajtvl,temp1,temp2,temp3
real(dprec)   Aj(nn),AjX(nn)
real(dprec)   t_J,t_X,t_JV,t_XJ,t_XV,t_JXV,t_XJV

n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
if (i==j) then
 do p=1,n
   Aj(p)=inv_tAkl(p,i)
 enddo
else
 do p=1,n
   Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
 enddo
endif

!Compute AjX'=Aj'*X
do p=1,n
  temp1=ZERO
  do q=1,n
    temp1=temp1+Aj(q)*X(q,p)
  enddo
  AjX(p)=temp1
enddo

!Compute t_J=tr[inv_tAkl*Jij]
if (i==j) then
  t_J=inv_tAkl(i,i)
else
  t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
endif

!Compute Ajtvl=Aj'*tvl
!        t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!        t_JV=tr[inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvk'*Aj*Ajtvl
!        t_JXV=tr[inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvk'*Aj*AjX'*inv_tAkltvl
Ajtvl=ZERO
t_XJ=ZERO
temp1=ZERO
temp2=ZERO
do p=1,n
  Ajtvl=Ajtvl+Aj(p)*tvl(p)
  t_XJ=t_XJ+AjX(p)*Aj(p)
  temp1=temp1+tvk(p)*Aj(p)
  temp2=temp2+AjX(p)*inv_tAkltvl(p)
enddo
t_JV=temp1*Ajtvl
t_JXV=temp1*temp2

!Compute t_X=tr[inv_tAkl*X]
!        t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*inv_tAkltvl
!        t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*Aj*Ajtvl
t_X=ZERO
t_XV=ZERO
temp2=ZERO
do p=1,n
  temp1=ZERO
  do q=1,n
    t_X=t_X+inv_tAkl(q,p)*X(q,p)
    temp1=temp1+tvkinv_tAkl(q)*X(q,p)
  enddo
  t_XV=t_XV+temp1*inv_tAkltvl(p)
  temp2=temp2+temp1*Aj(p)
enddo
t_XJV=temp2*Ajtvl

temp1=(TWO/SQRTPI)*Skl
temp2=1/t_V
temp3=1/t_J
ME_rXr_over_rij=temp1*temp2*temp3*sqrt(temp3)*( &
THREEHALF*t_J*t_V*t_X - ONEHALF*(t_X*t_JV + t_V*t_XJ) + t_J*t_XV + ONEHALF*temp3*t_JV*t_XJ  &
- ONETHIRD*(t_JXV + t_XJV) &
)
end function ME_rXr_over_rij


function ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,t_V,Skl)
!function ME_rXr_rYr_over_rij computes the following matrix element:
!<\tilde phi_k| (r' X r)(r' Y r)/r_ij |\tilde phi_l>
!Here X and Y are some real symmetric matrices. If matrices X or Y are not symmetric
!then user needs to symmetrize them before calling this function.
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   ME_rXr_rYr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn)
integer       i,j,tvk(nn),tvl(nn)
real(dprec)   inv_tAkltvl(nn),tvkinv_tAkl(nn),t_V,Skl
!Local variables:
integer       p,q,s,n
real(dprec)   temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9
real(dprec)   AX(nn,nn),AY(nn,nn)
real(dprec)   Aj(nn),AjX(nn),AjY(nn),AXAj(nn),AYAj(nn)
real(dprec)   AXinv_tAkltvl(nn),AYinv_tAkltvl(nn),tvkinv_tAklX(nn),tvkinv_tAklY(nn)
real(dprec)   Ajtvl,tvkAj
real(dprec)   t_J,t_X,t_Y
real(dprec)   t_JV,t_XJ,t_YJ,t_XV,t_YV,t_XY
real(dprec)   t_JXV,t_JYV,t_XJV,t_YJV,t_XYV,t_YXV,t_XYJ,t_YXJ
real(dprec)   t_YJXV,t_XJYV,t_XYJV,t_YXJV,t_JXYV,t_JYXV

n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
if (i==j) then
 do p=1,n
   Aj(p)=inv_tAkl(p,i)
 enddo
else
 do p=1,n
   Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
 enddo
endif

!Compute AjX'=Aj'*X
!    and AjY'=Aj'*Y
do p=1,n
  temp1=ZERO
  temp2=ZERO
  do q=1,n
    temp1=temp1+Aj(q)*X(q,p)
    temp2=temp2+Aj(q)*Y(q,p)
  enddo
  AjX(p)=temp1
  AjY(p)=temp2
enddo

!Compute AX=inv_tAkl*X  t_X=tr[inv_tAkl*X]
!        AY=inv_tAkl*Y  t_Y=tr[inv_tAkl*Y]
t_X=ZERO
t_Y=ZERO
do p=1,n
  do q=1,n
    temp1=ZERO
    temp2=ZERO
    do s=1,n
      temp1=temp1+inv_tAkl(s,q)*X(p,s)
      temp2=temp2+inv_tAkl(s,q)*Y(p,s)
    enddo
    AX(q,p)=temp1
    AY(q,p)=temp2
  enddo
  t_X=t_X+AX(p,p)
  t_Y=t_Y+AY(p,p)
enddo

!Compute t_J=tr[inv_tAkl*Jij]
if (i==j) then
  t_J=inv_tAkl(i,i)
else
  t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
endif

!Compute t_XY=tr[inv_tAkl*X*inv_tAkl*Y]=tr[AX*AY]
!        AXinv_tAkltvl=AX*inv_tAkltvl
!        AYinv_tAkltvl=AY*inv_tAkltvl
!        AXAj=AX*Aj
!        AYAj=AY*Aj
!        tvkinv_tAklX'=tvkinv_tAkl'*X
!        tvkinv_tAklY'=tvkinv_tAkl'*Y
t_XY=ZERO
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  temp5=ZERO
  temp6=ZERO
  do q=1,n
    t_XY=t_XY+AX(p,q)*AY(q,p)
    temp1=temp1+AX(p,q)*inv_tAkltvl(q)
    temp2=temp2+AY(p,q)*inv_tAkltvl(q)
    temp3=temp3+AX(p,q)*Aj(q)
    temp4=temp4+AY(p,q)*Aj(q)
    temp5=temp5+tvkinv_tAkl(q)*X(q,p)
    temp6=temp6+tvkinv_tAkl(q)*Y(q,p)
  enddo
  AXinv_tAkltvl(p)=temp1
  AYinv_tAkltvl(p)=temp2
  AXAj(p)=temp3
  AYAj(p)=temp4
  tvkinv_tAklX(p)=temp5
  tvkinv_tAklY(p)=temp6
enddo

!Compute
!Ajtvl=Aj'*tvl
!tvkAj=tvk'*Aj
!t_JV=tr[inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvk'*Aj*Ajtvl
!t_JXV=tr[inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvk'*Aj*AjX'*inv_tAkltvl
!t_JYV=tr[inv_tAkl*Jij*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvk'*Aj*AjY'*inv_tAkltvl
!t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_YJ=tr[inv_tAkl*Y*inv_tAkl*Jij]=AjY'*Aj
!t_XYJ=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij]=AjX'*AYAj
!t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*inv_tAkltvl
!t_YV=tr[inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*inv_tAkltvl
!t_XYV=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*AYinv_tAkltvl
!t_YXV=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*AXinv_tAkltvl
!t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*Aj*Ajtvl
!t_YJV=tr[inv_tAkl*Y*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*Aj*Ajtvl
!t_XYJV=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*AYAj*Ajtvl
!t_YXJV=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*AXAj*Ajtvl
!t_XJYV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*Aj*AjY'*inv_tAkltvl
!t_YJXV=tr[inv_tAkl*Y*inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*Aj*AjX'*inv_tAkltvl
!t_JXYV=tr[inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkAj*AjX'*AYinv_tAkltvl
!t_JYXV=tr[inv_tAkl*Jij*inv_tAkl*Y*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkAj*AjY'*AXinv_tAkltvl
Ajtvl=ZERO
tvkAj=ZERO
temp1=ZERO
temp2=ZERO
temp3=ZERO
t_XJ=ZERO
t_YJ=ZERO
t_XYJ=ZERO
t_XV=ZERO
t_YV=ZERO
t_XYV=ZERO
t_YXV=ZERO
temp4=ZERO
temp5=ZERO
temp6=ZERO
temp7=ZERO
temp8=ZERO
temp9=ZERO
do p=1,n
  Ajtvl=Ajtvl+Aj(p)*tvl(p)
  tvkAj=tvkAj+tvk(p)*Aj(p)
  temp1=temp1+tvk(p)*Aj(p)
  temp2=temp2+AjX(p)*inv_tAkltvl(p)
  temp3=temp3+AjY(p)*inv_tAkltvl(p)
  t_XJ=t_XJ+AjX(p)*Aj(p)
  t_YJ=t_YJ+AjY(p)*Aj(p)
  t_XYJ=t_XYJ+AjX(p)*AYAj(p)
  t_XV=t_XV+tvkinv_tAklX(p)*inv_tAkltvl(p)
  t_YV=t_YV+tvkinv_tAklY(p)*inv_tAkltvl(p)
  t_XYV=t_XYV+tvkinv_tAklX(p)*AYinv_tAkltvl(p)
  t_YXV=t_YXV+tvkinv_tAklY(p)*AXinv_tAkltvl(p)
  temp4=temp4+tvkinv_tAklX(p)*Aj(p)
  temp5=temp5+tvkinv_tAklY(p)*Aj(p)
  temp6=temp6+tvkinv_tAklX(p)*AYAj(p)
  temp7=temp7+tvkinv_tAklY(p)*AXAj(p)
  temp8=temp8+AjX(p)*AYinv_tAkltvl(p)
  temp9=temp9+AjY(p)*AXinv_tAkltvl(p)
enddo
t_JV=temp1*Ajtvl
t_JXV=temp1*temp2
t_JYV=temp1*temp3
t_XJV=temp4*Ajtvl
t_YJV=temp5*Ajtvl
t_XYJV=temp6*Ajtvl
t_YXJV=temp7*Ajtvl
t_XJYV=temp4*temp3
t_YJXV=temp5*temp2
t_JXYV=tvkAj*temp8
t_JYXV=tvkAj*temp9

!Compute t_YXJ=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij]
t_YXJ=t_XYJ

temp1=(TWO/SQRTPI)*Skl
temp2=1/t_V
temp3=1/t_J
ME_rXr_rYr_over_rij=temp1*temp2*temp3*sqrt(temp3)*(  &
NINE*ONEFOURTH*t_J*t_V*t_X*t_Y   &
-THREE*ONEFOURTH*(t_X*t_Y*t_JV + t_V*t_Y*t_XJ + t_V*t_X*t_YJ)  &
+THREEHALF*(t_J*t_Y*t_XV + t_J*t_X*t_YV + t_J*t_V*t_XY)  &
-ONEHALF*(t_JV*t_XY + t_XV*t_YJ + t_XJ*t_YV)  &
+THREE*ONEFOURTH*temp3*(t_Y*t_JV*t_XJ+t_X*t_JV*t_YJ+t_V*t_XJ*t_YJ)  &
-FIVE*ONEFOURTH*temp3*temp3*t_JV*t_XJ*t_YJ   &
-ONEHALF*(t_X*t_JYV + t_V*t_XYJ + t_Y*t_JXV + t_Y*t_XJV + t_X*t_YJV + t_V*t_YXJ) + t_J*t_XYV + t_J*t_YXV  &
+ONEHALF*temp3*(t_XJ*t_JYV + t_JV*t_XYJ + t_YJ*t_JXV + t_YJ*t_XJV + t_XJ*t_YJV + t_JV*t_YXJ) &
-ONETHIRD*(t_YJXV + t_XYJV + t_YXJV + t_JXYV + t_JYXV + t_XJYV)  &
)
!    t_J    t_X    t_Y  (t_V is passed)
!    t_XJ   t_YJ   t_JV   t_XV   t_YV   t_XY
!    t_JXV  t_JYV  t_XJV  t_YJV  t_XYJ  t_YXJ=t_XYJ  t_XYV  t_YXV

end function ME_rXr_rYr_over_rij




function ME_d_X_over_rij_d(i,j,X,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl, &
                                   tvkinv_tAkl,tr_AJ,tr_AV,tr_AJAV,Skl)
!function ME_d_X_over_rij_d computes the following matrix element:
!<\nabla_r \tilde phi_k| (1/r_ij) X | \nabla_r \tilde phi_l>
!Here X is a some real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!           i,j :: indicies of r_ij
!            X  :: n x n real matrix
!           tAk :: n x n real matrix of tAk
!           tAl :: n x n real matrix of tAl
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!!        tr_AJ :: scalar, tr_AJ = TrAJ(i,j) = tr[inv_tAkl*Jij]
!         tr_AV :: scalar, tr_AV = tau3 = tr[inv_tAkl*tvl*tvk']
!       tr_AJAV :: scalar, tr_AJAV = eta2(i,j) = tr[inv_tAkl*Jij*inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
real(dprec)   ME_d_X_over_rij_d
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn)
integer       i,j,tvk(nn),tvl(nn)
real(dprec)   inv_tAkltvl(nn),tvkinv_tAkl(nn),tr_AJ,tr_AV,tr_AJAV,Skl
!Local variables:
integer       k,n,s,c
real(dprec)   temp1,temp2
real(dprec)   tAlX(nn,nn),tAlXtAk(nn,nn),Aj(nn),tAkX(nn,nn),tAkXtAl(nn,nn)
real(dprec)   theta,kappa,lambda,omega,chi,theta1,kappa1,chi1,omega1
real(dprec)   h,m,Xtvl(nn),tvkX(nn)

!Doing multiplication Xtvl=X*tvl,tvkX=tvk*X
n=Glob_n
do s=1,n
  temp2=ZERO
  do c=1,n
	temp2=temp2+X(c,s)*tvl(c)
  enddo
  Xtvl(s)=temp2
enddo
do s=1,n
     temp2=ZERO
       do c=1,n
	temp2=temp2+tvk(c)*X(c,s)
  enddo
  tvkX(s)=temp2
enddo
!Doing multiplication tAlX=tAl*X,tAkX=tAk*X
do s=1,n
  do c=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+tAl(c,k)*X(k,s)
      temp2=temp2+tAk(c,k)*X(k,s)
    enddo
    tAlX(c,s)=temp1
    tAkX(c,s)=temp2
  enddo
enddo

!Doing multiplication tAlXtAk=tAlX*tAk,tAkXtAl=tAkX*tAl
do s=1,n
  do c=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+tAlX(c,k)*tAk(k,s)
      temp2=temp2+tAkX(c,k)*tAl(k,s)
    enddo
    tAlXtAk(c,s)=temp1
    tAkXtAl(c,s)=temp2
  enddo
enddo


!Compute theta=tr[inv_tAkl*(tAl*X*tAk+tAk*X*tAl)]/2
theta=ZERO
theta1=ZERO
do s=1,n
    temp1=ZERO
    temp2=ZERO
  do k=1,n
    temp1=temp1+inv_tAkl(s,k)*tAlXtAk(k,s)
    temp2=temp2+inv_tAkl(s,k)*tAkXtAl(k,s)
  enddo
  theta=theta+temp1
  theta1=theta1+temp2
enddo
theta=(theta+theta1)/2

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'

if (i==j) then
 do s=1,n
   Aj(s)=inv_tAkl(s,i)
 enddo
else
 do s=1,n
   Aj(s)=inv_tAkl(s,i)-inv_tAkl(s,j)
 enddo
endif

!Compute kappa=tr[inv_tAkl*Jij*inv_tAkl*(tAl*X*tAk+tAk*X*tAl)]/2 = Aj' (tAlXtAk +tAlXtAk)Aj/2
kappa=ZERO
kappa1=ZERO
do s=1,n
  temp1=ZERO
  temp2=ZERO
  do c=1,n
    temp1=temp1+Aj(c)*tAlXtAk(c,s)*Aj(s)
    temp2=temp2+Aj(c)*tAkXtAl(c,s)*Aj(s)
  enddo
  kappa=kappa+temp1
  kappa1=kappa1+temp2
enddo
kappa=(kappa+kappa1)/2

!Compute temp1 = j^{ij}' inv_tAkl tvl = j^{ij} ' inv_tAkltvl
!        temp2 = tvk' inv_tAkl j^{ij} = tvkinv_tAkl' j^{ij}
if (i==j) then
  temp1=inv_tAkltvl(i)
  temp2=tvkinv_tAkl(i)
else
  temp1=inv_tAkltvl(i)-inv_tAkltvl(j)
  temp2=tvkinv_tAkl(i)-tvkinv_tAkl(j)
endif

!Compute the following quantities
!lambda = tr[inv_tAkl*tAl*X*tAk*inv_tAkl*V] = tvkinv_tAkl' tAlXtAk inv_tAkltvl
! omega = tr[inv_tAkl*(tAl*X*tAk+tAk*X*tAl)*inv_tAkl*Jij*inv_tAkl*V]/2 = (tvkinv_tAkl' (tAlXtAk+tAkXtAl) Aj) * temp1/2
!   chi = tr[inv_tAkl*Jij*inv_tAkl*(tAl*X*tAk+tAk*X*tAl)*inv_tAkl*V]/2 = temp2 * (Aj' (tAlXtAk+tAkXtAl) inv_tAkl)/2
!h =  tr[inv_tAkl*Jij*inv_tAkl*tAk*X*tvl*tvk]
!m=   tr[X*tAl*inv_tAkl*Jij*inv_tAkl*tvl*tvk]
lambda=ZERO
omega=ZERO
omega1=ZERO
chi=ZERO
chi1=ZERO
h=ZERO
m=ZERO
do s=1,n
  do c=1,n
    lambda=lambda+tvkinv_tAkl(c)*tAlXtAk(c,s)*inv_tAkltvl(s)
    omega=omega+tvkinv_tAkl(c)*tAlXtAk(c,s)*Aj(s)
    omega1=omega1+tvkinv_tAkl(c)*tAkXtAl(c,s)*Aj(s)
    m=m+tvkX(c)*tAl(c,s)*Aj(s)
    chi=chi+Aj(c)*tAlXtAk(c,s)*inv_tAkltvl(s)
    chi1=chi1+Aj(c)*tAkXtAl(c,s)*inv_tAkltvl(s)
    h=h+Aj(c)*tAk(c,s)*Xtvl(s)
  enddo
enddo
omega=omega*temp1
omega1=omega1*temp1
omega=(omega+omega1)/2
m=m*temp1
chi=temp2*chi
chi1=temp2*chi1
chi=(chi+chi1)/2
h=temp2*h

temp1=(4*Skl)/(15*SQRTPI*tr_AV*tr_AJ*tr_AJ*sqrt(tr_AJ))
temp2=15*kappa*tr_AJAV + &
  5*tr_AJ*(6*lambda*tr_AJ+9*tr_AV*theta*tr_AJ-3*kappa*tr_AV-2*omega-2*chi+h+m-3*theta*tr_AJAV)
ME_d_X_over_rij_d=temp1*temp2

end function ME_d_X_over_rij_d



function ftransaux(x)
!This function evaluates
!f(x) = ( |x| - sqrt(1-x^2)arccos(sqrt(1-x^2)) ) / (x |x|)
!given a real valued -1<x<1 argument. |x| stands for absolute value.
!A series representation is employed for |x|<xmin
!Depending on the kind parameter (dprec=4,8,10,16 - double, extended, or quadruple
!precision) for real numbers, a different xmin is used.
!In all cases the accuracy is close to machine precision corresponding to
!that kind parameter (1-2 last significant figures might be inaccurate in
!the worst case).
real(dprec) ftransaux,x
real(dprec),parameter :: xmin_4=0.30_dprec !for single precision
real(dprec),parameter :: xmin_8=0.27_dprec !for double precision
real(dprec),parameter :: xmin_10=0.2_dprec !for extended precision
real(dprec),parameter :: xmin_16=0.065_dprec !for quadruple precision
!Local variables
real(dprec) x2,ax,t,xmin

ax=abs(x)
selectcase (dprec)
  case(0:4)
    xmin=xmin_4
  case(5:8)
    xmin=xmin_8
  case(9:10)
    xmin=xmin_10
  case(11:16)
    xmin=xmin_16
endselect
if (ax<xmin) then
  x2=x*x
  t=(524288.0_dprec/50702925.0_dprec)
  t=(262144.0_dprec/22309287.0_dprec)+x2*t
  t=(65536.0_dprec/4849845.0_dprec)+x2*t
  t=(32768.0_dprec/2078505.0_dprec)+x2*t
  t=(2048.0_dprec/109395.0_dprec)+x2*t
  t=(1024.0_dprec/45045.0_dprec)+x2*t
  t=(256.0_dprec/9009.0_dprec)+x2*t
  t=(128.0_dprec/3465.0_dprec)+x2*t
  t=(16.0_dprec/315.0_dprec)+x2*t
  t=(8.0_dprec/105.0_dprec)+x2*t
  t=(2.0_dprec/15.0_dprec)+x2*t
  t=(1.0_dprec/3.0_dprec)+x2*t
  ftransaux=x*t
else
  t=sqrt(1.0_dprec-x*x)
  ftransaux=(ax-t*acos(t))/(ax*x)
endif

end function ftransaux


subroutine spinPreCalc(n, nFactorial, SziME, parityFactor, SSFmassChargeCoefficient, SSNCmassChargeCoefficient, &
  SOmassChargeCoefficient, AMMmassChargeCoefficient, AnihMassChargeCoefficient, ketMatrix, spatialYoung, &
  positronPosition, numberOfSpinFunctions, spinFreeME, SiSjME, SSNCspinME)
  use spinStuff
  implicit none

  character(len = maxLen), intent(in) :: spatialYoung
  integer, intent(in) :: n, nFactorial

  real(dprec), dimension(nFactorial), intent(out) :: parityFactor
  real(dprec), dimension(n, n, 4), intent(out) :: SOmassChargeCoefficient, AMMmassChargeCoefficient
  real(dprec), dimension(n, n, nFactorial), intent(out) :: ketMatrix
  real(dprec), dimension(n, n), intent(out) :: SSFmassChargeCoefficient, AnihMassChargeCoefficient, SSNCmassChargeCoefficient
  
  integer, intent(out) :: positronPosition, numberOfSpinFunctions

  real(dprec), dimension(nFactorial), intent(out) :: spinFreeME
  real(dprec), dimension(n, 2, nFactorial), intent(out) :: SziME
  real(kind = dprec), dimension(n, n, 2, nFactorial), intent(out) :: SiSjME, SSNCspinME
  
  ! local variables
  integer :: i, j, k, l, m
  character(len = maxLen) :: mySpatialYoung
  integer, dimension(nFactorial) :: parities
  integer, dimension(n, n, nFactorial) :: allPermutations

  SSFmassChargeCoefficient = ZERO
  SSNCmassChargeCoefficient = ZERO
  do i = 1, n
    do j = 1, n
      SSFmassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
           (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * EIGHT * PI / THREE
      SSNCmassChargeCoefficient(i, j) = Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
           (Glob_Mass(i + 1) * Glob_Mass(j + 1))
      
    enddo
  enddo

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
  allPermutations, parities, spinFreeME, SziME, SiSjMe, SSNCspinME)

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


subroutine spinDependentMatrixElements(m_k, m_l, vechLk, vechLl, Pket, &
     SziME, SSNCspinME, SSNCmassChargeCoefficient, SOmassChargeCoefficient, &
     AMMmassChargeCoefficient, SSNCkl, SO1kl, SO2kl, &
     AMM1kl, AMM2kl, numberOfSpinFunctions)
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
  integer,intent(in)       :: m_k, m_l, numberOfSpinFunctions
  real(dprec),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
  real(dprec),intent(in)   :: Pket(Glob_n,Glob_n)

  real(dprec), dimension(numberOfSpinFunctions), intent(out)  :: SO1kl, SO2kl, AMM1kl, AMM2kl, SSNCkl
  !Parameters (These are needed to declare static arrays. Using static
  !arrays makes the function call a little faster in comparison with
  !the case when arrays are dynamically allocated in stack)
  integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
  integer,parameter :: nnp=nn*(nn+1)/2
  real(dprec),intent(in)   :: SSNCspinME(Glob_n, Glob_n, numberOfSpinFunctions), &
                              SziME(Glob_n, numberOfSpinFunctions), &
                              SOmassChargeCoefficient(Glob_n, Glob_n, 4), &
                              AMMmassChargeCoefficient(Glob_n, Glob_n, 4), &
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
  real(dprec) :: commonFactor, gamma, gamma_diag, jiVl, jiAlAklinvVk, jiAlAklinvVl, jiAklinvVk, jiAklinvVl, &
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
  commonFactor = TWO * Glob_Piraised3n2 / (SQRTPI * det_tAkl * sqrt(det_tAkl))

  SO1kl = ZERO
  SO2kl = ZERO

  AMM1kl = ZERO
  AMM2kl = ZERO

  do indexI = 1, n

    i = 0
    do k = 1, numberOfSpinFunctions
       if (abs(SziME(indexI, k)) < localEps) i = i + 1
    enddo
    if (i == numberOfSpinFunctions) cycle ! corresponding <Szi> = 0

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

    do k = 1, numberOfSpinFunctions
      SO1kl(k) = SO1kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
      AMM1kl(k) = AMM1kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
    enddo


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
      
      do k = 1, numberOfSpinFunctions
        SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
      enddo
      
      temp1 = &
      gamma**3 / THREE * (jjVl * (jjAklinvVk - jiAklinvVk) + &
      jjAlAklinvVk * (jjAklinvVl - jiAklinvVl) + &
      jjAlAklinvVl * (jiAklinvVk - jjAklinvVk)) 

      do k = 1, numberOfSpinFunctions
        SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
        AMM2kl(k) = AMM2kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
      enddo

      
      temp1 = &
      gamma**3 / THREE * (jiVl * (jiAklinvVk - jjAklinvVk) + &
      jiAlAklinvVk * (jiAklinvVl - jjAklinvVl) + &
      jiAlAklinvVl * (jjAklinvVk - jiAklinvVk)) 

      do k = 1, numberOfSpinFunctions
        SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
        AMM2kl(k) = AMM2kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
     enddo
     

    enddo ! indexJ cycle

  enddo ! indexI cycle


  ! SS term (separate loop for not to interfere with the case Siz == 0 at previous loop)
  SSNCkl = ZERO
  do indexI = 1,n
     do indexJ = indexI + 1,n !(we need only indexJ > indexI)
        
        gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
        inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))
  
        jjAklinvVk = inv_tAkl(indexJ, pm_k)
        jjAklinvVl = inv_tAkl(indexJ, pm_l)
        jiAklinvVl = inv_tAkl(indexI, pm_l)
        jiAklinvVk = inv_tAkl(indexI, pm_k)

        temp1 = &
           (gamma**5 / 15._dprec) * ( jiAklinvVk * (jiAklinvVl  - jjAklinvVl) + &
           jjAklinvVk * (jjAklinvVl - jiAklinvVl) ) !additional factor of 1/sqrt(6) is taken from spin part 

        do k = 1, NumberOfSpinFunctions  
            if (abs(SSNCspinME(indexI, indexJ, k)) < localEps) cycle
            SSNCkl(k) = SSNCkl(k) + SSNCspinME(indexI, indexJ, k) * SSNCmassChargeCoefficient(indexI, indexJ) * temp1
        enddo
        
     enddo !indexJ loop
  enddo !indexI loop
    
  SSNCkl = SSNCkl * commonFactor
  SO1kl = SO1kl * commonFactor
  SO2kl = SO2kl * commonFactor
  AMM1kl = AMM1kl * commonFactor
  AMM2kl = AMM2kl * commonFactor


end subroutine spinDependentMatrixElements


subroutine overlapMatrixElementsL1(m_k, vechLk, P, Skk)
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


end subroutine overlapMatrixElementsL1





end module matelem
