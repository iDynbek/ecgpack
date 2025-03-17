module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=1 Gaussians.
use globvars
implicit none

contains

subroutine MatrixElementsLD(m_k, mm_k, m_l, mm_l, vechLk, vechLl, P, &
               Hkl, Skl, Tkl,Vkl, Dk, Dl, grad_k, grad_l)
!This subroutine computes symmetry adapted matrix element with
!two real L=1 correlated Gaussians:
! 
!fk = [x_{m_k}x_{mm_k}+y_{mm_k}y_{m_k}-2*z_{mm_k}z_{m_k}] exp[-r'(Lk*Lk')r] 
!
!where m_k,mm_k is some integers between 1 and n (n is the number of 
!pseudoparticles). Symmetry adaption is applied to the ket using 
!permutation matrix P
!
!Input:     
!   m_k,m_l,mm_k, mm_l :: integers that determine which x or y-components is in the
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
integer,intent(in)          :: m_k,m_l,mm_k,mm_l
real(dprec),intent(in)      :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec),intent(in)      :: P(Glob_n,Glob_n)
real(dprec),intent(out)     :: Skl,Hkl,Tkl, Vkl
real(dprec),intent(out)     :: Dk(2*Glob_np),Dl(2*Glob_np)
logical,intent(in)          :: grad_k, grad_l

!Parameters (These are needed to declare static arrays. Using static 
!arrays makes the function call a little faster in comparison with 
!the case when arrays are dynamically allocated in stack)
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
integer           n, np
integer           vl(nn),bl(nn),bk(nn)
real(dprec)       dHkldvechLk(nnp), dHkldvechLl(nnp)
real(dprec)       dSkldvechLk(nnp), dSkldvechLl(nnp)
real(dprec)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
real(dprec)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
real(dprec)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn)
real(dprec)       inv_tAkltAl(nn,nn),inv_tAkltAlM(nn,nn)
real(dprec)       inv_tAklAk(nn,nn),inv_tAklAkM(nn,nn),MAkinv_tAkl(nn,nn),MAk(nn,nn)
real(dprec)       eta1(nn,nn),sqrt_eta1(nn,nn),eta2(nn,nn),Rkl(nn,nn)
real(dprec)       eta11(nn,nn),eta22(nn,nn),eta(nn,nn)
real(dprec)       eta221(nn,nn),eta222(nn,nn),eta223(nn,nn),eta224(nn,nn),etan(nn,nn)
real(dprec)       W1(nn,nn),W2(nn,nn),W3(nn,nn),W4(nn,nn)
real(dprec)       twosym_tFkl(nn,nn),two_Fkk(nn,nn),two_Fll(nn,nn),twosym_tGkl(nn,nn)
real(dprec)       tKkl1(nn,nn),tUkl(nn,nn),tWkl(nn,nn), Pkl(nn,nn), Mkl(nn,nn),tKkll(nn,nn)
real(dprec)       tKkl2(nn,nn),tKkl3(nn,nn),tKkl4(nn,nn),tKkl5(nn,nn),tKkl6(nn,nn)
real(dprec)       gkl(nn,nn)
real(dprec)       twosym_tQkl(nn,nn),twosym_tDkl(nn,nn)
real(dprec)       inv_tAklvl(nn),vkinv_tAkl(nn),vkinv_tAkltAlM(nn),inv_tAklbk(nn),vlinv_tAkl(nn)
real(dprec)       inv_tAklbl(nn),bkinv_tAkl(nn),bkinv_tAkltAlM(nn),vlinv_tAklAkM(nn),vlMAkinv_tAkl(nn)
real(dprec)       u1(nn),u2(nn),u3(nn),u331(nn),u332(nn),u333(nn),u334(nn),vlinv_tAkltAlM(nn)
real(dprec)       u11(nn),u22(nn),u33(nn),u111(nn),u112(nn),u113(nn),u114(nn)
real(dprec)       u221(nn),u222(nn),u223(nn),u224(nn)
real(dprec)       temp1, temp2, temp3, temp4, temp5, temp6,temp44,temp444,temp11,temp22,temp66,temp666
real(dprec)       temp661,temp662,temp663,temp664,h1
real(dprec)       det_Lk, det_Ll, det_tAkl
real(dprec)       tau1,tau2,tau3,tau11,tau22,tau33,m,kkl
!real(dprec)       Tkl, Vkl
real(dprec)       m1,m2,m3,tau331,tau332,tau333,tau334,tau221,tau222,tau223,tau224,tau221_new,tau222_new,tau221_ad,tau222_ad
real(dprec)       temp441,temp442,temp443,temp4440,temp4441,temp4442,term1,term2,h,temp_n,temp4_new
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
  vl(i)=P(m_l,i)
  bl(i)=P(mm_l,i)
  bk(i)=P(mm_k,i)
enddo

!Compute inv_tAklvl = inv_tAkl * vl, inv_tAklbl = inv_tAkl * bl
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do j=1,n
    temp1=temp1+inv_tAkl(j,i)*vl(j)
    temp2=temp2+inv_tAkl(j,i)*bl(j)
    temp3=temp3+inv_tAkl(j,i)*bk(j)
  enddo
  inv_tAklvl(i)=temp1
  inv_tAklbl(i)=temp2
  inv_tAklbk(i)=temp3
enddo

!Compute vkinv_tAkl=vk'*inv_tAkl, bkinv_tAkl=bk'*inv_tAkl
do i=1,n
  vkinv_tAkl(i)=inv_tAkl(m_k,i)
  vlinv_tAkl(i)=inv_tAkl(m_l,i)
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
m=m1+m3  !look formula 40 in document

!Evaluating overlap
!temp1=ZERO
temp1=FOUR*det_tAkl*sqrt(det_tAkl)
Skl=Glob_Piraised3n2*m/temp1

!Doing multiplication inv_tAkltAl=inv_tAkl*tAl
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
      temp2=temp2+inv_tAkl(j,k)*Ak(k,i) 
    enddo
    inv_tAkltAl(j,i)=temp1
    inv_tAklAk(j,i)=temp2
  enddo
enddo

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkltAl(j,k)*Glob_MassMatrix(k,i)
      temp2=temp2+inv_tAklAk(j,k)*Glob_MassMatrix(k,i)  
    enddo
    inv_tAkltAlM(j,i)=temp1
    inv_tAklAkM(j,i)=temp2
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
!Computing tau2 = vk'*inv_tAkltAlM*Ak*inv_tAklvl
!We do it by multiplying twice the row-vector on the left
!by a matrix on the right and computing a dot product in the end.

!vkinv_tAkltAlM'=vk'*inv_tAkltAlM ,  bkinv_tAkltAlM'=bk'*inv_tAkltAlM

do i=1,n
  vkinv_tAkltAlM(i)=inv_tAkltAlM(m_k,i)
  bkinv_tAkltAlM(i)=inv_tAkltAlM(mm_k,i)
  vlinv_tAkltAlM(i)=inv_tAkltAlM(m_l,i)
  vlinv_tAklAkM(i)=inv_tAklAkM(m_l,i)
enddo




!u1=vkinv_tAkltAlM'*Ak, u11=bkinv_tAkltAlM'*Ak
!tau2=u1'*inv_tAklvl (storage for u1 as such is not needed, we use temp1=u1(i))
!tau22=u11'*inv_tAklbl
tau2=ZERO
tau22=ZERO
tau223=ZERO
tau224=ZERO
temp4_new=ZERO
tau221_new=ZERO
tau222_new=ZERO
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+vkinv_tAkltAlM(j)*Ak(j,i)
    temp2=temp2+bkinv_tAkltAlM(j)*Ak(j,i)     
  enddo
  tau2=tau2+temp1*inv_tAklvl(i)
  tau22=tau22+temp2*inv_tAklbl(i)
  tau223=tau223+temp1*inv_tAklbl(i)
  tau224=tau224+temp2*inv_tAklvl(i)
enddo

h=tau3*tau22+tau33*tau2+tau333*tau224+tau334*tau223  

Tkl=Skl*(SIX*tau1+FOUR*h/m)
!Evaluating eta1(i,j), sqrt_eta1(i,j), eta2(i,j), Rkl(i,j),
!and the potential energy. Notice that only the lower triangles
!of eta1, sqrt_eta1, eta2, and Rkl are filled.
!temp1=ZERO
Vkl=ZERO
temp1=Skl*(TWO/SQRTPI)
!temp1=Glob_Piraised3n2/(TWO*SQRTPI*det_tAkl*sqrt(det_tAkl))
do i=1,n
  temp2=inv_tAkl(i,i)
  temp3=sqrt(temp2)
  eta1(i,i)=temp2
  sqrt_eta1(i,i)=temp3
  !Getting row m_k of matrix inv_tAkl*Jii*inv_tAkl
  !as only this row is needed to compute eta2(i,i)
  do k=1,n
    u1(k)=inv_tAkl(i,m_k)*inv_tAkl(k,i) 
    u11(k)=inv_tAkl(i,mm_k)*inv_tAkl(k,i) 
  enddo
  temp4=ZERO
  temp44=ZERO
  temp443=ZERO
  temp444=ZERO
  temp4440=ZERO
  temp4442=ZERO
  do k=1,n
    temp4=temp4+u1(k)*vl(k)
    temp44=temp44+u11(k)*bl(k)
    temp443=temp443+u1(k)*bl(k)
    temp444=temp444+u11(k)*vl(k)
  enddo
  eta2(i,i)=temp44 ! to make consistent notation in document
  eta22(i,i)=temp4
  eta223(i,i)= temp444
  eta224(i,i)= temp443  
  eta(i,i)=temp4*temp44+temp443*temp444
  Rkl(i,i)=temp1*(ONE-ONETHIRD*(tau3*temp44+tau33*temp4+tau333*temp444+tau334*temp443)/(m*temp2) &
  + ONEFIFTH*(temp4*temp44+temp443*temp444)/(m*temp2*temp2))/temp3
  Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*Rkl(i,i)
enddo
do i=1,n
  do j=i+1,n
    temp2=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    temp3=sqrt(temp2)
    eta1(j,i)=temp2
    sqrt_eta1(j,i)=temp3
    do k=1,n
      u1(k)=(inv_tAkl(i,m_k)-inv_tAkl(j,m_k))*(inv_tAkl(k,i)-inv_tAkl(k,j))
      u11(k)=(inv_tAkl(i,mm_k)-inv_tAkl(j,mm_k))*(inv_tAkl(k,i)-inv_tAkl(k,j))
    enddo
    temp4=ZERO
    temp44=ZERO
    temp443=ZERO
    temp444=ZERO
    temp4440=ZERO
    temp4442=ZERO
    do k=1,n
      temp4=temp4+u1(k)*vl(k)
      temp44=temp44+u11(k)*bl(k)   
      temp443=temp443+u1(k)*bl(k)
      temp444=temp444+u11(k)*vl(k)
    enddo
    eta2(j,i)=temp44
    eta22(j,i)=temp4
    eta223(j,i)= temp444
    eta224(j,i)= temp443  
    eta(j,i)=temp4*temp44+temp443*temp444
    Rkl(j,i)=temp1*(ONE-ONETHIRD*(tau3*temp44+tau33*temp4+tau333*temp444+tau334*temp443)/(m*temp2)+ &
                 ONEFIFTH*(temp4*temp44+temp443*temp444)/(m*temp2*temp2))/temp3
    Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*Rkl(j,i)
  enddo  
enddo
!Hkl=ZERO
Hkl=Tkl+Vkl

!Now we start computing the gradient of Skl  (in gradient I follow the notation in document))

if (grad_k.or.grad_l) then
  !Evaluating matrix tKkl = inv_tAklvl * vkinv_tAkl'
  !which will be used a lot below. 
  do i=1,n
    do j=1,n
      tKkl1(i,j)=inv_tAklbl(i)*bkinv_tAkl(j)
      tKkl2(i,j)=inv_tAklvl(i)*vkinv_tAkl(j)
      tKkl5(i,j)=inv_tAklvl(i)*bkinv_tAkl(j)
      tKkl6(i,j)=inv_tAklbl(i)*vkinv_tAkl(j)   
      tKkll(i,j)= tKkl1(i,j)*tau3+tKkl2(i,j)*tau33+tKkl5(i,j)*tau333+tKkl6(i,j)*tau334                                          
      gkl(i,j)= -tKkl1(i,j)*tau2-tKkl2(i,j)*tau22-tKkl5(i,j)*tau223-tKkl6(i,j)*tau224      
    enddo
  enddo
  
  
do i=1,n  
        do j=1,i     
      twosym_tFkl(i,j)=THREE*inv_tAkl(j,i)+(tKkll(i,j)+tKkll(j,i))/m
      twosym_tFkl(j,i)=twosym_tFkl(i,j)
    enddo
  enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (grad_k) then  
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
    temp2=ZERO
    temp3=ZERO
    do j=1,n
      temp1=temp1+vkinv_tAkltAlM(j)*inv_tAkltAl(i,j)
      temp2=temp2+bkinv_tAkltAlM(j)*inv_tAkltAl(i,j)
      temp3=temp3+vlinv_tAklAkM(j)*inv_tAkltAl(i,j)
    enddo
    u1(i)=temp1
    u11(i)=temp2
    u111(i)=temp3
  enddo
  !Computing u2=inv_tAkltAlM*Ak*inv_tAklvl
  do i=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    do j=1,n
      temp1=temp1+Ak(i,j)*inv_tAklvl(j)
      temp2=temp2+Ak(i,j)*inv_tAklbl(j)
      temp3=temp3+tAl(i,j)*inv_tAklbk(j)
    enddo
    u3(i)=temp1
    u33(i)=temp2
    u333(i)=temp3
  enddo
  
  do i=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    do j=1,n
      temp1=temp1+inv_tAkltAlM(i,j)*u3(j)
      temp2=temp2+inv_tAkltAlM(i,j)*u33(j)
      temp3=temp3+inv_tAkltAlM(i,j)*u333(j)
    enddo
    u2(i)=temp1
    u22(i)=temp2
    u222(i)=temp3
  enddo   
  temp1=FOUR/m
  temp2=temp1*h/m
  do i=1,n
    do j=1,n                        
       tUkl(j,i)= W1(j,i)+temp2*tKkll(j,i)+temp1*(gkl(j,i)&     
       +tau33*(inv_tAklvl(j)*u1(i)-u2(j)*vkinv_tAkl(i))+tau3*(inv_tAklbl(j)*u11(i)-u22(j)*bkinv_tAkl(i))&                   
       +tau333*(inv_tAklvl(j)*u11(i)-u2(j)*bkinv_tAkl(i))+tau334*(inv_tAklbl(j)*u1(i)-u22(j)*vkinv_tAkl(i)))
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
  !Computing u1=inv_tAklAkM*inv_tAklAk'*vl
  do i=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    do j=1,n
      temp1=temp1+inv_tAklAk(j,i)*vl(j)
      temp2=temp2+inv_tAklAk(j,i)*bl(j)
      temp3=temp3+inv_tAkltAl(j,i)*bk(j)
    enddo
    u3(i)=temp1
    u33(i)=temp2
    u333(i)=temp3
  enddo
  
  do i=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    do j=1,n
      temp1=temp1+inv_tAklAkM(i,j)*u3(j)
      temp2=temp2+inv_tAklAkM(i,j)*u33(j)
      temp3=temp3+inv_tAklAkM(i,j)*u333(j)
    enddo
    u1(i)=temp1
    u11(i)=temp2
    u111(i)=temp3
  enddo 
  
  !Computing u2'=vkinv_tAkltAlM'*inv_tAklAk'
  do i=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    do j=1,n
      temp1=temp1+vkinv_tAkltAlM(j)*inv_tAklAk(i,j)
      temp2=temp2+bkinv_tAkltAlM(j)*inv_tAklAk(i,j)
      temp3=temp3+vlinv_tAklAkM(j)*inv_tAklAk(i,j)
    enddo
    u2(i)=temp1
    u22(i)=temp2
    u222(i)=temp3
  enddo     

  temp1=FOUR/m
  temp2=temp1*h/m
  do i=1,n
    do j=1,n
      W3(j,i)= W1(j,i)+temp2*tKkll(j,i)+temp1*( gkl(j,i) &  
      +tau33*(u1(j)*vkinv_tAkl(i)-inv_tAklvl(j)*u2(i))+tau3*(u11(j)*bkinv_tAkl(i)-inv_tAklbl(j)*u22(i))&
      + tau333*(u1(j)*bkinv_tAkl(i)-inv_tAklvl(j)*u22(i))+tau334*(u11(j)*vkinv_tAkl(i)-inv_tAklbl(j)*u2(i)))
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
    temp_n= tau3*eta2(i,i)+tau33*eta22(i,i)+tau333*eta223(i,i)+tau334*eta224(i,i)
    temp1=(TWO/SQRTPI)/(eta1(i,i)*sqrt_eta1(i,i))
    temp2=ONE+(eta(i,i)/eta1(i,i)-temp_n)/(eta1(i,i)*m) !first term
    temp3=(ONEFIFTH*eta(i,i)/eta1(i,i)-ONETHIRD*temp_n)/(m*m) !second term
    temp4=ONETHIRD/m    !3rd term
    temp11=ONEFIFTH/(eta1(i,i)*m)  !4th term                                            
    do t=1,n
      do q=t,n
        temp5=inv_tAkl(q,i)*inv_tAkl(i,t)              
        temp6=inv_tAkl(q,i)*(tKkl1(i,t)+tKkl1(t,i))+inv_tAkl(t,i)*(tKkl1(i,q)+tKkl1(q,i)) !old
        temp66=inv_tAkl(q,i)*(tKkl2(i,t)+tKkl2(t,i))+inv_tAkl(t,i)*(tKkl2(i,q)+tKkl2(q,i))                 
        temp663=inv_tAkl(q,i)*(tKkl5(i,t)+tKkl5(t,i))+inv_tAkl(t,i)*(tKkl5(i,q)+tKkl5(q,i)) 
        temp664=inv_tAkl(q,i)*(tKkl6(i,t)+tKkl6(t,i))+inv_tAkl(t,i)*(tKkl6(i,q)+tKkl6(q,i))    
        
        temp44=eta2(i,i)*temp66+eta22(i,i)*temp6+eta223(i,i)*temp664+eta224(i,i)*temp663       
               
        temp666= eta2(i,i)*(tKkl2(q,t)+tKkl2(t,q))+eta22(i,i)*(tKkl1(q,t)+tKkl1(t,q))+&                  
               eta223(i,i)*(tKkl6(q,t)+tKkl6(t,q))+eta224(i,i)*(tKkl5(q,t)+tKkl5(t,q))+ &
               tau3*temp6+tau33*temp66+tau333*temp663+tau334*temp664  !3rd term              
        twosym_tQkl(q,t)=temp1*(temp2*temp5+temp3*(tKkll(q,t)+tKkll(t,q))+temp4*temp666-temp11*temp44)  
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
       temp_n= tau3*eta2(j,i)+tau33*eta22(j,i)+tau333*eta223(j,i)+tau334*eta224(j,i)
       temp1=(TWO/SQRTPI)/(eta1(j,i)*sqrt_eta1(j,i))
       temp2=ONE+(eta(j,i)/eta1(j,i)-temp_n)/(eta1(j,i)*m) !first term
       temp3=(ONEFIFTH*eta(j,i)/eta1(j,i)-ONETHIRD*temp_n)/(m*m) !second term
       temp4=ONETHIRD/m    !3rd term
       temp11=ONEFIFTH/(eta1(j,i)*m)  !4th term    
      do t=1,n
        do q=t,n
          temp5=(inv_tAkl(q,i)-inv_tAkl(q,j))*(inv_tAkl(i,t)-inv_tAkl(j,t))                
          temp6=(inv_tAkl(q,i)-inv_tAkl(q,j))*(tKkl1(i,t)-tKkl1(j,t)+tKkl1(t,i)-tKkl1(t,j))+ &
                (tKkl1(q,i)-tKkl1(q,j)+tKkl1(i,q)-tKkl1(j,q))*(inv_tAkl(t,i)-inv_tAkl(t,j))                   
          temp66=(inv_tAkl(q,i)-inv_tAkl(q,j))*(tKkl2(i,t)-tKkl2(j,t)+tKkl2(t,i)-tKkl2(t,j))+ &
                (tKkl2(q,i)-tKkl2(q,j)+tKkl2(i,q)-tKkl2(j,q))*(inv_tAkl(t,i)-inv_tAkl(t,j))                                 
          temp663=(inv_tAkl(q,i)-inv_tAkl(q,j))*(tKkl5(i,t)-tKkl5(j,t)+tKkl5(t,i)-tKkl5(t,j))+ &
                (tKkl5(q,i)-tKkl5(q,j)+tKkl5(i,q)-tKkl5(j,q))*(inv_tAkl(t,i)-inv_tAkl(t,j))                
          temp664=(inv_tAkl(q,i)-inv_tAkl(q,j))*(tKkl6(i,t)-tKkl6(j,t)+tKkl6(t,i)-tKkl6(t,j))+ &
                (tKkl6(q,i)-tKkl6(q,j)+tKkl6(i,q)-tKkl6(j,q))*(inv_tAkl(t,i)-inv_tAkl(t,j))      
          temp44=eta2(j,i)*temp66+eta22(j,i)*temp6+eta223(j,i)*temp664+eta224(j,i)*temp663  
          temp666= eta2(j,i)*(tKkl2(q,t)+tKkl2(t,q))+eta22(j,i)*(tKkl1(q,t)+tKkl1(t,q))+&                  
               eta223(j,i)*(tKkl6(q,t)+tKkl6(t,q))+eta224(j,i)*(tKkl5(q,t)+tKkl5(t,q))+ &
               tau3*temp6+tau33*temp66+tau333*temp663+tau334*temp664  !3rd term            
          twosym_tQkl(q,t)=temp1*(temp2*temp5+temp3*(tKkll(q,t)+tKkll(t,q))+temp4*temp666-temp11*temp44) 
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

end subroutine MatrixElementsLD

subroutine MatrixElementsL1ForExpcValsD(m_k, mm_k, m_l, mm_l, vechLk, vechLl, Pbra, Pket, &
           Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
           MVkl, drach_MVkl, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, prvalkl, &
           NumCFGridPoints, CFGrid, &
           CFkl, NumDensGridPoints, DensGrid, Denskl, AreCorrFuncNeeded, ArePartDensNeeded)

           
!Arguments
integer,intent(in)       :: m_k,m_l,mm_k,mm_l
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
logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded


integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
integer           n,np
integer           tvk(nn),tvl(nn),tbk(nn),tbl(nn)
real(dprec)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn),tvl8(nn)
real(dprec)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn),rmkl1(nn,nn),Vkl1,Hkl1
real(dprec)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn), inv_tAkltAl(nn,nn)
real(dprec)       eta2(nn,nn),inv_tAkltAlM(nn,nn),eta22(nn,nn)
real(dprec)       eta(nn,nn)
real(dprec)       W1(nn,nn),W2(nn,nn),W3(nn,nn),W4(nn,nn),W5(nn,nn),W6(nn,nn),W7(nn,nn)
real(dprec)       W44(nn,nn),W55(nn,nn),W77(nn,nn)
real(dprec)       W4b(nn,nn),W5b(nn,nn),W7b(nn,nn)
real(dprec)       W44b(nn,nn),W55b(nn,nn),W77b(nn,nn),temp444,temp4444,temp444b,temp4444b
real(dprec)       inv_tAkltvl(nn),tvkinv_tAkl(nn),tvkinv_tAkltAlM(nn),u1(nn) 
real(dprec)       inv_tAkltbl(nn),tbkinv_tAkl(nn),tbkinv_tAkltAlM(nn),u11(nn) 
real(dprec)       temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp44
real(dprec)       temp55,temp66,temp55b,temp66b,temp44b,temp4b,temp5b,temp6b,temp1b,temp11b
real(dprec)       temp10,temp11,temp12,temp13,temp14,threshold,tr1, tr2, tr3, tr4,tr4vkvl,tr4bkbl,tr4bkvl,tr4vkbl
real(dprec)       det_Lk, det_Ll, det_tAkl,tau1,tau2,tau3,inv_tau3 ,V2kl,tau22,tau33,m
integer           i,j,k,t,indx,p,q
real(dprec)       TrAJ(nn,nn),sqrtTrAJ(nn,nn),TrAJAJ(nn,nn,nn,nn)
real(dprec)       jAj(nn,nn,nn,nn),jAtvl(nn,nn),tvkAj(nn,nn),Mass_For_Darwin(0:nn)

real(dprec)   m1,m2,m3,tau331,tau332,tau333,tau334,tau221,tau222,tau223,tau224,temp38,temp338,temp36,temp336,temp37,temp337
real(dprec)   temp31,temp331,temp32,temp332,temp33,temp333,temp34,temp334,temp35,temp335,templast
real(dprec)   eta221(nn,nn),eta222(nn,nn),eta223(nn,nn),eta224(nn,nn),u111(nn) 
real(dprec)   temp441,temp442,temp443,temp4440,temp4441,temp4442,h,term1,term2
real(dprec)   inv_tAkltAk(nn,nn),inv_tAkltAkM(nn,nn)
real(dprec)   inv_tAkltbk(nn),tvlinv_tAkl(nn),tvlinv_tAkltAkM(nn),tvlinv_tAkltAlM(nn)
real(dprec)   tbltbk(nn,nn),tvltvk(nn,nn),tbltvk(nn,nn),tvltbk(nn,nn),tvktbk(nn,nn),tvltbl(nn,nn)


!Vars for calculating <1/rij 1/pq>
real(dprec) :: a, b, d, fij, fpq, tfij, tfpq, uij, upq, tuij, tupq, phi,  phi_sq, phi_cube, dsqab, &
acosphi, tau, ttau, myeta, myteta, &
commonFactor, arccosCommon, a1, a2, a3, a4, aone, atwo, arccosAns
real(dprec) :: R11, R12, R1, R21, R22, R23, R2, R31, R32, R33, R3, R4, R51, R52, R53, R5, R6, &
ROne, RTwo, radicalCommon, radicalAns, totalAns, commonArccosRadical, xx, &
RDZeroOne, RDZeroTwo, RDOneOne, RDOneTwo, RDTwoOne, RDTwoTwo, RDThreeOne, RDThreeTwo, &
RDFourOne, RDFourTwo
real(dprec)   local_eps_for_xx
!Vars to calculate delta-fucntions directly
real(dprec) :: myalpha, jijAVk, jijAVl, jijAWk, jijAWl
!Var to set if orbit-orbit correction is needed
logical :: isOOklNeeded

isOOklNeeded = .true.
local_eps_for_xx = 1.d-6
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

!Doing multiplication inv_tAkltAl=inv_tAkl*tAl, inv_tAkltAk=inv_tAkl*tAk
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
      temp2=temp2+inv_tAkl(j,k)*tAk(k,i)
    enddo
    inv_tAkltAl(j,i)=temp1
    inv_tAkltAk(j,i)=temp2
  enddo
enddo

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkltAl(j,k)*Glob_MassMatrix(k,i)
      temp2=temp2+inv_tAkltAk(j,k)*Glob_MassMatrix(k,i)
    enddo
    inv_tAkltAlM(j,i)=temp1
    inv_tAkltAkM(j,i)=temp2
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


!Computing tvk=Pbra'*vk and tvl=Pket'*vl
do i=1,n
  tvk(i)=Pbra(m_k,i)
  tvl(i)=Pket(m_l,i)
  tbk(i)=Pbra(mm_k,i)
  tbl(i)=Pket(mm_l,i)
enddo



!Compute inv_tAkltvl = inv_tAkl * tvl


do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do j=1,n
    temp1=temp1+inv_tAkl(j,i)*tvl(j)
    temp2=temp2+inv_tAkl(j,i)*tbl(j)
    temp3=temp3+inv_tAkl(j,i)*tbk(j)
  enddo
  inv_tAkltvl(i)=temp1
  inv_tAkltbl(i)=temp2
  inv_tAkltbk(i)=temp3
enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAkl(j,i)
    temp2=temp2+tbk(j)*inv_tAkl(j,i)
    temp3=temp3+tvl(j)*inv_tAkl(j,i)
  enddo
  tvkinv_tAkl(i)=temp1
  tbkinv_tAkl(i)=temp2
  tvlinv_tAkl(i)=temp3
enddo


!Compute tau3=tvkinv_tAkl*tvl
tau3=ZERO
tau33=ZERO
tau333=ZERO
tau334=ZERO
do i=1,n
  tau3=tau3+tvkinv_tAkl(i)*tvl(i)
  tau33=tau33+tbkinv_tAkl(i)*tbl(i)
  tau333=tau333+tvkinv_tAkl(i)*tbl(i)
  tau334=tau334+tbkinv_tAkl(i)*tvl(i)
enddo
m1=tau3*tau33
m3=tau333*tau334
m=m1+m3  !look formula 40 in document
!Evaluating overlap
!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
temp1=FOUR*det_tAkl*sqrt(det_tAkl)
Skl=Glob_Piraised3n2*m/temp1



do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAkltAlM(j,i)
    temp2=temp2+tbk(j)*inv_tAkltAlM(j,i)
    temp3=temp3+tvl(j)*inv_tAkltAlM(j,i)
    temp4=temp4+tvl(j)*inv_tAkltAkM(j,i)
  enddo
  tvkinv_tAkltAlM(i)=temp1
  tbkinv_tAkltAlM(i)=temp2
  tvlinv_tAkltAlM(i)=temp3
  tvlinv_tAkltAkM(i)=temp4
enddo

!u1=tvkinv_tAkltAlM'*tAk
!tau2=u1'*inv_tAkltvl (storage for u1 as such is not needed, we use temp1=u1(i))
tau2=ZERO
tau22=ZERO
tau223=ZERO
tau224=ZERO
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkltAlM(j)*tAk(j,i)
    temp2=temp2+tbkinv_tAkltAlM(j)*tAk(j,i)
  enddo
  tau2=tau2+temp1*inv_tAkltvl(i)
  tau22=tau22+temp2*inv_tAkltbl(i)
  tau223=tau223+temp1*inv_tAkltbl(i)
  tau224=tau224+temp2*inv_tAkltvl(i)
enddo
h=tau3*tau22+tau33*tau2+tau333*tau224+tau334*tau223 
!Evaluating the kinetic energy
Tkl=Skl*(SIX*tau1+FOUR*h/m)
Vkl=ZERO
Vkl1=ZERO
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
    temp44=ZERO
    temp444=ZERO
    do k=1,n
      temp4=temp4+tvk(k)*inv_tAkl(k,i)*inv_tAkl(q,i)
      temp44=temp44+tbk(k)*inv_tAkl(k,i)*inv_tAkl(q,i)
      temp444=temp444+tvl(k)*inv_tAkl(k,i)*inv_tAkl(q,i)
    enddo
    u1(q)=temp4
    u11(q)=temp44
    u111(q)=temp444
  enddo
  !eta2=u1'*tvl
  temp4=ZERO
  temp44=ZERO
  temp443=ZERO
  temp444=ZERO
  temp4440=ZERO
  temp4441=ZERO
  temp4442=ZERO
  do k=1,n
    temp4=temp4+u1(k)*tvl(k)
    temp44=temp44+u11(k)*tbl(k) 
    temp443=temp443+u1(k)*tbl(k)
    temp444=temp444+u11(k)*tvl(k)
  enddo
  eta2(i,i)=temp44
  eta22(i,i)=temp4
  eta223(i,i)= temp444
  eta224(i,i)= temp443  
  !temp444=temp4*temp44
  eta(i,i)=temp4*temp44+temp443*temp444
  term1=tau3*temp44+tau33*temp4+tau333*temp444+tau334*temp443
  term2=temp4*temp44+temp443*temp444
  rm2kl(i,i)=temp5*(ONE-TWO*ONETHIRD*term1/(m*temp2) + EIGHT*ONEFIFTH*term2/(THREE*m*temp2*temp2))/temp2
  rmkl(i,i)=temp1*(ONE-ONETHIRD*term1/(m*temp2) + ONEFIFTH*term2/(m*temp2*temp2))/temp3
  !rmkl(i,i)=ME_over_rij(i,i,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
  Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*rmkl(i,i)
  !Vkl1=Vkl1+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*rmkl1(i,i)
  rkl(i,i)= temp1*temp3*(ONE+ONETHIRD*term1/(m*temp2) - ONEFIFTH*term2/(THREE*m*temp2*temp2))
  r2kl(i,i)=Skl*THREEHALF*temp2*(ONE+TWO*ONETHIRD*term1/(m*temp2))
  temp10=temp8/(temp2*temp3)
  !deltarkl(i,i)=temp10*(ONE-term1/(m*temp2)+term2/(THREE*m*temp2*temp2))
  prvalkl(i,i)=PI*temp10*( TWO*(Glob_EulerConst+log(temp2))*(ONE-term1/(m*temp2)+term2/(m*temp2*temp2)) &
  + FOUR*(term1-TWO*term2/temp2)/(THREE*m*temp2)+EIGHT*term2/(15*m*temp2*temp2) )
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
      temp44=ZERO
      temp444=ZERO
      temp10=ZERO
      do k=1,n
        temp4=temp4+tvk(k)*(inv_tAkl(k,i)-inv_tAkl(k,j))*(inv_tAkl(q,i)-inv_tAkl(q,j))
        temp44=temp44+tbk(k)*(inv_tAkl(k,i)-inv_tAkl(k,j))*(inv_tAkl(q,i)-inv_tAkl(q,j))
        temp444=temp444+tvl(k)*(inv_tAkl(k,i)-inv_tAkl(k,j))*(inv_tAkl(q,i)-inv_tAkl(q,j))
      enddo
      u1(q)=temp4
      u11(q)=temp44
      u111(q)=temp444
    enddo
    temp4=ZERO
    temp44=ZERO
    temp443=ZERO
    temp444=ZERO
    temp4440=ZERO
    temp4441=ZERO
    temp4442=ZERO
    do k=1,n
      temp4=temp4+u1(k)*tvl(k)
      temp44=temp44+u11(k)*tbl(k) 
      temp443=temp443+u1(k)*tbl(k)
      temp444=temp444+u11(k)*tvl(k)
    enddo
    !temp444=temp4*temp44
    eta2(j,i)=temp44
    eta2(i,j)=temp44   
    eta22(j,i)=temp4
    eta22(i,j)=temp4  
    eta223(j,i)= temp444
    eta224(j,i)= temp443  
    eta(j,i)=temp4*temp44+temp443*temp444
    eta(i,j)=temp4*temp44+temp443*temp444  
    term1=tau3*temp44+tau33*temp4+tau333*temp444+tau334*temp443
    term2=temp4*temp44+temp443*temp444
    rm2kl(j,i)=temp5*(ONE-TWO*ONETHIRD*term1/(m*temp2) + EIGHT*ONEFIFTH*term2/(THREE*m*temp2*temp2))/temp2
    rm2kl(i,j)=rm2kl(j,i)
    rmkl(j,i)=temp1*(ONE-ONETHIRD*term1/(m*temp2)+ONEFIFTH*term2/(m*temp2*temp2))/temp3
    !rmkl(j,i)=ME_over_rij(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
    rmkl(i,j)=rmkl(j,i)
    !rmkl1(i,j)=rmkl1(j,i)
    Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*rmkl(j,i)
    !Vkl1=Vkl1+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*rmkl1(j,i)
    rkl(j,i)=temp1*temp3*(ONE+ONETHIRD*term1/(m*temp2) - ONEFIFTH*term2/(THREE*m*temp2*temp2))
    rkl(i,j)=rkl(j,i)
    r2kl(j,i)=Skl*THREEHALF*temp2*(ONE+TWO*ONETHIRD*term1/(m*temp2))
    r2kl(i,j)=r2kl(j,i)
    temp10=temp8/(temp2*temp3)
    !deltarkl(j,i)=temp10*(ONE-term1/(m*temp2)+term2/(THREE*m*temp2*temp2))
    !deltarkl(i,j)=deltarkl(j,i)
    prvalkl(j,i)=PI*temp10*( TWO*(Glob_EulerConst+log(temp2))*(ONE-term1/(m*temp2)+term2/(m*temp2*temp2)) &
      + FOUR*(term1-TWO*term2/temp2)/(THREE*m*temp2)+EIGHT*term2/(15*m*temp2*temp2) )
    prvalkl(i,j)=prvalkl(j,i)
  enddo  
enddo
Hkl=Tkl+Vkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1/(r_ij*r_pq) is not implemented yet

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

!Evaluation of (1/r_{ij}*1/r_{pq}))_kl is not implemented yet
tau = tau3                                                                                     !tau
ttau = tau33   
myeta =  tau333
myteta = tau334
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
          if (i==j .and. p/=q) then 
            a = inv_tAkl(i,i)
            b = inv_tAkl(p,p) + inv_tAkl(q,q) -  inv_tAkl(p,q) - inv_tAkl(q,p) 
            d = inv_tAkl(i, p) - inv_tAkl(i, q)
            fij = tvkinv_tAkl(i)
            tfij = tbkinv_tAkl(i)
            uij = inv_tAkltvl(i)
            tuij = inv_tAkltbl(i)
            fpq =  tvkinv_tAkl(p) - tvkinv_tAkl(q)                                                            !fij
            tfpq = tbkinv_tAkl(p) - tbkinv_tAkl(q)   
            upq = inv_tAkltvl(p) - inv_tAkltvl(q)                                                  !uij  
            tupq = inv_tAkltbl(p) - inv_tAkltbl(q) 
          else if (i/=j .and. p==q) then
            a = inv_tAkl(i,i) + inv_tAkl(j,j) -  inv_tAkl(i,j) - inv_tAkl(j,i)
            b = inv_tAkl(p,p)
            d = inv_tAkl(i, p) - inv_tAkl(j, p)
            fij =  tvkinv_tAkl(i) -   tvkinv_tAkl(j)                                                     !fij
            tfij = tbkinv_tAkl(i) - tbkinv_tAkl(j)
            uij = inv_tAkltvl(i) - inv_tAkltvl(j)                                                             !uij  
            tuij = inv_tAkltbl(i) - inv_tAkltbl(j)
            fpq = tvkinv_tAkl(p)
            tfpq = tbkinv_tAkl(p)
            upq = inv_tAkltvl(p)
            tupq =inv_tAkltbl(p)
          else if (i==j .and. p==q) then
            a = inv_tAkl(i,i)
            b = inv_tAkl(p,p)
            d = inv_tAkl(i,p)
            fij = tvkinv_tAkl(i)                                                          !fij
            tfij = tbkinv_tAkl(i)
            uij = inv_tAkltvl(i)                                                          !uij  
            tuij = inv_tAkltbl(i)
            fpq = tvkinv_tAkl(p)
            tfpq = tbkinv_tAkl(p)
            upq = inv_tAkltvl(p)
            tupq = inv_tAkltbl(p)
          else 
            a = inv_tAkl(i,i) + inv_tAkl(j,j) -  inv_tAkl(i,j) - inv_tAkl(j,i)
            b = inv_tAkl(p,p) + inv_tAkl(q,q) -  inv_tAkl(p,q) - inv_tAkl(q,p)     !a
            d = inv_tAkl(i, p) + inv_tAkl(j, q) - inv_tAkl(i, q) - inv_tAkl(j, p)   !d
            fij = tvkinv_tAkl(i) - tvkinv_tAkl(j)                                                            !fij
            tfij = tbkinv_tAkl(i) - tbkinv_tAkl(j)  
            fpq = tvkinv_tAkl(p) - tvkinv_tAkl(q)                                                             !fij
            tfpq = tbkinv_tAkl(p) - tbkinv_tAkl(q) 
            uij = inv_tAkltvl(i) - inv_tAkltvl(j)                                                             !uij  
            tuij = inv_tAkltbl(i) - inv_tAkltbl(j)   
            upq = inv_tAkltvl(p) - inv_tAkltvl(q)                                                             !uij  
            tupq = inv_tAkltbl(p) - inv_tAkltbl(q)
          endif
          dsqab = d/(sqrt(a*b))
          xx = dsqab*dsqab
          phi = sqrt(ONE-xx)
          phi_sq = ONE - xx
          phi_cube = phi_sq*phi
          commonFactor = (Glob_Piraised3n2/PI)/(abs(det_tAkl)*sqrt(abs(det_tAkl)))

          if (xx < local_eps_for_xx) then

            !Zero-order correction
            radicalCommon = 1._dprec/(45._dprec*(a*b)**(2)*sqrt(a*b))
            R1 = 15._dprec*a*b*(3._dprec*a*b*tau - b*fij*uij - a*fpq*upq)*ttau
            R2 = b*(-15._dprec*a*b*tau + 9._dprec*b*fij*uij + 5._dprec*a*fpq*upq)*tfij*tuij
            R3 = a*(-15._dprec*a*b*tau + 5._dprec*b*fij*uij + 9._dprec*a*fpq*upq)*tfpq*tupq
            RDZeroOne = R1 + R2 + R3

            R1 = 15._dprec*a*b*(3._dprec*a*b*myeta - b*fij*tuij - a*fpq*tupq)*myteta
            R2 = b*(-15._dprec*a*b*myeta + 9._dprec*b*fij*tuij + 5._dprec*a*fpq*tupq)*tfij*uij
            R3 = a*(-15._dprec*a*b*myeta + 5._dprec*b*fij*tuij + 9._dprec*a*fpq*tupq)*tfpq*upq
            RDZeroTwo = R1 + R2 + R3

            !First-order correction
            R1 = (FIVE*a*b*ttau-THREE*b*tfij*tuij-THREE*a*tfpq*tupq)*(fij*upq+fpq*uij)
            R2 = (FIVE*a*b*tau-THREE*b*fij*uij-THREE*a*fpq*upq)*(tfij*tupq+tfpq*tuij)
            RDOneOne = (R1 + R2)*d

            R1 = (FIVE*a*b*myteta-THREE*b*tfij*uij-THREE*a*tfpq*upq)*(fij*tupq+fpq*tuij)
            R2 = (FIVE*a*b*myeta-THREE*b*fij*tuij-THREE*a*fpq*tupq)*(tfij*upq+tfpq*uij)
            RDOneTwo = (R1 + R2)*d

            !Second-order correction
            R1 = 25._dprec*a*b*ttau*(a*b*tau - b*fij*uij - a*fpq*upq)
            R2 = b*tfij*tuij*(-25._dprec*a*b*tau+25._dprec*b*fij*uij+21._dprec*a*fpq*upq)
            R3 = 6._dprec*a*b*tfij*tupq*(fpq*uij + fij*upq)
            R4 = -25._dprec*(a**2)*b*tau*tfpq*tupq
            R5 = 3._dprec*a*b*tfpq*fij*(2._dprec*upq*tuij + 7._dprec*uij*tupq)
            R6 = a*tfpq*fpq*(6._dprec*b*uij*tuij + 25._dprec*a*upq*tupq)
            RDTwoOne = (R1 + R2 + R3 + R4 + R5 + R6)*d**2/(150._dprec*(a*b)**3*sqrt(a*b))

            R1 = 25._dprec*a*b*myteta*(a*b*myeta - b*fij*tuij - a*fpq*tupq)
            R2 = b*tfij*uij*(-25._dprec*a*b*myeta+25._dprec*b*fij*tuij+21._dprec*a*fpq*tupq)
            R3 = 6._dprec*a*b*tfij*upq*(fpq*tuij + fij*tupq)
            R4 = -25._dprec*(a**2)*b*myeta*tfpq*upq
            R5 = 3._dprec*a*b*tfpq*fij*(2._dprec*tupq*uij + 7._dprec*tuij*upq)
            R6 = a*tfpq*fpq*(6._dprec*b*tuij*uij + 25._dprec*a*tupq*upq)
            RDTwoTwo = (R1 + R2 + R3 + R4 + R5 + R6)*d**2/(150._dprec*(a*b)**3*sqrt(a*b))


            !Third-order correction
            R1 = (a*b*ttau-b*tfij*tuij-a*tfpq*tupq)*(fij*upq+fpq*uij)
            R2 = (a*b*tau-b*fij*uij-a*fpq*upq)*(tfij*tupq+tfpq*tuij)
            RDThreeOne = (R1 + R2)*d**3/(10._dprec*(a*b)**3*sqrt(a*b))

            R1 = (a*b*myteta-b*tfij*uij-a*tfpq*upq)*(fij*tupq+fpq*tuij)
            R2 = (a*b*myeta-b*fij*tuij-fpq*tupq)*(tfij*upq+tfpq*uij)
            RDThreeTwo = (R1 + R2)*d**3/(10._dprec*(a*b)**3*sqrt(a*b))

            !Fourth-order correction
            R1 = 7._dprec*a*b*ttau*(3._dprec*a*b*tau - 5._dprec*b*fij*uij - 5._dprec*a*fpq*upq)
            R2 = b*tfij*tuij*(-35._dprec*a*b*tau + 49._dprec*b*fij*uij + 45._dprec*a*fpq*upq)
            R3 = 20._dprec*a*b*tfij*tupq*(fpq*uij+fij*upq)
            R4 = a*tfpq*fpq*(20._dprec*b*uij*tuij + 49._dprec*a*upq*tupq)
            R5 = -35._dprec*(a**2)*b*tau*tupq*tfpq
            R6 = 5._dprec*a*b*tfpq*fij*(4._dprec*upq*tuij + 9._dprec*uij*tupq)
            RDFourOne = (R1 + R2 + R3 + R4 + R5 + R6)*(d**4)/(280._dprec*((a*b)**4)*sqrt(a*b))

            R1 = 7._dprec*a*b*myteta*(3._dprec*a*b*myeta - 5._dprec*b*fij*tuij - 5._dprec*a*fpq*tupq)
            R2 = b*tfij*uij*(-35._dprec*a*b*myeta + 49._dprec*b*fij*tuij + 45._dprec*a*fpq*tupq)
            R3 = 20._dprec*a*b*tfij*upq*(fpq*tuij+fij*tupq)
            R4 = a*tfpq*fpq*(20._dprec*b*tuij*uij + 49._dprec*a*tupq*upq)
            R5 = -35._dprec*(a**2)*b*myeta*upq*tfpq
            R6 = 5._dprec*a*b*tfpq*fij*(4._dprec*tupq*uij + 9._dprec*tuij*upq)
            RDFourTwo = (R1 + R2 + R3 + R4 + R5 + R6)*(d**4)/(280._dprec*((a*b)**4)*sqrt(a*b))

            
            totalAns = commonFactor*radicalCommon*(RDZeroOne + RDZeroTwo + RDOneOne + RDOneTwo) + &
            commonFactor*(RDTwoOne + RDTwoTwo + RDThreeOne + RDThreeTwo + RDFourOne + RDFourTwo)

          else
            acosphi=asin(abs(dsqab))

            commonArccosRadical = 1._dprec/(15._dprec*abs(d)**3*sqrt(a*b)*(a*b)**3*phi_cube)
            !!! Calculation of arccos part  !!!
            arccosCommon =-(a*b)**3*sqrt(a*b)*phi_cube

            a1 = 5._dprec*d*ttau*(fpq*uij + fij*upq - 3._dprec*d*tau)
            a2 = (5._dprec*d*tau - 3._dprec*fij*upq - 3._dprec*fpq*uij)*(tfij*tupq + tfpq*tuij)
            a3 = 2._dprec*fij*uij*tfpq*tupq
            a4 = 2._dprec*fpq*upq*tfij*tuij
            aone = a1 + a2 + a3 + a4

            !Second term
            a1 = 5._dprec*d*myteta*(fpq*tuij + fij*tupq - 3._dprec*d*myeta)
            a2 = (5._dprec*d*myeta - 3._dprec*fij*tupq - 3._dprec*fpq*tuij)*(tfij*upq + tfpq*uij)
            a3 = 2._dprec*fij*tuij*tfpq*upq
            a4 = 2._dprec*fpq*tupq*tfij*uij
            atwo = a1 + a2 + a3 + a4

            arccosAns = arccosCommon * (aone + atwo) * acosphi


            !!! Calculation of radical part  !!!
            radicalCommon=abs(d)

            !First term
            R11 = (a*b*upq - b*d*uij)*fij
            R12 = (a*b*uij - a*d*upq)*fpq
            R1 = 5._dprec*(a**2)*(b**2)*d*(phi_sq)*ttau*(R11+R12)
            R21 = (d**2)*((2._dprec*(d**2) - 3._dprec*a*b)*uij + a*d*upq)*fij
            R22 = 5._dprec*(a**2)*b*(d**2)*(phi_sq)*tau
            R23 = a*((d**3)*uij + a*((d**2)-2._dprec*a*b)*upq)*fpq
            R2 = -(b**2)*tfij*tuij*(R21 + R22 + R23)
            R31 = -b*((d**3)*uij + a*(3._dprec*a*b-4._dprec*(d**2))*upq)*fij
            R32 = 5._dprec*(a**2)*(b**2)*d*(phi_sq)*tau
            R33 = -a*(b*(3._dprec*a*b-4._dprec*(d**2))*uij + (d**3)*upq)*fpq
            R3 = a*b*tfij*tupq*(R31 + R32 + R33)
            R4 = a*b*tfpq*tuij*(R31 + R32 + R33)
            R51 = (d**2)*((2._dprec*(d**2) - 3._dprec*a*b)*upq + b*d*uij)*fpq
            R52 = 5._dprec*a*(b**2)*(d**2)*(phi_sq)*tau
            R53 = b*((d**3) * upq + b*((d**2) - 2._dprec*a*b)*uij)*fij
            R5 = -(a**2)*tfpq*tupq*(R51 + R52 + R53)
            ROne = R1 + R2 + R3 + R4 + R5

            !Second term
            R11 = (a*b*tupq - b*d*tuij)*fij
            R12 = (a*b*tuij - a*d*tupq)*fpq
            R1 = 5._dprec*(a**2)*(b**2)*d*(phi_sq)*myteta*(R11 + R12)
            R21 = (d**2)*((2._dprec*(d**2) - 3._dprec*a*b)*tuij + a*d*tupq)*fij
            R22 = 5._dprec*(a**2)*b*(d**2)*(phi_sq)*myeta
            R23 = a*((d**3)*tuij + a*((d**2)-2._dprec*a*b)*tupq)*fpq
            R2 = -(b**2)*tfij*uij*(R21 + R22 + R23)
            R31 = -b*((d**3)*tuij + a*(3._dprec*a*b-4._dprec*(d**2))*tupq)*fij
            R32 = 5._dprec*(a**2)*(b**2)*d*(phi_sq)*myeta
            R33 = -a*(b*(3._dprec*a*b-4._dprec*(d**2))*tuij + (d**3)*tupq)*fpq
            R3 = a*b*tfij*upq*(R31 + R32 + R33)
            R4 = a*b*tfpq*uij*(R31 + R32 + R33)
            R51 = (d**2)*((2._dprec*(d**2) - 3._dprec*a*b)*tupq + b*d*tuij)*fpq
            R52 = 5._dprec*a*(b**2)*(d**2)*(phi_sq)*myeta
            R53 = b*((d**3) * tupq + b*((d**2) - 2._dprec*a*b)*tuij)*fij
            R5 = -(a**2)*tfpq*upq*(R51 + R52 + R53)
            RTwo = R1 + R2 + R3 + R4 + R5

            radicalAns = radicalCommon * (ROne + RTwo)

           
            totalAns = commonFactor * commonArccosRadical * (arccosAns + radicalAns)

          endif
          rmrmkl(i,j,p,q)=totalAns
          rmrmkl(j,i,p,q)=totalAns
          rmrmkl(i,j,q,p)=totalAns
          rmrmkl(j,i,q,p)=totalAns
          rmrmkl(p,q,i,j)=totalAns
          rmrmkl(p,q,j,i)=totalAns
          rmrmkl(q,p,i,j)=totalAns
          rmrmkl(q,p,j,i)=totalAns
        endif  
      enddo
    enddo
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!evaluate delta-functions directly
!V ---- tau
!tV --- myeta
!W ---- ttau
!tW --- myteta
deltarkl = ZERO
do i = 1,n 
  do j = i,n
    if (i==j) then 
      myalpha = inv_tAkl(i,i)
      jijAvk = tvkinv_tAkl(i)
      jijAvl = inv_tAkltvl(i)
      jijAwk = tbkinv_tAkl(i)
      jijAwl = inv_tAkltbl(i)
    else
      myalpha = inv_tAkl(i,i) + inv_tAkl(j,j) - inv_tAkl(i,j) - inv_tAkl(j,i)
      jijAvk = tvkinv_tAkl(i) - tvkinv_tAkl(j)
      jijAvl = inv_tAkltvl(i) - inv_tAkltvl(j)
      jijAwk = tbkinv_tAkl(i) - tbkinv_tAkl(j)
      jijAwl = inv_tAkltbl(i) - inv_tAkltbl(j)
    endif
    R1 = (myalpha**2)*(tau*ttau + myeta*myteta)
    R2 = -myalpha*(jijAvk*jijAvl*ttau + tau*jijAwk*jijAwl+ &
    jijAvk*jijAwl*myteta + myeta*jijAwk*jijAvl)
    R3 = TWO*(jijAvk*jijAvl*jijAwk*jijAwl)

    deltarkl(i,j) = (R1 + R2 + R3)*&
    Glob_Piraised3n2/(FOUR*PI*SQRTPI*(myalpha**3)*sqrt(myalpha)*det_tAkl*sqrt(det_tAkl))

    deltarkl(j,i) = deltarkl(i,j)
  enddo
enddo

!evaluate drachmanized delta-function and V2kl operator
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
    !temp2=ZERO
    temp2 = ME_d_X_over_rij_d(p,q,Glob_dmvM,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl, &
    tvkinv_tAkl, tbkinv_tAkl, inv_tAkltvl, inv_tAkltbl)
    drach_deltarkl(p,q)=(Glob_CurrEnergy*rmkl(p,q)-temp1-temp2)/temp4
    !drach_deltarkl(p,q)=temp2
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
!MVkl=ME_dWd2(W1,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)/temp1
MVkl=dXddYd(W1,W1,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)/temp1
W1(1:n,1:n)=ZERO
do i=1,n 
  W1(i,i)=ONE
  temp1=Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1)
  !MVkl=MVkl+ME_dWd2(W1,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)/temp1
  MVkl=MVkl+dXddYd(W1,W1,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)/temp1
  W1(i,i)=ZERO
enddo
MVkl=-MVkl/8


drach_MVkl=ZERO		   
!drach_MVkl=ME_dWd21(Glob_dmvM,Glob_dmvMB,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl) &
drach_MVkl=dXddYd(Glob_dmvM,Glob_dmvMB,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)&
  -V2kl-Glob_CurrEnergy*Glob_CurrEnergy*Skl+2*Glob_CurrEnergy*Vkl+ &
  !Glob_CurrEnergy*ME_dXd(Glob_dmvB,tvk,tvl,inv_tAkltvl,inv_tAkl,tAk,tAl,inv_tAkltAl,Skl,tau3)!+&
  Glob_CurrEnergy*dXddYd(Glob_dmvB,Glob_dmvB,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334, &
                         inv_tAkltAl,inv_tAkltAk)
do i=1,n                                        
  drach_MVkl=drach_MVkl-ScaledChargeProd(Glob_PseudoCharge0,Glob_PseudoCharge(i)) &
                        *ME_d_X_over_rij_d(i,i,Glob_dmvB,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl,&
                        tvkinv_tAkl, tbkinv_tAkl, inv_tAkltvl, inv_tAkltbl)           
                                   
  do j=i+1,n                                                                  
  drach_MVkl=drach_MVkl-ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j)) &
                        *ME_d_X_over_rij_d(i,j,Glob_dmvB,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl,&
                        tvkinv_tAkl, tbkinv_tAkl, inv_tAkltvl, inv_tAkltbl)                            
  enddo                
enddo
drach_MVkl = drach_MVkl*Glob_dmva2 + MVkl
   
if (isOOklNeeded) then
  do i=1,n
    do j=1,n
      tbltbk(i,j)=tbl(j)*tbk(i)
      tvltvk(i,j)=tvl(j)*tvk(i)
      tbltvk(i,j)=tbl(j)*tvk(i)
      tvltbk(i,j)=tvl(j)*tbk(i)
      
      tvltbl(i,j)=tvl(j)*tbl(i)
      tvktbk(i,j)=tvk(j)*tbk(i)
    enddo
  enddo


  call symmetrize_matrix(tbltbk)
  call symmetrize_matrix(tvltvk)
  call symmetrize_matrix(tvltbk)
  call symmetrize_matrix(tbltvk)
  call symmetrize_matrix(tvltbl)
  call symmetrize_matrix(tvktbk)


  !Evaluating Orbit-Orbit (OO) matrix element (without the factor of alpha**2)
  OOkl=ZERO
  !First double loop for OO
  do i=1,n
    do j=1,n  
      tr1=tAl(j,i)
      tr2=tAl(i,j)
      tr3=3*tAl(j,j)
      do p=1,n
        do q=1,n
          W1(p,q)=tAl(p,i)*tAl(j,q)+tAl(p,j)*tAl(i,q)
        enddo
      enddo
      do p=1,n
        W1(p,j)=W1(p,j)+tAkl(p,j)*tAl(j,i)      
      enddo

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

      do p=1,n
        do q=1,n
          W4(p,q) = 2*tAl(p,j)*tvl(i)*tvk(q) + 2*tAl(p,i)*tvl(j)*tvk(q)
          W44(p,q) = 2*tAl(p,j)*tbk(i)*tvk(q) + 2*tAl(p,i)*tbk(j)*tvk(q)        
          W4b(p,q) = 2*tAl(p,j)*tbl(i)*tbk(q) + 2*tAl(p,i)*tbl(j)*tbk(q)
          W44b(p,q) = 2*tAl(p,j)*tvk(i)*tbk(q) + 2*tAl(p,i)*tvk(j)*tbk(q)
        enddo
      enddo
      !W4 = W4 + vl vk' Ejj Al Eij
      do p=1,n
        W4(p,j) = W4(p,j) + tvl(p)*tvk(j)*tAl(j,i)
        W44(p,j) = W44(p,j) + tbk(p)*tvk(j)*tAl(j,i)     
        W4b(p,j) = W4b(p,j) + tbl(p)*tbk(j)*tAl(j,i)
        W44b(p,j) = W44b(p,j) + tvk(p)*tbk(j)*tAl(j,i)
      enddo
      do q=1,n
        W4(j,q) = W4(j,q) + 2*tAl(i,j)*tvl(j)*tvk(q) + tr3*tvl(i)*tvk(q) + tvl(i)*tvk(j)*tAl(j,q)     
        W44(j,q) = W44(j,q) + 2*tAl(i,j)*tbk(j)*tvk(q) + tr3*tbk(i)*tvk(q) + tbk(i)*tvk(j)*tAl(j,q)      
        W4b(j,q) = W4b(j,q) + 2*tAl(i,j)*tbl(j)*tbk(q) + tr3*tbl(i)*tbk(q) + tbl(i)*tbk(j)*tAl(j,q)     
        W44b(j,q) = W44b(j,q) + 2*tAl(i,j)*tvk(j)*tbk(q) + tr3*tvk(i)*tbk(q) + tvk(i)*tbk(j)*tAl(j,q)
      enddo

      do p=1,n
        do q=1,n
          W5(p,q) = tAkl(p,j)*tvl(j)*tvk(q) + tAl(p,j)*tvl(j)*tvk(q) + tvl(p)*tvk(j)*tAl(j,q)
          W55(p,q) = tAkl(p,j)*tbk(j)*tvk(q) + tAl(p,j)*tbk(j)*tvk(q) + tbk(p)*tvk(j)*tAl(j,q)       
          W5b(p,q) = tAkl(p,j)*tbl(j)*tbk(q) + tAl(p,j)*tbl(j)*tbk(q) + tbl(p)*tbk(j)*tAl(j,q)
          W55b(p,q) = tAkl(p,j)*tvk(j)*tbk(q) + tAl(p,j)*tvk(j)*tbk(q) + tvk(p)*tbk(j)*tAl(j,q)
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
      W77(1:n,1:n)=ZERO
      W7b(1:n,1:n)=ZERO
      W77b(1:n,1:n)=ZERO
      do q=1,n
        W7(j,q) = tvl(i)*tvk(q)
        W77(j,q) = tbl(i)*tvl(q)
        W7b(j,q) = tbl(i)*tbk(q)
        W77b(j,q) = tvl(i)*tbl(q)
      enddo
      
      call symmetrize_matrix(W1)
      call symmetrize_matrix(W2)
      call symmetrize_matrix(W3)
      call symmetrize_matrix(W4)
      call symmetrize_matrix(W5)
      call symmetrize_matrix(W44)
      call symmetrize_matrix(W55)
      call symmetrize_matrix(W6)
      call symmetrize_matrix(W7)
      call symmetrize_matrix(W77)    
      call symmetrize_matrix(W4b)
      call symmetrize_matrix(W5b)
      call symmetrize_matrix(W44b)
      call symmetrize_matrix(W55b)
      call symmetrize_matrix(W7b)
      call symmetrize_matrix(W77b)
      !compute integrals   
                            temp1=ME_rXr_over_rij(j,j,W1,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
                        temp2=ME_rXr_rYr_over_rij(j,j,W2,W3,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)              
            temp4=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W4,tbltbk,inv_tAkl,det_tAkl)
          temp44=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W44,tvltbl,inv_tAkl,det_tAkl)
          temp4b=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W4b,tvltvk,inv_tAkl,det_tAkl)
          temp44b=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W44b,tvltbl,inv_tAkl,det_tAkl)    
          temp444=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W7,W5b,inv_tAkl,det_tAkl)    
        temp4444=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W77,W55b,inv_tAkl,det_tAkl)    
        temp444b=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W7b,W5,inv_tAkl,det_tAkl)
        temp4444b=ONETHIRD*SG_ME_rXr_rYr_over_rij(j,j,W77b,W55,inv_tAkl,det_tAkl)    
        temp5=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W3,W5,tbltbk,inv_tAkl,det_tAkl)
      temp55=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W3,W55,tvltbl,inv_tAkl,det_tAkl)  
      temp5b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W3,W5b,tvltvk,inv_tAkl,det_tAkl)
      temp55b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W3,W55b,tvltbl,inv_tAkl,det_tAkl) 
        temp6=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W6,W7,tbltbk,inv_tAkl,det_tAkl)
      temp66=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W6,W77,tvktbk,inv_tAkl,det_tAkl)
      temp6b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W6,W7b,tvltvk,inv_tAkl,det_tAkl)
      temp66b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(j,j,W6,W77b,tvktbk,inv_tAkl,det_tAkl)
      temp7=-6*(tr1+tr2)*rmkl(j,j)+4*temp1-8*temp2-2*(temp4-temp44)+4*(temp5-temp55)+4*(temp6-temp66)    
      temp7=temp7-2*(temp4b-temp44b)+4*(temp5b-temp55b)+4*(temp6b-temp66b)-2*(temp444-temp4444)-2*(temp444b-temp444b)
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

      do p=1,n
        do q=1,n
          W4(p,q) = 2*tAl(p,j)*tvl(i)*tvk(q) + 2*tAl(p,i)*tvl(j)*tvk(q)
          W44(p,q) = 2*tAl(p,j)*tbk(i)*tvk(q) + 2*tAl(p,i)*tbk(j)*tvk(q)        
          W4b(p,q) = 2*tAl(p,j)*tbl(i)*tbk(q) + 2*tAl(p,i)*tbl(j)*tbk(q)
          W44b(p,q) = 2*tAl(p,j)*tvk(i)*tbk(q) + 2*tAl(p,i)*tvk(j)*tbk(q)
        enddo
      enddo
      !W4 = W4 + vl vk' Ejj Al (Eij-Eii)
      do p=1,n
        W4(p,j) = W4(p,j) + tvl(p)*tvk(j)*tAl(j,i)
        W4(p,i) = W4(p,i) - tvl(p)*tvk(j)*tAl(j,i)
        W44(p,j) = W44(p,j) + tbk(p)*tvk(j)*tAl(j,i)
        W44(p,i) = W44(p,i) - tbk(p)*tvk(j)*tAl(j,i)      
        W4b(p,j) = W4b(p,j) + tbl(p)*tbk(j)*tAl(j,i)
        W4b(p,i) = W4b(p,i) - tbl(p)*tbk(j)*tAl(j,i)     
        W44b(p,j) = W44b(p,j) + tvk(p)*tbk(j)*tAl(j,i)
        W44b(p,i) = W44b(p,i) - tvk(p)*tbk(j)*tAl(j,i)      
      enddo

      do q=1,n
        temp1 = 2*tAl(i,j)*tvl(j)*tvk(q) + tr3*tvl(i)*tvk(q) + tvl(i)*tvk(j)*tAl(j,q)      
        temp11 = 2*tAl(i,j)*tbk(j)*tvk(q) + tr3*tbk(i)*tvk(q) + tbk(i)*tvk(j)*tAl(j,q)     
        temp1b = 2*tAl(i,j)*tbl(j)*tbk(q) + tr3*tbl(i)*tbk(q) + tbl(i)*tbk(j)*tAl(j,q)
        temp11b = 2*tAl(i,j)*tvk(j)*tbk(q) + tr3*tvk(i)*tbk(q) + tvk(i)*tbk(j)*tAl(j,q)
        W4(j,q) = W4(j,q) + temp1
        W4(i,q) = W4(i,q) - temp1
        W44(j,q) = W44(j,q) + temp11
        W44(i,q) = W44(i,q) - temp11      
        W4b(j,q) = W4b(j,q) + temp1b
        W4b(i,q) = W4b(i,q) - temp1b
        W44b(j,q) = W44b(j,q) + temp11b
        W44b(i,q) = W44b(i,q) - temp11b
      enddo
      
      !W5 = Akl Ejj vl vk' + Al Ejj vl vk' + vl vk' Ejj Al
      do p=1,n
        do q=1,n
          W5(p,q) = tAkl(p,j)*tvl(j)*tvk(q) + tAl(p,j)*tvl(j)*tvk(q) + tvl(p)*tvk(j)*tAl(j,q)        
          W55(p,q) = tAkl(p,j)*tbk(j)*tvk(q) + tAl(p,j)*tbk(j)*tvk(q) + tbk(p)*tvk(j)*tAl(j,q)       
          W5b(p,q) = tAkl(p,j)*tbl(j)*tbk(q) + tAl(p,j)*tbl(j)*tbk(q) + tbl(p)*tbk(j)*tAl(j,q)
          W55b(p,q) = tAkl(p,j)*tvk(j)*tbk(q) + tAl(p,j)*tvk(j)*tbk(q) + tvk(p)*tbk(j)*tAl(j,q)
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
      W77(1:n,1:n)=ZERO
      W7b(1:n,1:n)=ZERO
      W77b(1:n,1:n)=ZERO
      do q=1,n
        W7(j,q) = tvl(i)*tvk(q)
        W7(i,q) = -tvl(i)*tvk(q)      
        W77(j,q) = tbl(i)*tvl(q)
        W77(i,q) = -tbl(i)*tvl(q)
        W7b(j,q) = tbl(i)*tbk(q)
        W7b(i,q) = -tbl(i)*tbk(q)
        W77b(j,q) = tvl(i)*tbl(q)
        W77b(i,q) = -tvl(i)*tbl(q)
      enddo
      
      call symmetrize_matrix(W1)
      call symmetrize_matrix(W2)
      call symmetrize_matrix(W3)
      call symmetrize_matrix(W4)
      call symmetrize_matrix(W44)
      call symmetrize_matrix(W5)
      call symmetrize_matrix(W55)
      call symmetrize_matrix(W6)
      call symmetrize_matrix(W7)
      call symmetrize_matrix(W77)
      call symmetrize_matrix(W4b)
      call symmetrize_matrix(W5b)
      call symmetrize_matrix(W44b)
      call symmetrize_matrix(W55b)
      call symmetrize_matrix(W7b)
      call symmetrize_matrix(W77b)
      !compute integrals (remember to take the difference)
      temp1=ME_rXr_over_rij(i,j,W1,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
      temp2=ME_rXr_rYr_over_rij(i,j,W2,W3,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
        temp4=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W4,tbltbk,inv_tAkl,det_tAkl)
        temp44=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W44,tvltbl,inv_tAkl,det_tAkl)     
        temp4b=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W4b,tvltvk,inv_tAkl,det_tAkl)     
      temp44b=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W44b,tvltbl,inv_tAkl,det_tAkl)      
      temp444=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W7,W5b,inv_tAkl,det_tAkl)    
      temp4444=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W77,W55b,inv_tAkl,det_tAkl)    
      temp444b=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W7b,W5,inv_tAkl,det_tAkl)
    temp4444b=ONETHIRD*SG_ME_rXr_rYr_over_rij(i,j,W77b,W55,inv_tAkl,det_tAkl)    
    temp5=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W3,W5,tbltbk,inv_tAkl,det_tAkl)
    temp55=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W3,W55,tvltbl,inv_tAkl,det_tAkl)     
    temp5b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W3,W5b,tvltvk,inv_tAkl,det_tAkl)
  temp55b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W3,W55b,tvltbl,inv_tAkl,det_tAkl)     
    temp6=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W6,W7,tbltbk,inv_tAkl,det_tAkl)
    temp66=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W6,W77,tvktbk,inv_tAkl,det_tAkl)    
    temp6b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W6,W7b,tvltvk,inv_tAkl,det_tAkl)
  temp66b=ONETHIRD*SG_ME_rXr_rYr_rZr_over_rij(i,j,W6,W77b,tvktbk,inv_tAkl,det_tAkl)
      temp7=-6*(tr1+tr2)*rmkl(i,j)+4*temp1-8*temp2-2*(temp4-temp44)+4*(temp5-temp55)+4*(temp6-temp66)
      temp7=temp7-2*(temp4b-temp44b)+4*(temp5b-temp55b)+4*(temp6b-temp66b)-2*(temp444-temp4444)-2*(temp444b-temp4444b)
      OOkl=OOkl+&
        temp7*ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))/(Glob_Mass(i+1)*Glob_Mass(j+1))
    enddo
  enddo
  OOkl=OOkl/2
else
  OOkl = ZERO
endif
!!

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

end subroutine MatrixElementsL1ForExpcValsD

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


function SG_ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,det_tAkl)
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
real(dprec)   t_XYJ,t_YXJ,det_tAkl

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

temp1=Glob_Piraised3n2/(SQRTPI*det_tAkl**(THREEHALF))
temp3=1/t_J
SG_ME_rXr_rYr_over_rij=THREE*temp1*temp3*sqrt(temp3)*(  &
THREEHALF*t_J*t_X*t_Y - ONEHALF*(t_Y*t_XJ + t_X*t_YJ) + &
t_J*t_XY  - ONETHIRD*(t_XYJ + t_YXJ) + ONEHALF*temp3*t_XJ*t_YJ&
)

end function SG_ME_rXr_rYr_over_rij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function SG_ME_rXr_rYr_rZr_over_rij(i,j,X,Y,Z,inv_tAkl,det_tAkl)
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
real(dprec)   SG_ME_rXr_rYr_rZr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),Y(nn,nn),Z(nn,nn),inv_tAkl(nn,nn)
integer       i,j
!Local variables:
integer       p,q,s,n
real(dprec)   temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12
real(dprec)   AX(nn,nn),AY(nn,nn),AZ(nn,nn)
real(dprec)   Aj(nn),AjX(nn),AjY(nn),AjZ(nn),AXAj(nn),AYAj(nn),AZAj(nn)
real(dprec)   AZAY(nn,nn),AYAZ(nn,nn),AZAX(nn,nn),AXAZ(nn,nn),AXAY(nn,nn),AYAX(nn,nn)
real(dprec)   AZAYAj(nn),AYAZAj(nn),AZAXAj(nn),AXAZAj(nn),AXAYAj(nn),AYAXAj(nn)
real(dprec)   t_J,t_X,t_Y,t_Z
real(dprec)   t_XJ,t_YJ,t_ZJ
real(dprec)   t_XY,t_ZY,t_ZX,t_YX,t_ZYX,t_YZX
real(dprec)   t_XYJ,t_YXJ,t_ZYJ,t_YZJ,t_XZJ,t_ZXJ
real(dprec)   t_ZYXJ,t_YZXJ,t_YXZJ,t_XYZJ,t_ZXYJ,t_XZYJ
real(dprec)   det_tAkl,term1,term2,term3,term4,term5,term6,term7,term8

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
  temp3=ZERO
  do q=1,n
    temp1=temp1+Aj(q)*X(q,p)
    temp2=temp2+Aj(q)*Y(q,p)
    temp3=temp3+Aj(q)*Z(q,p)
  enddo
  AjX(p)=temp1
  AjY(p)=temp2
  AjZ(p)=temp3
enddo

!Compute AX=inv_tAkl*X  t_X=tr[inv_tAkl*X]
!        AY=inv_tAkl*Y  t_Y=tr[inv_tAkl*Y]
t_X=ZERO
t_Y=ZERO
t_Z=ZERO
do p=1,n
  do q=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    do s=1,n
      temp1=temp1+inv_tAkl(s,q)*X(p,s)
      temp2=temp2+inv_tAkl(s,q)*Y(p,s)
      temp3=temp3+inv_tAkl(s,q)*Z(p,s)
    enddo
    AX(q,p)=temp1
    AY(q,p)=temp2
    AZ(q,p)=temp3
  enddo
  t_X=t_X+AX(p,p)
  t_Y=t_Y+AY(p,p)
  t_Z=t_Z+AZ(p,p)
enddo

do p=1,n
  do q=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    temp4=ZERO
    temp5=ZERO
    temp6=ZERO
    do s=1,n
      temp1=temp1+AZ(s,q)*AY(p,s)
      temp2=temp2+AY(s,q)*AZ(p,s)
      temp3=temp3+AZ(s,q)*AX(p,s)
      temp4=temp4+AX(s,q)*AZ(p,s)
      temp5=temp5+AX(s,q)*AY(p,s)
      temp6=temp6+AY(s,q)*AX(p,s)
    enddo
    AZAY(q,p)=temp1
    AYAZ(q,p)=temp2
    AZAX(q,p)=temp3
    AXAZ(q,p)=temp4
    AXAY(q,p)=temp5
    AYAX(q,p)=temp6
  enddo
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
t_ZX=ZERO
t_ZY=ZERO
t_ZYX=ZERO
t_YZX=ZERO
do p=1,n
  temp3=ZERO
  temp4=ZERO
  temp5=ZERO  
  temp6=ZERO
  temp7=ZERO
  temp8=ZERO
  temp9=ZERO
  temp10=ZERO
  temp11=ZERO
  do q=1,n
    t_XY=t_XY+AX(p,q)*AY(q,p)
    t_ZX=t_ZX+AZ(p,q)*AX(q,p)
    t_ZY=t_ZY+AZ(p,q)*AY(q,p)
    t_ZYX=t_ZYX+AZAY(p,q)*AX(q,p)
    t_YZX=t_YZX+AYAZ(p,q)*AX(q,p)
    temp3=temp3+AX(p,q)*Aj(q)
    temp4=temp4+AY(p,q)*Aj(q)
    temp5=temp5+AZ(p,q)*Aj(q)  
    temp6=temp6+AZAY(p,q)*Aj(q)
    temp7=temp7+AYAZ(p,q)*Aj(q)
    temp8=temp8+AZAX(p,q)*Aj(q)
    temp9=temp9+AXAZ(p,q)*Aj(q)
    temp10=temp10+AXAY(p,q)*Aj(q)
    temp11=temp11+AYAX(p,q)*Aj(q)
  enddo
  AXAj(p)=temp3
  AYAj(p)=temp4
  AZAj(p)=temp5
  AZAYAj(p)=temp6
  AYAZAj(p)=temp7
  AZAXAj(p)=temp8
  AXAZAj(p)=temp9
  AXAYAj(p)=temp10
  AYAXAj(p)=temp11
enddo

!Compute 
!t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_YJ=tr[inv_tAkl*Y*inv_tAkl*Jij]=AjY'*Aj 
!t_XYJ=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij]=AjX'*AYAj  

t_XJ=ZERO
t_YJ=ZERO
t_ZJ=ZERO
t_XYJ=ZERO
t_ZYJ=ZERO
t_ZXJ=ZERO
t_ZYXJ=ZERO
t_YZXJ=ZERO
t_ZXYJ=ZERO
t_XZYJ=ZERO
t_YXZJ=ZERO
t_XYZJ=ZERO
do p=1,n
  t_XJ=t_XJ+AjX(p)*Aj(p)
  t_YJ=t_YJ+AjY(p)*Aj(p)
  t_ZJ=t_ZJ+AjZ(p)*Aj(p)  
  t_XYJ=t_XYJ+AjX(p)*AYAj(p)
  t_ZYJ=t_ZYJ+AjZ(p)*AYAj(p)
  t_ZXJ=t_ZXJ+AjZ(p)*AXAj(p) 
  t_ZYXJ=t_ZYXJ+AjZ(p)*AYAXAj(p)
  t_YZXJ=t_YZXJ+AjY(p)*AZAXAj(p)
  t_ZXYJ=t_ZXYJ+AjZ(p)*AXAYAj(p)
  t_XZYJ=t_XZYJ+AjX(p)*AZAYAj(p)
  t_YXZJ=t_YXZJ+AjY(p)*AXAZAj(p)
  t_XYZJ=t_XYZJ+AjX(p)*AYAZAj(p)
enddo

!Compute t_YXJ=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij]
t_YXJ=t_XYJ
t_YZJ=t_ZYJ
t_XZJ=t_ZXJ



temp1=Glob_Piraised3n2/(SQRTPI*det_tAkl**(THREEHALF))
temp3=1/t_J
term1=THREE*temp1*temp3*sqrt(temp3)*THREEHALF*t_Z*(  &
THREEHALF*t_J*t_X*t_Y - ONEHALF*(t_Y*t_XJ + t_X*t_YJ) + &
t_J*t_XY  - ONETHIRD*(t_XYJ + t_YXJ) + ONEHALF*temp3*t_XJ*t_YJ&
)
term2=-THREEHALF*t_ZY*(3*t_X/sqrt(t_J)-t_XJ/(t_J*sqrt(t_J)))
term3=THREEHALF*t_Y*(THREEHALF*t_X*t_ZJ/(t_J*sqrt(t_J))-3*t_ZX/sqrt(t_J)&
     -THREEHALF*t_XJ*t_ZJ/(sqrt(t_J)*t_J*t_J)+(t_ZXJ+t_XZJ)/(t_J*sqrt(t_J)))
term4= THREEHALF*t_YX*t_ZJ/(t_J*sqrt(t_J))
term5=-3*(t_ZYX+t_YZX)/sqrt(t_J)-9*t_X*t_YJ*t_ZJ/(4*t_J*t_J*sqrt(t_J))
term6=THREEHALF*(t_ZX*t_YJ+t_X*(t_ZYJ+t_YZJ))/(t_J*sqrt(t_J))+15*t_XJ*t_YJ*t_ZJ/(4*t_J*t_J*t_J*sqrt(t_J))
term7=-THREEHALF*(t_YJ*(t_ZXJ+t_XZJ)+t_XJ*(t_ZYJ+t_YZJ)+t_ZJ*(t_YXJ+t_XYJ))/(t_J*t_J*sqrt(t_J))
term8=(t_ZYXJ+t_YZXJ+t_YXZJ+t_XYZJ+t_ZXYJ+t_XZYJ)/(t_J*sqrt(t_J))

SG_ME_rXr_rYr_rZr_over_rij=term1-temp1*(term2+term3+term4+term5+term6+term7+term8)
end function SG_ME_rXr_rYr_rZr_over_rij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ME_rXr_over_rij(i,j,X,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
real(dprec)   X(nn,nn),inv_tAkl(nn,nn),det_tAkl
integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

real(dprec)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
real(dprec)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
real(dprec)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
integer       p,q,n
real(dprec)   Aj(nn),AjX(nn)
real(dprec)   t_J,t_X,t_XJ
real(dprec)   t_JV1,t_XV1,t_JXV1,t_XJV1
real(dprec)   t_JV2,t_XV2,t_JXV2,t_XJV2
real(dprec)   t_JV5,t_XV5,t_JXV5,t_XJV5
real(dprec)   t_JV6,t_XV6,t_JXV6,t_XJV6
real(dprec)   Ajtvl,Ajtbl
real(dprec)    temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
real(dprec)    temp11,temp22,temp33,mu,mX,mXJ,u,Xmu,muXJ

n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do q=1,n
    temp1=temp1+inv_tAkl(q,p)*tvl(q)
    temp2=temp2+inv_tAkl(q,p)*tbl(q)
    temp3=temp3+inv_tAkl(q,p)*tbk(q)
  enddo
  inv_tAkltvl(p)=temp1
  inv_tAkltbl(p)=temp2
  inv_tAkltbk(p)=temp3
enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do q=1,n
    temp1=temp1+tvk(q)*inv_tAkl(q,p)
    temp2=temp2+tbk(q)*inv_tAkl(q,p)
    temp3=temp3+tvl(q)*inv_tAkl(q,p)
  enddo
  tvkinv_tAkl(p)=temp1
  tbkinv_tAkl(p)=temp2
  tvlinv_tAkl(p)=temp3
enddo


!Compute tau3=tvkinv_tAkl*tvl
tau3=ZERO
tau33=ZERO
tau333=ZERO
tau334=ZERO
do p=1,n
  tau3=tau3+tvkinv_tAkl(p)*tvl(p)
  tau33=tau33+tbkinv_tAkl(p)*tbl(p)
  tau333=tau333+tvkinv_tAkl(p)*tbl(p)
  tau334=tau334+tbkinv_tAkl(p)*tvl(p)
enddo
m1=tau3*tau33
m3=tau333*tau334
m=m1+m3  !look formula 40 in document

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
Ajtbl=ZERO
t_XJ=ZERO
temp1=ZERO
temp2=ZERO
temp11=ZERO
temp22=ZERO
do p=1,n
  Ajtvl=Ajtvl+Aj(p)*tvl(p)
  Ajtbl=Ajtbl+Aj(p)*tbl(p)  
  t_XJ=t_XJ+AjX(p)*Aj(p) 
  temp1=temp1+tvk(p)*Aj(p)
  temp2=temp2+AjX(p)*inv_tAkltvl(p) 
  temp11=temp11+tbk(p)*Aj(p)
  temp22=temp22+AjX(p)*inv_tAkltbl(p)
enddo
t_JV1=temp11*Ajtbl
t_JXV1=temp11*temp22
t_JV2=temp1*Ajtvl
t_JXV2=temp1*temp2
t_JV5=temp11*Ajtvl
t_JXV5=temp11*temp2
t_JV6=temp1*Ajtbl
t_JXV6=temp1*temp22
!Compute t_X=tr[inv_tAkl*X]
!        t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*inv_tAkltvl
!        t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*Aj*Ajtvl
t_X=ZERO
t_XV1=ZERO
t_XV2=ZERO
t_XV5=ZERO
t_XV6=ZERO
temp2=ZERO
temp22=ZERO
do p=1,n
  temp1=ZERO
  temp11=ZERO
  do q=1,n
    t_X=t_X+inv_tAkl(q,p)*X(q,p)    
    temp1=temp1+tvkinv_tAkl(q)*X(q,p)
    temp11=temp11+tbkinv_tAkl(q)*X(q,p)
  enddo
  t_XV1=t_XV1+temp11*inv_tAkltbl(p)
  t_XV2=t_XV2+temp1*inv_tAkltvl(p)
  t_XV5=t_XV5+temp11*inv_tAkltvl(p)
  t_XV6=t_XV6+temp1*inv_tAkltbl(p)
  temp2=temp2+temp1*Aj(p)
  temp22=temp22+temp11*Aj(p)
enddo
t_XJV1=temp22*Ajtbl
t_XJV2=temp2*Ajtvl
t_XJV5=temp22*Ajtvl
t_XJV6=temp2*Ajtbl

mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
mX=tau3*t_XV1+tau33*t_XV2+tau333*t_XV5+tau334*t_XV6
u=t_JV1*t_JV2+t_JV5*t_JV6
mXJ=tau3*(t_XJV1+t_JXV1)+tau33*(t_XJV2+t_JXV2)+tau333*(t_XJV5+t_JXV5)+tau334*(t_XJV6+t_JXV6)
muXJ=t_JV2*(t_XJV1+t_JXV1)+t_JV1*(t_XJV2+t_JXV2)+t_JV6*(t_XJV5+t_JXV5)+t_JV5*(t_XJV6+t_JXV6)
Xmu=t_XV2*t_JV1+t_XV1*t_JV2+t_XV6*t_JV5+t_XV5*t_JV6


temp1=Glob_Piraised3n2/(TWO*SQRTPI*det_tAkl**(THREEHALF))
temp2=+THREEHALF*t_X*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)
temp3=-m*t_XJ/(TWO*t_J*sqrt(t_J))
temp4=+mX/sqrt(t_J)
temp5=+mu*t_XJ/(TWO*t_J*t_J*sqrt(t_J))
temp6=-(mXJ+Xmu)/(THREE*t_J*sqrt(t_J))
temp7=-u*t_XJ/(TWO*t_J*t_J*t_J*sqrt(t_J))
temp8=+ONEFIFTH*muXJ/(t_J*t_J*sqrt(t_J))
ME_rXr_over_rij=temp1*(temp2+temp3+temp4+temp5+temp6+temp7+temp8)

end function ME_rXr_over_rij


function ME_over_rij(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
real(dprec)   ME_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   inv_tAkl(nn,nn),det_tAkl
integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

real(dprec)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
real(dprec)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
real(dprec)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
integer       p,q,n
real(dprec)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl    
real(dprec)    temp1,temp2,temp3,temp11,temp22,temp33,mu,u

n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do q=1,n
    temp1=temp1+inv_tAkl(q,p)*tvl(q)
    temp2=temp2+inv_tAkl(q,p)*tbl(q)
    temp3=temp3+inv_tAkl(q,p)*tbk(q)
  enddo
  inv_tAkltvl(p)=temp1
  inv_tAkltbl(p)=temp2
  inv_tAkltbk(p)=temp3
enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do q=1,n
    temp1=temp1+tvk(q)*inv_tAkl(q,p)
    temp2=temp2+tbk(q)*inv_tAkl(q,p)
    temp3=temp3+tvl(q)*inv_tAkl(q,p)
  enddo
  tvkinv_tAkl(p)=temp1
  tbkinv_tAkl(p)=temp2
  tvlinv_tAkl(p)=temp3
enddo


!Compute tau3=tvkinv_tAkl*tvl
tau3=ZERO
tau33=ZERO
tau333=ZERO
tau334=ZERO
do p=1,n
  tau3=tau3+tvkinv_tAkl(p)*tvl(p)
  tau33=tau33+tbkinv_tAkl(p)*tbl(p)
  tau333=tau333+tvkinv_tAkl(p)*tbl(p)
  tau334=tau334+tbkinv_tAkl(p)*tvl(p)
enddo
m1=tau3*tau33
m3=tau333*tau334
m=m1+m3  !look formula 40 in document

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
Ajtbl=ZERO
temp1=ZERO
temp11=ZERO
do p=1,n
  Ajtvl=Ajtvl+Aj(p)*tvl(p)
  Ajtbl=Ajtbl+Aj(p)*tbl(p)  
  temp1=temp1+tvk(p)*Aj(p)
  temp11=temp11+tbk(p)*Aj(p)
enddo
t_JV1=temp11*Ajtbl
t_JV2=temp1*Ajtvl
t_JV5=temp11*Ajtvl
t_JV6=temp1*Ajtbl

mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
u=t_JV1*t_JV2+t_JV5*t_JV6

temp1=Glob_Piraised3n2/(TWO*SQRTPI*det_tAkl**(THREEHALF))

ME_over_rij=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

end function ME_over_rij


function ME_d_X_over_rij_d(p,q,X,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,twk,twl, &
  tvkinv_tAkl, twkinv_tAkl, inv_tAkltvl, inv_tAkltwl)
    
  real(dprec)   ME_d_X_over_rij_d
  integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
  !Arguments:
  real(dprec)   X(nn,nn),tAl(nn,nn),tAk(nn,nn),inv_tAkl(nn,nn),det_tAkl, &
  tvkinv_tAkl(nn), twkinv_tAkl(nn), inv_tAkltvl(nn), inv_tAkltwl(nn)
  integer       p,q,tvk(nn),tvl(nn),twk(nn),twl(nn)

  
  !Local variables:
  integer       c,s,n,k,i,j
  real(dprec)   tvkXtAl(nn),tbkXtAl(nn)
  real(dprec)   tAkX(nn,nn),tAkXtAl(nn,nn),tAkXtvl(nn),tAkXtbl(nn),XtAl(nn,nn)


  real(dprec) :: commonFactor, gamma, temp, temp1, temp2, temp3

  !Vars for Q-part
  real(dprec) :: trXs, Xij, V, Vij, VX, VijX, VXij, VijXij, &
  W, Wij, WX, WijX, WXij, WijXij, &
  jijAvk, jijAvl, jijAwk, jijAwl, &
  jijAXsAvk, jijAXsAvl, jijAXsAwk, jijAXsAwl 
  real(dprec) :: tV, tVij, tVX, tVijX, tVXij, tVijXij, &
  tW, tWij, tWX, tWijX, tWXij, tWijXij
  real(dprec) :: I11, I12, I13, I1, &
  I21, I22, I23, I2, &
  I31, I32, I33, I3, &
  I41, I42, I43, I4
  real(dprec) :: XAl(nn, nn), AXsA(nn, nn), Xs(nn, nn), XsA(nn, nn)
  real(dprec) :: AXsA_Vl(nn), AXsA_Wl(nn), Vk_AXsA(nn), Wk_AXsA(nn)
  real(dprec) :: Qans
  !Vars for RVk part
  real(dprec) :: AlA(nn,nn), XAlA(nn,nn), XAlA_Vl(nn), XAlA_Wl(nn), &
  Vk_XAlA(nn), Wk_XAlA(nn), VkXAlAjij, VkXAlAVl, VkXAlAWl
  real(dprec) :: RVk, RVk1, RVk2, RVk3
  !Vars for RWk part
  real(dprec) :: WkXAlAVl, WkXAlAWl, WkXAlAjij, RWk, RWk1, RWk2, RWk3
  !Vars for RVl part
  real(dprec) :: AkX(nn,nn), AAkX(nn,nn),  AAkX_Vl(nn), AAkX_Wl(nn), &
  VkAAkXVl, WkAAkXWl, VkAAkXWl, WkAAkXVl, &
  jijAAkXVl, jijAAkXWl
  real(dprec) :: RVl, RVl1, RVl2, RVl3
  !Vars for RWl part
  real(dprec) :: RWl, RWl1, RWl2, RWl3
  !Vars for D2 terms
  real(dprec) :: X_Vl(nn), X_Wl(nn), &
  VkXVl, VkXWl, WkXWl, WkXVl, DTwo1, DTwo2, DTwo3, DTwo4, DTwo

  n = Glob_n

  !!! Q-part  !!!
  !Build Xs matrix
  XAl = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n 
        temp = temp + X(i,k)*tAl(k,j)
      enddo
      XAl(i,j) = temp
    enddo
  enddo
  Xs = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n 
        temp = temp + tAk(i,k) * XAl(k,j)
      enddo
      Xs(i,j) = temp
    enddo
  enddo

  !Symmetrize MS
  do i = 1,n
    do j = i+1,n
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
    enddo
  enddo

  !Find trA
  trXs = ZERO
  do i=1,n
    do j=1,n
      trXs = trXs + inv_tAkl(i,j)*Xs(j,i) 
    enddo
  enddo

  XsA = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + Xs(i,k)*inv_tAkl(k,j)
      enddo
      XsA(i,j) = temp
    enddo
  enddo

  AXsA = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + inv_tAkl(i,k)*XsA(k,j)
      enddo
      AXsA(i,j) = temp
    enddo
  enddo

  V=ZERO
  W=ZERO
  tV=ZERO
  tW=ZERO
  do i=1,n
    V=V+tvkinv_tAkl(i)*tvl(i)
    W=W+twkinv_tAkl(i)*twl(i)
    tV=tV+tvkinv_tAkl(i)*twl(i)
    tW=tW+twkinv_tAkl(i)*tvl(i)
  enddo

  
  do i=1,n
    temp = ZERO
    temp1 = ZERO
    temp2 = ZERO
    temp3 = ZERO
    do j=1,n 
      temp = temp + AXsA(i,j)*tvl(j)
      temp1 = temp1 + AXsA(i,j)*twl(j)
      temp2 = temp2 + tvk(j)*AXsA(j,i)
      temp3 = temp3 + twk(j)*AXsA(j,i)
    enddo
    AXsA_Vl(i) = temp
    AXsA_Wl(i) = temp1
    Vk_AXsA(i) = temp2
    Wk_AXsA(i) = temp3
  enddo

  VX = ZERO
  WX = ZERO
  tVX = ZERO
  tWX = ZERO
  do i=1,n 
    VX = VX + tvk(i)*AXsA_vl(i)
    WX = WX + twk(i)*AXsA_wl(i)
    tVX = tVX + tvk(i)*AXsA_wl(i)
    tWX = tWX + twk(i)*AXsA_vl(i)
  enddo
  !!! END Q part !!!!

  !!!!! RVk & RWk part of M-matelem !!!!
  AlA = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + tAl(i,k)*inv_tAkl(k,j)
      enddo
      AlA(i,j) = temp
    enddo
  enddo

  XAlA = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + X(i,k)*AlA(k,j)
      enddo
      XAlA(i,j) = temp
    enddo
  enddo

  do i=1,n
    temp = ZERO
    temp1 = ZERO
    temp2 = ZERO
    temp3 = ZERO
    do j=1,n 
      temp = temp + XAlA(i,j)*tvl(j)
      temp1 = temp1 + XAlA(i,j)*twl(j)
      temp2 = temp2 + tvk(j)*XAlA(j,i)
      temp3 = temp3 + twk(j)*XAlA(j,i)
    enddo
    XAlA_Vl(i) = temp
    XAlA_Wl(i) = temp1
    Vk_XAlA(i) = temp2
    Wk_XAlA(i) = temp3
  enddo

  VkXAlAVl = ZERO
  VkXAlAWl = ZERO
  WkXAlAWl = ZERO
  WkXAlAVl = ZERO
  do i=1,n 
    VkXAlAVl = VkXAlAVl + tvk(i)*XAlA_Vl(i)
    VkXAlAWl = VkXAlAWl + tvk(i)*XAlA_Wl(i)
    WkXAlAWl = WkXAlAWl + twk(i)*XAlA_Wl(i)
    WkXAlAVl = WkXAlAVl + twk(i)*XAlA_Vl(i)
  enddo
  !!!!! END RVk & RWk part of M-matelem !!!!

  !!!!! RVl part of M-matelem !!!!
  AkX = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + tAk(i,k)*X(k,j)
      enddo
      AkX(i,j) = temp
    enddo
  enddo

  AAkX = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + inv_tAkl(i,k)*AkX(k,j)
      enddo
      AAkX(i,j) = temp
    enddo
  enddo

  do i=1,n
    temp = ZERO
    temp1 = ZERO
    do j=1,n 
      temp = temp + AAkX(i,j)*tvl(j)
      temp1 = temp1 + AAkX(i,j)*twl(j)
    enddo
    AAkX_Vl(i) = temp
    AAkX_Wl(i) = temp1
  enddo

  VkAAkXVl = ZERO
  VkAAkXWl = ZERO
  WkAAkXWl = ZERO
  WkAAkXVl = ZERO
  do i=1,n 
    VkAAkXVl =  VkAAkXVl + tvk(i)*AAkX_Vl(i)
    VkAAkXWl = VkAAkXWl + tvk(i)*AAkX_Wl(i)
    WkAAkXWl = WkAAkXWl + twk(i)*AAkX_Wl(i)
    WkAAkXVl = WkAAkXVl + twk(i)*AAkX_Vl(i)
  enddo
  !!!!! END RVl part of M-matelem !!!!

  !!!!! D2 terms !!!! 
  do i=1,n
    temp = ZERO
    temp1 = ZERO
    do j=1,n 
      temp = temp + X(i,j)*tvl(j)
      temp1 = temp1 + X(i,j)*twl(j)
    enddo
    X_Vl(i) = temp
    X_Wl(i) = temp1
  enddo

  VkXVl = ZERO
  VkXWl = ZERO
  WkXWl = ZERO 
  WkXVl = ZERO
  do i=1,n 
    VkXVl = VkXVl + tvk(i)*X_Vl(i) 
    VkXWl = VkXWl + tvk(i)*X_Wl(i) 
    WkXWl = WkXWl + tWk(i)*X_Wl(i)
    WkXVl = WkXVl + tWk(i)*X_Vl(i)
  enddo

  !!!!! END D2 terms !!!!


  if (p==q) then
    !Common part
    gamma = inv_tAkl(p,p)
    gamma = ONE/sqrt(gamma)
    !Q part
    Xij = AXsA(p, p)
    jijAvk = tvkinv_tAkl(p)
    jijAvl = inv_tAkltvl(p)
    jijAwk = twkinv_tAkl(p)
    jijAwl = inv_tAkltwl(p)
    jijAXsAvl = AXsA_Vl(p)
    jijAXsAvk = Vk_AXsA(p)
    jijAXsAwl = AXsA_Wl(p)
    jijAXsAwk = Wk_AXsA(p)
    !RVk part
    VkXAlAjij = Vk_XAlA(p)
    !RWk part
    WkXAlAjij = Wk_XAlA(p) 
    ! RVl part 
    jijAAkXvl = AAkX_Vl(p)
    jijAAkXWl = AAkX_Wl(p)
  else
    !Common part
    gamma = inv_tAkl(p,p) + inv_tAkl(q,q) - inv_tAkl(p,q) - inv_tAkl(q,p)
    gamma = ONE/sqrt(gamma)
    !Q part
    Xij = AXsA(p,p) + AXsA(q,q) - AXsA(p,q) - AXsA(q,p)   
    jijAvk = tvkinv_tAkl(p) - tvkinv_tAkl(q)
    jijAvl = inv_tAkltvl(p) - inv_tAkltvl(q)
    jijAwk = twkinv_tAkl(p) - twkinv_tAkl(q)
    jijAwl = inv_tAkltwl(p) - inv_tAkltwl(q)
    jijAXsAvl = AXsA_Vl(p) - AXsA_Vl(q)
    jijAXsAvk = Vk_AXsA(p) - Vk_AXsA(q)
    jijAXsAwl = AXsA_Wl(p) - AXsA_Wl(q)
    jijAXsAwk = Wk_AXsA(p) - Wk_AXsA(q)
    !RVk part
    VkXAlAjij = Vk_XAlA(p) - Vk_XAlA(q) 
    !RWk part
    WkXAlAjij = Wk_XAlA(p) - Wk_XAlA(q)
    ! RVl part 
    jijAAkXvl = AAkX_Vl(p) - AAkX_Vl(q)
    jijAAkXWl = AAkX_Wl(p) - AAkX_Wl(q)
  endif

  commonFactor = Glob_Piraised3n2/(sqrt(PI)*det_tAkl*sqrt(det_tAkl))

  !Q-part
  Vij = jijAvk * jijAvl
  VijX = jijAvk * jijAXsAvl 
  VXij = jijAvl * jijAXsAvk
  VijXij = jijAvk * Xij * jijAvl

  Wij = jijAwk * jijAwl
  WijX = jijAwk * jijAXsAwl 
  WXij = jijAwl * jijAXsAwk
  WijXij = jijAwk * Xij * jijAwl

  tVij = jijAvk * jijAwl
  tVijX = jijAvk * jijAXsAwl 
  tVXij = jijAwl * jijAXsAvk
  tVijXij = jijAvk * Xij * jijAwl

  tWij = jijAwk * jijAvl
  tWijX = jijAwk * jijAXsAvl 
  tWXij = jijAvl * jijAXsAwk
  tWijXij = jijAwk * Xij * jijAvl

  I11 = trXs*V*W + trXs*tV*tW
  I12 = VX*W + tVX*tW
  I13 = V*WX + tV*tWX
  I1 = (THREEHALF*I11 + I12 + I13)*gamma

  I21 = (Xij*V*W + trXs*Vij*W + trXs*V*Wij) + (Xij*tV*tW + trXs*tVij*tW + trXs*tV*tWij)
  I22 = (VijX*W + VXij*W + VX*Wij) + (tVijX*tW + tVXij*tW + tVX*tWij) 
  I23 = (V*WijX + V*WXij + Vij*WX) + (tV*tWijX + tV*tWXij + tVij*tWX)
  I2 = -(THREEHALF*I21 + I22 + I23)*gamma**3/THREE

  I31 = (Xij*Vij*W + Xij*V*Wij + trXs*Vij*Wij) + (Xij*tVij*tW + Xij*tV*tWij + trXs*tVij*tWij)
  I32 = (VijXij*W + VijX*Wij + VXij*Wij) + (tVijXij*tW + tVijX*tWij + tVXij*tWij) 
  I33 = (V*WijXij + Vij*WijX + Vij*WXij) + (tV*tWijXij + tVij*tWijX + tVij*tWXij)
  I3 = (THREEHALF*I31 + I32 + I33)*gamma**5/FIVE

  I41 = Xij*Vij*Wij + Xij*tVij*tWij 
  I42 = VijXij*Wij + tVijXij*tWij
  I43 = Vij*WijXij +  tVij*tWijXij
  I4 = -(THREEHALF*I41 + I42 + I43)*gamma**7/SEVEN

  Qans = (I1 + I2 + I3 + I4)*TWO*commonFactor

  !RVk-part of M-matelem
  RVk1 = gamma*(VkXAlAVl*W +VkXAlAWl*tW)
  RVk2 = -gamma**3/THREE*(VkXAlAjij*jijAvl*W + VkXAlAVl*Wij + &
  VkXAlAjij*jijAWl*tW + VkXAlAWl*tWij)
  RVk3 = gamma**5/FIVE*(VkXAlAjij*jijAvl*Wij + VkXAlAjij*jijAWl*tWij)
  RVk = -(RVk1 + RVk2 + RVk3)*commonFactor

  !RWk-part of M-matelem
  RWk1 = gamma*(V*WkXAlAWl + tV*WkXAlAVl)
  RWk2 = -gamma**3/THREE*(Vij*WkXAlAWl + V*WkXAlAjij*jijAWl + &
  tVij*WkXAlAVl + tV*WkXAlAjij*jijAVl)
  RWk3 = gamma**5/FIVE*(Vij*WkXAlAjij*jijAWl + tVij*WkXAlAjij*jijAVl)
  RWk = -(RWk1 + RWk2 + RWk3)*commonFactor

  !RVl-part of M-matelem
  RVl1 = gamma*(VkAAkXVl*W + tV*WkAAkXVl)
  RVl2 = -gamma**3/THREE*(jijAVk*jijAAkXVl*W + VkAAkXVl*Wij + &
  jijAVk*jijAWl*WkAAkXVl +  tV*jijAWk*jijAAkXVl)
  RVl3 = gamma**5/FIVE*(jijAVk*jijAAkXVl*Wij + &
  tVij*jijAWk*jijAAkXVl)
  
  RVl = -(RVl1 + RVl2 + RVl3)*Glob_Piraised3n2/(sqrt(PI)*det_tAkl*sqrt(det_tAkl))

  !RWl-part of M-matelem 
  RWl1 = gamma*(V*WkAAkXWl + VkAAkXWl*tW)
  RWl2 = -gamma**3/THREE*(Vij*WkAAkXWl + V*jijAWk*jijAAkXWl + &
  jijAVk*jijAAkXWl*tW + VkAAkXWl*tWij)
  RWl3 = gamma**5/FIVE*(Vij*jijAWk*jijAAkXWl + jijAVk*jijAAkXWl*tWij)
  RWl = -(RWl1 + RWl2 + RWl3)*Glob_Piraised3n2/(sqrt(PI)*det_tAkl*sqrt(det_tAkl))
  !END RWl-part of M-matelem !

  
  !D2 part
  DTwo1 = VkXWl*(gamma*tW - ONE/THREE*(gamma**3) * jijAwk * jijAVl)
  DTwo2 = WkXVl*(gamma*tV - ONE/THREE*(gamma**3) * jijAVk * jijAWl)
  DTwo3 = VkXVl*(gamma*W - ONE/THREE*(gamma**3) * jijAWk * jijAWl)
  DTwo4 = WkXWl*(gamma*V - ONE/THREE*(gamma**3) * jijAVk * jijAVl)
  DTwo = (DTwo1 + DTwo2 + DTwo3 + DTwo4)*commonFactor

  ME_d_X_over_rij_d = Qans + RVk + RVl + RWk + RWl + DTwo

  end function ME_d_X_over_rij_d



function ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
real(dprec)   ME_rXr_rYr_over_rij
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
real(dprec)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),det_tAkl
integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

real(dprec)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
real(dprec)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
real(dprec)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
integer       p,q,n,s,k
real(dprec)   Aj(nn),AjX(nn),AjY(nn)
real(dprec)   t_J,t_X,t_XJ,t_Y,t_YJ,t_XYJ
real(dprec)   t_JV1,t_XV1,t_JXV1,t_XJV1,t_YV1,t_JYV1,t_YJV1,t_XYV1,t_YXV1,t_XYJV1,t_YXJV1,t_XJYV1,t_YJXV1,t_JXYV1,t_JYXV1
real(dprec)   t_JV2,t_XV2,t_JXV2,t_XJV2,t_YV2,t_JYV2,t_YJV2,t_XYV2,t_YXV2,t_XYJV2,t_YXJV2,t_XJYV2,t_YJXV2,t_JXYV2,t_JYXV2
real(dprec)   t_JV5,t_XV5,t_JXV5,t_XJV5,t_YV5,t_JYV5,t_YJV5,t_XYV5,t_YXV5,t_XYJV5,t_YXJV5,t_XJYV5,t_YJXV5,t_JXYV5,t_JYXV5
real(dprec)   t_JV6,t_XV6,t_JXV6,t_XJV6,t_YV6,t_JYV6,t_YJV6,t_XYV6,t_YXV6,t_XYJV6,t_YXJV6,t_XJYV6,t_YJXV6,t_JXYV6,t_JYXV6
real(dprec)    Ajtvl,Ajtbl,tvkAj,tbkAj,prod
real(dprec)    temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp44,temp77,temp88,temp99
real(dprec)    temp11,temp22,temp33,mu,mX,mXJ,u,Xmu,muXJ,mY,mYJ,Ymu,muYJ,t_XYJV11,t_XYJV22,t_XYJV55,t_XYJV66,YX,mXYJ
real(dprec)    temp2Y,temp22Y,t_XY,temp55,temp66,AXAj(nn),AYAj(nn),muYX,muXYJ,XYJ,YXJ,t_YXJ,XJYJ,mYX
real(dprec)    tvkinv_tAklX(nn),tbkinv_tAklX(nn),tvkinv_tAklY(nn),tbkinv_tAklY(nn),AY(nn,nn),AX(nn,nn)
real(dprec)    AXinv_tAkltvl(nn),AXinv_tAkltbl(nn),AYinv_tAkltvl(nn),AYinv_tAkltbl(nn)
real(dprec)    term1,term2,term3,term4,term5,term6,term7,term8,term9,term10,term11,term12,term13,temp1Y,temp11Y


n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do q=1,n
    temp1=temp1+inv_tAkl(q,p)*tvl(q)
    temp2=temp2+inv_tAkl(q,p)*tbl(q)
    temp3=temp3+inv_tAkl(q,p)*tbk(q)
  enddo
  inv_tAkltvl(p)=temp1
  inv_tAkltbl(p)=temp2
  inv_tAkltbk(p)=temp3
enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do q=1,n
    temp1=temp1+tvk(q)*inv_tAkl(q,p)
    temp2=temp2+tbk(q)*inv_tAkl(q,p)
    temp3=temp3+tvl(q)*inv_tAkl(q,p)
  enddo
  tvkinv_tAkl(p)=temp1
  tbkinv_tAkl(p)=temp2
  tvlinv_tAkl(p)=temp3
enddo


!Compute tau3=tvkinv_tAkl*tvl
tau3=ZERO
tau33=ZERO
tau333=ZERO
tau334=ZERO
do p=1,n
  tau3=tau3+tvkinv_tAkl(p)*tvl(p)
  tau33=tau33+tbkinv_tAkl(p)*tbl(p)
  tau333=tau333+tvkinv_tAkl(p)*tbl(p)
  tau334=tau334+tbkinv_tAkl(p)*tvl(p)
enddo
m1=tau3*tau33
m3=tau333*tau334
m=m1+m3  !look formula 40 in document

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
  temp2=ZERO
  do q=1,n
    temp1=temp1+Aj(q)*X(q,p)
    temp2=temp2+Aj(q)*Y(q,p)
  enddo
  AjX(p)=temp1
  AjY(p)=temp2
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
Ajtbl=ZERO
t_XJ=ZERO
t_YJ=ZERO
temp1=ZERO
temp2=ZERO
temp2Y=ZERO
temp11=ZERO
temp22=ZERO
temp22Y=ZERO
do p=1,n
  Ajtvl=Ajtvl+Aj(p)*tvl(p)
  Ajtbl=Ajtbl+Aj(p)*tbl(p)  
  t_XJ=t_XJ+AjX(p)*Aj(p)
  t_YJ=t_YJ+AjY(p)*Aj(p)
  temp1=temp1+tvk(p)*Aj(p)
  temp2=temp2+AjX(p)*inv_tAkltvl(p)
  temp2Y=temp2Y+AjY(p)*inv_tAkltvl(p)
  temp11=temp11+tbk(p)*Aj(p)
  temp22=temp22+AjX(p)*inv_tAkltbl(p)
  temp22Y=temp22Y+AjY(p)*inv_tAkltbl(p)
enddo
t_JV1=temp11*Ajtbl
t_JXV1=temp11*temp22
t_JYV1=temp11*temp22Y
t_JV2=temp1*Ajtvl
t_JXV2=temp1*temp2
t_JYV2=temp1*temp2Y
t_JV5=temp11*Ajtvl
t_JXV5=temp11*temp2
t_JYV5=temp11*temp2Y
t_JV6=temp1*Ajtbl
t_JXV6=temp1*temp22
t_JYV6=temp1*temp22Y
!Compute t_X=tr[inv_tAkl*X]
!        t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*inv_tAkltvl
!        t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*Aj*Ajtvl
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
enddo

t_X=ZERO
t_Y=ZERO
t_XV1=ZERO
t_XV2=ZERO
t_XV5=ZERO
t_XV6=ZERO
t_YV1=ZERO
t_YV2=ZERO
t_YV5=ZERO
t_YV6=ZERO
temp2=ZERO
temp22=ZERO
temp2Y=ZERO
temp22Y=ZERO
do p=1,n
  temp1=ZERO
  temp11=ZERO
  temp1Y=ZERO
  temp11Y=ZERO
  do q=1,n
    t_X=t_X+inv_tAkl(q,p)*X(q,p)
    t_Y=t_Y+inv_tAkl(q,p)*Y(q,p)    
    temp1=temp1+tvkinv_tAkl(q)*X(q,p)
    temp11=temp11+tbkinv_tAkl(q)*X(q,p)
    temp1Y=temp1Y+tvkinv_tAkl(q)*Y(q,p)
    temp11Y=temp11Y+tbkinv_tAkl(q)*Y(q,p)
  enddo
  t_XV1=t_XV1+temp11*inv_tAkltbl(p)
  t_XV2=t_XV2+temp1*inv_tAkltvl(p)
  t_XV5=t_XV5+temp11*inv_tAkltvl(p)
  t_XV6=t_XV6+temp1*inv_tAkltbl(p)
  t_YV1=t_YV1+temp11Y*inv_tAkltbl(p)
  t_YV2=t_YV2+temp1Y*inv_tAkltvl(p)
  t_YV5=t_YV5+temp11Y*inv_tAkltvl(p)
  t_YV6=t_YV6+temp1Y*inv_tAkltbl(p)
  temp2=temp2+temp1*Aj(p)
  temp22=temp22+temp11*Aj(p)
  temp2Y=temp2Y+temp1Y*Aj(p)
  temp22Y=temp22Y+temp11Y*Aj(p)
enddo
t_XJV1=temp22*Ajtbl
t_XJV2=temp2*Ajtvl
t_XJV5=temp22*Ajtvl
t_XJV6=temp2*Ajtbl
t_YJV1=temp22Y*Ajtbl
t_YJV2=temp2Y*Ajtvl
t_YJV5=temp22Y*Ajtvl
t_YJV6=temp2Y*Ajtbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
t_XY=ZERO
do p=1,n
  temp1=ZERO
  temp2=ZERO
  temp11=ZERO
  temp22=ZERO
  temp3=ZERO
  temp4=ZERO
  temp5=ZERO
  temp6=ZERO 
  temp55=ZERO
  temp66=ZERO 
  do q=1,n
    t_XY=t_XY+AX(p,q)*AY(q,p)
    temp1=temp1+AX(p,q)*inv_tAkltvl(q)
    temp2=temp2+AY(p,q)*inv_tAkltvl(q)
    temp11=temp11+AX(p,q)*inv_tAkltbl(q)
    temp22=temp22+AY(p,q)*inv_tAkltbl(q)
    temp3=temp3+AX(p,q)*Aj(q)
    temp4=temp4+AY(p,q)*Aj(q)
    temp5=temp5+tvkinv_tAkl(q)*X(q,p)
    temp6=temp6+tvkinv_tAkl(q)*Y(q,p)  
    temp55=temp55+tbkinv_tAkl(q)*X(q,p)
    temp66=temp66+tbkinv_tAkl(q)*Y(q,p)  
  enddo
  AXinv_tAkltvl(p)=temp1
  AYinv_tAkltvl(p)=temp2
  AXinv_tAkltbl(p)=temp11
  AYinv_tAkltbl(p)=temp22
  AXAj(p)=temp3
  AYAj(p)=temp4
  tvkinv_tAklX(p)=temp5
  tvkinv_tAklY(p)=temp6 
  tbkinv_tAklX(p)=temp55
  tbkinv_tAklY(p)=temp66 
enddo

tvkAj=ZERO
tbkAj=ZERO
temp1=ZERO
temp11=ZERO
temp2=ZERO
temp22=ZERO
temp3=ZERO
temp33=ZERO
t_XYJ=ZERO
t_YXJ=ZERO
t_XYV1=ZERO
t_YXV1=ZERO
t_XYV2=ZERO
t_YXV2=ZERO
t_XYV5=ZERO
t_YXV5=ZERO
t_XYV6=ZERO
t_YXV6=ZERO
temp4=ZERO
temp5=ZERO
temp44=ZERO
temp55=ZERO
temp6=ZERO
temp7=ZERO
temp8=ZERO
temp9=ZERO
temp66=ZERO
temp77=ZERO
temp88=ZERO
temp99=ZERO
do p=1,n
  tvkAj=tvkAj+tvk(p)*Aj(p)
  tbkAj=tbkAj+tbk(p)*Aj(p)
!  temp1=temp1+tvk(p)*Aj(p)
!  temp11=temp11+tbk(p)*Aj(p)
  temp2=temp2+AjX(p)*inv_tAkltvl(p)
  temp22=temp22+AjX(p)*inv_tAkltbl(p)
  temp3=temp3+AjY(p)*inv_tAkltvl(p)
  temp33=temp33+AjY(p)*inv_tAkltbl(p)
  t_XYJ=t_XYJ+AjX(p)*AYAj(p)
  t_YXJ=t_YXJ+AjY(p)*AXAj(p)
  t_XYV1=t_XYV1+tbkinv_tAklX(p)*AYinv_tAkltbl(p)
  t_XYV2=t_XYV2+tvkinv_tAklX(p)*AYinv_tAkltvl(p)
  t_XYV5=t_XYV5+tbkinv_tAklX(p)*AYinv_tAkltvl(p)
  t_XYV6=t_XYV6+tvkinv_tAklX(p)*AYinv_tAkltbl(p) 
  t_YXV1=t_YXV1+tbkinv_tAklY(p)*AXinv_tAkltbl(p)
  t_YXV2=t_YXV2+tvkinv_tAklY(p)*AXinv_tAkltvl(p)
  t_YXV5=t_YXV5+tbkinv_tAklY(p)*AXinv_tAkltvl(p)
  t_YXV6=t_YXV6+tvkinv_tAklY(p)*AXinv_tAkltbl(p)  
  temp4=temp4+tvkinv_tAklX(p)*Aj(p)
  temp5=temp5+tvkinv_tAklY(p)*Aj(p)
  temp44=temp44+tbkinv_tAklX(p)*Aj(p)
  temp55=temp55+tbkinv_tAklY(p)*Aj(p)  
  temp6=temp6+tvkinv_tAklX(p)*AYAj(p)
  temp7=temp7+tvkinv_tAklY(p)*AXAj(p) 
  temp8=temp8+AjX(p)*AYinv_tAkltvl(p) 
  temp9=temp9+AjY(p)*AXinv_tAkltvl(p)
  temp66=temp66+tbkinv_tAklX(p)*AYAj(p)
  temp77=temp77+tbkinv_tAklY(p)*AXAj(p) 
  temp88=temp88+AjX(p)*AYinv_tAkltbl(p) 
  temp99=temp99+AjY(p)*AXinv_tAkltbl(p) 
enddo
t_XYJV1=temp66*Ajtbl
t_XYJV2=temp6*Ajtvl
t_XYJV5=temp66*Ajtvl
t_XYJV6=temp6*Ajtbl
t_YXJV1=temp77*Ajtbl
t_YXJV2=temp7*Ajtvl
t_YXJV5=temp77*Ajtvl
t_YXJV6=temp7*Ajtbl
t_XJYV1=temp44*temp33
t_XJYV2=temp4*temp3
t_XJYV5=temp44*temp3
t_XJYV6=temp4*temp33
t_YJXV1=temp55*temp22
t_YJXV2=temp5*temp2
t_YJXV5=temp55*temp2
t_YJXV6=temp5*temp22
t_JXYV1=tbkAj*temp88
t_JXYV2=tvkAj*temp8
t_JXYV5=tbkAj*temp8
t_JXYV6=tvkAj*temp88
t_JYXV1=tbkAj*temp99
t_JYXV2=tvkAj*temp9
t_JYXV5=tbkAj*temp9
t_JYXV6=tvkAj*temp99
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
t_XYJV11=(t_XYJV1+t_YXJV1+t_XJYV1+t_YJXV1+t_JXYV1+t_JYXV1)
t_XYJV22=(t_XYJV2+t_YXJV2+t_XJYV2+t_YJXV2+t_JXYV2+t_JYXV2)
t_XYJV55=(t_XYJV5+t_YXJV5+t_XJYV5+t_YJXV5+t_JXYV5+t_JYXV5)
t_XYJV66=(t_XYJV6+t_YXJV6+t_XJYV6+t_YJXV6+t_JXYV6+t_JYXV6)
!!!!!!!!!!!!!!!!!!!!!!
mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
mX=tau3*t_XV1+tau33*t_XV2+tau333*t_XV5+tau334*t_XV6
mY=tau3*t_YV1+tau33*t_YV2+tau333*t_YV5+tau334*t_YV6
mYX=tau3*t_YXV1+tau33*t_YXV2+tau333*t_YXV5+tau334*t_YXV6
u=t_JV1*t_JV2+t_JV5*t_JV6
mXJ=tau3*(t_XJV1+t_JXV1)+tau33*(t_XJV2+t_JXV2)+tau333*(t_XJV5+t_JXV5)+tau334*(t_XJV6+t_JXV6)
mYJ=tau3*(t_YJV1+t_JYV1)+tau33*(t_YJV2+t_JYV2)+tau333*(t_YJV5+t_JYV5)+tau334*(t_YJV6+t_JYV6)
muXJ=t_JV2*(t_XJV1+t_JXV1)+t_JV1*(t_XJV2+t_JXV2)+t_JV6*(t_XJV5+t_JXV5)+t_JV5*(t_XJV6+t_JXV6)
muYJ=t_JV2*(t_YJV1+t_JYV1)+t_JV1*(t_YJV2+t_JYV2)+t_JV6*(t_YJV5+t_JYV5)+t_JV5*(t_YJV6+t_JYV6)
Xmu=t_XV2*t_JV1+t_XV1*t_JV2+t_XV6*t_JV5+t_XV5*t_JV6
Ymu=t_YV2*t_JV1+t_YV1*t_JV2+t_YV6*t_JV5+t_YV5*t_JV6

mXYJ=tau3*t_XYJV11+tau33*t_XYJV22+tau333*t_XYJV55+tau334*t_XYJV66
YX=t_YV2*t_XV1+t_YV1*t_XV2+t_YV6*t_XV5+t_YV5*t_XV6
muYX=t_JV1*t_YXV2+t_JV2*t_YXV1+t_JV6*t_YXV5+t_JV5*t_YXV6
muXYJ=t_JV1*t_XYJV22+t_JV2*t_XYJV11+t_JV6*t_XYJV55+t_JV5*t_XYJV66
XYJ=t_XV1*t_YJV2+t_XV2*t_YJV1+t_XV5*t_YJV6+t_XV6*t_YJV5
YXJ=t_YV1*t_XJV2+t_YV2*t_XJV1+t_YV5*t_XJV6+t_YV6*t_XJV5
XJYJ=(t_XJV1+t_JXV1)*(t_YJV2+t_JYV2)+(t_XJV2+t_JXV2)*(t_YJV1+t_JYV1)+(t_XJV5+t_JXV5)*(t_YJV6+t_JYV6)+(t_XJV6+t_JXV6)*(t_YJV5+t_JYV5)

temp1=Glob_Piraised3n2/(TWO*SQRTPI*det_tAkl**(THREEHALF))
temp2=+THREEHALF*t_X*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)
temp3=-m*t_XJ/(TWO*t_J*sqrt(t_J))
temp4=+mX/sqrt(t_J)
temp5=+mu*t_XJ/(TWO*t_J*t_J*sqrt(t_J))
temp6=-(mXJ+Xmu)/(THREE*t_J*sqrt(t_J))
temp7=-u*t_XJ/(TWO*t_J*t_J*t_J*sqrt(t_J))
temp8=+ONEFIFTH*muXJ/(t_J*t_J*sqrt(t_J))

term1=THREEHALF*temp1*t_Y*(temp2+temp3+temp4+temp5+temp6+temp7+temp8)
term2=-THREEHALF*t_XY*(m-ONETHIRD*mu/t_J+ONEFIFTH*u/(t_J*t_J))/sqrt(t_J)
term3=+THREEHALF*t_X*(t_YJ*m/(TWO*t_J)-mY-t_YJ*mu/(TWO*t_J*t_J) &
     +ONETHIRD*(Ymu+mYJ)/t_J+t_YJ*u/(TWO*t_J**3)-ONEFIFTH*muYJ/(t_J*t_J))/sqrt(t_J)
term4=-3*m*t_XJ*t_YJ/(4*t_J*t_J*sqrt(t_J))
term5=+(mY*t_XJ+m*(t_YXJ+t_XYJ)+mX*t_YJ)/(2*t_J*sqrt(t_J))
term6=-(YX+mYX)/sqrt(t_J)
term7=+5*mu*t_XJ*t_YJ/(4*sqrt(t_J)*t_J**3)
term8=-(mu*(t_YXJ+t_XYJ)+t_XJ*(Ymu+mYJ)+t_YJ*(Xmu+mXJ))/(2*sqrt(t_J)*t_J**2)
term9=+ONETHIRD*(muYX+XYJ+YXJ+mXYJ)/(sqrt(t_J)*t_J)
term10=-7*t_YJ*t_XJ*u/(4*sqrt(t_J)*t_J**4)
term11=+(u*(t_YXJ+t_XYJ)+t_XJ*muYJ+t_YJ*muXJ)/(2*sqrt(t_J)*t_J**3)
term12=-ONEFIFTH*(muXYJ+XJYJ)/(sqrt(t_J)*t_J**2)

term13=temp1*(term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12)

ME_rXr_rYr_over_rij=term1-term13

end function ME_rXr_rYr_over_rij


function rPr_rQr(P,Q,tvk,tbk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,tvkinv_tAkl,tbkinv_tAkl,inv_tAkltvl,inv_tAkltbl)
real(dprec)   rPr_rQr
!arguments
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer       tvk(nn),tbk(nn)
real(dprec)   P(nn,nn),Q(nn,nn),inv_tAkl(nn,nn),tau3,tau33,tau333,tau334,det_tAkl
real(dprec)   inv_tAkltvl(nn),inv_tAkltbl(nn),tvkinv_tAkl(nn),tbkinv_tAkl(nn)

!decleration
integer       i,j,n,k
real(dprec)   temp1,temp2,temp3,temp4,P1,P2,P5,P6,Q1,Q2,Q5,Q6,trP,trQ,trQP
real(dprec)   PQ1,PQ2,PQ5,PQ6,QP1,QP2,QP5,QP6
real(dprec)   inv_tAklP(nn,nn),inv_tAklQ(nn,nn),inv_tAklQP(nn,nn)
real(dprec)   tvkinv_tAklP(nn),tbkinv_tAklP(nn),tvkinv_tAklQ(nn),tbkinv_tAklQ(nn)
real(dprec)   Pgamma,Qgamma,PQgamma,QPgamma
n=Glob_n

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*P(k,i)
      temp2=temp2+inv_tAkl(j,k)*Q(k,i)     
    enddo
    inv_tAklP(j,i)=temp1
    inv_tAklQ(j,i)=temp2
  enddo
enddo

do i=1,n
  do j=1,n
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAklQ(j,k)*inv_tAklP(k,i)    
    enddo
    inv_tAklQP(j,i)=temp1
  enddo
enddo

trQP=ZERO
trP=ZERO
trQ=ZERO
do i=1,n
  trQP=trQP+inv_tAklQP(i,i)
enddo

do i=1,n
  do j=1,n
    trP=trP+inv_tAkl(i,j)*P(j,i)
    trQ=trQ+inv_tAkl(i,j)*Q(j,i)
  enddo
enddo

Q1=ZERO
Q2=ZERO
Q5=ZERO
Q6=ZERO
P1=ZERO
P2=ZERO
P5=ZERO
P6=ZERO
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  do j=1,n
    temp1=temp1+tbkinv_tAkl(j)*P(j,i)
    temp2=temp2+tvkinv_tAkl(j)*P(j,i)   
    temp3=temp3+tbkinv_tAkl(j)*Q(j,i)
    temp4=temp4+tvkinv_tAkl(j)*Q(j,i)
  enddo
  P1=P1+temp1*inv_tAkltbl(i)
  P2=P2+temp2*inv_tAkltvl(i) 
  P5=P5+temp1*inv_tAkltvl(i)
  P6=P6+temp2*inv_tAkltbl(i)
  Q1=Q1+temp3*inv_tAkltbl(i)
  Q2=Q2+temp4*inv_tAkltvl(i) 
  Q5=Q5+temp3*inv_tAkltvl(i)
  Q6=Q6+temp4*inv_tAkltbl(i)
enddo

do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAklP(j,i)
    temp2=temp2+tbk(j)*inv_tAklP(j,i)
    temp3=temp3+tvk(j)*inv_tAklQ(j,i)
    temp4=temp4+tbk(j)*inv_tAklQ(j,i)
  enddo
  tvkinv_tAklP(i)=temp1
  tbkinv_tAklP(i)=temp2
  tvkinv_tAklQ(i)=temp3
  tbkinv_tAklQ(i)=temp4
enddo

PQ1=ZERO
PQ2=ZERO
PQ5=ZERO
PQ6=ZERO
QP1=ZERO
QP2=ZERO
QP5=ZERO
QP6=ZERO
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  do j=1,n
    temp1=temp1+tbkinv_tAklP(j)*inv_tAklQ(j,i)
    temp2=temp2+tvkinv_tAklP(j)*inv_tAklQ(j,i)   
    temp3=temp3+tbkinv_tAklQ(j)*inv_tAklP(j,i)
    temp4=temp4+tvkinv_tAklQ(j)*inv_tAklP(j,i)
  enddo
  PQ1=PQ1+temp1*inv_tAkltbl(i)
  PQ2=PQ2+temp2*inv_tAkltvl(i)
  PQ5=PQ5+temp1*inv_tAkltvl(i)
  PQ6=PQ6+temp2*inv_tAkltbl(i) 
  QP1=QP1+temp3*inv_tAkltbl(i)
  QP2=QP2+temp4*inv_tAkltvl(i)
  QP5=QP5+temp3*inv_tAkltvl(i)
  QP6=QP6+temp4*inv_tAkltbl(i)
enddo

temp1=Glob_Piraised3n2/(FOUR*det_tAkl**(THREEHALF))
temp2=THREEHALF*(tau3*tau33+tau333*tau334)*(trQP+THREEHALF*trQ*trP)
temp3=Q1*P2+Q2*P1+Q5*P6+Q6*P5
Pgamma=P1*tau3+P2*tau33+P5*tau333+P6*tau334
Qgamma=Q1*tau3+Q2*tau33+Q5*tau333+Q6*tau334
PQgamma=PQ1*tau3+PQ2*tau33+PQ5*tau333+PQ6*tau334
QPgamma=QP1*tau3+QP2*tau33+QP5*tau333+QP6*tau334

rPr_rQr= temp1*(THREEHALF*trQ*Pgamma+THREEHALF*trP*Qgamma+PQgamma+QPgamma+temp2+temp3)

end function rPr_rQr


function dXddYd(X,Y,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)
real(dprec)   dXddYd
!arguments
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
integer       tvk(nn),tbk(nn),tvl(nn),tbl(nn)
real(dprec)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),tau3,tau33,tau333,tau334,det_tAkl
real(dprec)   tAk(nn,nn),tAl(nn,nn)
real(dprec)   inv_tAkltvl(nn),inv_tAkltbl(nn),tvkinv_tAkl(nn),tbkinv_tAkl(nn)
real(dprec)   inv_tAkltAl(nn,nn),inv_tAkltAk(nn,nn)
!decleration
integer       i,j,n,k
real(dprec)   P(nn,nn),Q(nn,nn),trAkX,trAlY
real(dprec)   temp1,temp2,temp3,temp4,temp33,temp44,P1,P2,P5,P6,Q1,Q2,Q5,Q6,trP,trQ,trQP
real(dprec)   PQ1,PQ2,PQ5,PQ6,QP1,QP2,QP5,QP6
real(dprec)   inv_tAklP(nn,nn),inv_tAklQ(nn,nn),inv_tAklQP(nn,nn),tvkXtAk(nn),tbkXtAk(nn)
real(dprec)   tvkinv_tAklP(nn),tbkinv_tAklP(nn),tvkinv_tAklQ(nn),tbkinv_tAklQ(nn)
real(dprec)   XtAk(nn,nn),YtAl(nn,nn)
real(dprec)   Pgamma,Qgamma,PQgamma,QPgamma,big_term1,big_term2,big_term3,big_term4
real(dprec)   term1,term2,term3,term4,term5,term6,term7,term8,prod,gamma
real(dprec)   term9,term10,term11,term12,term13,term14,term15,term16
real(dprec)   term17,term18,term19,term20,term21,term22,term23,term24,term25
real(dprec)   tAkX(nn,nn),tAlY(nn,nn),Ytvl(nn),Ytbl(nn),Xtbk(nn),tvkX(nn),tbkX(nn),tvlY(nn)
real(dprec)   term3_1,term3_2,term4_1,term4_2,term8_1,term8_2,term9_1,term9_2,term11_1,term11_2,term16_1,term16_2
real(dprec)   term12_1,term12_2,term17_1,term17_2,term13_1,term14_2,term18_2,term19_1,temp5,temp6

n=Glob_n
prod=Glob_Piraised3n2/(FOUR*det_tAkl**(THREEHALF))
gamma=tau3*tau33+tau333*tau334


do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+inv_tAkl(i,j)*tvl(j)
    temp2=temp2+inv_tAkl(i,j)*tbl(j)
  enddo
  inv_tAkltvl(i)=temp1
  inv_tAkltbl(i)=temp2
enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAkl(j,i)
    temp2=temp2+tbk(j)*inv_tAkl(j,i)
  enddo
  tvkinv_tAkl(i)=temp1
  tbkinv_tAkl(i)=temp2
enddo
 !carefully look at the i,j indices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do j=1,n
    temp1=temp1+Y(i,j)*tvl(j)
    temp2=temp2+Y(i,j)*tbl(j)
    temp3=temp3+X(i,j)*tbk(j)
  enddo
  Ytvl(i)=temp1
  Ytbl(i)=temp2
  Xtbk(i)=temp3
enddo

do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*X(j,i)
    temp2=temp2+tbk(j)*X(j,i)
    temp3=temp3+tvl(j)*Y(j,i)
  enddo
  tvkX(i)=temp1
  tbkX(i)=temp2
  tvlY(i)=temp3
enddo

trAkX=ZERO
trAlY=ZERO
do i=1,n
  do j=1,n
    trAkX=trAkX+tAk(i,j)*X(j,i)
    trAlY=trAlY+tAl(i,j)*Y(j,i)
  enddo
enddo

do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    temp4=ZERO
    do k=1,n
      temp1=temp1+tAk(j,k)*X(k,i)
      temp2=temp2+tAl(j,k)*Y(k,i)  
      temp3=temp3+X(j,k)*tAk(k,i)
      temp4=temp4+Y(j,k)*tAl(k,i)
    enddo
    tAkX(j,i)=temp1
    tAlY(j,i)=temp2
    XtAk(j,i)=temp3
    YtAl(j,i)=temp4
  enddo
enddo

do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+tAkX(j,k)*tAk(k,i)
      temp2=temp2+tAlY(j,k)*tAl(k,i)     
    enddo
    P(j,i)=temp1
    Q(j,i)=temp2
  enddo
enddo

do i=1,n
  do j=1,n
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(j,k)*P(k,i)
      temp2=temp2+inv_tAkl(j,k)*Q(k,i)     
    enddo
    inv_tAklP(j,i)=temp1
    inv_tAklQ(j,i)=temp2
  enddo
enddo



trP=ZERO
trQ=ZERO
do i=1,n
  do j=1,n
    trP=trP+inv_tAkl(i,j)*P(j,i)
    trQ=trQ+inv_tAkl(i,j)*Q(j,i)
  enddo
enddo

Q1=ZERO
Q2=ZERO
Q5=ZERO
Q6=ZERO
P1=ZERO
P2=ZERO
P5=ZERO
P6=ZERO
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  do j=1,n
    temp1=temp1+tbkinv_tAkl(j)*P(j,i)
    temp2=temp2+tvkinv_tAkl(j)*P(j,i)   
    temp3=temp3+tbkinv_tAkl(j)*Q(j,i)
    temp4=temp4+tvkinv_tAkl(j)*Q(j,i)
  enddo
  P1=P1+temp1*inv_tAkltbl(i)
  P2=P2+temp2*inv_tAkltvl(i) 
  P5=P5+temp1*inv_tAkltvl(i)
  P6=P6+temp2*inv_tAkltbl(i)
  Q1=Q1+temp3*inv_tAkltbl(i)
  Q2=Q2+temp4*inv_tAkltvl(i) 
  Q5=Q5+temp3*inv_tAkltvl(i)
  Q6=Q6+temp4*inv_tAkltbl(i)
enddo
Qgamma=Q1*tau3+Q2*tau33+Q5*tau333+Q6*tau334
Pgamma=P1*tau3+P2*tau33+P5*tau333+P6*tau334

term3_1=ZERO
term3_2=ZERO
term4_1=ZERO
term4_2=ZERO
term11_1=ZERO
term11_2=ZERO
term16_1=ZERO
term16_2=ZERO
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  temp33=ZERO
  temp44=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAkl(j)*tAl(j,i)
    temp2=temp2+tbkinv_tAkl(j)*tAl(j,i)
    temp3=temp3+tbkinv_tAkl(j)*tAl(j,i)
    temp4=temp4+tvkinv_tAkl(j)*tAl(j,i)    
    temp33=temp33+tvkX(j)*tAk(j,i)
    temp44=temp44+tbkX(j)*tAk(j,i)
  enddo
  term3_1=term3_1+temp1*Ytvl(i)
  term3_2=term3_2+temp2*Ytvl(i)
  term4_1=term4_1+temp3*Ytbl(i)
  term4_2=term4_2+temp4*Ytbl(i)
  term11_1=term11_1+temp33*inv_tAkltvl(i)
  term11_2=term11_2+temp33*inv_tAkltbl(i)
  term16_1=term16_1+temp44*inv_tAkltbl(i)
  term16_2=term16_2+temp44*inv_tAkltvl(i)
enddo

do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  do j=1,n
    temp1=temp1+tvk(j)*inv_tAklP(j,i)
    temp2=temp2+tbk(j)*inv_tAklP(j,i)
    temp3=temp3+tvk(j)*inv_tAklQ(j,i)
    temp4=temp4+tbk(j)*inv_tAklQ(j,i)
  enddo
  tvkinv_tAklP(i)=temp1
  tbkinv_tAklP(i)=temp2
  tvkinv_tAklQ(i)=temp3
  tbkinv_tAklQ(i)=temp4
enddo
do i=1,n
  temp1=ZERO
  temp2=ZERO
  do j=1,n
    temp1=temp1+tvkX(j)*tAk(j,i)
    temp2=temp2+tbkX(j)*tAk(j,i)
  enddo
  tvkXtAk(i)=temp1
  tbkXtAk(i)=temp2
enddo
term8_1=ZERO
term8_2=ZERO
term9_1=ZERO
term9_2=ZERO
term12_1=ZERO
term12_2=ZERO
term17_1=ZERO
term17_2=ZERO
term13_1=ZERO
term14_2=ZERO
term18_2=ZERO
term19_1=ZERO
do i=1,n
  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO
  temp5=ZERO
  temp6=ZERO
  do j=1,n
    temp1=temp1+tvkinv_tAklP(j)*inv_tAkltAl(j,i)
    temp2=temp2+tbkinv_tAklP(j)*inv_tAkltAl(j,i)  
    temp3=temp3+tvkXtAk(j)*inv_tAklQ(j,i)
    temp4=temp4+tbkXtAk(j)*inv_tAklQ(j,i)
    temp5=temp5+tvkXtAk(j)*inv_tAkltAl(j,i)
    temp6=temp6+tbkXtAk(j)*inv_tAkltAl(j,i)
  enddo
  term8_1=term8_1+temp1*Ytvl(i)
  term8_2=term8_2+temp2*Ytvl(i)
  term9_1=term9_1+temp2*Ytbl(i)
  term9_2=term9_2+temp1*Ytbl(i)
  term12_1=term12_1+temp3*inv_tAkltvl(i)
  term12_2=term12_2+temp3*inv_tAkltbl(i)
  term17_1=term17_1+temp4*inv_tAkltbl(i)
  term17_2=term17_2+temp4*inv_tAkltvl(i)
  term13_1=term13_1+temp5*Ytvl(i)
  term14_2=term14_2+temp5*Ytbl(i)
  term18_2=term18_2+temp6*Ytvl(i)
  term19_1=term19_1+temp6*Ytbl(i)
enddo
term1=36*trAkX*trAlY*prod*gamma
term2=-24*trAkX*prod*(THREEHALF*trQ*gamma + Qgamma)
term3=24*trAkX*prod*(tau33*term3_1+tau333*term3_2)
term4=24*trAkX*prod*(tau3*term4_1+tau334*term4_2)
term6=-24*trAlY*prod*(THREEHALF*trP*gamma+Pgamma)
term7=16*rPr_rQr(P,Q,tvk,tbk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,tvkinv_tAkl,tbkinv_tAkl,inv_tAkltvl,inv_tAkltbl)
term8=-24*trP*prod*(tau33*term3_1+tau333*term3_2)-16*prod*(tau33*term8_1+P1*term3_1+P6*term3_2+tau333*term8_2)
term9=-24*trP*prod*(tau3*term4_1+tau334*term4_2)-16*prod*(P2*term4_1+tau3*term9_1+tau334*term9_2+P5*term4_2)
term11=24*trAlY*prod*(tau33*term11_1+tau334*term11_2)
term12=-24*trQ*prod*(tau33*term11_1+tau334*term11_2)-16*prod*(tau33*term12_1+Q1*term11_1+tau334*term12_2+Q5*term11_2)
term13=16*prod*(tau33*term13_1+term11_2*term3_2)
term14=16*prod*(term11_1*term4_1+tau334*term14_2)
term16=24*trAlY*prod*(tau3*term16_1+tau333*term16_2)
term17=-24*trQ*prod*(tau3*term16_1+tau333*term16_2)&
-16*prod*(Q2*term16_1+tau3*term17_1+Q6*term16_2+tau333*term17_2)
term18=16*prod*(term3_1*term16_1+tau333*term18_2)
term19=16*prod*(tau3*term19_1+term4_2*term16_1)

big_term1=term1+term2+term3+term4+term6
big_term2=term7+term8+term9+term11
big_term3=term12+term13+term14
big_term4=term16+term17+term18+term19
dXddYd= big_term1+big_term2+big_term3+big_term4

end function dXddYd

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

        ketMatrix(k, l, i) = real(allPermutations(l, k, i), kind=dprec)
        ! note the transposition here

      enddo
    enddo

  enddo

  do i = 1, nFactorial
    parityFactor(i) = real(parities(i), kind=dprec)
  enddo

end subroutine spinPreCalc



subroutine spinDependentMatrixElements(m_k, m_l, mm_k, mm_l, vechLk, vechLl, Pket, &
  SziME, SSNCspinME, SSNCmassChargeCoefficient, SOmassChargeCoefficient, &
  AMMmassChargeCoefficient, SSNCkl, SO1kl, SO2kl, &
  AMM1kl, AMM2kl, numberOfSpinFunctions)
!This subroutine computes symmetry adapted matrix element
!with two real L=1 correlated Gaussians. These matrix element
!is used in calculations of expectation values.

!Input:
!   m_k,m_l,mm_k, mm_l :: integers that determine which x or y-components is in the
!                premultiplier of the Gaussian
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.

!Output (all matrix elements are computed with normilized functions):

!   SSNCkl :: Non-contact spin-spin term (without the factor of alpha**2)
!   SO1kl, SO2kl  :: Spin-Orbit corrections (without the factor of alpha**2)
!         1 and 2 stay for spin-same orbit and spin-another orbit contributions
!   AMM1kl, AMM2kl  :: AMM corrections (without the factor of alpha**2)
!         1 and 2 stay for spin-same orbit and spin-another orbit contributions

!Arguments
integer,intent(in)       :: m_k, m_l, mm_k, mm_l, numberOfSpinFunctions
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
real(dprec)       temp1, temp2, temp3, temp1010, temp1001, temp0110, temp0101, t1, t2, det_tAkl
integer :: i, j, k, indx


integer :: pm_k, pm_l, pmm_k, pmm_l ! new non-zero components of v_k and v_l
real(dprec) :: commonFactor, gamma, gamma_diag, jiAlAklinvVl, &
              jjAlAklinvVl, localEps

! V-quantities 
real(dprec) :: VkAklinvVl, jiAkAklinvVl, jiAlAklinvVk, jiAklinvVk, jiAklinvVl, &
jjAkAklinvVl, jjAlAklinvVk, jjAklinvVk, jjAklinvVl

! W-quantities
real(dprec) :: WkAklinvWl, jiAklinvWk, &
jiAkAklinvWl, jiAlAklinvWk, jiAklinvWl, &
jjAkAklinvWl, jjAlAklinvWk, jjAklinvWk, jjAklinvWl

!mixed quantities
real(dprec) :: VkAklinvWl, WkAklinvVl

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
! new v_l = (P TRANSPOSED) * v_l
pm_k = m_k
pmm_k = mm_k

pm_l = m_l
pmm_l = mm_l
do i = 1, n
 if (abs(Pket(m_l, i) - 1.d0) < 1.d-13) pm_l = i
 if (abs(Pket(mm_l, i) - 1.d0) < 1.d-13) pmm_l = i
enddo


commonFactor = Glob_Piraised3n2 / (SQRTPI * det_tAkl * sqrt(det_tAkl))

SO1kl = ZERO
SO2kl = ZERO

AMM1kl = ZERO
AMM2kl = ZERO

SSNCkl = ZERO

do indexI = 1, n
 
 ! gamma diagonal coefficient
 gamma_diag = ONE / sqrt(inv_tAkl(indexI, indexI))

 ! calculating all the traces we need
 ! tr(Axy') is computed as (y, Ax) everywhere
 ! variable names: jiAlAklinvVk = (j^i, A_l A_{kl}^(-1) v_k) (names doesnt account for permutations)


 jiAkAklinvVl = ZERO
 do i = 1, n
   jiAkAklinvVl = jiAkAklinvVl + tAk(indexI, i) * inv_tAkl(i, pm_l)
 enddo
 jiAkAklinvWl = ZERO
 do i = 1, n
   jiAkAklinvWl = jiAkAklinvWl + tAk(indexI, i) * inv_tAkl(i, pmm_l)
 enddo

 jiAlAklinvVk = ZERO
 do i = 1, n
   jiAlAklinvVk = jiAlAklinvVk + tAl(indexI, i) * inv_tAkl(i, pm_k)
 enddo
 jiAlAklinvWk = ZERO
 do i = 1, n
   jiAlAklinvWk = jiAlAklinvWk + tAl(indexI, i) * inv_tAkl(i, pmm_k)
 enddo

 jiAklinvVl = inv_tAkl(indexI, pm_l)
 jiAklinvVk = inv_tAkl(indexI, pm_k)
 jiAklinvWl = inv_tAkl(indexI, pmm_l)
 jiAklinvWk = inv_tAkl(indexI, pmm_k)

 WkAklinvWl = inv_tAkl(pmm_k, pmm_l)
 WkAklinvVl = inv_tAkl(pmm_k, pm_l)
 VkAklinvWl = inv_tAkl(pm_k, pmm_l)
 VkAklinvVl = inv_tAkl(pm_k, pm_l)


 ! I term -> diagonal (spin-same orbit) matrix element (f[ii, ii])
 t1 = jiAklinvVk * jiAkAklinvVl + jiAlAklinvVk * jiAklinvVl
 temp1010 =  gamma_diag**3 / THREE * t1 *  WkAklinvWl - &
 gamma_diag**5 / FIVE *  t1 * (jiAklinvWk * jiAklinvWl)
 
 ! Vl <-> Wl
 t1 = jiAklinvVk * jiAkAklinvWl + jiAlAklinvVk * jiAklinvWl
 temp1001 = gamma_diag**3 / THREE * t1 * WkAklinvVl - &
 gamma_diag**5 / FIVE *t1 * (jiAklinvWk * jiAklinvVl)

 ! Vk <-> Wk
 t1 = jiAklinvWk * jiAkAklinvVl + jiAlAklinvWk * jiAklinvVl
 temp0110 = gamma_diag**3 / THREE * t1 * VkAklinvWl - &
 gamma_diag**5 / FIVE * t1* (jiAklinvVk * jiAklinvWl)

 ! Vl <-> Wl, Vk <-> Wk
 t1 = jiAklinvWk * jiAkAklinvWl + jiAlAklinvWk * jiAklinvWl
 temp0101 =  gamma_diag**3 / THREE * t1 * VkAklinvVl - &
 gamma_diag**5 / FIVE * t1 * (jiAklinvVk * jiAklinvVl)

 temp1 = temp1010 + temp0101 + temp1001 + temp0110
 
 do k = 1, numberOfSpinFunctions
   !if (abs(SziME(indexI, k)) < localEps) cycle
   SO1kl(k) = SO1kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
   AMM1kl(k) = AMM1kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
 enddo


 ! these traces are needed for spin-other orbit contribution and SSNC
 do indexJ = 1, n
   if (indexI == indexJ) cycle

   gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
   inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))

   jjAkAklinvVl = ZERO
   do i = 1, n
    jjAkAklinvVl = jjAkAklinvVl + tAk(indexJ, i) * inv_tAkl(i, pm_l)
   enddo
   jjAkAklinvWl = ZERO
   do i = 1, n
    jjAkAklinvWl = jjAkAklinvWl + tAk(indexJ, i) * inv_tAkl(i, pmm_l)
   enddo

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
   jjAklinvVl = inv_tAkl(indexJ, pm_l)
   jjAklinvWl = inv_tAkl(indexJ, pmm_l)


   !! II term -> f[ii, ij]
   t1 = jiAklinvVk * jjAkAklinvVl + jjAlAklinvVk * jiAklinvVl
   temp1010 = gamma_diag**3 / THREE * t1 * WkAklinvWl - &
   gamma_diag**5 / FIVE * t1 * (jiAklinvWk * jiAklinvWl)

   !Vl <-> Wl
   t1 = jiAklinvVk * jjAkAklinvWl + jjAlAklinvVk * jiAklinvWl
   temp1001 = gamma_diag**3 / THREE * t1 * WkAklinvVl - &
   gamma_diag**5 / FIVE * t1 * (jiAklinvWk * jiAklinvVl)

   !Vk <-> Wk
   t1 = jiAklinvWk * jjAkAklinvVl + jjAlAklinvWk * jiAklinvVl
   temp0110 = gamma_diag**3 / THREE * t1 * VkAklinvWl - &
   gamma_diag**5 / FIVE * t1 * (jiAklinvVk * jiAklinvWl)

   !Vk <-> Wk, Vl <-> Wl
   t1 = jiAklinvWk * jjAkAklinvWl + jjAlAklinvWk * jiAklinvWl
   temp0101 = gamma_diag**3 / THREE * t1 * VkAklinvVl - &
   gamma_diag**5 / FIVE * t1 * (jiAklinvVk * jiAklinvVl)

   temp1 = temp1010 + temp0101 + temp1001 + temp0110

   do k = 1, numberOfSpinFunctions
     !if (abs(SziME(indexI, k)) < localEps) cycle
     SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
   enddo
   

   !! III term -> f[ij, jj]
   t1 = jjAkAklinvVl * (jjAklinvVk - jiAklinvVk) + jjAlAklinvVk * (jjAklinvVl - jiAklinvVl)
   t2 = jjAklinvWk * jjAklinvWl + jiAklinvWk * jiAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
   temp1010 = gamma**3 / THREE * t1 * WkAklinvWl - &
   gamma**5 / FIVE * t1 * t2

   !Vl <-> Wl
   t1 = jjAkAklinvWl * (jjAklinvVk - jiAklinvVk) + jjAlAklinvVk * (jjAklinvWl - jiAklinvWl)
   t2 = jjAklinvWk * jjAklinvVl + jiAklinvWk * jiAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
   temp1001 = gamma**3 / THREE * t1 * WkAklinvVl - &
   gamma**5 / FIVE * t1 * t2

   !Vk <-> Wk
   t1 = jjAkAklinvVl * (jjAklinvWk - jiAklinvWk) + jjAlAklinvWk * (jjAklinvVl - jiAklinvVl)
   t2 = jjAklinvVk * jjAklinvWl + jiAklinvVk * jiAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
   temp0110 = gamma**3 / THREE * t1 * VkAklinvWl - &
   gamma**5 / FIVE * t1 * t2

   !!Vl <-> Wl, Vk <-> Wk
   t1 = jjAkAklinvWl * (jjAklinvWk - jiAklinvWk) + jjAlAklinvWk * (jjAklinvWl - jiAklinvWl)
   t2 = jjAklinvVk * jjAklinvVl + jiAklinvVk * jiAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
   temp0101 = gamma**3 / THREE * t1 * VkAklinvVl - &
   gamma**5 / FIVE * t1 * t2

   temp1 = temp1010 + temp0101 + temp1001 + temp0110

   do k = 1, numberOfSpinFunctions
     !if (abs(SziME(indexI, k)) < localEps) cycle
     SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
     AMM2kl(k) = AMM2kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
   enddo


   !! IV term -> f[ij, ii] (i <->j of the III term)
   t1 = jiAkAklinvVl * (jiAklinvVk - jjAklinvVk) + jiAlAklinvVk * (jiAklinvVl - jjAklinvVl)
   t2 = jjAklinvWk * jjAklinvWl + jiAklinvWk * jiAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
   temp1010 = gamma**3 / THREE * t1 * WkAklinvWl - &
   gamma**5 / FIVE * t1 * t2

   !Vl <-> Wl
   t1 = jiAkAklinvWl * (jiAklinvVk - jjAklinvVk) + jiAlAklinvVk * (jiAklinvWl - jjAklinvWl)
   t2 = jjAklinvWk * jjAklinvVl + jiAklinvWk * jiAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
   temp1001 = gamma**3 / THREE * t1 * WkAklinvVl - &
   gamma**5 / FIVE * t1 * t2

   !Vk <-> Wk
   t1 = jiAkAklinvVl * (jiAklinvWk - jjAklinvWk) + jiAlAklinvWk * (jiAklinvVl - jjAklinvVl)
   t2 = jjAklinvVk * jjAklinvWl + jiAklinvVk * jiAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
   temp0110 = gamma**3 / THREE * t1 * VkAklinvWl - &
   gamma**5 / FIVE * t1 * t2

   !!Vl <-> Wl, Vk <-> Wk
   t1 = jiAkAklinvWl * (jiAklinvWk - jjAklinvWk) + jiAlAklinvWk * (jiAklinvWl - jjAklinvWl)
   t2 = jjAklinvVk * jjAklinvVl + jiAklinvVk * jiAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
   temp0101 = gamma**3 / THREE * t1 * VkAklinvVl - &
   gamma**5 / FIVE * t1 * t2

   temp1 = temp1010 + temp0101 + temp1001 + temp0110

   do k = 1, numberOfSpinFunctions
     !if (abs(SziME(indexI, k)) < localEps) cycle
     SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
     AMM2kl(k) = AMM2kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
  enddo
  

  if (indexJ <= indexI) cycle !we need only indexJ > indexI for this term
  !SSNC term
  !g1010
  t1 = jiAklinvVk * jiAklinvVl + jjAklinvVk * jjAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
  t2 = jiAklinvWk * jiAklinvWl + jjAklinvWk * jjAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
  temp1010 = gamma**5 / FIVE * (t1 * WkAklinvWl) - &
  gamma**7 / SEVEN * (t1 * t2)

  !Vl <-> Wl
  t1 = jiAklinvVk * jiAklinvWl + jjAklinvVk * jjAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
  t2 = jiAklinvWk * jiAklinvVl + jjAklinvWk * jjAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
  temp1001 = gamma**5 / FIVE * (t1 * WkAklinvVl) - &
  gamma**7 / SEVEN * (t1 * t2)

  !Vk <-> Wk
  t1 = jiAklinvWk * jiAklinvVl + jjAklinvWk * jjAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
  t2 = jiAklinvVk * jiAklinvWl + jjAklinvVk * jjAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
  temp0110 = gamma**5 / FIVE * (t1 * VkAklinvWl) - &
  gamma**7 / SEVEN * (t1 * t2)

  !Vl <-> Wl, Vk <-> Wk
  t1 = jiAklinvVk * jiAklinvVl + jjAklinvVk * jjAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
  t2 = jiAklinvWk * jiAklinvWl + jjAklinvWk * jjAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
  temp0101 = gamma**5 / FIVE * (t2 * VkAklinvVl) - &
  gamma**7 / SEVEN * (t1 * t2)

  temp1 = ONETHIRD * (temp1010 + temp1001 + temp0110 + temp0101) !ONETHIRD - from common factor and spin part (1/sqrt(6)) 

  
  do k = 1, NumberOfSpinFunctions  
    !if (abs(SSNCspinME(indexI, indexJ, k)) < localEps) cycle
    SSNCkl(k) = SSNCkl(k) + SSNCspinME(indexI, indexJ, k) * SSNCmassChargeCoefficient(indexI, indexJ) * temp1
  enddo


 enddo ! indexJ cycle
enddo ! indexI cycle

SSNCkl = SSNCkl * commonFactor
SO1kl = SO1kl * commonFactor 
SO2kl = SO2kl * commonFactor
AMM1kl = AMM1kl * commonFactor
AMM2kl = AMM2kl * commonFactor


end subroutine spinDependentMatrixElements



subroutine overlapMatrixElementsLD(m_k, mm_k, vechLk, P, Skk)
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
  m=m1+m3  !look formula 40 in document

  !Evaluating overlap
  !temp1=ZERO
  temp1=FOUR*det_tAkl*sqrt(det_tAkl)
  Skk=Glob_Piraised3n2*m/temp1
  
  end subroutine overlapMatrixElementsLD


!function SG_ME_rXr_over_rij(i,j,X,inv_tAkl,det_tAkl)
!!function SG_ME_rXr_over_rij computes the following matrix element:
!!<\tilde psi_k| (r' X r)/r_ij |\tilde psi_l>
!!where psi_k = exp[-(r' Ak r)] is a Simple Gaussian wavefunction
!!Here X is a some real symmetric matrix. If matrix X is not symmetric
!!then user needs to symmetrize it before calling this function.
!!Input:
!!            X  :: n x n real matrix
!!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
!real(dprec)   SG_ME_rXr_over_rij
!integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!!Arguments:
!real(dprec)   X(nn,nn),inv_tAkl(nn,nn),det_tAkl
!integer       i,j
!!Local variables:
!integer       p,q,n
!real(dprec)   temp1,temp2,temp3
!real(dprec)   Aj(nn),AjX(nn)
!real(dprec)   t_J,t_X,t_XJ
!
!n=Glob_n
!!Form Aj=inv_tAkl*ji        j/=i 
!!     Aj=inv_tAkl*(ji-jj)   j/=i 
!!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
!if (i==j) then
! do p=1,n
!   Aj(p)=inv_tAkl(p,i)
! enddo
!else
! do p=1,n
!   Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
! enddo
!endif 
!
!!Compute AjX'=Aj'*X
!do p=1,n
!  temp1=ZERO
!  do q=1,n
!    temp1=temp1+Aj(q)*X(q,p)
!  enddo
!  AjX(p)=temp1
!enddo
!
!!Compute t_J=tr[inv_tAkl*Jij] 
!if (i==j) then
!  t_J=inv_tAkl(i,i)
!else
!  t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
!endif
!
!!Compute t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_XJ=ZERO
!do p=1,n
!  t_XJ=t_XJ+AjX(p)*Aj(p)
!enddo
!
!!Compute t_X=tr[inv_tAkl*X]
!t_X=ZERO
!do p=1,n
!  do q=1,n
!    t_X=t_X+inv_tAkl(q,p)*X(q,p)
!  enddo
!enddo
!
!temp1=Glob_Piraised3n2/(SQRTPI*det_tAkl**(THREEHALF))
!temp3=1/t_J
!SG_ME_rXr_over_rij=temp1*temp3*sqrt(temp3)*(THREE*t_J*t_X - t_XJ)
!
!end function SG_ME_rXr_over_rij
!function ME_d_X_over_rij_d1(i,j,X,tAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
!    
!real(dprec)   ME_d_X_over_rij_d1
!integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!!Arguments:
!real(dprec)   X(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn),det_tAkl
!integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)
!
!!Local variables:
!integer       c,s,p,q,n,k
!real(dprec)   Xtvl(nn),Xtbl(nn),tAlX(nn,nn),tAlXtAl(nn,nn),tAlXtvl(nn),tAlXtbl(nn)
!real(dprec)   temp1,temp2,temp3,temp4,t_tAlX
!
!
!!Doing multiplication Xtvl=X*tvl,tvkX=tvk*X
!n=Glob_n
!do s=1,n
!  temp1=ZERO
!  temp2=ZERO
!  do c=1,n
!	temp1=temp1+X(c,s)*tvl(c)
!        temp2=temp2+X(c,s)*tbl(c)
!  enddo
!  Xtvl(s)=temp1
!  Xtbl(s)=temp2
!enddo
!
!do s=1,n
!  do c=1,n
!    temp1=ZERO
!    do k=1,n
!      temp1=temp1+tAl(c,k)*X(k,s) 
!    enddo
!    tAlX(c,s)=temp1
!  enddo
!enddo
!
!!Doing multiplication tAlXtAk=tAlX*tAk,tAkXtAl=tAkX*tAl
!do s=1,n
!  do c=1,n
!    temp1=ZERO
!    do k=1,n
!      temp1=temp1+tAlX(c,k)*tAl(k,s)
!    enddo
!    tAlXtAl(c,s)=temp1
!  enddo
!enddo
!
!do s=1,n
!  temp1=ZERO
!  temp2=ZERO
!  do c=1,n
!	temp1=temp1+tAlX(c,s)*tvl(c)
!        temp2=temp2+tAlX(c,s)*tbl(c)
!  enddo
!  tAlXtvl(s)=temp1
!  tAlXtbl(s)=temp2
!enddo
!
!t_tAlX=ZERO
!do p=1,n
!  do q=1,n
!    t_tAlX=t_tAlX+tAl(q,p)*X(q,p)  
!  enddo
!enddo
!
!temp1=-6*t_tAlX*ME_over_rij(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
!temp2=4*ME_rXr_over_rij(i,j,tAlXtAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
!temp3=-4*ME_over_rij_tvl(i,j,inv_tAkl,det_tAkl,tvk,tAlXtvl,tbk,tbl)
!temp4=-4*ME_over_rij_tbl(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tAlXtbl)
!ME_d_X_over_rij_d1=-(temp1+temp2+temp3+temp4)
!end function ME_d_X_over_rij_d1

end module matelem
