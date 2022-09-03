module matelem
!Module matelem contains subroutines for computing 
!matrix elements with real L=1 (l_1=1, l_2=1) P-state Gaussians.
use globvars
implicit none

contains

subroutine MatrixElementsL1(jj_k, jj_l, vechLk, vechLl, P, &
               Hkl, Skl, Dk, Dl, grad_k, grad_l)

!(x_ki*y_kj-x_kj*y_ki)(x_li*y_lj-x_lj*y_li)


!**********
!
! Arguments
!
!**********

integer,intent(in)	:: jj_k, jj_l    !j_k, j_l
real(dprec),intent(in)  :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec),intent(in)  :: P(Glob_n,Glob_n)
real(dprec),intent(out) :: Skl, Hkl
real(dprec),intent(out) :: Dk(2*Glob_np),Dl(2*Glob_np)
logical,intent(in)      :: grad_k, grad_l
integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles

!****************
!
! Local Variables
!
!****************

integer         n
real(dprec)	Lk(nn,nn), LkI(3*nn,3*nn), Ak(nn,nn), AkI(3*nn,3*nn) 
real(dprec)	Ll(nn,nn), pLl(nn,nn), LlI(3*nn,3*nn), tAl(nn,nn), tAlI(3*nn,3*nn)
real(dprec)	tAkl(nn,nn), inv_tAkl(nn,nn), inv_tAklI(3*nn,3*nn), inv_tAklp(nn,nn)
real(dprec)     Jm(nn,nn)
real(dprec)	Pp(3*nn,3*nn)
real(dprec)	Wk(3*nn,3*nn), Wl(3*nn,3*nn), Wlp(3*nn,3*nn)
real(dprec)     beta1(3*nn,3*nn), beta2(3*nn,3*nn), beta3(3*nn,3*nn), beta33(3*nn,3*nn)
real(dprec)	beta4(nn,nn), beta4I(3*nn,3*nn), beta5(nn,nn), beta5I(3*nn,3*nn)
real(dprec)     beta6(3*nn,3*nn), beta7(3*nn,3*nn), beta8(3*nn,3*nn), beta9(3*nn,3*nn)
real(dprec)     beta10(3*nn,3*nn), beta11(3*nn,3*nn)
real(dprec)     beta12(3*nn,3*nn), beta13(3*nn,3*nn), beta14(3*nn,3*nn)
real(dprec)     betaJm(nn,nn), betaJmI(3*nn,3*nn), MI(3*nn,3*nn)
real(dprec)	beta16(3*nn,3*nn), beta17(3*nn,3*nn), beta18(3*nn,3*nn), beta19(3*nn,3*nn)
real(dprec)	beta20(3*nn,3*nn), beta21(3*nn,3*nn), beta22(3*nn,3*nn), beta23(3*nn,3*nn)
real(dprec)	beta24(3*nn,3*nn), beta25(3*nn,3*nn), beta26(3*nn,3*nn), beta35(3*nn,3*nn)
real(dprec)	beta27(3*nn,3*nn), beta28(3*nn,3*nn), beta29(3*nn,3*nn)
real(dprec)     beta30(3*nn,3*nn), beta31(3*nn,3*nn), beta32(3*nn,3*nn), beta34(3*nn,3*nn)
real(dprec)	beta37(3*nn,3*nn), beta38(3*nn,3*nn), beta15(3*nn,3*nn)
real(dprec)     X1(nn,nn), X1I(3*nn,3*nn)
real(dprec)	nu1, nu2, nu3, nu6, nu7
real(dprec)	tau1, tau6, tau7, tau8, tau9, tau10, Xi
real(dprec)     det_tAkl
real(dprec)	temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16
real(dprec)	temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27
real(dprec)     Z1(nn,nn), Z2(nn,nn), Z3(3*nn,3*nn), Z4(3*nn,3*nn), Z5(3*nn,3*nn), Z6(3*nn,3*nn), Z7(3*nn,3*nn), Z8(3*nn,3*nn)
real(dprec)     Z9(3*nn,3*nn), Z10(3*nn,3*nn), Z11(3*nn,3*nn), Z12(3*nn,3*nn), Z13(3*nn,3*nn), Z14(3*nn,3*nn), Z15(3*nn,3*nn)
real(dprec)     Z16(3*nn,3*nn)
real(dprec)	alpha1(3*nn,3*nn), alpha1p(3*nn,3*nn), alpha2(3*nn,3*nn)
real(dprec)	alpha5(nn,nn), alpha6(nn,nn), alpha6p(nn,nn), alpha2p(3*nn,3*nn)
real(dprec)	alpha7(3*nn,3*nn), alpha8(3*nn,3*nn), alpha9(3*nn,3*nn), alpha10(3*nn,3*nn)
real(dprec)     alpha11(3*nn,3*nn), alpha12(3*nn,3*nn), alpha13(3*nn,3*nn), alpha14(3*nn,3*nn)
real(dprec)	alpha15(3*nn,3*nn), alpha16(3*nn,3*nn), alpha17(3*nn,3*nn), alpha18(3*nn,3*nn) 
real(dprec)	alpha19(3*nn,3*nn), alpha20(3*nn,3*nn)
real(dprec)     alpha7p(3*nn,3*nn), alpha8p(3*nn,3*nn), alpha9p(3*nn,3*nn), alpha10p(3*nn,3*nn)
real(dprec)     alpha11p(3*nn,3*nn), alpha12p(3*nn,3*nn), alpha13p(3*nn,3*nn), alpha14p(3*nn,3*nn)
real(dprec)     alpha15p(3*nn,3*nn), alpha16p(3*nn,3*nn), alpha17p(3*nn,3*nn), alpha18p(3*nn,3*nn)
real(dprec)     alpha19p(3*nn,3*nn), alpha20p(3*nn,3*nn)
real(dprec)     alpha21(nn,nn), alpha22(3*nn,3*nn), alpha22T(3*nn,3*nn), alpha23(3*nn,3*nn), alpha23T(3*nn,3*nn)
real(dprec)	alpha24(3*nn,3*nn), alpha24T(3*nn,3*nn)
real(dprec)     alpha21p(nn,nn), alpha22p(3*nn,3*nn), alpha22Tp(3*nn,3*nn), alpha23p(3*nn,3*nn), alpha23Tp(3*nn,3*nn)
real(dprec)     alpha24p(3*nn,3*nn), alpha24Tp(3*nn,3*nn), alpha28(3*nn,3*nn), alpha28p(3*nn,3*nn) 
real(dprec)     alpha25(nn,nn), alpha25p(nn,nn), alpha26(3*nn,3*nn), alpha26p(3*nn,3*nn), alpha27(3*nn,3*nn), alpha27p(3*nn,3*nn)
real(dprec)	alpha29(3*nn,3*nn), alpha29p(3*nn,3*nn), alpha30(3*nn,3*nn), alpha30p(3*nn,3*nn)
real(dprec)     alpha31(3*nn,3*nn), alpha31p(3*nn,3*nn)
real(dprec)     Tkl, Rkl, Vkl, pk(Glob_np), pl(Glob_np)
real(dprec)	trans(3*Glob_n*(3*Glob_n+1)/2,Glob_np), trans1(3*Glob_n*(3*Glob_n+1)/2,Glob_np)
real(dprec)     result1k(Glob_n*(Glob_n+1)/2), result2k(3*GLob_n*(3*Glob_n+1)/2)
real(dprec)     result3k(Glob_n*(Glob_n+1)/2), result4k(Glob_n*(Glob_n+1)/2), result5k(3*GLob_n*(3*Glob_n+1)/2)
real(dprec)	result1l(Glob_n*(Glob_n+1)/2), result2l(3*GLob_n*(3*Glob_n+1)/2)
real(dprec)     result3l(Glob_n*(Glob_n+1)/2), result4l(Glob_n*(Glob_n+1)/2), result5l(3*GLob_n*(3*Glob_n+1)/2)
real(dprec)	resultJk(Glob_n*(Glob_n+1)/2), resultJl(Glob_n*(Glob_n+1)/2)
real(dprec)	result6k(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     result6l(3*Glob_n*(3*Glob_n+1)/2)
integer		MMsum, Msum, M1, M2, M3
integer         i, j, k, r, t, indx, v, coef, mk, nk, rk, sk, nl, ml, rl, sl 
integer         i_k, i_l, j_k, j_l
!integer         i_k, i_l
real(dprec)     SSkl

9292 format(1x,18f10.7)
9293 format(1x,21f4.1)
n=Glob_n


i_k = 1
i_l = 1
j_k = jj_k
j_l = jj_l
if (jj_k.eq.1) j_k = 2
if (jj_l.eq.1) j_l = 2

!write (*,*) ' i_k, j_k, i_l, j_l',  i_k, j_k, i_l, j_l

!********************
!
! END OF DECLARATIONS
!
!********************

!write(*,*) 'Glob_np, Glob_n', Glob_np, Glob_n

!***************************
!
! CONSTRUCTION OF W MATRICES
!
!***************************

!(x_ki*y_kj-x_kj*y_ki)(x_li*y_lj-x_lj*y_li)

Wk=ZERO
Wl=ZERO
Hkl=ZERO
Dk=ZERO
Dl=ZERO

   Wk((i_k-1)*3+1,(j_k-1)*3+2)=ONEHALF
   Wk((j_k-1)*3+2,(i_k-1)*3+1)=ONEHALF
   Wk((j_k-1)*3+1,(i_k-1)*3+2)=-ONEHALF
   Wk((i_k-1)*3+2,(j_k-1)*3+1)=-ONEHALF
!nk=ONE
!mk=FIVE
!rk=TWO
!sk=FOUR
nk=(i_k-1)*3+1
mk=(j_k-1)*3+2
rk=(i_k-1)*3+2
sk=(j_k-1)*3+1

   Wl((i_l-1)*3+1,(j_l-1)*3+2)=ONEHALF
   Wl((j_l-1)*3+2,(i_l-1)*3+1)=ONEHALF
   Wl((j_l-1)*3+1,(i_l-1)*3+2)=-ONEHALF
   Wl((i_l-1)*3+2,(j_l-1)*3+1)=-ONEHALF
!nl=ONE
!ml=FIVE
!rl=TWO
!sl=FOUR
nl=(i_l-1)*3+1
ml=(j_l-1)*3+2
rl=(i_l-1)*3+2
sl=(j_l-1)*3+1

!write(*,*)'n,nn',n,nn

!write(*,*)'i_k,j_k,i_l,j_l',i_k,j_k,i_l,j_l
!write(*,*)'nk,mk,rk,sk',nk,mk,rk,sk

!write(*,*)'Wk'
!do i=1,3*n
!  write(*,9292)(Wk(i,j),j=1,3*n)
!enddo

!write(*,*)'Wl'
!do i=1,3*n
!  write(*,9292)(Wl(i,j),j=1,3*n)
!enddo


!***************************
!
! CONSTRUCTION OF A MATRICES
! 
!***************************

!First built are matrices Lk, LkI and  Ll, LlI from vechL and vechLl.

indx=ZERO

do i=1,n
  do j=i,n
      indx=indx+1

      Lk(i,j)=ZERO
      Lk(j,i)=vechLk(indx)
     
      Ll(i,j)=ZERO
      Ll(j,i)=vechLl(indx)

      if (grad_k.or.grad_l) then 
        LkI(3*(i-1)+1,3*(j-1)+1)=ZERO
        LkI(3*(j-1)+1,3*(i-1)+1)=vechLk(indx)
     
        LkI(3*(i-1)+2,3*(j-1)+2)=ZERO
        LkI(3*(j-1)+2,3*(i-1)+2)=vechLk(indx)
     
        LkI(3*(i-1)+3,3*(j-1)+3)=ZERO
        LkI(3*(j-1)+3,3*(i-1)+3)=vechLk(indx)

        LlI(3*(i-1)+1,3*(j-1)+1)=ZERO
        LlI(3*(j-1)+1,3*(i-1)+1)=vechLl(indx)

        LlI(3*(i-1)+2,3*(j-1)+2)=ZERO
        LlI(3*(j-1)+2,3*(i-1)+2)=vechLl(indx)
     
        LlI(3*(i-1)+3,3*(j-1)+3)=ZERO
        LlI(3*(j-1)+3,3*(i-1)+3)=vechLl(indx)

      endif

  enddo 
enddo

Pp=ZERO
MI=ZERO
do i=1,n
  do j=1,n

      Pp(3*(i-1)+1,3*(j-1)+1)=P(i,j)
      Pp(3*(i-1)+2,3*(j-1)+2)=P(i,j)
      Pp(3*(i-1)+3,3*(j-1)+3)=P(i,j)

      MI(3*(i-1)+1,3*(j-1)+1)=Glob_MassMatrix(i,j)
      MI(3*(i-1)+2,3*(j-1)+2)=Glob_MassMatrix(i,j)
      MI(3*(i-1)+3,3*(j-1)+3)=Glob_MassMatrix(i,j)


  enddo
enddo

!write(*,*)'Lk'
!do i=1,n
!  write(*,9292)(Lk(i,j),j=1,n)
!enddo

!write(*,*)'Ll'
!do i=1,n
!  write(*,9292)(Ll(i,j),j=1,n)
!enddo

!write(*,*)'LkI'
!do i=1,3*n
!  write(*,9292)(LkI(i,j),j=1,3*n)
!enddo

!write(*,*)'Pp'
!do i=1,3*n
!  write(*,9292)(Pp(i,j),j=1,3*n)
!enddo

!write(*,*)'MI'
!do i=1,3*n
!  write(*,9292)(MI(i,j),j=1,3*n)
!enddo

!Then we build matrices Ak and tAl

do i=1,n
  do j=i,n

      temp1=ZERO
      temp2=ZERO
      
      do k=1,j

        temp1=temp1+Lk(i,k)*Lk(j,k)
	temp2=temp2+Ll(i,k)*Ll(j,k)

      enddo

      Ak(i,j)=temp1
      Ak(j,i)=temp1

      tAl(i,j)=temp2
      tAl(j,i)=temp2
  
  enddo
enddo

!write(*,*)'Ak'
!do i=1,n
!  write(*,9292)(Ak(i,j),j=1,n)
!enddo

!write(*,*)'tAl'
!do i=1,n
!  write(*,9292)(tAl(i,j),j=1,n)
!enddo

!Next elements of Al and Wl are permuted to 
!account for the action of the permutation
!matrix tAl=P'*Al*P.
!Also formed is matrix tAkl=Ak+tAl.

do i=1,n
  do j=1,n
    
    temp1=ZERO
    temp2=ZERO

    do k=1,n
      temp1=temp1+P(k,j)*tAl(k,i)
      temp2=temp2+P(k,j)*Ll(k,i)
    enddo
      Z1(j,i)=temp1
      Z2(j,i)=temp2
  
  enddo
enddo

do i=1,n
  do j=i,n
    
    temp1=ZERO

    do k=1,n

      temp1=temp1+Z1(i,k)*P(k,j)

    enddo

    tAl(i,j)=temp1
    tAl(j,i)=temp1

    tAkl(i,j)=Ak(i,j)+temp1
    tAkl(j,i)=tAkl(i,j)

 enddo
enddo

do i=1,n
  do j=1,n

    temp1=ZERO
    do k=1,n
      temp1=temp1+Z2(i,k)*P(k,j)
    enddo
     pLl(i,j)=temp1
 enddo
enddo

!write(*,*)'tAl'
!do i=1,n
!  write(*,9292)(tAl(i,j),j=1,n)
!enddo

!write(*,*)'tAkl'
!do i=1,n
!  write(*,9292)(tAkl(i,j),j=1,n)
!enddo

!write(*,*)'pLl'
!do i=1,n
!  write(*,9292)(pLl(i,j),j=1,n)
!enddo

Z3=ZERO
do i=1,3*n
  do j=1,3*n

    temp1=ZERO
    do k=1,3*n
      temp1=temp1+Pp(k,j)*Wl(k,i)
    enddo
    Z3(j,i)=temp1

  enddo
enddo

   nl=ZERO
   ml=ZERO
   rl=ZERO
   sl=ZERO

do i=1,3*n
  do j=i,3*n

    temp1=ZERO
    do k=1,3*n
      temp1=temp1+Z3(i,k)*Pp(k,j)
    enddo

    Wlp(i,j)=temp1
    Wlp(j,i)=temp1

    if(temp1.gt.ONEFOURTH)then
      nl=i
      ml=j
    endif

    if(temp1.lt.-ONEFOURTH)then
      rl=i
      sl=j
    endif

  enddo
enddo

!write(*,*)'Wlp'
!do i=1,3*n
!  write(*,9292)(Wlp(i,j),j=1,3*n)
!enddo

!After this we can do the Cholesky factorization of tAkl.
!The Cholesky factor will be temporarily stored in the 
!lower triangle of W1

det_tAkl=ONE

do i=1,n
  do j=i,n
    
    temp1=tAkl(i,j)
    do k=i-1,1,-1
      temp1=temp1-Z1(i,k)*Z1(j,k)
    enddo
    
    if (i==j) then
      Z1(i,i)=sqrt(temp1)
      det_tAkl=det_tAkl*temp1
    else
      Z1(j,i)=temp1/Z1(i,i)
      Z1(i,j)=ZERO
    endif
  
  enddo
enddo

!Inverting tAkl and using its Cholesky factor
!(stored in Z1) and placing into inv_tAkl

do i=1,n
  
  Z1(i,i)=ONE/Z1(i,i)
  do j=i+1,n
    
    temp1=ZERO
    do k=i,j-1
      temp1=temp1-Z1(j,k)*Z1(k,i)
    enddo
      Z1(j,i)=temp1/Z1(j,j)
  
  enddo
enddo

do i=1,n
  do j=i,n
    
    temp1=ZERO
    do k=j,n
      temp1=temp1+Z1(k,i)*Z1(k,j)
    enddo
      inv_tAkl(i,j)=temp1
      inv_tAkl(j,i)=temp1
  
  enddo
enddo

!write(*,*)'inv_tAkl'
!do i=1,n
!  write(*,9292)(inv_tAkl(i,j),j=1,n)
!enddo

!Performing Kronecker Product of inv_tAkl with I3

inv_tAklI=ZERO

do i=1,n
  do j=1,n
    
    inv_tAklI(3*(i-1)+1,3*(j-1)+1)=inv_tAkl(i,j)
    inv_tAklI(3*(i-1)+2,3*(j-1)+2)=inv_tAkl(i,j)
    inv_tAklI(3*(i-1)+3,3*(j-1)+3)=inv_tAkl(i,j)

  enddo
enddo

!write(*,*)'inv_tAklI'
!do i=1,3*n
!  write(*,9292)(inv_tAkl(i,j),j=1,3*n)
!enddo

!****************************
!
! DOING MATRIX MULTIPLICATION
! 
!****************************

!beta1=inv_tAklI*Wk
!beta2=inv_tAklI*Wlp
!beta6=MI*Wlp

beta1=ZERO
beta2=ZERO
beta6=ZERO

do i=1,3*n
    
  beta1(i,nk)=inv_tAklI(i,mk)*Wk(mk,nk)
  beta1(i,mk)=inv_tAklI(i,nk)*Wk(nk,mk)
  beta1(i,rk)=inv_tAklI(i,sk)*Wk(sk,rk)
  beta1(i,sk)=inv_tAklI(i,rk)*Wk(rk,sk)

  beta2(i,nl)=inv_tAklI(i,ml)*Wlp(ml,nl)
  beta2(i,ml)=inv_tAklI(i,nl)*Wlp(nl,ml)  
  beta2(i,rl)=inv_tAklI(i,sl)*Wlp(sl,rl)
  beta2(i,sl)=inv_tAklI(i,rl)*Wlp(rl,sl)
 
  beta6(i,nl)=MI(i,ml)*Wlp(ml,nl)
  beta6(i,ml)=MI(i,nl)*Wlp(nl,ml)
  beta6(i,rl)=MI(i,sl)*Wlp(sl,rl)
  beta6(i,sl)=MI(i,rl)*Wlp(rl,sl)

enddo

!write(*,*)'beta6'
!do i=1,3*n
!  write(*,9292)(beta6(i,j),j=1,3*n)
!enddo

!beta3=beta1*beta2

beta3=ZERO
beta33=ZERO

do i=1,3*n

  beta3(i,nl)=beta1(i,nk)*beta2(nk,nl)+beta1(i,mk)*beta2(mk,nl)+beta1(i,rk)*beta2(rk,nl)+beta1(i,sk)*beta2(sk,nl)
  beta3(i,ml)=beta1(i,nk)*beta2(nk,ml)+beta1(i,mk)*beta2(mk,ml)+beta1(i,rk)*beta2(rk,ml)+beta1(i,sk)*beta2(sk,ml)
  beta3(i,rl)=beta1(i,nk)*beta2(nk,rl)+beta1(i,mk)*beta2(mk,rl)+beta1(i,rk)*beta2(rk,rl)+beta1(i,sk)*beta2(sk,rl)
  beta3(i,sl)=beta1(i,nk)*beta2(nk,sl)+beta1(i,mk)*beta2(mk,sl)+beta1(i,rk)*beta2(rk,sl)+beta1(i,sk)*beta2(sk,sl)

  beta33(i,nk)=beta2(i,nl)*beta1(nl,nk)+beta2(i,ml)*beta1(ml,nk)+beta2(i,rl)*beta1(rl,nk)+beta2(i,sl)*beta1(sl,nk)
  beta33(i,mk)=beta2(i,nl)*beta1(nl,mk)+beta2(i,ml)*beta1(ml,mk)+beta2(i,rl)*beta1(rl,mk)+beta2(i,sl)*beta1(sl,mk)
  beta33(i,rk)=beta2(i,nl)*beta1(nl,rk)+beta2(i,ml)*beta1(ml,rk)+beta2(i,rl)*beta1(rl,rk)+beta2(i,sl)*beta1(sl,rk)
  beta33(i,sk)=beta2(i,nl)*beta1(nl,sk)+beta2(i,ml)*beta1(ml,sk)+beta2(i,rl)*beta1(rl,sk)+beta2(i,sl)*beta1(sl,sk)

enddo

!write(*,*)'beta3'
!do i=1,3*n
!  write(*,9292)(beta3(i,j),j=1,3*n)
!enddo

!write(*,*)'beta33'
!do i=1,3*n
!  write(*,9292)(beta33(i,j),j=1,3*n)
!enddo


!Finding trace
!nu3=tr[beta3]

nu3=ZERO

do i=1,3*n
   
   nu3=nu3+beta3(i,i)

enddo

!******************
!
! COMPUTING OVERLAP
! (non-normalized)
!******************

Skl=ONEHALF*(PI**((THREE*n)/TWO))&
*(det_tAkl**(-THREE/TWO))*nu3

!write(*,*) 'Skl', Skl

!*********************************************************
!
!Derivative of Overlap with respect to vechL_k and vechL_l
!
!*********************************************************

!*******************
!
! Constructing Trans
!
!*******************

if (grad_k.or.grad_l) then

trans=ZERO

!trans(1,1)=ONE
!trans(4,2)=ONE
!trans(7,1)=ONE
!trans(10,2)=ONE
!trans(12,1)=ONE
!trans(15,2)=ONE
!trans(16,3)=ONE
!trans(19,3)=ONE
!trans(21,3)=ONE


!trans(1,1)=ONE
!trans(19,1)=ONE
!trans(36,1)=ONE
!trans(4,2)=ONE
!trans(22,2)=ONE
!trans(39,2)=ONE
!trans(7,3)=ONE
!trans(25,3)=ONE
!trans(42,3)=ONE
!trans(10,4)=ONE
!trans(28,4)=ONE
!trans(45,4)=ONE
!trans(13,5)=ONE
!trans(31,5)=ONE
!trans(48,5)=ONE
!trans(16,6)=ONE
!trans(34,6)=ONE
!trans(51,6)=ONE
!trans(52,7)=ONE
!trans(67,7)=ONE
!trans(81,7)=ONE
!trans(55,8)=ONE
!trans(70,8)=ONE
!trans(84,8)=ONE
!trans(58,9)=ONE
!trans(73,9)=ONE
!trans(87,9)=ONE
!trans(61,10)=ONE
!trans(76,10)=ONE
!trans(90,10)=ONE
!trans(64,11)=ONE
!trans(79,11)=ONE
!trans(93,11)=ONE
!trans(94,12)=ONE
!trans(106,12)=ONE
!trans(117,12)=ONE
!trans(97,13)=ONE
!trans(109,13)=ONE
!trans(120,13)=ONE
!trans(100,14)=ONE
!trans(112,14)=ONE
!trans(123,14)=ONE
!trans(103,15)=ONE
!trans(115,15)=ONE
!trans(126,15)=ONE
!trans(127,16)=ONE
!trans(136,16)=ONE
!trans(144,16)=ONE
!trans(130,17)=ONE
!trans(139,17)=ONE
!trans(147,17)=ONE
!trans(133,18)=ONE
!trans(142,18)=ONE
!trans(150,18)=ONE
!trans(151,19)=ONE
!trans(157,19)=ONE
!trans(162,19)=ONE
!trans(154,20)=ONE
!trans(160,20)=ONE
!trans(165,20)=ONE
!trans(166,21)=ONE
!trans(169,21)=ONE
!trans(171,21)=ONE

!write(*,*)'Trans'
!do i=1,3*n*(3*n+1)/2
!write (*,9293)(trans(i,j),j=1,n*(n+1)/2) 
!enddo

trans=ZERO
MMsum=ZERO
do j=1,n
  do i=j,n
    MMsum=MMsum+1
    Msum=ZERO
   do k=1,j-1
     Msum=Msum+(3*n-3*k+3)+(3*n-3*k+2)+(3*n-3*k+1)
   enddo
   M1=Msum+(3*(i-j)+1)
   M2=Msum+(3*n-3*(j-1))+(3*(i-j)+1)
   M3=Msum+(3*n-3*(j-1))+(3*n-3*(j-1)-1)+(3*(i-j)+1)
  trans(M1,MMsum)=ONE
  trans(M2,MMsum)=ONE
  trans(M3,MMsum)=ONE

!if (i.eq.2.and.j.eq.2) then
!write(*,*)'i,j,M1,M2,M3,MMsum',i,j,M1,M2,M3,MMsum
!endif  

  enddo
enddo

!write(*,*)'Trans1'
!do i=1,3*n*(3*n+1)/2
!write(*,9293)(trans1(i,j),j=1,n*(n+1)/2)
!enddo

!write(*,*)'comparing'
!do i=1,3*n*(3*n+1)/2
!  do j=1,n*(n+1)/2
!    if (abs(trans1(i,j)-trans(i,j)).gt.0.1D+00)write(*,*)i,j
!  enddo
!enddo

endif
 
!result1k=inv_tAkl*Lk

if (grad_k) then 
  result1k=ZERO
  indx=ZERO
  do j=1,n
    do i=j,n
      indx=indx+1
      temp1=ZERO
      do k=j,n

        temp1=temp1+(inv_tAkl(k,i)+inv_tAkl(i,k))*Lk(k,j)

      enddo

      result1k(indx)=temp1

    enddo
  enddo
endif

!alpha1=beta3*inv_tAklI

if (grad_k.or.grad_l) then

  alpha1=ZERO
  alpha2=ZERO

  do i=1,3*n
    do j=1,3*n

      alpha1(i,j)=beta3(i,nl)*inv_tAklI(nl,j)&
                 +beta3(i,ml)*inv_tAklI(ml,j)&
                 +beta3(i,rl)*inv_tAklI(rl,j)&
                 +beta3(i,sl)*inv_tAklI(sl,j)
  
      alpha2(i,j)=beta33(i,nk)*inv_tAklI(nk,j)&
                 +beta33(i,mk)*inv_tAklI(mk,j)&
                 +beta33(i,rk)*inv_tAklI(rk,j)&
                 +beta33(i,sk)*inv_tAklI(sk,j)

    enddo
  enddo
endif

!result2k=alpha1*LkI

if (grad_k) then

  result2k=ZERO
  indx=ZERO
  
  do j=1,3*n
    do i=j,3*n
    
      indx=indx+1
      temp1=ZERO
      temp2=ZERO
      do k=j,3*n

        temp1=temp1+(alpha1(i,k)+alpha1(k,i))*LkI(k,j)
        temp2=temp2+(alpha2(i,k)+alpha2(k,i))*LkI(k,j)

      enddo

      result2k(indx)=temp1+temp2

    enddo
  enddo

  do i=1,n*(n+1)/2
    
    temp1=ZERO
    do k=1,3*n*(3*n+1)/2
    
       temp1=temp1+result2k(k)*trans(k,i)
    
    enddo
    
    Dk(Glob_np+i)=Dk(Glob_np+i)-ONEHALF*(PI**(THREEHALF*n))*(det_tAkl**(-THREEHALF))*(THREEHALF*result1k(i)*nu3+temp1)
  
  enddo
endif

!write(*,*) 'Dk',Dk

!result1l=inv_tAkl*pLl

if (grad_l) then 
  Z1=ZERO
  do i=1,n
    do j=1,n
      temp1=ZERO
      do k=1,n  
        temp1=temp1+P(i,k)*inv_tAkl(k,j)
      enddo  
      Z1(i,j)=temp1
    enddo
  enddo  
  do i=1,n
    do j=1,n
      temp1=ZERO
      do k=1,n  
        temp1=temp1+Z1(i,k)*P(j,k)
      enddo  
      inv_tAklp(i,j)=temp1
    enddo
  enddo  

  result1l=ZERO
  indx=ZERO
  do j=1,n
    do i=j,n
      
      indx=indx+1
      temp1=ZERO
      do k=j,n

        temp1=temp1+(inv_tAklp(i,k)+inv_tAklp(k,i))*Ll(k,j)

      enddo

      result1l(indx)=temp1

    enddo
  enddo

  Z3=ZERO
  Z4=ZERO
  do i=1,3*n
    do j=1,3*n
      temp1=ZERO
      temp2=ZERO
      do k=1,3*n  
        temp1=temp1+pP(i,k)*alpha1(k,j)
        temp2=temp2+pP(i,k)*alpha2(k,j)
      enddo  
      Z3(i,j)=temp1
      Z4(i,j)=temp2
    enddo
  enddo  
  do i=1,3*n
    do j=1,3*n
      temp1=ZERO
      temp2=ZERO
      do k=1,3*n  
        temp1=temp1+Z3(i,k)*pP(j,k)
        temp2=temp2+Z3(i,k)*pP(j,k)
      enddo  
      alpha1p(i,j)=temp1
      alpha2p(i,j)=temp2
    enddo
  enddo  
 
!result2l=alpha1*LlI
 
  result2l=ZERO
  indx=ZERO
  do j=1,3*n
    do i=j,3*n
      
      indx=indx+1
      temp1=ZERO
      temp2=ZERO
      do k=j,3*n

        temp1=temp1+(alpha1p(i,k)+alpha1p(k,i))*LlI(k,j)
        temp2=temp2+(alpha2p(i,k)+alpha2p(k,i))*LlI(k,j)

      enddo

      result2l(indx)=temp1+temp2

    enddo
  enddo

  do i=1,n*(n+1)/2
    temp1=ZERO

    do k=1,3*n*(3*n+1)/2

      temp1=temp1+result2l(k)*trans(k,i)

    enddo

    Dl(Glob_np+i)=Dl(Glob_np+i)-ONEHALF*(PI**(THREEHALF*n))*&
(det_tAkl**(-THREEHALF))*(THREEHALF*result1l(i)*nu3+temp1)

  enddo
endif

!write(*,*) 'Dl',Dl

!go to 1111

!***********************************
!
! EVALUATING KINETIC ENERGY INTEGRAL
!
!***********************************

!****************************
!
! DOING MATRIX MULTIPLICATION
!
!****************************

!beta4=inv_tAkl*Ak
!beta5=M*tAl

do i=1,n
  do j=1,n
    
    temp1=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(i,k)*Ak(k,j)
    enddo
      beta4(i,j)=temp1
   
    temp1=ZERO
    do k=1,n
      temp1=temp1+Glob_MassMatrix(i,k)*tAl(k,j)
    enddo
      beta5(i,j)=temp1

  enddo
enddo 

!write(*,*)'beta4'
!do i=1,n
!  write(*,9292)(beta4(i,j),j=1,n)
!enddo

!write(*,*)'beta5'
!do i=1,n
!  write(*,9292)(beta5(i,j),j=1,n)
!enddo

!X1=beta4*beta5

do i=1,n
  do j=1,n
   
    temp1=ZERO
    do k=1,n
      temp1=temp1+beta4(i,k)*beta5(k,j)
    enddo
      X1(i,j)=temp1

  enddo
enddo 

!write(*,*)'X1'
!do i=1,n
!  write(*,9292)(X1(i,j),j=1,n)
!enddo

!Performing Kroneker Product of X1, beta4, beta5  with I3

X1I=ZERO
beta4I=ZERO
beta5I=ZERO

do i=1,n  
  do j=1,n
  
    X1I(3*(i-1)+1,3*(j-1)+1)=X1(i,j)
    X1I(3*(i-1)+2,3*(j-1)+2)=X1(i,j)
    X1I(3*(i-1)+3,3*(j-1)+3)=X1(i,j)

    beta4I(3*(i-1)+1,3*(j-1)+1)=beta4(i,j)
    beta4I(3*(i-1)+2,3*(j-1)+2)=beta4(i,j)
    beta4I(3*(i-1)+3,3*(j-1)+3)=beta4(i,j)

    beta5I(3*(i-1)+1,3*(j-1)+1)=beta5(i,j)
    beta5I(3*(i-1)+2,3*(j-1)+2)=beta5(i,j)
    beta5I(3*(i-1)+3,3*(j-1)+3)=beta5(i,j)


  enddo
enddo

!beta7=beta1*X1I
!beta10=beta1*beta4I
!beta12=beta1*beta5I
!beta14=beta1*beta6

beta7=ZERO
beta10=ZERO
beta12=ZERO
beta14=ZERO

do i=1,3*n
  do j=1,3*n

    beta7(i,j)=beta1(i,nk)*X1I(nk,j)+beta1(i,mk)*X1I(mk,j)+beta1(i,rk)*X1I(rk,j)+beta1(i,sk)*X1I(sk,j)
    beta10(i,j)=beta1(i,nk)*beta4I(nk,j)+beta1(i,mk)*beta4I(mk,j)+beta1(i,rk)*beta4I(rk,j)+beta1(i,sk)*beta4I(sk,j)
    beta12(i,j)=beta1(i,nk)*beta5I(nk,j)+beta1(i,mk)*beta5I(mk,j)+beta1(i,rk)*beta5I(rk,j)+beta1(i,sk)*beta5I(sk,j)
    beta14(i,j)=beta1(i,nk)*beta6(nk,j)+beta1(i,mk)*beta6(mk,j)+beta1(i,rk)*beta6(rk,j)+beta1(i,sk)*beta6(sk,j)

  enddo
enddo

!beta8=beta7*beta2
!beta9=X1I*beta3
!beta13=beta12*beta2
!beta11=beta10*beta6

beta8=ZERO
beta9=ZERO
beta13=ZERO
beta11=ZERO

do i=1,3*n

  temp1=ZERO
  temp2=ZERO
  temp3=ZERO
  temp4=ZERO

  temp5=ZERO
  temp6=ZERO
  temp7=ZERO
  temp8=ZERO

  temp9=ZERO
  temp10=ZERO
  temp11=ZERO
  temp12=ZERO

  temp13=ZERO
  temp14=ZERO
  temp15=ZERO
  temp16=ZERO

  do k=1,3*n

    temp1=temp1+beta7(i,k)*beta2(k,nl)
    temp2=temp2+beta7(i,k)*beta2(k,ml)
    temp3=temp3+beta7(i,k)*beta2(k,rl)
    temp4=temp4+beta7(i,k)*beta2(k,sl)
    
    temp5=temp5+beta12(i,k)*beta2(k,nl)
    temp6=temp6+beta12(i,k)*beta2(k,ml)
    temp7=temp7+beta12(i,k)*beta2(k,rl)
    temp8=temp8+beta12(i,k)*beta2(k,sl)

    temp9=temp9+beta10(i,k)*beta6(k,nl)
    temp10=temp10+beta10(i,k)*beta6(k,ml)
    temp11=temp11+beta10(i,k)*beta6(k,rl)
    temp12=temp12+beta10(i,k)*beta6(k,sl)

    temp13=temp13+X1I(i,k)*beta3(k,nl)
    temp14=temp14+X1I(i,k)*beta3(k,ml)
    temp15=temp15+X1I(i,k)*beta3(k,rl)
    temp16=temp16+X1I(i,k)*beta3(k,sl)
 
  enddo

    beta8(i,nl)=temp1
    beta8(i,ml)=temp2
    beta8(i,rl)=temp3
    beta8(i,sl)=temp4
 
    beta13(i,nl)=temp5
    beta13(i,ml)=temp6
    beta13(i,rl)=temp7
    beta13(i,sl)=temp8

    beta11(i,nl)=temp9
    beta11(i,ml)=temp10
    beta11(i,rl)=temp11
    beta11(i,sl)=temp12

    beta9(i,nl)=temp13
    beta9(i,ml)=temp14
    beta9(i,rl)=temp15
    beta9(i,sl)=temp16

enddo

!**************************
!
! FINDING TRACE OF MATRICES
!
!**************************

!tau1=tr[X1]

tau1=ZERO

do i=1,n

  tau1=tau1+X1(i,i)

enddo

!tau6=tr[beta8]
!tau7=tr[beta9]
!tau8=tr[beta11]
!tau9=tr[beta13]
!tau10=tr[beta14]

tau6=ZERO
tau7=ZERO
tau8=ZERO
tau9=ZERO
tau10=ZERO

do i=1,3*n
  
  tau6=tau6+beta8(i,i)
  tau7=tau7+beta9(i,i)
  tau8=tau8+beta11(i,i)
  tau9=tau9+beta13(i,i)
  tau10=tau10+beta14(i,i)
 
enddo

!*************************
!
! COMPUTING KINETIC ENERGY
!     (non-normalized)
!*************************

Tkl=(PI**((THREE*n)/TWO))&
*(det_tAkl**(-THREE/TWO))&
*(THREE*tau1*nu3+TWO*(tau6+tau7-tau8-tau9+tau10))

!write(*,*) 'Tkl',Tkl

!************************
!
! Kinetic Energy Gradient
!
!************************

!go to 1111

if (grad_k.or.grad_l) then

!beta18=X1I*beta33
!beta19=beta33*X1I
!beta20=beta5I*beta33
!beta21=beta3*X1I
!beta22=beta2*X1I

 beta18=ZERO
 beta19=ZERO
 beta20=ZERO
 beta21=ZERO
 beta22=ZERO
 beta24=ZERO
 beta25=ZERO
 beta27=ZERO
 beta28=ZERO

 do i=1,3*n

   temp1=ZERO
   temp2=ZERO
   temp3=ZERO
   temp4=ZERO

   temp5=ZERO
   temp6=ZERO
   temp7=ZERO
   temp8=ZERO

   temp9=ZERO
   temp10=ZERO
   temp11=ZERO
   temp12=ZERO

   temp13=ZERO
   temp14=ZERO
   temp15=ZERO
   temp16=ZERO

   do j=1,3*n

  
     temp1=temp1+X1I(i,j)*beta33(j,nk)
     temp2=temp2+X1I(i,j)*beta33(j,mk)
     temp3=temp3+X1I(i,j)*beta33(j,rk)
     temp4=temp4+X1I(i,j)*beta33(j,sk)

     temp5=temp5+beta5I(i,j)*beta33(j,nk)
     temp6=temp6+beta5I(i,j)*beta33(j,mk)
     temp7=temp7+beta5I(i,j)*beta33(j,rk)
     temp8=temp8+beta5I(i,j)*beta33(j,sk)

     temp9=temp9+beta5I(i,j)*beta3(j,nl)
     temp10=temp10+beta5I(i,j)*beta3(j,ml)
     temp11=temp11+beta5I(i,j)*beta3(j,rl)
     temp12=temp12+beta5I(i,j)*beta3(j,sl)

     temp13=temp13+beta4I(i,j)*beta6(j,nl)
     temp14=temp14+beta4I(i,j)*beta6(j,ml)
     temp15=temp15+beta4I(i,j)*beta6(j,rl)
     temp16=temp16+beta4I(i,j)*beta6(j,sl)

     beta19(i,j)=beta33(i,nk)*X1I(nk,j)&
                +beta33(i,mk)*X1I(mk,j)&
                +beta33(i,rk)*X1I(rk,j)&
                +beta33(i,sk)*X1I(sk,j)

     beta21(i,j)=beta3(i,nl)*X1I(nl,j)&
                +beta3(i,ml)*X1I(ml,j)&
                +beta3(i,rl)*X1I(rl,j)&
                +beta3(i,sl)*X1I(sl,j)

     beta22(i,j)=beta2(i,nl)*X1I(nl,j)&
                +beta2(i,ml)*X1I(ml,j)&
                +beta2(i,rl)*X1I(rl,j)&
                +beta2(i,sl)*X1I(sl,j)

     beta28(i,j)=beta33(i,nk)*beta5I(nk,j)&  
                +beta33(i,mk)*beta5I(mk,j)&  
                +beta33(i,rk)*beta5I(rk,j)&  
                +beta33(i,sk)*beta5I(sk,j)   

   enddo
  
   beta18(i,nk)=temp1
   beta18(i,mk)=temp2
   beta18(i,rk)=temp3
   beta18(i,sk)=temp4

   beta20(i,nk)=temp5
   beta20(i,mk)=temp6
   beta20(i,rk)=temp7
   beta20(i,sk)=temp8

   beta24(i,nl)=temp9
   beta24(i,ml)=temp10
   beta24(i,rl)=temp11
   beta24(i,sl)=temp12

   beta25(i,nl)=temp13
   beta25(i,ml)=temp14
   beta25(i,rl)=temp15
   beta25(i,sl)=temp16
   
   beta27(i,nk)=beta6(i,nl)*beta1(nl,nk)&
               +beta6(i,ml)*beta1(ml,nk)&
               +beta6(i,rl)*beta1(rl,nk)&
               +beta6(i,sl)*beta1(sl,nk)

   beta27(i,mk)=beta6(i,nl)*beta1(nl,mk)&
               +beta6(i,ml)*beta1(ml,mk)&
               +beta6(i,rl)*beta1(rl,mk)&
               +beta6(i,sl)*beta1(sl,mk)

   beta27(i,rk)=beta6(i,nl)*beta1(nl,rk)&
               +beta6(i,ml)*beta1(ml,rk)&
               +beta6(i,rl)*beta1(rl,rk)&
               +beta6(i,sl)*beta1(sl,rk)

   beta27(i,sk)=beta6(i,nl)*beta1(nl,sk)&
               +beta6(i,ml)*beta1(ml,sk)&
               +beta6(i,rl)*beta1(rl,sk)&
               +beta6(i,sl)*beta1(sl,sk)

 enddo

! write(*,*)'beta18'
! do i=1,3*n
!   write(*,9292)(beta18(i,j),j=1,3*n)
! enddo
 
! write(*,*)'beta19'
! do i=1,3*n
!   write(*,9292)(beta19(i,j),j=1,3*n)
! enddo

! write(*,*)'beta20'
! do i=1,3*n
!   write(*,9292)(beta20(i,j),j=1,3*n)
! enddo

 !beta23=beta22*beta1

 beta23=ZERO
 beta26=ZERO

 do i=1,3*n
   
   temp1=ZERO
   temp2=ZERO
   temp3=ZERO
   temp4=ZERO

   do j=1,3*n
     
     temp1=temp1+beta22(i,j)*beta1(j,nk)
     temp2=temp2+beta22(i,j)*beta1(j,mk)
     temp3=temp3+beta22(i,j)*beta1(j,rk)
     temp4=temp4+beta22(i,j)*beta1(j,sk)

   enddo
  
   beta23(i,nk)=temp1
   beta23(i,mk)=temp2
   beta23(i,rk)=temp3
   beta23(i,sk)=temp4
   
   beta26(i,nk)=beta25(i,nl)*beta1(nl,nk)&
               +beta25(i,ml)*beta1(ml,nk)&
               +beta25(i,rl)*beta1(rl,nk)&
               +beta25(i,sl)*beta1(sl,nk)

   beta26(i,mk)=beta25(i,nl)*beta1(nl,mk)&
               +beta25(i,ml)*beta1(ml,mk)&
               +beta25(i,rl)*beta1(rl,mk)&
               +beta25(i,sl)*beta1(sl,mk)

   beta26(i,rk)=beta25(i,nl)*beta1(nl,rk)&
               +beta25(i,ml)*beta1(ml,rk)&
               +beta25(i,rl)*beta1(rl,rk)&
               +beta25(i,sl)*beta1(sl,rk)

   beta26(i,sk)=beta25(i,nl)*beta1(nl,sk)&
               +beta25(i,ml)*beta1(ml,sk)&
               +beta25(i,rl)*beta1(rl,sk)&
               +beta25(i,sl)*beta1(sl,sk)

 enddo

 
 !alpha5=beta5*inv_tAkl
 !alpha6=X1*inv_tAkl

 alpha5=ZERO
 alpha6=ZERO

 do i=1,n
   do j=1,n
     temp1=ZERO
     temp2=ZERO
     do k=1,n

       temp1=temp1+beta5(i,k)*inv_tAkl(k,j)
       temp2=temp2+X1(i,k)*inv_tAkl(k,j)

     enddo

     alpha5(i,j)=temp1
     alpha6(i,j)=temp2

   enddo
 enddo

 !alpha7=beta8*inv_tAklI
 !alpha8=beta18*inv_tAklI
 !alpha9=beta19*inv_tAklI
 !alpha10=beta20*inv_tAklI
 !alpha11=beta21*inv_tAklI
 !alpha20=beta14*inv_tAklI

 alpha7=ZERO
 alpha8=ZERO
 alpha9=ZERO
 alpha10=ZERO
 alpha11=ZERO
 alpha12=ZERO
 alpha13=ZERO
 alpha14=ZERO
 alpha15=ZERO
 alpha16=ZERO
 alpha17=ZERO
 alpha18=ZERO 
 alpha19=ZERO 
 alpha20=ZERO

 do i=1,3*n
   do j=1,3*n
      
     temp1=ZERO
     temp2=ZERO
     temp3=ZERO
     temp4=ZERO

     do k=1,3*n
        
       temp1=temp1+beta14(i,k)*inv_tAklI(k,j)
       temp2=temp2+beta19(i,k)*inv_tAklI(k,j)
       temp3=temp3+beta21(i,k)*inv_tAklI(k,j)
       temp4=temp4+beta28(i,k)*inv_tAklI(k,j)

     enddo
      
     alpha20(i,j)=temp1

     alpha9(i,j)=temp2

     alpha12(i,j)=temp3

     alpha19(i,j)=temp4

     alpha7(i,j)=beta8(i,nl)*inv_tAklI(nl,j)&
                +beta8(i,ml)*inv_tAklI(ml,j)&
                +beta8(i,rl)*inv_tAklI(rl,j)&
                +beta8(i,sl)*inv_tAklI(sl,j)
     
     alpha8(i,j)=beta18(i,nk)*inv_tAklI(nk,j)&
                +beta18(i,mk)*inv_tAklI(mk,j)&
                +beta18(i,rk)*inv_tAklI(rk,j)&
                +beta18(i,sk)*inv_tAklI(sk,j)

     alpha10(i,j)=beta20(i,nk)*inv_tAklI(nk,j)&
                 +beta20(i,mk)*inv_tAklI(mk,j)&
                 +beta20(i,rk)*inv_tAklI(rk,j)&
                 +beta20(i,sk)*inv_tAklI(sk,j)
    
     alpha11(i,j)=beta9(i,nl)*inv_tAklI(nl,j)&
                 +beta9(i,ml)*inv_tAklI(ml,j)&
                 +beta9(i,rl)*inv_tAklI(rl,j)&
                 +beta9(i,sl)*inv_tAklI(sl,j)

     alpha13(i,j)=beta23(i,nk)*inv_tAklI(nk,j)&
                 +beta23(i,mk)*inv_tAklI(mk,j)&
                 +beta23(i,rk)*inv_tAklI(rk,j)&
                 +beta23(i,sk)*inv_tAklI(sk,j)

     alpha14(i,j)=beta24(i,nl)*inv_tAklI(nl,j)&
                 +beta24(i,ml)*inv_tAklI(ml,j)&
                 +beta24(i,rl)*inv_tAklI(rl,j)&
                 +beta24(i,sl)*inv_tAklI(sl,j)

     alpha15(i,j)=beta11(i,nl)*inv_tAklI(nl,j)&
                 +beta11(i,ml)*inv_tAklI(ml,j)&
                 +beta11(i,rl)*inv_tAklI(rl,j)&
                 +beta11(i,sl)*inv_tAklI(sl,j)

     alpha16(i,j)=beta26(i,nk)*inv_tAklI(nk,j)&
                 +beta26(i,mk)*inv_tAklI(mk,j)&
                 +beta26(i,rk)*inv_tAklI(rk,j)&
                 +beta26(i,sk)*inv_tAklI(sk,j)

     alpha17(i,j)=beta27(i,nk)*inv_tAklI(nk,j)&
                 +beta27(i,mk)*inv_tAklI(mk,j)&
                 +beta27(i,rk)*inv_tAklI(rk,j)&
                 +beta27(i,sk)*inv_tAklI(sk,j)

     alpha18(i,j)=beta13(i,nl)*inv_tAklI(nl,j)& 
                 +beta13(i,ml)*inv_tAklI(ml,j)& 
                 +beta13(i,rl)*inv_tAklI(rl,j)& 
                 +beta13(i,sl)*inv_tAklI(sl,j) 

   enddo
 enddo

endif

if (grad_k) then

 !result3k=alpha5*Lk
 !result4k=alpha6*Lk


  result3k=ZERO
  result4k=ZERO
  indx=ZERO
  do j=1,n
    do i=j,n
      indx=indx+1
      temp1=ZERO
      temp2=ZERO
      do k=1,n

        temp1=temp1+(alpha5(i,k)+alpha5(k,i))*Lk(k,j)
        temp2=temp2+(alpha6(i,k)+alpha6(k,i))*Lk(k,j)

      enddo

      result3k(indx)=temp1
      result4k(indx)=temp2

    enddo
  enddo

  result5k=ZERO
  indx=ZERO
  do j=1,3*n
    do i=j,3*n

      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      temp7=ZERO
      temp8=ZERO
      temp9=ZERO
      temp10=ZERO
      temp11=ZERO
      temp12=ZERO
      temp13=ZERO
      temp14=ZERO

      do k=1,3*n
          
        temp1=temp1+(alpha7(i,k)+alpha7(k,i))*LkI(k,j)
        temp2=temp2+(alpha8(i,k)+alpha8(k,i))*LkI(k,j)
        temp3=temp3+(alpha9(i,k)+alpha9(k,i))*LkI(k,j)
        temp4=temp4+(alpha10(i,k)+alpha10(k,i))*LkI(k,j)
        temp5=temp5+(alpha11(i,k)+alpha11(k,i))*LkI(k,j)
        temp6=temp6+(alpha12(i,k)+alpha12(k,i))*LkI(k,j) 
        temp7=temp7+(alpha13(i,k)+alpha13(k,i))*LkI(k,j)       
        temp8=temp8+(alpha14(i,k)+alpha14(k,i))*LkI(k,j)
        temp9=temp9+(alpha15(i,k)+alpha15(k,i))*LkI(k,j)
        temp10=temp10+(alpha16(i,k)+alpha16(k,i))*LkI(k,j)
        temp11=temp11+(alpha17(i,k)+alpha17(k,i))*LkI(k,j)
        temp12=temp12+(alpha18(i,k)+alpha18(k,i))*LkI(k,j)
        temp13=temp13+(alpha19(i,k)+alpha19(k,i))*LkI(k,j) 
        temp14=temp14+(alpha20(i,k)+alpha20(k,i))*LkI(k,j)
 
      enddo
      indx=indx+1

      result5k(indx)=-THREE*tau1*result2k(indx)&
+TWO*(-temp1-temp2-temp3+temp4-temp5-temp6-temp7&
+temp8+temp9+temp10-temp11+temp12+temp13-temp14)

    enddo
  enddo

  do i=1,n*(n+1)/2
    temp1=ZERO
    do k=1,3*n*(3*n+1)/2
      temp1=temp1+result5k(k)*trans(k,i)
    enddo

    Dk(i)=Dk(i)+(PI**(THREEHALF*n))*(det_tAkl**(-THREEHALF))&
*(-THREEHALF*result1k(i)*(THREE*nu3*tau1+TWO*(tau6+tau7-tau8-tau9+tau10))&
+THREE*nu3*(result3k(i)-result4k(i))+temp1)

  enddo
endif

!write(*,*) 'Dk',Dk

!alpha21=beta4*M

if (grad_l) then

  alpha21=ZERO

  do i=1,n
    do j=1,n

      temp1=ZERO
      do k=1,n

        temp1=temp1+beta4(i,k)*Glob_MassMatrix(k,j)

      enddo

      alpha21(i,j)=temp1

    enddo
  enddo
 
!beta31=beta33*beta4I

  beta29=ZERO
  beta30=ZERO

  do i=1,3*n
    do j=1,3*n

      beta29(i,j)=beta33(i,nk)*beta4I(nk,j)&
                 +beta33(i,mk)*beta4I(mk,j)&
                 +beta33(i,rk)*beta4I(rk,j)&
                 +beta33(i,sk)*beta4I(sk,j)
     
      beta30(i,j)=beta3(i,nl)*beta4I(nl,j)&
                 +beta3(i,ml)*beta4I(ml,j)&
                 +beta3(i,rl)*beta4I(rl,j)&
                 +beta3(i,sl)*beta4I(sl,j) 
    enddo
  enddo

!alpha22=beta29*MI
!alpha24=beta30*MI

  alpha22=ZERO
  alpha23=ZERO
  alpha24=ZERO

  do i=1,3*n
    do j=1,3*n

      temp1=ZERO
      temp2=ZERO

      do k=1,3*n

        temp1=temp1+beta29(i,k)*MI(k,j)
        temp2=temp2+beta30(i,k)*MI(k,j)

      enddo

      alpha22(i,j)=temp1

      alpha23(i,j)=temp2

      alpha24(i,j)=beta33(i,nk)*MI(nk,j)&   
                  +beta33(i,mk)*MI(mk,j)&  
                  +beta33(i,rk)*MI(rk,j)&
                  +beta33(i,sk)*MI(sk,j)

    enddo
  enddo

  Z3=ZERO
  Z4=ZERO         
  Z5=ZERO
  Z6=ZERO
  Z7=ZERO
  Z8=ZERO
  Z9=ZERO
  Z10=ZERO
  Z11=ZERO
  Z12=ZERO
  Z13=ZERO
  Z14=ZERO
  Z15=ZERO
  Z16=ZERO
  do i=1,3*n      
    do j=1,3*n    
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      temp7=ZERO
      temp8=ZERO
      temp9=ZERO
      temp10=ZERO
      temp11=ZERO
      temp12=ZERO
      temp13=ZERO
      temp14=ZERO
      do k=1,3*n  
        temp1=temp1+pP(i,k)*alpha7(k,j)
        temp2=temp2+pP(i,k)*alpha8(k,j)
        temp3=temp3+pP(i,k)*alpha9(k,j)
        temp4=temp4+pP(i,k)*alpha22(k,j)
        temp5=temp5+pP(i,k)*alpha11(k,j)
        temp6=temp6+pP(i,k)*alpha12(k,j)
        temp7=temp7+pP(i,k)*alpha13(k,j)
        temp8=temp8+pP(i,k)*alpha23(k,j)
        temp9=temp9+pP(i,k)*alpha15(k,j)
        temp10=temp10+pP(i,k)*alpha16(k,j)
        temp11=temp11+pP(i,k)*alpha24(k,j)
        temp12=temp12+pP(i,k)*alpha18(k,j)
        temp13=temp13+pP(i,k)*alpha19(k,j)
        temp14=temp14+pP(i,k)*alpha20(k,j)
      enddo  
      Z3(i,j)=temp1
      Z4(i,j)=temp2
      Z5(i,j)=temp3
      Z6(i,j)=temp4
      Z7(i,j)=temp5
      Z8(i,j)=temp6
      Z9(i,j)=temp7
      Z10(i,j)=temp8
      Z11(i,j)=temp9
      Z12(i,j)=temp10
      Z13(i,j)=temp11
      Z14(i,j)=temp12
      Z15(i,j)=temp13
      Z16(i,j)=temp14
   enddo         
  enddo           
  do i=1,3*n      
    do j=1,3*n    
      temp1=ZERO
      temp2=ZERO  
      temp3=ZERO  
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      temp7=ZERO 
      temp8=ZERO
      temp9=ZERO
      temp10=ZERO
      temp11=ZERO
      temp12=ZERO
      temp13=ZERO
      temp14=ZERO
      do k=1,3*n  
        temp1=temp1+Z3(i,k)*pP(j,k)
        temp2=temp2+Z4(i,k)*pP(j,k)
        temp3=temp3+Z5(i,k)*pP(j,k)
        temp4=temp4+Z6(i,k)*pP(j,k)
        temp5=temp5+Z7(i,k)*pP(j,k)
        temp6=temp6+Z8(i,k)*pP(j,k)
        temp7=temp7+Z9(i,k)*pP(j,k)
        temp8=temp8+Z10(i,k)*pP(j,k)
        temp9=temp9+Z11(i,k)*pP(j,k)
        temp10=temp10+Z12(i,k)*pP(j,k)
        temp11=temp11+Z13(i,k)*pP(j,k)
        temp12=temp12+Z14(i,k)*pP(j,k)
        temp13=temp13+Z15(i,k)*pP(j,k)
        temp14=temp14+Z16(i,k)*pP(j,k)
      enddo  
      alpha7p(i,j)=temp1
      alpha8p(i,j)=temp2
      alpha9p(i,j)=temp3
      alpha22p(i,j)=temp4
      alpha11p(i,j)=temp5
      alpha12p(i,j)=temp6
      alpha13p(i,j)=temp7
      alpha23p(i,j)=temp8
      alpha15p(i,j)=temp9
      alpha16p(i,j)=temp10
      alpha24p(i,j)=temp11
      alpha18p(i,j)=temp12
      alpha19p(i,j)=temp13
      alpha20p(i,j)=temp14
    enddo
  enddo

  Z1=ZERO
  Z2=ZERO
  do i=1,n
    do j=1,n
      temp1=ZERO
      temp2=ZERO
      do k=1,n
        temp1=temp1+P(i,k)*alpha6(k,j)
        temp2=temp2+P(i,k)*alpha21(k,j)
      enddo  
      Z1(i,j)=temp1
      Z2(i,j)=temp2
    enddo
  enddo  
  do i=1,n
    do j=1,n
      temp1=ZERO
      temp2=ZERO
      do k=1,n
        temp1=temp1+Z1(i,k)*P(j,k)
        temp2=temp2+Z2(i,k)*P(j,k)
      enddo  
      alpha6p(i,j)=temp1
      alpha21p(i,j)=temp2
    enddo
  enddo  

!result3l=alpha21*pLl
!result4l=alpha6*pLl

  result3l=ZERO
  result4l=ZERO
  indx=ZERO

  do j=1,n
    do i=j,n

      indx=indx+1
      temp1=ZERO
      temp2=ZERO

      do k=j,n

        temp1=temp1+(alpha21p(i,k)+alpha21p(k,i))*Ll(k,j)
        temp2=temp2+(alpha6p(i,k)+alpha6p(k,i))*Ll(k,j)

      enddo

      result3l(indx)=temp1
      result4l(indx)=temp2

    enddo
  enddo

  result5l=ZERO
  indx=ZERO
  do j=1,3*n
    do i=j,3*n

      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      temp7=ZERO
      temp8=ZERO
      temp9=ZERO 
      temp10=ZERO
      temp11=ZERO
      temp12=ZERO
      temp13=ZERO
      temp14=ZERO

      do k=1,3*n  
          
        temp1=temp1+(alpha7p(i,k)+alpha7p(k,i))*LlI(k,j)
        temp2=temp2+(alpha8p(i,k)+alpha8p(k,i))*LlI(k,j)
        temp3=temp3+(alpha9p(i,k)+alpha9p(k,i))*LlI(k,j)
        temp4=temp4+(alpha22p(i,k)+alpha22p(k,i))*LlI(k,j)
        temp5=temp5+(alpha11p(i,k)+alpha11p(k,i))*LlI(k,j)
        temp6=temp6+(alpha12p(i,k)+alpha12p(k,i))*LlI(k,j)
        temp7=temp7+(alpha13p(i,k)+alpha13p(k,i))*LlI(k,j)
        temp8=temp8+(alpha23p(i,k)+alpha23p(k,i))*LlI(k,j)
        temp9=temp9+(alpha15p(i,k)+alpha15p(k,i))*LlI(k,j)
        temp10=temp10+(alpha16p(i,k)+alpha16p(k,i))*LlI(k,j)
        temp11=temp11+(alpha24p(i,k)+alpha24p(k,i))*LlI(k,j)
        temp12=temp12+(alpha18p(i,k)+alpha18p(k,i))*LlI(k,j)
        temp13=temp13+(alpha19p(i,k)+alpha19p(k,i))*LlI(k,j)
        temp14=temp14+(alpha20p(i,k)+alpha20p(k,i))*LlI(k,j)
 
      enddo
      indx=indx+1

      result5l(indx)=-THREE*tau1*result2l(indx)&
+TWO*(-temp1-temp2-temp3+temp4-temp5-temp6-temp7&
+temp8+temp9+temp10-temp11+temp12+temp13-temp14)

    enddo
  enddo  

  do i=1,n*(n+1)/2
    temp1=ZERO
    do k=1,3*n*(3*n+1)/2
      temp1=temp1+result5l(k)*trans(k,i)
    enddo

    Dl(i)=Dl(i)+(PI**(THREEHALF*n))*(det_tAkl**(-THREEHALF))&
*(-THREEHALF*result1l(i)*(THREE*nu3*tau1+TWO*(tau6+tau7-tau8-tau9+tau10))&
+THREE*nu3*(result3l(i)-result4l(i))+temp1)

  enddo
endif

!write(*,*) 'Dl',Dl

!1111 continue

!go to 1111

!****************************
!
! EVALUATING POTENTIAL ENERGY 
!      (non-normalized)
!****************************

!*************************
!
! CONSTRUCTION OF J MATRIX
!
!*************************

Vkl=ZERO
do i=1,n !POTENTIAL LOOP
  do j=i,n !POTENTIAL LOOP

    Jm=ZERO
    Jm(i,j)=-ONE
    Jm(j,i)=-ONE
    Jm(i,i)=ONE
    Jm(j,j)=ONE
    if (i.eq.j) Jm(i,i) = ONE

!****************************
!
! DOING MATRIX MULTIPLICATION
!
!****************************

!betaJm=inv_tAkl*Jm

betaJm=ZERO

do r=1,n
  
   betaJm(r,i)=inv_tAkl(r,i)*Jm(i,i)+inv_tAkl(r,j)*Jm(j,i)
   betaJm(r,j)=inv_tAkl(r,i)*Jm(i,j)+inv_tAkl(r,j)*Jm(j,j)
 
  if (i.eq.j) then 

    betaJm(r,i)=inv_tAkl(r,i)*Jm(i,i)

  endif

enddo   

! Preforming Kronecker Product of betaJm with I3

betaJmI=ZERO

do r=1,n
  do t=1,n

    betaJmI(3*(r-1)+1,3*(t-1)+1)=betaJm(r,t)
    betaJmI(3*(r-1)+2,3*(t-1)+2)=betaJm(r,t)
    betaJmI(3*(r-1)+3,3*(t-1)+3)=betaJm(r,t)

  enddo
enddo

!beta16=beta2*betaJmI
!beta17=beta3*betaJmI

beta15=ZERO
beta16=ZERO
beta17=ZERO
beta18=ZERO

do r=1,3*n
  do t=1,3*n

    if (grad_k.or.grad_l) then

      beta15(r,t)=beta1(r,nk)*betaJmI(nk,t)+beta1(r,mk)*betaJmI(mk,t)&
                 +beta1(r,rk)*betaJmI(rk,t)+beta1(r,sk)*betaJmI(sk,t)

      beta16(r,t)=beta2(r,nl)*betaJmI(nl,t)+beta2(r,ml)*betaJmI(ml,t)&
                 +beta2(r,rl)*betaJmI(rl,t)+beta2(r,sl)*betaJmI(sl,t)

    endif 

  beta17(r,t)=beta3(r,nl)*betaJmI(nl,t)+beta3(r,ml)*betaJmI(ml,t)&
             +beta3(r,rl)*betaJmI(rl,t)+beta3(r,sl)*betaJmI(sl,t)

  beta18(r,t)=beta33(r,nk)*betaJmI(nk,t)+beta33(r,mk)*betaJmI(mk,t)&
             +beta33(r,rk)*betaJmI(rk,t)+beta33(r,sk)*betaJmI(sk,t)
 
  enddo	
enddo

!**************************
!
! FINDING TRACE OF MATRICES
!
!**************************

!Xi=tr[betaJm]

Xi=ZERO

do r=1,n
  
  Xi=Xi+betaJm(r,r)

enddo


nu6=ZERO
nu7=ZERO

do r=1,3*n
  
  nu6=nu6+beta17(r,r)
  nu7=nu7+beta18(r,r)

enddo

!***************************
!
! COMPUTING POTENTIAL ENERGY
!
!***************************

Rkl=TWO*(PI**((THREE*n-ONE)/TWO))&
*(Xi**(-ONE/TWO))&
*(det_tAkl**(-THREEHALF))&  
*(ONEHALF*nu3& 
-(ONE/SIX)*(Xi**(-ONE))*(nu6+nu7))

if (i.ne.j) Vkl=Vkl+Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*Rkl

if (i.eq.j) Vkl=Vkl+Glob_PseudoCharge0*Glob_PseudoCharge(i)*Rkl

!**************************
!
! Potential Energy Gradient
!
!**************************

if (grad_k.or. grad_l) then

!alpha25=betaJm*inv_tAkl

  alpha25=ZERO

  do r=1,n
    do t=1,n

      temp1=ZERO

      do k=1,n

        temp1=temp1+betaJm(r,k)*inv_tAkl(k,t)

      enddo

      alpha25(r,t)=temp1

    enddo
  enddo

!write(*,*)'alpha25'
!do r=1,n
!  write(*,9292)(alpha25(r,t),t=1,n)
!enddo

!beta34=beta16*beta1

  beta34=ZERO
  beta35=ZERO
  beta37=ZERO
  beta38=ZERO 

  do r=1,3*n
    
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    temp4=ZERO
    temp5=ZERO
    temp6=ZERO
    temp7=ZERO
    temp8=ZERO
    temp9=ZERO
    temp10=ZERO
    temp11=ZERO
    temp12=ZERO
    temp13=ZERO
    temp14=ZERO
    temp15=ZERO
    temp16=ZERO

    do t=1,3*n

        temp1=temp1+beta16(r,t)*beta1(t,nk)
        temp2=temp2+beta16(r,t)*beta1(t,mk)
        temp3=temp3+beta16(r,t)*beta1(t,rk)
        temp4=temp4+beta16(r,t)*beta1(t,sk)

        temp5=temp5+betaJmI(r,t)*beta33(t,nk)
        temp6=temp6+betaJmI(r,t)*beta33(t,mk)
        temp7=temp7+betaJmI(r,t)*beta33(t,rk)
        temp8=temp8+betaJmI(r,t)*beta33(t,sk)

        temp9=temp9+beta15(r,t)*beta2(t,nl)
        temp10=temp10+beta15(r,t)*beta2(t,ml)
        temp11=temp11+beta15(r,t)*beta2(t,rl)
        temp12=temp12+beta15(r,t)*beta2(t,sl)

        temp13=temp13+betaJmI(r,t)*beta3(t,nl)
        temp14=temp14+betaJmI(r,t)*beta3(t,ml)
        temp15=temp15+betaJmI(r,t)*beta3(t,rl)
        temp16=temp16+betaJmI(r,t)*beta3(t,sl)

    enddo
 
    beta34(r,nk)=temp1
    beta34(r,mk)=temp2
    beta34(r,rk)=temp3
    beta34(r,sk)=temp4

    beta35(r,nk)=temp5
    beta35(r,mk)=temp6
    beta35(r,rk)=temp7
    beta35(r,sk)=temp8

    beta37(r,nl)=temp9
    beta37(r,ml)=temp10
    beta37(r,rl)=temp11
    beta37(r,sl)=temp12

    beta38(r,nl)=temp13
    beta38(r,ml)=temp14
    beta38(r,rl)=temp15
    beta38(r,sl)=temp16

  enddo

!write(*,*)'beta34'
!do r=1,3*n
!  write(*,9292)(beta34(r,t),t=1,3*n)
!enddo

!alpha26=beta17*inv_tAklI
!alpha27=beta34*inv_tAklI

  alpha26=ZERO
  alpha27=ZERO
  alpha28=ZERO
  alpha29=ZERO
  alpha30=ZERO
  alpha31=ZERO

  do r=1,3*n
   do t=1,3*n

      temp1=ZERO
      temp2=ZERO

      do k=1,3*n

        temp1=temp1+beta17(r,k)*inv_tAklI(k,t) 
        temp2=temp2+beta18(r,k)*inv_tAklI(k,t)
      
      enddo

      alpha26(r,t)=temp1 
      alpha29(r,t)=temp2
 
      alpha27(r,t)=beta34(r,nk)*inv_tAklI(nk,t)&
                  +beta34(r,mk)*inv_tAklI(mk,t)&
                  +beta34(r,rk)*inv_tAklI(rk,t)&
                  +beta34(r,sk)*inv_tAklI(sk,t)

      alpha28(r,t)=beta35(r,nk)*inv_tAklI(nk,t)&
                  +beta35(r,mk)*inv_tAklI(mk,t)&
                  +beta35(r,rk)*inv_tAklI(rk,t)&
                  +beta35(r,sk)*inv_tAklI(sk,t) 

      alpha30(r,t)=beta37(r,nl)*inv_tAklI(nl,t)&
                  +beta37(r,ml)*inv_tAklI(ml,t)&
                  +beta37(r,rl)*inv_tAklI(rl,t)&
                  +beta37(r,sl)*inv_tAklI(sl,t)

      alpha31(r,t)=beta38(r,nl)*inv_tAklI(nl,t)&
                  +beta38(r,ml)*inv_tAklI(ml,t)&
                  +beta38(r,rl)*inv_tAklI(rl,t)&
                  +beta38(r,sl)*inv_tAklI(sl,t)

    enddo
  enddo
endif

!write(*,*)'alpha26'
!do r=1,3*n
!  write(*,9292)(alpha26(r,t),t=1,3*n)
!enddo

!write(*,*)'alpha27'
!do r=1,3*n
!  write(*,9292)(alpha27(r,t),t=1,3*n)
!enddo

if (grad_k) then

  result6k=ZERO
  indx=ZERO   
  do t=1,3*n 
    do r=t,3*n

      indx=indx+1

      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO

      do k=1,3*n

        temp1=temp1+(alpha26(r,k)+alpha26(k,r))*LkI(k,t)
        temp2=temp2+(alpha27(r,k)+alpha27(k,r))*LkI(k,t)
        temp3=temp3+(alpha28(r,k)+alpha28(k,r))*LkI(k,t)
        temp4=temp4+(alpha29(r,k)+alpha29(k,r))*LkI(k,t)
        temp5=temp5+(alpha30(r,k)+alpha30(k,r))*LkI(k,t)
        temp6=temp6+(alpha31(r,k)+alpha31(k,r))*LkI(k,t)

      enddo
      result6k(indx)=temp1+temp2+temp3+temp4+temp5+temp6
 
    enddo  
  enddo

!write(*,*)'result6k',result6k


  resultJk=ZERO  
  indx=ZERO
  do t=1,n  
    do r=t,n  
      indx=indx+1
      temp1=ZERO  
      do k=1,n  
        temp1=temp1+(alpha25(r,k)+alpha25(k,r))*Lk(k,t)
      enddo
      resultJk(indx)=temp1
    enddo
  enddo

!write(*,*)'resultJk',resultJk

  pk=ZERO

  do r=1,n*(n+1)/2

    temp1=ZERO
    temp2=ZERO

    do k=1,3*n*(3*n+1)/2

      temp1=temp1+result6k(k)*trans(k,r)
      temp2=temp2+result2k(k)*trans(k,r)

    enddo

    pk(r)=pk(r)+TWO*(PI**((THREE*n-ONE)/TWO))&
*(det_tAkl**(-THREEHALF))&
*(Xi**(-ONEHALF))&
*((-THREEHALF*result1k(r)+ONEHALF*(Xi**(-ONE))*resultJk(r))*&
(ONEHALF*nu3-(ONE/SIX)*(Xi**(-ONE))*(nu6+nu7))&
-ONEHALF*temp2&
+(ONE/SIX)*(-(Xi**(-TWO))*resultJk(r)*(nu6+nu7)&
+(Xi**(-ONE))*temp1))

  enddo

  do r=1,n*(n+1)/2

!    if (i.ne.j) write(*,*)'i,j,pk',i,j,Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*pk(r)
!    if (i.eq.j) write(*,*)'i,j,pk',i,j,Glob_PseudoCharge0*Glob_PseudoCharge(i)*pk(r)
 
    if (i.ne.j) Dk(r)=Dk(r)+Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*pk(r)

    if (i.eq.j) Dk(r)=Dk(r)+Glob_PseudoCharge0*Glob_PseudoCharge(i)*pk(r)

  enddo
endif

if (grad_l) then 

  Z1=ZERO
  do r=1,n
    do t=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+P(r,k)*alpha25(k,t)
      enddo
      Z1(r,t)=temp1
    enddo
  enddo
  do r=1,n
    do t=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+Z1(r,k)*P(t,k)
      enddo
      alpha25p(r,t)=temp1
    enddo
  enddo

  resultJl=ZERO
  indx=ZERO
  do t=1,n
    do r=t,n
      indx=indx+1  
      temp1=ZERO
      do k=1,n    

        temp1=temp1+(alpha25p(r,k)+alpha25p(k,r))*Ll(k,t)

      enddo

      resultJl(indx)=temp1

    enddo
  enddo

  Z3=ZERO
  Z4=ZERO
  Z5=ZERO
  Z6=ZERO
  Z7=ZERO
  Z8=ZERO
  do r=1,3*n
    do t=1,3*n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      do k=1,3*n
        temp1=temp1+pP(r,k)*alpha26(k,t)
        temp2=temp2+pP(r,k)*alpha27(k,t)
        temp3=temp3+pP(r,k)*alpha28(k,t)
        temp4=temp4+pP(r,k)*alpha29(k,t)
        temp5=temp5+pP(r,k)*alpha30(k,t)
        temp6=temp6+pP(r,k)*alpha31(k,t)
      enddo
      Z3(r,t)=temp1
      Z4(r,t)=temp2
      Z5(r,t)=temp3
      Z6(r,t)=temp4
      Z7(r,t)=temp5
      Z8(r,t)=temp6
    enddo
  enddo
  do r=1,3*n
    do t=1,3*n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      do k=1,3*n
        temp1=temp1+Z3(r,k)*pP(t,k)
        temp2=temp2+Z4(r,k)*pP(t,k)
        temp3=temp3+Z5(r,k)*pP(t,k)
        temp4=temp4+Z6(r,k)*pP(t,k)
        temp5=temp5+Z7(r,k)*pP(t,k)
        temp6=temp6+Z8(r,k)*pp(t,k)
      enddo
      alpha26p(r,t)=temp1
      alpha27p(r,t)=temp2
      alpha28p(r,t)=temp3
      alpha29p(r,t)=temp4
      alpha30p(r,t)=temp5
      alpha31p(r,t)=temp6
    enddo
  enddo

  result6l=ZERO
  indx=ZERO
  do t=1,3*n
    do r=t,3*n
      indx=indx+1

      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      do k=1,3*n
        temp1=temp1+(alpha26p(r,k)+alpha26p(k,r))*LlI(k,t)
        temp2=temp2+(alpha27p(r,k)+alpha27p(k,r))*LlI(k,t)
        temp3=temp3+(alpha28p(r,k)+alpha28p(k,r))*LlI(k,t)
        temp4=temp4+(alpha29p(r,k)+alpha29p(k,r))*LlI(k,t)
        temp5=temp5+(alpha30p(r,k)+alpha30p(k,r))*LlI(k,t)
        temp6=temp6+(alpha31p(r,k)+alpha31p(k,r))*LlI(k,t)
      enddo
      result6l(indx)=temp1+temp2+temp3+temp4+temp5+temp6

    enddo
  enddo

  pl=ZERO

  do r=1,n*(n+1)/2 

    temp1=ZERO
    temp2=ZERO

    do k=1,3*n*(3*n+1)/2

      temp1=temp1+result6l(k)*trans(k,r)
      temp2=temp2+result2l(k)*trans(k,r)

    enddo

    pl(r)=pl(r)+TWO*(PI**((THREE*n-ONE)/TWO))&
*(det_tAkl**(-THREEHALF))&
*(Xi**(-ONEHALF))&
*((-THREEHALF*result1l(r)+ONEHALF*(Xi**(-ONE))*resultJl(r))*&
(ONEHALF*nu3-(ONE/SIX)*(Xi**(-ONE))*(nu6+nu7))& 
-ONEHALF*temp2&
+(ONE/SIX)*(-(Xi**(-TWO))*resultJl(r)*(nu6+nu7)&
+(Xi**(-ONE))*temp1))
 
  enddo

  do r=1,n*(n+1)/2
!    if (i.ne.j) write(*,*)'i,j,pl',i,j,Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*pl(r)
!    if (i.eq.j) write(*,*)'i,j,pl',i,j,Glob_PseudoCharge0*Glob_PseudoCharge(i)*pl(r)

    if (i.ne.j) Dl(r)=Dl(r)+Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*pl(r)

    if (i.eq.j) Dl(r)=Dl(r)+Glob_PseudoCharge0*Glob_PseudoCharge(i)*pl(r)
  enddo 
endif

  enddo !POTENTIAL LOOP
enddo  !POTENTIAL LOOP

!1111 continue

!write(*,*)'Vkl',Vkl

!write(*,*)'Dk',Dk
!write(*,*)'Dl',Dl

Hkl=Tkl+Vkl

!write(*,*)'Hkl',Hkl

end subroutine

subroutine MatrixElementsL1ForExpcVals(m_k, m_l, vechLk, vechLl, Pbra, Pket, Hkl, Skl, Tkl, Vkl, &
               rm2kl, rmkl, rkl, r2kl, deltarkl, MVkl, Darwinkl, OOkl)
!This subroutine computes symmetry adapted matrix elements 
!with two real L=1 correlated Gaussians. These matrix elements
!are used in calculations of expectation values.
!Symmetry adaption is applied to the bra and ket using permutation matrices Pbra and Pket
!
!Input:     
!   m_k, m_l :: integers that determine which z-component is in the
!       premultiplier of the Gaussian
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters. 
!   Pbra :: The symmetry permutation matrix of size n x n that is applied to bra
!   Pket :: The symmetry permutation matrix of size n x n that is applied to ket
!Output (all matrix elements are computed with normilized functions):
!   Hkl	 ::	Hamiltonian
!   Skl	 ::	Overlap
!   Tkl  :: Kinetic energy
!   Vkl  :: Potential energy
!   rm2kl :: 1/r_i^2, 1/r_{ij}^2
!   rmkl :: 1/r_i, 1/r_{ij}
!   rkl  :: r_i, r_{ij}
!   r2kl :: r_i^2, r_{ij}^2
!deltarkl:: delta(r_i), delta(r_{ij})
!   MVkl :: Mass-velocity correction (without the factor of alpha**2)
!Darwinkl:: Darwin correction (without the factor of alpha**2)
!   OOkl :: Orbit-Orbit correction (without the factor of alpha**2)

!Arguments
integer,intent(in)       :: m_k,m_l
real(dprec),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
real(dprec),intent(in)   :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
real(dprec),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,Darwinkl,OOkl
real(dprec),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
real(dprec),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
real(dprec),intent(out)  :: deltarkl(Glob_n,Glob_n)

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
real(dprec)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn)
real(dprec)       inv_tAkltAl(nn,nn),inv_tAkltAlM(nn,nn)
real(dprec)       inv_tAklAk(nn,nn),inv_tAklAkM(nn,nn)
real(dprec)       eta1(nn,nn),sqrt_eta1(nn,nn),eta2(nn,nn)
real(dprec)       W1(nn,nn),W2(nn,nn),W3(nn,nn),W4(nn,nn)
real(dprec)       inv_tAkltvl(nn),tvkinv_tAkl(nn),tvkinv_tAkltAlM(nn)
real(dprec)       u1(nn),u2(nn),u3(nn)
real(dprec)       temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
real(dprec)       det_Lk, det_Ll, det_tAkl
real(dprec)       tau1,tau2,tau3
integer           i,j,k,q,t,indx
real(dprec)       Mass_For_Darwin(0:nn)

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
	do k=1,j
	  temp1=temp1+Lk(i,k)*Lk(j,k)
	enddo 
	tAk(i,j)=temp1
	tAk(j,i)=temp1
    temp1=ZERO
	do k=1,j
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

!Evaluating overlap
temp1=abs(det_Ll*det_Lk)/det_tAkl
Skl=Glob_2raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))

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

!Evaluating eta1(i,j)=Tr[inv_tAkl*Jij], sqrt_eta1(i,j),
!           eta2(i,j)=tvk'*inv_tAkl*Jij*inv_tAkl*tvl,
!Vkl, (1/r_{ij}^2)_kl, (1/r_{ij})_kl, (r_{ij})_kl, (r_{ij}^2)_kl
!and delta(r_{ij})_kl
Vkl=ZERO
temp5=Skl*TWO
temp1=temp5/SQRTPI
temp8=Skl/(PI*SQRTPI)
do i=1,n
  temp2=inv_tAkl(i,i)
  temp3=sqrt(temp2)
  eta1(i,i)=temp2
  sqrt_eta1(i,i)=temp3  
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
  rm2kl(i,i)=temp5*(ONE-2*temp6)/temp2
  rmkl(i,i)=temp1*(ONE-temp6)/temp3
  Vkl=Vkl+Glob_PseudoCharge(i)*rmkl(i,i)
  rkl(i,i)=temp1*(ONE+temp6)*temp3
  r2kl(i,i)=Skl*THREEHALF*(ONE+2*temp6)*temp2
  deltarkl(i,i)=temp8*(ONE-temp7)/(temp2*temp3)
enddo
Vkl=Vkl*Glob_PseudoCharge0
do i=1,n
  do j=i+1,n
    temp2=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    temp3=sqrt(temp2)
    eta1(j,i)=temp2
    sqrt_eta1(j,i)=temp3
    !u1'=tvk'*inv_tAkl*Jij*inv_tAkl
    do q=1,n
      temp4=ZERO
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
    temp7=temp4/(temp2*tau3)
    temp6=temp7/THREE
    rm2kl(j,i)=temp5*(ONE-2*temp6)/temp2
    rm2kl(i,j)=rm2kl(j,i)
    rmkl(j,i)=temp1*(ONE-temp6)/temp3
    rmkl(i,j)=rmkl(j,i)
    Vkl=Vkl+Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*rmkl(j,i)
    rkl(j,i)=temp1*(ONE+temp6)*temp3
    rkl(i,j)=rkl(j,i)
    r2kl(j,i)=Skl*THREEHALF*(ONE+2*temp6)*temp2
    r2kl(i,j)=r2kl(j,i)
    deltarkl(j,i)=temp8*(ONE-temp7)/(temp2*temp3)
    deltarkl(i,j)=deltarkl(j,i)
  enddo  
enddo

Hkl=Tkl+Vkl

!Evaluation of the Darwin correction
Mass_For_Darwin(0)=Glob_Mass(1)
Mass_For_Darwin(1:n)=Glob_Mass(2:n+1)
!Mass_For_Darwin(0)=10.0D20
!Mass_For_Darwin(1)=10.0D20
!Mass_For_Darwin(2)=10.0D20
Darwinkl=ZERO
do i=1,n
  Darwinkl=Darwinkl+(   &
     ONE/(Mass_For_Darwin(0)*Mass_For_Darwin(0)) &
    +ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
    )*Glob_PseudoCharge0*Glob_PseudoCharge(i)*deltarkl(i,i)
enddo
do i=1,n
  do j=1,n
    if(j/=i) then
      Darwinkl=Darwinkl+   &
        ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
       *Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*deltarkl(i,j)
    endif   
  enddo  
enddo
Darwinkl=-Darwinkl*PI/2

!Mass-velocity and Orbit-Orbit are yet to be implemented
MVkl=Skl
OOkl=Skl

end subroutine MatrixElementsL1ForExpcVals



end module matelem
