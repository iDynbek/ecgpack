module matelem
!Module matelem contains subroutines for computing 
!matrix elements with real L=2 Gaussians.
use globvars
implicit none

!********************************************************************
!                            D states 
!
!
!      Module matelem contains subroutines for computing
!            matrix elements with real Gaussians.
!
!         Also contains unpacking procedure of i_k, i_l 
!           parameters to i_kk, j_kk, and i_ll, j_ll.
!
!            xi*xj+yi*yj-2xi*xj and xi^2+yi^2-2zi^2
!  < r' W_k r exp [ -r' A_k r ] |^O| r' W_l r exp [ -r' A_l r ] > 
!
!********************************************************************

contains

subroutine MatrixElementsL1(i_k, i_l, vechLk, vechLl, P, &
               Hkl, Skl, Dk, Dl, grad_k, grad_l)
  
!**********
!
! Arguments
!
!**********

integer,intent(in)	:: i_k, i_l
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

integer         n ! n = number pseudo-electrons
integer		nk, mk, rk, sk, hk, gk, nl, ml, rl, sl, hl, gl, i_kk, j_kk, i_ll, j_ll ! particle indices
integer         i, j, k, r, s, t, indx, coef
real(dprec)	temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
real(dprec)     temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19
real(dprec)     temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27, temp28, temp29 
real(dprec)     temp30, temp31, temp32, temp33, temp34, temp35, temp36, temp37
! matrices
real(dprec)	Lk(nn,nn), Ll(nn,nn), pLl(nn,nn), LkI(3*nn,3*nn), LlI(3*nn,3*nn)
real(dprec)	Ak(nn,nn), tAl(nn,nn), tAkl(nn,nn), AkI(3*nn,3*nn), tAlI(3*nn,3*nn)
real(dprec)	inv_tAkl(nn,nn), inv_tAklI(3*nn,3*nn), inv_tAklp(nn,nn)
real(dprec)	Pp(3*nn,3*nn), Jm(nn,nn), betaJm(nn,nn), betaJmI(3*nn,3*nn), MI(3*nn,3*nn)
real(dprec)	Wk(3*nn,3*nn), Wl(3*nn,3*nn), Wlp(3*nn,3*nn)
real(dprec)     beta1(3*nn,3*nn), beta2(3*nn,3*nn), beta3(3*nn,3*nn), beta33(3*nn,3*nn)
real(dprec)	beta4(nn,nn), beta4I(3*nn,3*nn), beta5(nn,nn), beta5I(3*nn,3*nn)
real(dprec)     beta6(3*nn,3*nn), beta7(3*nn,3*nn), beta8(3*nn,3*nn), beta9(3*nn,3*nn)
real(dprec)     beta10(3*nn,3*nn), beta11(3*nn,3*nn), beta12(3*nn,3*nn), beta13(3*nn,3*nn)
real(dprec)     beta14(3*nn,3*nn), beta15(3*nn,3*nn), beta16(3*nn,3*nn), beta17(3*nn,3*nn)
real(dprec)     beta18(3*nn,3*nn), beta19(3*nn,3*nn), beta20(3*nn,3*nn), beta21(3*nn,3*nn)
real(dprec)     beta22(3*nn,3*nn), beta23(3*nn,3*nn), beta24(3*nn,3*nn), beta25(3*nn,3*nn)
real(dprec)     beta26(3*nn,3*nn), beta27(3*nn,3*nn), beta28(3*nn,3*nn), beta29(3*nn,3*nn)
real(dprec)     beta30(3*nn,3*nn), beta31(3*nn,3*nn), beta32(3*nn,3*nn), beta34(3*nn,3*nn)
real(dprec)     beta36(3*nn,3*nn), beta35(3*nn,3*nn), beta37(3*nn,3*nn), beta38(3*nn,3*nn)
real(dprec)     beta39(3*nn,3*nn), beta40(3*nn,3*nn), beta41(3*nn,3*nn), beta42(3*nn,3*nn)
real(dprec)     beta43(3*nn,3*nn), beta44(3*nn,3*nn), beta45(3*nn,3*nn), beta46(3*nn,3*nn)
reaL(dprec)     beta47(3*nn,3*nn), X1(nn,nn), X1I(3*nn,3*nn)
real(dprec)     Z1(nn,nn), Z2(nn,nn), Z3(3*nn,3*nn), Z4(3*nn,3*nn), Z5(3*nn,3*nn), Z6(3*nn,3*nn)
real(dprec)     Z7(3*nn,3*nn), Z8(3*nn,3*nn), Z9(3*nn,3*nn), Z10(3*nn,3*nn), Z11(3*nn,3*nn)
real(dprec)     Z12(3*nn,3*nn), Z13(3*nn,3*nn), Z14(3*nn,3*nn), Z15(3*nn,3*nn), Z16(3*nn,3*nn)
! traces and determinant
real(dprec)	nu3, nu6, nu7, nu8, tau1, tau6, tau7, tau8, tau9, tau10, Xi
real(dprec)     det_tAkl
!derivative components
integer		MMsum, Msum, M1, M2, M3
!derivative matrices
real(dprec)     dnu3(3*nn,3*nn), dnu3p(3*nn,3*nn), dnu33(3*nn,3*nn), dnu33p(3*nn,3*nn)
real(dprec)     dtau6(3*nn,3*nn), dtau6p(3*nn,3*nn), dtau66(3*nn,3*nn), dtau66p(3*nn,3*nn)
real(dprec)     dtau666(3*nn,3*nn), dtau666p(3*nn,3*nn), dtau6666(3*nn,3*nn), ddtau6666p(3*nn,3*nn)
real(dprec)     dtau6666M(3*nn,3*nn), dtau6666Mp(3*nn,3*nn)
real(dprec)     dtau7(3*nn,3*nn), dtau7p(3*nn,3*nn), dtau77(3*nn,3*nn), dtau77p(3*nn,3*nn)
real(dprec)     dtau777(3*nn,3*nn), dtau777p(3*nn,3*nn), dtau7777(3*nn,3*nn), ddtau7777p(3*nn,3*nn)
real(dprec)     dtau7777M(3*nn,3*nn), dtau7777Mp(3*nn,3*nn)
real(dprec)     dtau8(3*nn,3*nn), dtau8p(3*nn,3*nn), dtau88(3*nn,3*nn), dtau88p(3*nn,3*nn)
real(dprec)     dtau888(3*nn,3*nn), dtau888p(3*nn,3*nn)
real(dprec)     dtau9(3*nn,3*nn), dtau9p(3*nn,3*nn), dtau99(3*nn,3*nn), dtau99p(3*nn,3*nn)
real(dprec)     dtau999(3*nn,3*nn), dtau999p(3*nn,3*nn), dtau10(3*nn,3*nn), dtau10p(3*nn,3*nn)
real(dprec)     dnu6(3*nn,3*nn), dnu6p(3*nn,3*nn), dnu66(3*nn,3*nn), dnu66p(3*nn,3*nn)
real(dprec)     dnu666(3*nn,3*nn), dnu666p(3*nn,3*nn)
real(dprec)     dnu7(3*nn,3*nn), dnu7p(3*nn,3*nn), dnu77(3*nn,3*nn), dnu77p(3*nn,3*nn)
real(dprec)     dnu777(3*nn,3*nn), dnu777p(3*nn,3*nn)
real(dprec)     dnu8(3*nn,3*nn), dnu8p(3*nn,3*nn), dnu88(3*nn,3*nn), dnu88p(3*nn,3*nn)
real(dprec)     dnu888(3*nn,3*nn), dnu888p(3*nn,3*nn), dnu8888(3*nn,3*nn), dnu8888p(3*nn,3*nn)
real(dprec)     alpha1(nn,nn), alpha1p(nn,nn), alpha2(nn,nn), alpha2p(nn,nn)
real(dprec)	mu(3*nn,3*nn), mup(3*nn,3*nn), mu2(3*nn,3*nn), mu2p(3*nn,3*nn)
real(dprec)	alpha25(nn,nn), alpha25p(nn,nn), pk(Glob_np), pl(Glob_np)
!derivative vech
real(dprec)	trans(3*Glob_n*(3*Glob_n+1)/2,Glob_np), trans1(3*Glob_n*(3*Glob_n+1)/2,Glob_np)
real(dprec)     resultAk(Glob_n*(Glob_n+1)/2), resultAl(Glob_n*(Glob_n+1)/2)
real(dprec)	resultX1k(Glob_n*(Glob_n+1)/2), resultX1l(Glob_n*(Glob_n+1)/2)
real(dprec)	result4k(Glob_n*(Glob_n+1)/2), result4l(Glob_n*(Glob_n+1)/2)
real(dprec)	resultJk(Glob_n*(Glob_n+1)/2), resultJl(Glob_n*(Glob_n+1)/2)
real(dprec)     dtau6k(3*Glob_n*(3*Glob_n+1)/2), dtau6l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dtau7k(3*Glob_n*(3*Glob_n+1)/2), dtau7l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dtau8k(3*Glob_n*(3*Glob_n+1)/2), dtau8l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dtau9k(3*Glob_n*(3*Glob_n+1)/2), dtau9l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dtau10k(3*Glob_n*(3*Glob_n+1)/2), dtau10l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dnu3k(3*Glob_n*(3*Glob_n+1)/2), dnu3l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dnu6k(3*Glob_n*(3*Glob_n+1)/2), dnu6l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dnu7k(3*Glob_n*(3*Glob_n+1)/2), dnu7l(3*Glob_n*(3*Glob_n+1)/2)
real(dprec)     dnu8k(3*Glob_n*(3*Glob_n+1)/2), dnu8l(3*Glob_n*(3*Glob_n+1)/2)
!quantities
real(dprec)     Tkl, Rkl, Vkl, SSkl

!write formating
9292 format(1x,10f14.7)
9293 format(1x,21f4.1)

n=Glob_n

!write(*,*) 'Glob_np, Glob_n', Glob_np, Glob_n
!write(*,*)'n,nn',n,nn

!write(*,*)'i_k,i_l',i_k,i_l

!********************
!
! END OF DECLARATIONS
!
!********************

!*******************************
!
! ZEROING ALL DESIRED QUANTITIES 
!
!*******************************

Skl=ZERO
Tkl=ZERO
Vkl=ZERO
Hkl=ZERO
Dk=ZERO
Dl=ZERO

!*******************************
!
! UNPACKING INTEGERS i_k AND i_l
!
! that come from the head
! of the subroutine
! i_k correspounds to basis k
! i_l correspounds to basis l
!
! to make i_kk, j_kk and
! to make i_ll, j_ll
!
!*******************************

i_kk=ZERO
j_kk=ZERO

if (i_k.eq.1) i_kk= 1
if (i_k.eq.1) j_kk= 1

if (i_k.eq.2) i_kk= 1
if (i_k.eq.2) j_kk= 2

if (i_k.eq.3) i_kk= 2
if (i_k.eq.3) j_kk= 2

if (i_k.eq.4) i_kk= 1
if (i_k.eq.4) j_kk= 3

if (i_k.eq.5) i_kk= 2 
if (i_k.eq.5) j_kk= 3 

if (i_k.eq.6) i_kk= 3 
if (i_k.eq.6) j_kk= 3 

if (i_k.eq.7) i_kk= 1 
if (i_k.eq.7) j_kk= 4 

if (i_k.eq.8) i_kk= 2 
if (i_k.eq.8) j_kk= 4 

if (i_k.eq.9) i_kk= 3 
if (i_k.eq.9) j_kk= 4 

if (i_k.eq.10) i_kk= 4 
if (i_k.eq.10) j_kk= 4 

if (i_k.eq.11) i_kk= 1 
if (i_k.eq.11) j_kk= 5 

if (i_k.eq.12) i_kk= 2 
if (i_k.eq.12) j_kk= 5 

if (i_k.eq.13) i_kk= 3 
if (i_k.eq.13) j_kk= 5 

if (i_k.eq.14) i_kk= 4 
if (i_k.eq.14) j_kk= 5 

if (i_k.eq.15) i_kk= 5 
if (i_k.eq.15) j_kk= 5 


i_ll=ZERO
j_ll=ZERO

if (i_l.eq.1) i_ll= 1 
if (i_l.eq.1) j_ll= 1 

if (i_l.eq.2) i_ll= 1 
if (i_l.eq.2) j_ll= 2 	

if (i_l.eq.3) i_ll= 2 
if (i_l.eq.3) j_ll= 2 	

if (i_l.eq.4) i_ll= 1 
if (i_l.eq.4) j_ll= 3 

if (i_l.eq.5) i_ll= 2 
if (i_l.eq.5) j_ll= 3 

if (i_l.eq.6) i_ll= 3 
if (i_l.eq.6) j_ll= 3 

if (i_l.eq.7) i_ll= 1 
if (i_l.eq.7) j_ll= 4 

if (i_l.eq.8) i_ll= 2 
if (i_l.eq.8) j_ll= 4 

if (i_l.eq.9) i_ll= 3 
if (i_l.eq.9) j_ll= 4 

if (i_l.eq.10) i_ll= 4 
if (i_l.eq.10) j_ll= 4

if (i_l.eq.11) i_ll= 1 
if (i_l.eq.11) j_ll= 5 

if (i_l.eq.12) i_ll= 2
if (i_l.eq.12) j_ll= 5 

if (i_l.eq.13) i_ll= 3 
if (i_l.eq.13) j_ll= 5 

if (i_l.eq.14) i_ll= 4 
if (i_l.eq.14) j_ll= 5 

if (i_l.eq.15) i_ll= 5 
if (i_l.eq.15) j_ll= 5 

!write(*,*)'i_kk,j_kk,i_ll,j_ll',i_kk,j_kk,i_ll,j_ll

!******************
!
! Kronecker Product
!
! of Mass Matrix
! and permutation 
! operator 
!
!******************

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

!*******************
!
! BUILDING Wk AND Wl
!
!*******************

Wk=ZERO

if (i_kk.ne.j_kk) then

   Wk((i_kk-1)*3+1,(j_kk-1)*3+1)=ONEHALF
   Wk((i_kk-1)*3+2,(j_kk-1)*3+2)=ONEHALF
   Wk((i_kk-1)*3+3,(j_kk-1)*3+3)=-ONE

   Wk((j_kk-1)*3+1,(i_kk-1)*3+1)=ONEHALF
   Wk((j_kk-1)*3+2,(i_kk-1)*3+2)=ONEHALF
   Wk((j_kk-1)*3+3,(i_kk-1)*3+3)=-ONE

nk=3*(i_kk-1)+1
mk=3*(i_kk-1)+2
rk=3*(i_kk-1)+3
sk=3*(j_kk-1)+1
gk=3*(j_kk-1)+2
hk=3*(j_kk-1)+3

endif

if (i_kk.eq.j_kk) then

   Wk((i_kk-1)*3+1,(i_kk-1)*3+1)=ONE
   Wk((i_kk-1)*3+2,(i_kk-1)*3+2)=ONE
   Wk((i_kk-1)*3+3,(i_kk-1)*3+3)=-TWO

nk=3*(i_kk-1)+1
mk=3*(i_kk-1)+2
rk=3*(i_kk-1)+3

endif

Wl=ZERO

if (i_ll.ne.j_ll) then

   Wl((i_ll-1)*3+1,(j_ll-1)*3+1)=ONEHALF
   Wl((i_ll-1)*3+2,(j_ll-1)*3+2)=ONEHALF
   Wl((i_ll-1)*3+3,(j_ll-1)*3+3)=-ONE
        
   Wl((j_ll-1)*3+1,(i_ll-1)*3+1)=ONEHALF
   Wl((j_ll-1)*3+2,(i_ll-1)*3+2)=ONEHALF
   Wl((j_ll-1)*3+3,(i_ll-1)*3+3)=-ONE

nl=3*(i_ll-1)+1
ml=3*(i_ll-1)+2
rl=3*(i_ll-1)+3
sl=3*(j_ll-1)+1
gl=3*(j_ll-1)+2
hl=3*(j_ll-1)+3

endif

if (i_ll.eq.j_ll) then

   Wl((i_ll-1)*3+1,(i_ll-1)*3+1)=ONE
   Wl((i_ll-1)*3+2,(i_ll-1)*3+2)=ONE
   Wl((i_ll-1)*3+3,(i_ll-1)*3+3)=-TWO

nl=3*(i_ll-1)+1
ml=3*(i_ll-1)+2
rl=3*(i_ll-1)+3

endif

!write(*,*)'nl,ml,rl,sl,gl,hl',nl,ml,rl,sl,gl,hl

!write(*,*)'Wk'
!do i=1,3*n
!  write(*,9293)(Wk(i,j),j=1,3*n)
!enddo

!write(*,*)'Wl'
!do i=1,3*n
!  write(*,9293)(Wl(i,j),j=1,3*n)
!enddo


!*************
!
! Permuting Wl
!
!*************

Z3=ZERO

do i=1,n

if (i_ll.ne.j_ll) then 

  Z3(3*(i-1)+1,nl)=Pp(sl,3*(i-1)+1)*Wl(sl,nl)
  Z3(3*(i-1)+2,ml)=Pp(gl,3*(i-1)+2)*Wl(gl,ml)
  Z3(3*(i-1)+3,rl)=Pp(hl,3*(i-1)+3)*Wl(hl,rl)
  Z3(3*(i-1)+1,sl)=Pp(nl,3*(i-1)+1)*Wl(nl,sl)
  Z3(3*(i-1)+2,gl)=Pp(ml,3*(i-1)+2)*Wl(ml,gl)
  Z3(3*(i-1)+3,hl)=Pp(rl,3*(i-1)+3)*Wl(rl,hl)

endif

if (i_ll.eq.j_ll) then

  Z3(3*(i-1)+1,nl)=Pp(nl,3*(i-1)+1)*Wl(nl,nl)
  Z3(3*(i-1)+2,ml)=Pp(ml,3*(i-1)+2)*Wl(ml,ml)
  Z3(3*(i-1)+3,rl)=Pp(rl,3*(i-1)+3)*Wl(rl,rl)

endif

enddo

!write(*,*)'Z3'
!do i=1,3*n
!  write(*,9292)(Z3(i,j),j=1,3*n)
!enddo

do i=1,n
  do j=1,n
    
if (i_ll.ne.j_ll) then

    Wlp(3*(i-1)+1,3*(j-1)+1)=Z3(3*(i-1)+1,nl)*Pp(nl,3*(j-1)+1)&
                            +Z3(3*(i-1)+1,sl)*Pp(sl,3*(j-1)+1)
    Wlp(3*(i-1)+2,3*(j-1)+2)=Z3(3*(i-1)+2,ml)*Pp(ml,3*(j-1)+2)&
                            +Z3(3*(i-1)+2,gl)*Pp(gl,3*(j-1)+2)
    Wlp(3*(i-1)+3,3*(j-1)+3)=Z3(3*(i-1)+3,rl)*Pp(rl,3*(j-1)+3)&
                            +Z3(3*(i-1)+3,hl)*Pp(hl,3*(j-1)+3)
endif

if (i_ll.eq.j_ll) then

    Wlp(3*(i-1)+1,3*(j-1)+1)=Z3(3*(i-1)+1,nl)*Pp(nl,3*(j-1)+1)
    Wlp(3*(i-1)+2,3*(j-1)+2)=Z3(3*(i-1)+2,ml)*Pp(ml,3*(j-1)+2)
    Wlp(3*(i-1)+3,3*(j-1)+3)=Z3(3*(i-1)+3,rl)*Pp(rl,3*(j-1)+3)

endif

  enddo
enddo

do i=1,n
  do j=1,n

if (i_ll.ne.j_ll) then
  
  if (Wlp(3*(i-1)+1,3*(j-1)+1).gt.ZERO) then
     sl=3*(i-1)+1
     nl=3*(j-1)+1
   endif
   
   if (Wlp(3*(i-1)+2,3*(j-1)+2).gt.ZERO) then
     gl=3*(i-1)+2
     ml=3*(j-1)+2
   endif

   if (Wlp(3*(i-1)+3,3*(j-1)+3).lt.ZERO) then
     hl=3*(i-1)+3
     rl=3*(j-1)+3 
   endif
endif

if (i_ll.eq.j_ll) then

  if (Wlp(3*(i-1)+1,3*(j-1)+1).eq.ONE) then
    nl=3*(i-1)+1
  endif

  if (Wlp(3*(i-1)+2,3*(j-1)+2).eq.ONE) then
   ml=3*(i-1)+2
  endif

  if (Wlp(3*(i-1)+3,3*(j-1)+3).eq.-TWO) then
    rl=3*(i-1)+3
  endif

endif

  enddo
enddo

!write(*,*)'nl,ml,rl,sl,gl,hl',nl,ml,rl,sl,gl,hl

!write(*,*)'Wlp'
!do i=1,3*n
!  write(*,9293)(Wlp(i,j),j=1,3*n)
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

!Building Pp and at the same time building MI.

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

!beta1=inv_tAklI*Wl
!beta2=inv_tAklI*Wlp

beta1=ZERO
beta2=ZERO
beta6=ZERO

do i=1,n
    
if (i_kk.ne.j_kk) then

  beta1(3*(i-1)+1,nk)=inv_tAklI(3*(i-1)+1,sk)*Wk(sk,nk)
  beta1(3*(i-1)+2,mk)=inv_tAklI(3*(i-1)+2,gk)*Wk(gk,mk)
  beta1(3*(i-1)+3,rk)=inv_tAklI(3*(i-1)+3,hk)*Wk(hk,rk)
  beta1(3*(i-1)+1,sk)=inv_tAklI(3*(i-1)+1,nk)*Wk(nk,sk)
  beta1(3*(i-1)+2,gk)=inv_tAklI(3*(i-1)+2,mk)*Wk(mk,gk)
  beta1(3*(i-1)+3,hk)=inv_tAklI(3*(i-1)+3,rk)*Wk(rk,hk)

endif

if (i_kk.eq.j_kk) then

  beta1(3*(i-1)+1,nk)=inv_tAklI(3*(i-1)+1,nk)*Wk(nk,nk)
  beta1(3*(i-1)+2,mk)=inv_tAklI(3*(i-1)+2,mk)*Wk(mk,mk)
  beta1(3*(i-1)+3,rk)=inv_tAklI(3*(i-1)+3,rk)*Wk(rk,rk)

endif

if (i_ll.ne.j_ll) then 

  beta2(3*(i-1)+1,nl)=inv_tAklI(3*(i-1)+1,sl)*Wlp(sl,nl) 
  beta2(3*(i-1)+2,ml)=inv_tAklI(3*(i-1)+2,gl)*Wlp(gl,ml) 
  beta2(3*(i-1)+3,rl)=inv_tAklI(3*(i-1)+3,hl)*Wlp(hl,rl) 
  beta2(3*(i-1)+1,sl)=inv_tAklI(3*(i-1)+1,nl)*Wlp(nl,sl)
  beta2(3*(i-1)+2,gl)=inv_tAklI(3*(i-1)+2,ml)*Wlp(ml,gl)
  beta2(3*(i-1)+3,hl)=inv_tAklI(3*(i-1)+3,rl)*Wlp(rl,hl)

  beta6(3*(i-1)+1,nl)=MI(3*(i-1)+1,sl)*Wlp(sl,nl)
  beta6(3*(i-1)+2,ml)=MI(3*(i-1)+2,gl)*Wlp(gl,ml)
  beta6(3*(i-1)+3,rl)=MI(3*(i-1)+3,hl)*Wlp(hl,rl)
  beta6(3*(i-1)+1,sl)=MI(3*(i-1)+1,nl)*Wlp(nl,sl)
  beta6(3*(i-1)+2,gl)=MI(3*(i-1)+2,ml)*Wlp(ml,gl)
  beta6(3*(i-1)+3,hl)=MI(3*(i-1)+3,rl)*Wlp(rl,hl)

endif

if (i_ll.eq.j_ll) then

  beta2(3*(i-1)+1,nl)=inv_tAklI(3*(i-1)+1,nl)*Wlp(nl,nl)
  beta2(3*(i-1)+2,ml)=inv_tAklI(3*(i-1)+2,ml)*Wlp(ml,ml)
  beta2(3*(i-1)+3,rl)=inv_tAklI(3*(i-1)+3,rl)*Wlp(rl,rl)

  beta6(3*(i-1)+1,nl)=MI(3*(i-1)+1,nl)*Wlp(nl,nl)
  beta6(3*(i-1)+2,ml)=MI(3*(i-1)+2,ml)*Wlp(ml,ml)
  beta6(3*(i-1)+3,rl)=MI(3*(i-1)+3,rl)*Wlp(rl,rl)

endif

enddo 

!write(*,*)'beta1'
!do i=1,3*n
!  write(*,9292)(beta1(i,j),j=1,3*n)
!enddo

!write(*,*)'beta2'
!do i=1,3*n
!  write(*,9292)(beta2(i,j),j=1,3*n)
!enddo

!beta3=beta1*beta2
!beta33=beta2*beta1

beta3=ZERO
beta33=ZERO

do i=1,n
 
if (i_kk.ne.j_kk.and.i_ll.ne.j_ll) then

  beta3(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta2(nk,nl)&
                     +beta1(3*(i-1)+1,sk)*beta2(sk,nl) 
  beta3(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta2(mk,ml)&  
                     +beta1(3*(i-1)+2,gk)*beta2(gk,ml)    
  beta3(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta2(rk,rl)&  
                     +beta1(3*(i-1)+3,hk)*beta2(hk,rl)    
  beta3(3*(i-1)+1,sl)=beta1(3*(i-1)+1,nk)*beta2(nk,sl)&  
                     +beta1(3*(i-1)+1,sk)*beta2(sk,sl)    
  beta3(3*(i-1)+2,gl)=beta1(3*(i-1)+2,mk)*beta2(mk,gl)&  
                     +beta1(3*(i-1)+2,gk)*beta2(gk,gl)    
  beta3(3*(i-1)+3,hl)=beta1(3*(i-1)+3,rk)*beta2(rk,hl)&  
                     +beta1(3*(i-1)+3,hk)*beta2(hk,hl)    

  beta33(3*(i-1)+1,nk)=beta2(3*(i-1)+1,nl)*beta1(nl,nk)&
                      +beta2(3*(i-1)+1,sl)*beta1(sl,nk)
  beta33(3*(i-1)+2,mk)=beta2(3*(i-1)+2,ml)*beta1(ml,mk)&
                      +beta2(3*(i-1)+2,gl)*beta1(gl,mk)
  beta33(3*(i-1)+3,rk)=beta2(3*(i-1)+3,rl)*beta1(rl,rk)&
                      +beta2(3*(i-1)+3,hl)*beta1(hl,rk)
  beta33(3*(i-1)+1,sk)=beta2(3*(i-1)+1,nl)*beta1(nl,sk)&
                      +beta2(3*(i-1)+1,sl)*beta1(sl,sk)
  beta33(3*(i-1)+2,gk)=beta2(3*(i-1)+2,ml)*beta1(ml,gk)&
                      +beta2(3*(i-1)+2,gl)*beta1(gl,gk)
  beta33(3*(i-1)+3,hk)=beta2(3*(i-1)+3,rl)*beta1(rl,hk)&
                      +beta2(3*(i-1)+3,hl)*beta1(hl,hk)

endif

if (i_kk.eq.j_kk.and.i_ll.eq.j_ll) then

  beta3(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta2(nk,nl)
  beta3(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta2(mk,ml)
  beta3(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta2(rk,rl)

  beta33(3*(i-1)+1,nk)=beta2(3*(i-1)+1,nl)*beta1(nl,nk)
  beta33(3*(i-1)+2,mk)=beta2(3*(i-1)+2,ml)*beta1(ml,mk)
  beta33(3*(i-1)+3,rk)=beta2(3*(i-1)+3,rl)*beta1(rl,rk)

endif

if (i_kk.eq.j_kk.and.i_ll.ne.j_ll) then

  beta3(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta2(nk,nl)
  beta3(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta2(mk,ml)
  beta3(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta2(rk,rl)
  beta3(3*(i-1)+1,sl)=beta1(3*(i-1)+1,nk)*beta2(nk,sl)
  beta3(3*(i-1)+2,gl)=beta1(3*(i-1)+2,mk)*beta2(mk,gl)
  beta3(3*(i-1)+3,hl)=beta1(3*(i-1)+3,rk)*beta2(rk,hl)

  beta33(3*(i-1)+1,nk)=beta2(3*(i-1)+1,nl)*beta1(nl,nk)&
                      +beta2(3*(i-1)+1,sl)*beta1(sl,nk)
  beta33(3*(i-1)+2,mk)=beta2(3*(i-1)+2,ml)*beta1(ml,mk)&
                      +beta2(3*(i-1)+2,gl)*beta1(gl,mk)
  beta33(3*(i-1)+3,rk)=beta2(3*(i-1)+3,rl)*beta1(rl,rk)&
                      +beta2(3*(i-1)+3,hl)*beta1(hl,rk)

endif

if (i_kk.ne.j_kk.and.i_ll.eq.j_ll) then

  beta3(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta2(nk,nl)&
                     +beta1(3*(i-1)+1,sk)*beta2(sk,nl)
  beta3(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta2(mk,ml)&
                     +beta1(3*(i-1)+2,gk)*beta2(gk,ml)
  beta3(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta2(rk,rl)&
                     +beta1(3*(i-1)+3,hk)*beta2(hk,rl)

  beta33(3*(i-1)+1,nk)=beta2(3*(i-1)+1,nl)*beta1(nl,nk)
  beta33(3*(i-1)+2,mk)=beta2(3*(i-1)+2,ml)*beta1(ml,mk)
  beta33(3*(i-1)+3,rk)=beta2(3*(i-1)+3,rl)*beta1(rl,rk)
  beta33(3*(i-1)+1,sk)=beta2(3*(i-1)+1,nl)*beta1(nl,sk)
  beta33(3*(i-1)+2,gk)=beta2(3*(i-1)+2,ml)*beta1(ml,gk)
  beta33(3*(i-1)+3,hk)=beta2(3*(i-1)+3,rl)*beta1(rl,hk)

endif

enddo

!write(*,*)'beta3'
!do i=1,3*n
!  write(*,9292)(beta3(i,j),j=1,3*n)
!enddo

!write(*,*)'beta33'
!do i=1,3*n
!  write(*,9292)(beta33(i,j),j=1,3*n)
!enddo

!***********************
!
! Taking trace of Betas
!
!***********************

!nu3=tr[beta3]

nu3=ZERO

do  i=1,n

  nu3=nu3+beta3(3*(i-1)+1,3*(i-1)+1)&
         +beta3(3*(i-1)+2,3*(i-1)+2)&
         +beta3(3*(i-1)+3,3*(i-1)+3)

enddo

!write(*,*)'nu3,det_tAkl',nu3,det_tAkl

!**************************************************************
!
! Derivative of Determinant with respect to vechL_k and vechL_l
!
!**************************************************************

if (grad_k) then 

  resultAk=ZERO
  indx=ZERO
  do j=1,n
    do i=j,n
      indx=indx+1
      temp1=ZERO
      do k=j,n

        temp1=temp1+(inv_tAkl(k,i)+inv_tAkl(i,k))*Lk(k,j)

      enddo

      resultAk(indx)=temp1

    enddo
  enddo
endif

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

  resultAl=ZERO
  indx=ZERO
  do j=1,n
    do i=j,n
      
      indx=indx+1
      temp1=ZERO
      do k=j,n

        temp1=temp1+(inv_tAklp(i,k)+inv_tAklp(k,i))*Ll(k,j)

      enddo

      resultAl(indx)=temp1

    enddo
  enddo
endif

!*******************
!
! Constructing Trans
!
!*******************

if (grad_k.or.grad_l) then

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

    enddo
  enddo

!write(*,*)'Trans'
!do i=1,3*n*(3*n+1)/2
!  write(*,9293)(trans(i,j),j=1,n*(n+1)/2)
!enddo

!dnu3=beta3*inv_tAklI
!dnu33=beta33*inv_tAklI

  dnu3=ZERO
  dnu33=ZERO

  do i=1,n
    do j=1,n

if (i_kk.ne.j_kk) then

      dnu33(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)&
                                +beta33(3*(i-1)+1,sk)*inv_tAklI(sk,3*(j-1)+1) 
      dnu33(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)&
                                +beta33(3*(i-1)+2,gk)*inv_tAklI(gk,3*(j-1)+2) 
      dnu33(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)&
                                +beta33(3*(i-1)+3,hk)*inv_tAklI(hk,3*(j-1)+3) 

endif

if (i_kk.eq.j_kk) then

      dnu33(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)
      dnu33(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)
      dnu33(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)

endif

if (i_ll.ne.j_ll) then

      dnu3(3*(i-1)+1,3*(j-1)+1)=beta3(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)&
                               +beta3(3*(i-1)+1,sl)*inv_tAklI(sl,3*(j-1)+1)
      dnu3(3*(i-1)+2,3*(j-1)+2)=beta3(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)&
                               +beta3(3*(i-1)+2,gl)*inv_tAklI(gl,3*(j-1)+2)
      dnu3(3*(i-1)+3,3*(j-1)+3)=beta3(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)&
                               +beta3(3*(i-1)+3,hl)*inv_tAklI(hl,3*(j-1)+3)

endif

if (i_ll.eq.j_ll) then

      dnu3(3*(i-1)+1,3*(j-1)+1)=beta3(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)
      dnu3(3*(i-1)+2,3*(j-1)+2)=beta3(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)
      dnu3(3*(i-1)+3,3*(j-1)+3)=beta3(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)

endif

    enddo
  enddo

endif

if (grad_k) then

  dnu3k=ZERO
  indx=ZERO
  
  do j=1,3*n

    do i=j,3*n
    
      indx=indx+1

      temp1=ZERO
      temp2=ZERO

      do k=j,3*n

        temp1=temp1+(dnu3(i,k)+dnu3(k,i))*LkI(k,j)
        temp2=temp2+(dnu33(i,k)+dnu33(k,i))*LkI(k,j)

      enddo

      dnu3k(indx)=temp1+temp2

    enddo
  enddo

endif

if (grad_l) then

  Z5=ZERO
  Z6=ZERO

  do i=1,3*n
    do j=1,3*n

      temp1=ZERO
      temp2=ZERO 

      do k=1,3*n  

        temp1=temp1+pP(i,k)*dnu3(k,j)
        temp2=temp2+pP(i,k)*dnu33(k,j)

      enddo  

      Z5(i,j)=temp1
      Z6(i,j)=temp2

    enddo
  enddo  

  do i=1,3*n
    do j=1,3*n

      temp1=ZERO
      temp2=ZERO
 
      do k=1,3*n  

        temp1=temp1+Z5(i,k)*pP(j,k)
        temp2=temp2+Z6(i,k)*pP(j,k)

      enddo  

      dnu3p(i,j)=temp1
      dnu33p(i,j)=temp2

    enddo
  enddo  
 
  dnu3l=ZERO
  indx=ZERO

  do j=1,3*n
    do i=j,3*n
      
      indx=indx+1
      
      temp1=ZERO
      temp2=ZERO
      
      do k=j,3*n

        temp1=temp1+(dnu3p(i,k)+dnu3p(k,i))*LlI(k,j)
        temp2=temp2+(dnu33p(i,k)+dnu33p(k,i))*LlI(k,j)

      enddo

      dnu3l(indx)=temp1+temp2   

    enddo
  enddo
endif

!*****************
!
! OVERLAP FORMULA
!
!*****************

!write(*,*)'i_k,i_l,i_kk,j_kk,i_ll,j_ll',i_k,i_l,i_kk,j_kk,i_ll,j_ll

Skl=ONEHALF*(PI**((THREE*n)/TWO))*(det_tAkl**(-THREE/TWO))*nu3

!write(*,*) 'Skl',Skl

!*****************
!
! OVERLAP GRADIENT
!
!*****************

  if (grad_k) then

    do i=1,n*(n+1)/2
    
      temp1=ZERO

      do k=1,3*n*(3*n+1)/2
    
        temp1=temp1+dnu3k(k)*trans(k,i)
    
      enddo
    
      Dk(Glob_np+i)=Dk(Glob_np+i)-ONEHALF*(PI**(THREEHALF*n))&
                   *(det_tAkl**(-THREEHALF))&
                   *(THREEHALF*resultAk(i)*nu3+temp1)

  enddo

! write(*,*) 'Dk',Dk

  endif

  if (grad_l) then

    do i=1,n*(n+1)/2
    
      temp1=ZERO

      do k=1,3*n*(3*n+1)/2

        temp1=temp1+dnu3l(k)*trans(k,i)

      enddo

      Dl(Glob_np+i)=Dl(Glob_np+i)-ONEHALF*(PI**(THREEHALF*n))&
                   *(det_tAkl**(-THREEHALF))&
                   *(THREEHALF*resultAl(i)*nu3+temp1)

    enddo

!    write(*,*) 'Dl',Dl

  endif  
 
!go to 1111

!***********************************
!
! KINETIC ENERGY INTEGRAL
!
!***********************************

!****************************
!
! DOING MATRIX MULTIPLICATION
!
!****************************

!beta4=inv_tAkl*Ak
!beta5=M*tAl

beta4=ZERO
beta5=ZERO

do i=1,n
  do j=1,n
    
    temp1=ZERO
    temp2=ZERO
    do k=1,n
      temp1=temp1+inv_tAkl(i,k)*Ak(k,j)
      temp2=temp2+Glob_MassMatrix(i,k)*tAl(k,j)
    enddo
      beta4(i,j)=temp1
      beta5(i,j)=temp2

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

X1=ZERO

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

    beta4I(3*(i-1)+1,3*(j-1)+1)=beta4(i,j)
    beta4I(3*(i-1)+2,3*(j-1)+2)=beta4(i,j)
    beta4I(3*(i-1)+3,3*(j-1)+3)=beta4(i,j)

    beta5I(3*(i-1)+1,3*(j-1)+1)=beta5(i,j)
    beta5I(3*(i-1)+2,3*(j-1)+2)=beta5(i,j)
    beta5I(3*(i-1)+3,3*(j-1)+3)=beta5(i,j)
 
    X1I(3*(i-1)+1,3*(j-1)+1)=X1(i,j)
    X1I(3*(i-1)+2,3*(j-1)+2)=X1(i,j)
    X1I(3*(i-1)+3,3*(j-1)+3)=X1(i,j)

  enddo
enddo

!write(*,*)'beta4I'
!do i=1,3*n
!  write(*,9292)(beta4I(i,j),j=1,3*n)
!enddo

!write(*,*)'beta5I'
!do i=1,3*n
!  write(*,9292)(beta5I(i,j),j=1,3*n)
!enddo

!tau6
beta7=ZERO
beta8=ZERO
beta10=ZERO
beta11=ZERO
beta12=ZERO

do i=1,n

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
  temp17=ZERO
  temp18=ZERO

  do j=1,n

if (i_ll.ne.j_ll) then

  beta7(3*(i-1)+1,3*(j-1)+1)=beta3(3*(i-1)+1,nl)*X1I(nl,3*(j-1)+1)&
                            +beta3(3*(i-1)+1,sl)*X1I(sl,3*(j-1)+1)
  beta7(3*(i-1)+2,3*(j-1)+2)=beta3(3*(i-1)+2,ml)*X1I(ml,3*(j-1)+2)&
                            +beta3(3*(i-1)+2,gl)*X1I(gl,3*(j-1)+2) 
  beta7(3*(i-1)+3,3*(j-1)+3)=beta3(3*(i-1)+3,rl)*X1I(rl,3*(j-1)+3)&
                            +beta3(3*(i-1)+3,hl)*X1I(hl,3*(j-1)+3)

  beta8(3*(i-1)+1,3*(j-1)+1)=beta2(3*(i-1)+1,nl)*X1I(nl,3*(j-1)+1)&
                            +beta2(3*(i-1)+1,sl)*X1I(sl,3*(j-1)+1)
  beta8(3*(i-1)+2,3*(j-1)+2)=beta2(3*(i-1)+2,ml)*X1I(ml,3*(j-1)+2)&
                            +beta2(3*(i-1)+2,gl)*X1I(gl,3*(j-1)+2)
  beta8(3*(i-1)+3,3*(j-1)+3)=beta2(3*(i-1)+3,rl)*X1I(rl,3*(j-1)+3)&
                            +beta2(3*(i-1)+3,hl)*X1I(hl,3*(j-1)+3)

  temp1=temp1+X1I(3*(i-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,nl)
  temp2=temp2+X1I(3*(i-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,ml)
  temp3=temp3+X1I(3*(i-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,rl)
  temp4=temp4+X1I(3*(i-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,sl)
  temp5=temp5+X1I(3*(i-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,gl)
  temp6=temp6+X1I(3*(i-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,hl)

  temp10=temp10+beta5I(3*(i-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,nl)
  temp11=temp11+beta5I(3*(i-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,ml)
  temp12=temp12+beta5I(3*(i-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,rl)
  temp13=temp13+beta5I(3*(i-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,sl)
  temp14=temp14+beta5I(3*(i-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,gl)
  temp15=temp15+beta5I(3*(i-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,hl)

  beta12(3*(i-1)+1,3*(j-1)+1)=beta3(3*(i-1)+1,nl)*beta4I(nl,3*(j-1)+1)&
                             +beta3(3*(i-1)+1,sl)*beta4I(sl,3*(j-1)+1)
  beta12(3*(i-1)+2,3*(j-1)+2)=beta3(3*(i-1)+2,ml)*beta4I(ml,3*(j-1)+2)&
                             +beta3(3*(i-1)+2,gl)*beta4I(gl,3*(j-1)+2)
  beta12(3*(i-1)+3,3*(j-1)+3)=beta3(3*(i-1)+3,rl)*beta4I(rl,3*(j-1)+3)&
                             +beta3(3*(i-1)+3,hl)*beta4I(hl,3*(j-1)+3)

endif

if (i_ll.eq.j_ll) then

  beta7(3*(i-1)+1,3*(j-1)+1)=beta3(3*(i-1)+1,nl)*X1I(nl,3*(j-1)+1)
  beta7(3*(i-1)+2,3*(j-1)+2)=beta3(3*(i-1)+2,ml)*X1I(ml,3*(j-1)+2)
  beta7(3*(i-1)+3,3*(j-1)+3)=beta3(3*(i-1)+3,rl)*X1I(rl,3*(j-1)+3)

  beta8(3*(i-1)+1,3*(j-1)+1)=beta2(3*(i-1)+1,nl)*X1I(nl,3*(j-1)+1)
  beta8(3*(i-1)+2,3*(j-1)+2)=beta2(3*(i-1)+2,ml)*X1I(ml,3*(j-1)+2)
  beta8(3*(i-1)+3,3*(j-1)+3)=beta2(3*(i-1)+3,rl)*X1I(rl,3*(j-1)+3)

  temp7=temp7+X1I(3*(i-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,nl)
  temp8=temp8+X1I(3*(i-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,ml)
  temp9=temp9+X1I(3*(i-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,rl)

  temp16=temp16+beta5I(3*(i-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,nl)
  temp17=temp17+beta5I(3*(i-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,ml)
  temp18=temp18+beta5I(3*(i-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,rl)

  beta12(3*(i-1)+1,3*(j-1)+1)=beta3(3*(i-1)+1,nl)*beta4I(nl,3*(j-1)+1)
  beta12(3*(i-1)+2,3*(j-1)+2)=beta3(3*(i-1)+2,ml)*beta4I(ml,3*(j-1)+2)
  beta12(3*(i-1)+3,3*(j-1)+3)=beta3(3*(i-1)+3,rl)*beta4I(rl,3*(j-1)+3)

endif

  enddo

if (i_ll.ne.j_ll) then 

  beta10(3*(i-1)+1,nl)=temp1
  beta10(3*(i-1)+2,ml)=temp2
  beta10(3*(i-1)+3,rl)=temp3
  beta10(3*(i-1)+1,sl)=temp4
  beta10(3*(i-1)+2,gl)=temp5
  beta10(3*(i-1)+3,hl)=temp6

  beta11(3*(i-1)+1,nl)=temp10
  beta11(3*(i-1)+2,ml)=temp11
  beta11(3*(i-1)+3,rl)=temp12
  beta11(3*(i-1)+1,sl)=temp13
  beta11(3*(i-1)+2,gl)=temp14
  beta11(3*(i-1)+3,hl)=temp15

endif

if (i_ll.eq.j_ll) then

  beta10(3*(i-1)+1,nl)=temp7
  beta10(3*(i-1)+2,ml)=temp8
  beta10(3*(i-1)+3,rl)=temp9

  beta11(3*(i-1)+1,nl)=temp16
  beta11(3*(i-1)+2,ml)=temp17
  beta11(3*(i-1)+3,rl)=temp18

endif

enddo

!tau6
beta9=ZERO
!tau7
beta13=ZERO
beta14=ZERO
beta16=ZERO
beta17=ZERO
beta18=ZERO
!tau8
beta19=ZERO
!tau9
beta24=ZERO
beta25=ZERO

do i=1,n

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
  temp17=ZERO
  temp18=ZERO
  temp19=ZERO
  temp20=ZERO
  temp21=ZERO
  temp22=ZERO
  temp23=ZERO
  temp24=ZERO
  temp25=ZERO
  temp26=ZERO
  temp27=ZERO

  do j=1,n

if (i_kk.eq.j_kk) then

    temp1=temp1+beta8(3*(i-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,nk)
    temp2=temp2+beta8(3*(i-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,mk)
    temp3=temp3+beta8(3*(i-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,rk)

    temp10=temp10+X1I(3*(i-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,nk)
    temp11=temp11+X1I(3*(i-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,mk)
    temp12=temp12+X1I(3*(i-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,rk)

    temp19=temp19+beta5I(3*(i-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,nk)
    temp20=temp20+beta5I(3*(i-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,mk)
    temp21=temp21+beta5I(3*(i-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,rk)

    beta13(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*X1I(nk,3*(j-1)+1)
    beta13(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*X1I(mk,3*(j-1)+2)
    beta13(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*X1I(rk,3*(j-1)+3)

    beta14(3*(i-1)+1,3*(j-1)+1)=beta1(3*(i-1)+1,nk)*X1I(nk,3*(j-1)+1)
    beta14(3*(i-1)+2,3*(j-1)+2)=beta1(3*(i-1)+2,mk)*X1I(mk,3*(j-1)+2)
    beta14(3*(i-1)+3,3*(j-1)+3)=beta1(3*(i-1)+3,rk)*X1I(rk,3*(j-1)+3)

    beta18(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*beta4I(nk,3*(j-1)+1)
    beta18(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*beta4I(mk,3*(j-1)+2)
    beta18(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*beta4I(rk,3*(j-1)+3)

    beta19(3*(i-1)+1,3*(j-1)+1)=beta1(3*(i-1)+1,nk)*beta4I(nk,3*(j-1)+1)
    beta19(3*(i-1)+2,3*(j-1)+2)=beta1(3*(i-1)+2,mk)*beta4I(mk,3*(j-1)+2)
    beta19(3*(i-1)+3,3*(j-1)+3)=beta1(3*(i-1)+3,rk)*beta4I(rk,3*(j-1)+3)

    beta24(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*beta5I(nk,3*(j-1)+1)
    beta24(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*beta5I(mk,3*(j-1)+2)
    beta24(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*beta5I(rk,3*(j-1)+3)

    beta25(3*(i-1)+1,3*(j-1)+1)=beta1(3*(i-1)+1,nk)*beta5I(nk,3*(j-1)+1)
    beta25(3*(i-1)+2,3*(j-1)+2)=beta1(3*(i-1)+2,mk)*beta5I(mk,3*(j-1)+2)
    beta25(3*(i-1)+3,3*(j-1)+3)=beta1(3*(i-1)+3,rk)*beta5I(rk,3*(j-1)+3)

endif

if (i_kk.ne.j_kk) then

    temp4=temp4+beta8(3*(i-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,nk)
    temp5=temp5+beta8(3*(i-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,mk)
    temp6=temp6+beta8(3*(i-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,rk)
    temp7=temp7+beta8(3*(i-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,sk)
    temp8=temp8+beta8(3*(i-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,gk)
    temp9=temp9+beta8(3*(i-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,hk)

    temp13=temp13+X1I(3*(i-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,nk)
    temp14=temp14+X1I(3*(i-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,mk)
    temp15=temp15+X1I(3*(i-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,rk)
    temp16=temp16+X1I(3*(i-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,sk)
    temp17=temp17+X1I(3*(i-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,gk)
    temp18=temp18+X1I(3*(i-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,hk)
 
    temp22=temp22+beta5I(3*(i-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,nk)
    temp23=temp23+beta5I(3*(i-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,mk)
    temp24=temp24+beta5I(3*(i-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,rk)
    temp25=temp25+beta5I(3*(i-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,sk)
    temp26=temp26+beta5I(3*(i-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,gk)
    temp27=temp27+beta5I(3*(i-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,hk)

    beta13(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*X1I(nk,3*(j-1)+1)&
                               +beta33(3*(i-1)+1,sk)*X1I(sk,3*(j-1)+1)
    beta13(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*X1I(mk,3*(j-1)+2)&
                               +beta33(3*(i-1)+2,gk)*X1I(gk,3*(j-1)+2)
    beta13(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*X1I(rk,3*(j-1)+3)&
                               +beta33(3*(i-1)+3,hk)*X1I(hk,3*(j-1)+3)

    beta14(3*(i-1)+1,3*(j-1)+1)=beta1(3*(i-1)+1,nk)*X1I(nk,3*(j-1)+1)&
                               +beta1(3*(i-1)+1,sk)*X1I(sk,3*(j-1)+1)
    beta14(3*(i-1)+2,3*(j-1)+2)=beta1(3*(i-1)+2,mk)*X1I(mk,3*(j-1)+2)&
                               +beta1(3*(i-1)+2,gk)*X1I(gk,3*(j-1)+2)
    beta14(3*(i-1)+3,3*(j-1)+3)=beta1(3*(i-1)+3,rk)*X1I(rk,3*(j-1)+3)&
                               +beta1(3*(i-1)+3,hk)*X1I(hk,3*(j-1)+3)

    beta18(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*beta4I(nk,3*(j-1)+1)&
                               +beta33(3*(i-1)+1,sk)*beta4I(sk,3*(j-1)+1)
    beta18(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*beta4I(mk,3*(j-1)+2)&
                               +beta33(3*(i-1)+2,gk)*beta4I(gk,3*(j-1)+2)
    beta18(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*beta4I(rk,3*(j-1)+3)&
                               +beta33(3*(i-1)+3,hk)*beta4I(hk,3*(j-1)+3)

    beta19(3*(i-1)+1,3*(j-1)+1)=beta1(3*(i-1)+1,nk)*beta4I(nk,3*(j-1)+1)&
                               +beta1(3*(i-1)+1,sk)*beta4I(sk,3*(j-1)+1)
    beta19(3*(i-1)+2,3*(j-1)+2)=beta1(3*(i-1)+2,mk)*beta4I(mk,3*(j-1)+2)&
                               +beta1(3*(i-1)+2,gk)*beta4I(gk,3*(j-1)+2)
    beta19(3*(i-1)+3,3*(j-1)+3)=beta1(3*(i-1)+3,rk)*beta4I(rk,3*(j-1)+3)&
                               +beta1(3*(i-1)+3,hk)*beta4I(hk,3*(j-1)+3)

    beta24(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*beta5I(nk,3*(j-1)+1)&
                               +beta33(3*(i-1)+1,sk)*beta5I(sk,3*(j-1)+1)
    beta24(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*beta5I(mk,3*(j-1)+2)&
                               +beta33(3*(i-1)+2,gk)*beta5I(gk,3*(j-1)+2)
    beta24(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*beta5I(rk,3*(j-1)+3)&
                               +beta33(3*(i-1)+3,hk)*beta5I(hk,3*(j-1)+3)

    beta25(3*(i-1)+1,3*(j-1)+1)=beta1(3*(i-1)+1,nk)*beta5I(nk,3*(j-1)+1)&
                               +beta1(3*(i-1)+1,sk)*beta5I(sk,3*(j-1)+1)
    beta25(3*(i-1)+2,3*(j-1)+2)=beta1(3*(i-1)+2,mk)*beta5I(mk,3*(j-1)+2)&
                               +beta1(3*(i-1)+2,gk)*beta5I(gk,3*(j-1)+2)
    beta25(3*(i-1)+3,3*(j-1)+3)=beta1(3*(i-1)+3,rk)*beta5I(rk,3*(j-1)+3)&
                               +beta1(3*(i-1)+3,hk)*beta5I(hk,3*(j-1)+3)

endif

  enddo

if (i_kk.eq.j_kk) then

  beta9(3*(i-1)+1,nk)=temp1
  beta9(3*(i-1)+2,mk)=temp2 
  beta9(3*(i-1)+3,rk)=temp3 
 
  beta16(3*(i-1)+1,nk)=temp10
  beta16(3*(i-1)+2,mk)=temp11
  beta16(3*(i-1)+3,rk)=temp12

  beta17(3*(i-1)+1,nk)=temp19
  beta17(3*(i-1)+2,mk)=temp20
  beta17(3*(i-1)+3,rk)=temp21

endif

if (i_kk.ne.j_kk) then

  beta9(3*(i-1)+1,nk)=temp4
  beta9(3*(i-1)+2,mk)=temp5
  beta9(3*(i-1)+3,rk)=temp6
  beta9(3*(i-1)+1,sk)=temp7
  beta9(3*(i-1)+2,gk)=temp8
  beta9(3*(i-1)+3,hk)=temp9

  beta16(3*(i-1)+1,nk)=temp13
  beta16(3*(i-1)+2,mk)=temp14
  beta16(3*(i-1)+3,rk)=temp15
  beta16(3*(i-1)+1,sk)=temp16
  beta16(3*(i-1)+2,gk)=temp17
  beta16(3*(i-1)+3,hk)=temp18

  beta17(3*(i-1)+1,nk)=temp22
  beta17(3*(i-1)+2,mk)=temp23
  beta17(3*(i-1)+3,rk)=temp24
  beta17(3*(i-1)+1,sk)=temp25
  beta17(3*(i-1)+2,gk)=temp26
  beta17(3*(i-1)+3,hk)=temp27

endif

enddo

!tau17
beta15=ZERO
!tau8
beta20=ZERO
beta21=ZERO
beta26=ZERO

do i=1,n

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
  temp17=ZERO
  temp18=ZERO
  temp19=ZERO
  temp20=ZERO
  temp21=ZERO
  temp22=ZERO
  temp23=ZERO
  temp24=ZERO
  temp25=ZERO
  temp26=ZERO
  temp27=ZERO
  temp28=ZERO
  temp29=ZERO
  temp30=ZERO
  temp31=ZERO
  temp32=ZERO
  temp33=ZERO
  temp34=ZERO
  temp35=ZERO
  temp36=ZERO

  do j=1,n
    
if (i_ll.eq.j_ll) then

    temp1=temp1+beta14(3*(i-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl)
    temp2=temp2+beta14(3*(i-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml)
    temp3=temp3+beta14(3*(i-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl)

    temp10=temp10+beta19(3*(i-1)+1,3*(j-1)+1)*beta6(3*(j-1)+1,nl)
    temp11=temp11+beta19(3*(i-1)+2,3*(j-1)+2)*beta6(3*(j-1)+2,ml)
    temp12=temp12+beta19(3*(i-1)+3,3*(j-1)+3)*beta6(3*(j-1)+3,rl)

    temp19=temp19+beta4I(3*(i-1)+1,3*(j-1)+1)*beta6(3*(j-1)+1,nl)
    temp20=temp20+beta4I(3*(i-1)+2,3*(j-1)+2)*beta6(3*(j-1)+2,ml)
    temp21=temp21+beta4I(3*(i-1)+3,3*(j-1)+3)*beta6(3*(j-1)+3,rl)

    temp28=temp28+beta25(3*(i-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl)
    temp29=temp29+beta25(3*(i-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml)
    temp30=temp30+beta25(3*(i-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl)

endif

if (i_ll.ne.j_ll) then

    temp4=temp4+beta14(3*(i-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl)
    temp5=temp5+beta14(3*(i-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml)
    temp6=temp6+beta14(3*(i-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl)
    temp7=temp7+beta14(3*(i-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,sl)
    temp8=temp8+beta14(3*(i-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,gl)
    temp9=temp9+beta14(3*(i-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,hl)

    temp13=temp13+beta19(3*(i-1)+1,3*(j-1)+1)*beta6(3*(j-1)+1,nl)
    temp14=temp14+beta19(3*(i-1)+2,3*(j-1)+2)*beta6(3*(j-1)+2,ml)
    temp15=temp15+beta19(3*(i-1)+3,3*(j-1)+3)*beta6(3*(j-1)+3,rl)
    temp16=temp16+beta19(3*(i-1)+1,3*(j-1)+1)*beta6(3*(j-1)+1,sl)
    temp17=temp17+beta19(3*(i-1)+2,3*(j-1)+2)*beta6(3*(j-1)+2,gl)
    temp18=temp18+beta19(3*(i-1)+3,3*(j-1)+3)*beta6(3*(j-1)+3,hl)

    temp22=temp22+beta4I(3*(i-1)+1,3*(j-1)+1)*beta6(3*(j-1)+1,nl)
    temp23=temp23+beta4I(3*(i-1)+2,3*(j-1)+2)*beta6(3*(j-1)+2,ml)
    temp24=temp24+beta4I(3*(i-1)+3,3*(j-1)+3)*beta6(3*(j-1)+3,rl)
    temp25=temp25+beta4I(3*(i-1)+1,3*(j-1)+1)*beta6(3*(j-1)+1,sl)
    temp26=temp26+beta4I(3*(i-1)+2,3*(j-1)+2)*beta6(3*(j-1)+2,gl)
    temp27=temp27+beta4I(3*(i-1)+3,3*(j-1)+3)*beta6(3*(j-1)+3,hl)

    temp31=temp31+beta25(3*(i-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl)
    temp32=temp32+beta25(3*(i-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml)
    temp33=temp33+beta25(3*(i-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl)
    temp34=temp34+beta25(3*(i-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,sl)
    temp35=temp35+beta25(3*(i-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,gl)
    temp36=temp36+beta25(3*(i-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,hl)

endif

  enddo

if (i_ll.eq.j_ll) then

  beta15(3*(i-1)+1,nl)=temp1
  beta15(3*(i-1)+2,ml)=temp2
  beta15(3*(i-1)+3,rl)=temp3

  beta20(3*(i-1)+1,nl)=temp10
  beta20(3*(i-1)+2,ml)=temp11
  beta20(3*(i-1)+3,rl)=temp12

  beta21(3*(i-1)+1,nl)=temp19
  beta21(3*(i-1)+2,ml)=temp20
  beta21(3*(i-1)+3,rl)=temp21

  beta26(3*(i-1)+1,nl)=temp28
  beta26(3*(i-1)+2,ml)=temp29
  beta26(3*(i-1)+3,rl)=temp30

endif

if (i_ll.ne.j_ll) then

  beta15(3*(i-1)+1,nl)=temp4
  beta15(3*(i-1)+2,ml)=temp5
  beta15(3*(i-1)+3,rl)=temp6
  beta15(3*(i-1)+1,sl)=temp7
  beta15(3*(i-1)+2,gl)=temp8
  beta15(3*(i-1)+3,hl)=temp9

  beta20(3*(i-1)+1,nl)=temp13
  beta20(3*(i-1)+2,ml)=temp14
  beta20(3*(i-1)+3,rl)=temp15
  beta20(3*(i-1)+1,sl)=temp16
  beta20(3*(i-1)+2,gl)=temp17
  beta20(3*(i-1)+3,hl)=temp18

  beta21(3*(i-1)+1,nl)=temp22
  beta21(3*(i-1)+2,ml)=temp23
  beta21(3*(i-1)+3,rl)=temp24
  beta21(3*(i-1)+1,sl)=temp25
  beta21(3*(i-1)+2,gl)=temp26
  beta21(3*(i-1)+3,hl)=temp27

  beta26(3*(i-1)+1,nl)=temp31
  beta26(3*(i-1)+2,ml)=temp32
  beta26(3*(i-1)+3,rl)=temp33
  beta26(3*(i-1)+1,sl)=temp34
  beta26(3*(i-1)+2,gl)=temp35
  beta26(3*(i-1)+3,hl)=temp36

endif

enddo

!tau8
beta22=ZERO
beta23=ZERO
!tau10
beta27=ZERO

do i=1,n

  if (i_kk.eq.j_kk.and.i_ll.eq.j_ll) then

    beta22(3*(i-1)+1,nk)=beta21(3*(i-1)+1,nl)*beta1(nl,nk)
    beta22(3*(i-1)+2,mk)=beta21(3*(i-1)+2,ml)*beta1(ml,mk)
    beta22(3*(i-1)+3,rk)=beta21(3*(i-1)+3,rl)*beta1(rl,rk)

    beta23(3*(i-1)+1,nk)=beta6(3*(i-1)+1,nl)*beta1(nl,nk)
    beta23(3*(i-1)+2,mk)=beta6(3*(i-1)+2,ml)*beta1(ml,mk)
    beta23(3*(i-1)+3,rk)=beta6(3*(i-1)+3,rl)*beta1(rl,rk)

    beta27(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta6(nk,nl)
    beta27(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta6(mk,ml)
    beta27(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta6(rk,rl)

  endif

  if (i_kk.ne.j_kk.and.i_ll.ne.j_ll) then

    beta22(3*(i-1)+1,nk)=beta21(3*(i-1)+1,nl)*beta1(nl,nk)&
                        +beta21(3*(i-1)+1,sl)*beta1(sl,nk)
    beta22(3*(i-1)+2,mk)=beta21(3*(i-1)+2,ml)*beta1(ml,mk)&
                        +beta21(3*(i-1)+2,gl)*beta1(gl,mk)
    beta22(3*(i-1)+3,rk)=beta21(3*(i-1)+3,rl)*beta1(rl,rk)&
                        +beta21(3*(i-1)+3,hl)*beta1(hl,rk)
    beta22(3*(i-1)+1,sk)=beta21(3*(i-1)+1,nl)*beta1(nl,sk)&
                        +beta21(3*(i-1)+1,sl)*beta1(sl,sk)
    beta22(3*(i-1)+2,gk)=beta21(3*(i-1)+2,ml)*beta1(ml,gk)&
                        +beta21(3*(i-1)+2,gl)*beta1(gl,gk)
    beta22(3*(i-1)+3,hk)=beta21(3*(i-1)+3,rl)*beta1(rl,hk)&
                        +beta21(3*(i-1)+3,hl)*beta1(hl,hk)

    beta23(3*(i-1)+1,nk)=beta6(3*(i-1)+1,nl)*beta1(nl,nk)&
                        +beta6(3*(i-1)+1,sl)*beta1(sl,nk)
    beta23(3*(i-1)+2,mk)=beta6(3*(i-1)+2,ml)*beta1(ml,mk)&
                        +beta6(3*(i-1)+2,gl)*beta1(gl,mk)
    beta23(3*(i-1)+3,rk)=beta6(3*(i-1)+3,rl)*beta1(rl,rk)&
                        +beta6(3*(i-1)+3,hl)*beta1(hl,rk)
    beta23(3*(i-1)+1,sk)=beta6(3*(i-1)+1,nl)*beta1(nl,sk)&
                        +beta6(3*(i-1)+1,sl)*beta1(sl,sk)
    beta23(3*(i-1)+2,gk)=beta6(3*(i-1)+2,ml)*beta1(ml,gk)&
                        +beta6(3*(i-1)+2,gl)*beta1(gl,gk)
    beta23(3*(i-1)+3,hk)=beta6(3*(i-1)+3,rl)*beta1(rl,hk)&
                        +beta6(3*(i-1)+3,hl)*beta1(hl,hk)

    beta27(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta6(nk,nl)&
                        +beta1(3*(i-1)+1,sk)*beta6(sk,nl)
    beta27(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta6(mk,ml)&
                        +beta1(3*(i-1)+2,gk)*beta6(gk,ml)
    beta27(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta6(rk,rl)&
                        +beta1(3*(i-1)+3,hk)*beta6(hk,rl)
    beta27(3*(i-1)+1,sl)=beta1(3*(i-1)+1,nk)*beta6(nk,sl)&
                        +beta1(3*(i-1)+1,sk)*beta6(sk,sl)
    beta27(3*(i-1)+2,gl)=beta1(3*(i-1)+2,mk)*beta6(mk,gl)&
                        +beta1(3*(i-1)+2,gk)*beta6(gk,gl)
    beta27(3*(i-1)+3,hl)=beta1(3*(i-1)+3,rk)*beta6(rk,hl)&
                        +beta1(3*(i-1)+3,hk)*beta6(hk,hl)

  endif

  if (i_kk.eq.j_kk.and.i_ll.ne.j_ll) then

    beta22(3*(i-1)+1,nk)=beta21(3*(i-1)+1,nl)*beta1(nl,nk)&
                        +beta21(3*(i-1)+1,sl)*beta1(sl,nk)
    beta22(3*(i-1)+2,mk)=beta21(3*(i-1)+2,ml)*beta1(ml,mk)&
                        +beta21(3*(i-1)+2,gl)*beta1(gl,mk)
    beta22(3*(i-1)+3,rk)=beta21(3*(i-1)+3,rl)*beta1(rl,rk)&
                        +beta21(3*(i-1)+3,hl)*beta1(hl,rk)

    beta23(3*(i-1)+1,nk)=beta6(3*(i-1)+1,nl)*beta1(nl,nk)&
                        +beta6(3*(i-1)+1,sl)*beta1(sl,nk)
    beta23(3*(i-1)+2,mk)=beta6(3*(i-1)+2,ml)*beta1(ml,mk)&
                        +beta6(3*(i-1)+2,gl)*beta1(gl,mk)
    beta23(3*(i-1)+3,rk)=beta6(3*(i-1)+3,rl)*beta1(rl,rk)&
                        +beta6(3*(i-1)+3,hl)*beta1(hl,rk)

    beta27(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta6(nk,nl)
    beta27(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta6(mk,ml)
    beta27(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta6(rk,rl)
    beta27(3*(i-1)+1,sl)=beta1(3*(i-1)+1,nk)*beta6(nk,sl)
    beta27(3*(i-1)+2,gl)=beta1(3*(i-1)+2,mk)*beta6(mk,gl)
    beta27(3*(i-1)+3,hl)=beta1(3*(i-1)+3,rk)*beta6(rk,hl)

  endif

  if (i_kk.ne.j_kk.and.i_ll.eq.j_ll) then

    beta22(3*(i-1)+1,nk)=beta21(3*(i-1)+1,nl)*beta1(nl,nk)
    beta22(3*(i-1)+2,mk)=beta21(3*(i-1)+2,ml)*beta1(ml,mk)
    beta22(3*(i-1)+3,rk)=beta21(3*(i-1)+3,rl)*beta1(rl,rk)
    beta22(3*(i-1)+1,sk)=beta21(3*(i-1)+1,nl)*beta1(nl,sk)
    beta22(3*(i-1)+2,gk)=beta21(3*(i-1)+2,ml)*beta1(ml,gk)
    beta22(3*(i-1)+3,hk)=beta21(3*(i-1)+3,rl)*beta1(rl,hk)

    beta23(3*(i-1)+1,nk)=beta6(3*(i-1)+1,nl)*beta1(nl,nk)
    beta23(3*(i-1)+2,mk)=beta6(3*(i-1)+2,ml)*beta1(ml,mk)
    beta23(3*(i-1)+3,rk)=beta6(3*(i-1)+3,rl)*beta1(rl,rk)
    beta23(3*(i-1)+1,sk)=beta6(3*(i-1)+1,nl)*beta1(nl,sk)
    beta23(3*(i-1)+2,gk)=beta6(3*(i-1)+2,ml)*beta1(ml,gk)
    beta23(3*(i-1)+3,hk)=beta6(3*(i-1)+3,rl)*beta1(rl,hk)

    beta27(3*(i-1)+1,nl)=beta1(3*(i-1)+1,nk)*beta6(nk,nl)&
                        +beta1(3*(i-1)+1,sk)*beta6(sk,nl)
    beta27(3*(i-1)+2,ml)=beta1(3*(i-1)+2,mk)*beta6(mk,ml)&
                        +beta1(3*(i-1)+2,gk)*beta6(gk,ml)
    beta27(3*(i-1)+3,rl)=beta1(3*(i-1)+3,rk)*beta6(rk,rl)&
                        +beta1(3*(i-1)+3,hk)*beta6(hk,rl)

  endif
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

tau6=ZERO
tau7=ZERO
tau8=ZERO
tau9=ZERO
tau10=ZERO

do i=1,3*n

  tau6=tau6+beta7(i,i)
  tau7=tau7+beta13(i,i)
  tau8=tau8+beta20(i,i)
  tau9=tau9+beta24(i,i)
  tau10=tau10+beta27(i,i)

enddo

!write(*,*)'tau1,tau6,tau7,tau8,tau9,tau10',tau6,tau7,tau8,tau9,tau10

!************************
!
! Kinetic Energy Gradient
!
!************************

if (grad_k.or.grad_l) then

alpha1=ZERO
alpha2=ZERO
mu2=ZERO
dtau6=ZERO
dtau66=ZERO
dtau666=ZERO
dtau6666=ZERO
dtau6666M=ZERO
dtau7=ZERO
dtau77=ZERO
dtau777=ZERO
dtau7777=ZERO
dtau7777M=ZERO
dtau8=ZERO
dtau88=ZERO
dtau888=ZERO
dtau9=ZERO
dtau99=ZERO
dtau999=ZERO
dtau10=ZERO

do i=1,n
  do j=1,n

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
     temp17=ZERO
     temp18=ZERO

     do k=1,n
    
       temp16=temp16+X1(i,k)*inv_tAkl(k,j)
       temp17=temp17+beta5(i,k)*inv_tAkl(k,j)
       temp18=temp18+beta4(i,k)*Glob_MassMatrix(k,j)

       temp1=temp1+beta7(3*(i-1)+1,3*(k-1)+1)*inv_tAklI(3*(k-1)+1,3*(j-1)+1)
       temp2=temp2+beta7(3*(i-1)+2,3*(k-1)+2)*inv_tAklI(3*(k-1)+2,3*(j-1)+2)
       temp3=temp3+beta7(3*(i-1)+3,3*(k-1)+3)*inv_tAklI(3*(k-1)+3,3*(j-1)+3)

       temp4=temp4+beta12(3*(i-1)+1,3*(k-1)+1)*MI(3*(k-1)+1,3*(j-1)+1)
       temp5=temp5+beta12(3*(i-1)+2,3*(k-1)+2)*MI(3*(k-1)+2,3*(j-1)+2)
       temp6=temp6+beta12(3*(i-1)+3,3*(k-1)+3)*MI(3*(k-1)+3,3*(j-1)+3)

       temp7=temp7+beta13(3*(i-1)+1,3*(k-1)+1)*inv_tAklI(3*(k-1)+1,3*(j-1)+1)
       temp8=temp8+beta13(3*(i-1)+2,3*(k-1)+2)*inv_tAklI(3*(k-1)+2,3*(j-1)+2)
       temp9=temp9+beta13(3*(i-1)+3,3*(k-1)+3)*inv_tAklI(3*(k-1)+3,3*(j-1)+3)

       temp10=temp10+beta18(3*(i-1)+1,3*(k-1)+1)*MI(3*(k-1)+1,3*(j-1)+1)
       temp11=temp11+beta18(3*(i-1)+2,3*(k-1)+2)*MI(3*(k-1)+2,3*(j-1)+2)
       temp12=temp12+beta18(3*(i-1)+3,3*(k-1)+3)*MI(3*(k-1)+3,3*(j-1)+3)

       temp13=temp13+beta24(3*(i-1)+1,3*(k-1)+1)*inv_tAklI(3*(k-1)+1,3*(j-1)+1)
       temp14=temp14+beta24(3*(i-1)+2,3*(k-1)+2)*inv_tAklI(3*(k-1)+2,3*(j-1)+2)
       temp15=temp15+beta24(3*(i-1)+3,3*(k-1)+3)*inv_tAklI(3*(k-1)+3,3*(j-1)+3)

     enddo 

     alpha1(i,j)=temp16
     alpha2(i,j)=temp17
     mu2(i,j)=temp18

     dtau6(3*(i-1)+1,3*(j-1)+1)=temp1
     dtau6(3*(i-1)+2,3*(j-1)+2)=temp2
     dtau6(3*(i-1)+3,3*(j-1)+3)=temp3

     dtau6666M(3*(i-1)+1,3*(j-1)+1)=temp4
     dtau6666M(3*(i-1)+2,3*(j-1)+2)=temp5
     dtau6666M(3*(i-1)+3,3*(j-1)+3)=temp6

     dtau7(3*(i-1)+1,3*(j-1)+1)=temp7
     dtau7(3*(i-1)+2,3*(j-1)+2)=temp8
     dtau7(3*(i-1)+3,3*(j-1)+3)=temp9

     dtau7777M(3*(i-1)+1,3*(j-1)+1)=temp10
     dtau7777M(3*(i-1)+2,3*(j-1)+2)=temp11
     dtau7777M(3*(i-1)+3,3*(j-1)+3)=temp12

     dtau9(3*(i-1)+1,3*(j-1)+1)=temp13
     dtau9(3*(i-1)+2,3*(j-1)+2)=temp14
     dtau9(3*(i-1)+3,3*(j-1)+3)=temp15

if (i_kk.ne.j_kk) then

     dtau66(3*(i-1)+1,3*(j-1)+1)=beta9(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)&
                                +beta9(3*(i-1)+1,sk)*inv_tAklI(sk,3*(j-1)+1)
     dtau66(3*(i-1)+2,3*(j-1)+2)=beta9(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)&
                                +beta9(3*(i-1)+2,gk)*inv_tAklI(gk,3*(j-1)+2)
     dtau66(3*(i-1)+3,3*(j-1)+3)=beta9(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)&
                                +beta9(3*(i-1)+3,hk)*inv_tAKlI(hk,3*(j-1)+3)

     dtau777(3*(i-1)+1,3*(j-1)+1)=beta16(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)&
                                 +beta16(3*(i-1)+1,sk)*inv_tAklI(sk,3*(j-1)+1)
     dtau777(3*(i-1)+2,3*(j-1)+2)=beta16(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)&
                                 +beta16(3*(i-1)+2,gk)*inv_tAklI(gk,3*(j-1)+2)
     dtau777(3*(i-1)+3,3*(j-1)+3)=beta16(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)&
                                 +beta16(3*(i-1)+3,hk)*inv_tAKlI(hk,3*(j-1)+3)

     dtau7777(3*(i-1)+1,3*(j-1)+1)=beta17(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)&
                                  +beta17(3*(i-1)+1,sk)*inv_tAklI(sk,3*(j-1)+1)
     dtau7777(3*(i-1)+2,3*(j-1)+2)=beta17(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)&
                                  +beta17(3*(i-1)+2,gk)*inv_tAklI(gk,3*(j-1)+2)
     dtau7777(3*(i-1)+3,3*(j-1)+3)=beta17(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)&
                                  +beta17(3*(i-1)+3,hk)*inv_tAKlI(hk,3*(j-1)+3)

     dtau88(3*(i-1)+1,3*(j-1)+1)=beta22(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)&
                                +beta22(3*(i-1)+1,sk)*inv_tAklI(sk,3*(j-1)+1)
     dtau88(3*(i-1)+2,3*(j-1)+2)=beta22(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)&
                                +beta22(3*(i-1)+2,gk)*inv_tAklI(gk,3*(j-1)+2)
     dtau88(3*(i-1)+3,3*(j-1)+3)=beta22(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)&
                                +beta22(3*(i-1)+3,hk)*inv_tAKlI(hk,3*(j-1)+3)

     dtau888(3*(i-1)+1,3*(j-1)+1)=beta23(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)&
                                 +beta23(3*(i-1)+1,sk)*inv_tAklI(sk,3*(j-1)+1)
     dtau888(3*(i-1)+2,3*(j-1)+2)=beta23(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)&
                                 +beta23(3*(i-1)+2,gk)*inv_tAklI(gk,3*(j-1)+2)
     dtau888(3*(i-1)+3,3*(j-1)+3)=beta23(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)&
                                 +beta23(3*(i-1)+3,hk)*inv_tAKlI(hk,3*(j-1)+3)

     dtau999(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*MI(nk,3*(j-1)+1)&
                                 +beta33(3*(i-1)+1,sk)*MI(sk,3*(j-1)+1)
     dtau999(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*MI(mk,3*(j-1)+2)&
                                 +beta33(3*(i-1)+2,gk)*MI(gk,3*(j-1)+2)
     dtau999(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*MI(rk,3*(j-1)+3)&
                                 +beta33(3*(i-1)+3,hk)*MI(hk,3*(j-1)+3)

endif

if (i_kk.eq.j_kk) then

     dtau66(3*(i-1)+1,3*(j-1)+1)=beta9(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)
     dtau66(3*(i-1)+2,3*(j-1)+2)=beta9(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)
     dtau66(3*(i-1)+3,3*(j-1)+3)=beta9(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)

     dtau777(3*(i-1)+1,3*(j-1)+1)=beta16(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)
     dtau777(3*(i-1)+2,3*(j-1)+2)=beta16(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)
     dtau777(3*(i-1)+3,3*(j-1)+3)=beta16(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)

     dtau7777(3*(i-1)+1,3*(j-1)+1)=beta17(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)
     dtau7777(3*(i-1)+2,3*(j-1)+2)=beta17(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)
     dtau7777(3*(i-1)+3,3*(j-1)+3)=beta17(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)

     dtau88(3*(i-1)+1,3*(j-1)+1)=beta22(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)
     dtau88(3*(i-1)+2,3*(j-1)+2)=beta22(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)
     dtau88(3*(i-1)+3,3*(j-1)+3)=beta22(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)

     dtau888(3*(i-1)+1,3*(j-1)+1)=beta23(3*(i-1)+1,nk)*inv_tAklI(nk,3*(j-1)+1)
     dtau888(3*(i-1)+2,3*(j-1)+2)=beta23(3*(i-1)+2,mk)*inv_tAklI(mk,3*(j-1)+2)
     dtau888(3*(i-1)+3,3*(j-1)+3)=beta23(3*(i-1)+3,rk)*inv_tAklI(rk,3*(j-1)+3)

     dtau999(3*(i-1)+1,3*(j-1)+1)=beta33(3*(i-1)+1,nk)*MI(nk,3*(j-1)+1)
     dtau999(3*(i-1)+2,3*(j-1)+2)=beta33(3*(i-1)+2,mk)*MI(mk,3*(j-1)+2)
     dtau999(3*(i-1)+3,3*(j-1)+3)=beta33(3*(i-1)+3,rk)*MI(rk,3*(j-1)+3)

endif

if (i_ll.ne.j_ll) then

     dtau666(3*(i-1)+1,3*(j-1)+1)=beta10(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)&
                                 +beta10(3*(i-1)+1,sl)*inv_tAklI(sl,3*(j-1)+1)
     dtau666(3*(i-1)+2,3*(j-1)+2)=beta10(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)&
                                 +beta10(3*(i-1)+2,gl)*inv_tAklI(gl,3*(j-1)+2)
     dtau666(3*(i-1)+3,3*(j-1)+3)=beta10(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)&
                                 +beta10(3*(i-1)+3,hl)*inv_tAKlI(hl,3*(j-1)+3)

     dtau6666(3*(i-1)+1,3*(j-1)+1)=beta11(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)&
                                  +beta11(3*(i-1)+1,sl)*inv_tAklI(sl,3*(j-1)+1)
     dtau6666(3*(i-1)+2,3*(j-1)+2)=beta11(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)&
                                  +beta11(3*(i-1)+2,gl)*inv_tAklI(gl,3*(j-1)+2)
     dtau6666(3*(i-1)+3,3*(j-1)+3)=beta11(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)&
                                  +beta11(3*(i-1)+3,hl)*inv_tAKlI(hl,3*(j-1)+3)

     dtau77(3*(i-1)+1,3*(j-1)+1)=beta15(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)&
                                +beta15(3*(i-1)+1,sl)*inv_tAklI(sl,3*(j-1)+1)
     dtau77(3*(i-1)+2,3*(j-1)+2)=beta15(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)&
                                +beta15(3*(i-1)+2,gl)*inv_tAklI(gl,3*(j-1)+2)
     dtau77(3*(i-1)+3,3*(j-1)+3)=beta15(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)&
                                +beta15(3*(i-1)+3,hl)*inv_tAKlI(hl,3*(j-1)+3)

     dtau8(3*(i-1)+1,3*(j-1)+1)=beta20(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)&
                               +beta20(3*(i-1)+1,sl)*inv_tAklI(sl,3*(j-1)+1)
     dtau8(3*(i-1)+2,3*(j-1)+2)=beta20(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)&
                               +beta20(3*(i-1)+2,gl)*inv_tAklI(gl,3*(j-1)+2)
     dtau8(3*(i-1)+3,3*(j-1)+3)=beta20(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)&
                               +beta20(3*(i-1)+3,hl)*inv_tAKlI(hl,3*(j-1)+3)

     dtau99(3*(i-1)+1,3*(j-1)+1)=beta26(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)&
                                +beta26(3*(i-1)+1,sl)*inv_tAklI(sl,3*(j-1)+1)
     dtau99(3*(i-1)+2,3*(j-1)+2)=beta26(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)&
                                +beta26(3*(i-1)+2,gl)*inv_tAklI(gl,3*(j-1)+2)
     dtau99(3*(i-1)+3,3*(j-1)+3)=beta26(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)&
                                +beta26(3*(i-1)+3,hl)*inv_tAKlI(hl,3*(j-1)+3)

     dtau10(3*(i-1)+1,3*(j-1)+1)=beta27(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)&
                                +beta27(3*(i-1)+1,sl)*inv_tAklI(sl,3*(j-1)+1)
     dtau10(3*(i-1)+2,3*(j-1)+2)=beta27(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)&
                                +beta27(3*(i-1)+2,gl)*inv_tAklI(gl,3*(j-1)+2)
     dtau10(3*(i-1)+3,3*(j-1)+3)=beta27(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)&
                                +beta27(3*(i-1)+3,hl)*inv_tAKlI(hl,3*(j-1)+3)

endif

if (i_ll.eq.j_ll) then

     dtau666(3*(i-1)+1,3*(j-1)+1)=beta10(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)
     dtau666(3*(i-1)+2,3*(j-1)+2)=beta10(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)
     dtau666(3*(i-1)+3,3*(j-1)+3)=beta10(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)

     dtau6666(3*(i-1)+1,3*(j-1)+1)=beta11(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)
     dtau6666(3*(i-1)+2,3*(j-1)+2)=beta11(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)
     dtau6666(3*(i-1)+3,3*(j-1)+3)=beta11(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)

     dtau77(3*(i-1)+1,3*(j-1)+1)=beta15(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)
     dtau77(3*(i-1)+2,3*(j-1)+2)=beta15(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)
     dtau77(3*(i-1)+3,3*(j-1)+3)=beta15(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)

     dtau8(3*(i-1)+1,3*(j-1)+1)=beta20(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)
     dtau8(3*(i-1)+2,3*(j-1)+2)=beta20(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)
     dtau8(3*(i-1)+3,3*(j-1)+3)=beta20(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)

     dtau99(3*(i-1)+1,3*(j-1)+1)=beta26(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)
     dtau99(3*(i-1)+2,3*(j-1)+2)=beta26(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)
     dtau99(3*(i-1)+3,3*(j-1)+3)=beta26(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)

     dtau10(3*(i-1)+1,3*(j-1)+1)=beta27(3*(i-1)+1,nl)*inv_tAklI(nl,3*(j-1)+1)
     dtau10(3*(i-1)+2,3*(j-1)+2)=beta27(3*(i-1)+2,ml)*inv_tAklI(ml,3*(j-1)+2)
     dtau10(3*(i-1)+3,3*(j-1)+3)=beta27(3*(i-1)+3,rl)*inv_tAklI(rl,3*(j-1)+3)

endif

  enddo
enddo

endif

if (grad_k) then

 !resultX1k=alpha1*Lk
 !result4k=alpha2*Lk


  resultX1k=ZERO
  indx=ZERO

  do j=1,n
    do i=j,n
      indx=indx+1

      temp1=ZERO
      temp2=ZERO

      do k=1,n

        temp1=temp1+(alpha1(i,k)+alpha1(k,i))*Lk(k,j)
        temp2=temp2+(alpha2(i,k)+alpha2(k,i))*Lk(k,j)

      enddo

      resultX1k(indx)=-temp1+temp2

    enddo
  enddo

  dtau6k=ZERO
  dtau7k=ZERO
  dtau8k=ZERO
  dtau9k=ZERO
  dtau10k=ZERO

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
          
        temp1=temp1+(dtau6(i,k)+dtau6(k,i))*LkI(k,j)
        temp2=temp2+(dtau66(i,k)+dtau66(k,i))*LkI(k,j)
        temp3=temp3+(dtau666(i,k)+dtau666(k,i))*LkI(k,j)
        temp4=temp4+(dtau6666(i,k)+dtau6666(k,i))*LkI(k,j)

        temp5=temp5+(dtau7(i,k)+dtau7(k,i))*LkI(k,j)
        temp6=temp6+(dtau77(i,k)+dtau77(k,i))*LkI(k,j)
        temp7=temp7+(dtau777(i,k)+dtau777(k,i))*LkI(k,j)
        temp8=temp8+(dtau7777(i,k)+dtau7777(k,i))*LkI(k,j)

        temp9=temp9+(dtau8(i,k)+dtau8(k,i))*LkI(k,j)
        temp10=temp10+(dtau88(i,k)+dtau88(k,i))*LkI(k,j)
        temp11=temp11+(dtau888(i,k)+dtau888(k,i))*LkI(k,j)

        temp12=temp12+(dtau9(i,k)+dtau9(k,i))*LkI(k,j)
        temp13=temp13+(dtau99(i,k)+dtau99(k,i))*LkI(k,j)

        temp14=temp14+(dtau10(i,k)+dtau10(k,i))*LkI(k,j)

      enddo
      indx=indx+1

      dtau6k(indx)=-temp1-temp2-temp3+temp4
      dtau7k(indx)=-temp5-temp6-temp7+temp8
      dtau8k(indx)=-temp9-temp10+temp11
      dtau9k(indx)=-temp12-temp13
      dtau10k(indx)=-temp14

    enddo
  enddo
endif

if (grad_l) then

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
       
        temp1=temp1+pP(i,k)*dtau6(k,j)
        temp2=temp2+pP(i,k)*dtau66(k,j)
        temp3=temp3+pP(i,k)*dtau666(k,j)
        temp4=temp4+pP(i,k)*dtau6666M(k,j)
        temp5=temp5+pP(i,k)*dtau7(k,j)
        temp6=temp6+pP(i,k)*dtau77(k,j)
        temp7=temp7+pP(i,k)*dtau777(k,j)
        temp8=temp8+pP(i,k)*dtau7777M(k,j)
        temp9=temp9+pP(i,k)*dtau8(k,j)
        temp10=temp10+pP(i,k)*dtau88(k,j)
        temp11=temp11+pP(i,k)*dtau9(k,j)
        temp12=temp12+pP(i,k)*dtau99(k,j)
        temp13=temp13+pP(i,k)*dtau999(k,j)
        temp14=temp14+pP(i,k)*dtau10(k,j)

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

      dtau6p(i,j)=temp1
      dtau66p(i,j)=temp2
      dtau666p(i,j)=temp3
      dtau6666Mp(i,j)=temp4

      dtau7p(i,j)=temp5
      dtau77p(i,j)=temp6
      dtau777p(i,j)=temp7
      dtau7777Mp(i,j)=temp8

      dtau8p(i,j)=temp9
      dtau88p(i,j)=temp10

      dtau9p(i,j)=temp11
      dtau99p(i,j)=temp12
      dtau999p(i,j)=temp13
 
      dtau10p(i,j)=temp14
 
    enddo
  enddo

  Z1=ZERO
  Z2=ZERO

  do i=1,n
    do j=1,n
      temp1=ZERO
      temp2=ZERO
      do k=1,n
        temp1=temp1+P(i,k)*mu2(k,j)
        temp2=temp2+P(i,k)*alpha1(k,j)
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
      mu2p(i,j)=temp1
      alpha1p(i,j)=temp2
    enddo
  enddo  

  resultX1l=ZERO
  indx=ZERO

  do j=1,n
    do i=j,n

      indx=indx+1
      temp1=ZERO
      temp2=ZERO

      do k=j,n

        temp1=temp1+(alpha1p(i,k)+alpha1p(k,i))*Ll(k,j)
        temp2=temp2+(mu2p(i,k)+mu2p(k,i))*Ll(k,j)

      enddo

      resultX1l(indx)=-temp1+temp2

    enddo
  enddo

  dtau6l=ZERO
  dtau7l=ZERO
  dtau8l=ZERO
  dtau9l=ZERO
  dtau10l=ZERO

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
          
        temp1=temp1+(dtau6p(i,k)+dtau6p(k,i))*LlI(k,j)
        temp2=temp2+(dtau66p(i,k)+dtau66p(k,i))*LlI(k,j)
        temp3=temp3+(dtau666p(i,k)+dtau666p(k,i))*LlI(k,j)
        temp4=temp4+(dtau6666Mp(i,k)+dtau6666Mp(k,i))*LlI(k,j)
        temp5=temp5+(dtau7p(i,k)+dtau7p(k,i))*LlI(k,j)
        temp6=temp6+(dtau77p(i,k)+dtau77p(k,i))*LlI(k,j)
        temp7=temp7+(dtau777p(i,k)+dtau777p(k,i))*LlI(k,j)
        temp8=temp8+(dtau7777Mp(i,k)+dtau7777Mp(k,i))*LlI(k,j)
        temp9=temp9+(dtau8p(i,k)+dtau8p(k,i))*LlI(k,j)
        temp10=temp10+(dtau88p(i,k)+dtau88p(k,i))*LlI(k,j)
        temp11=temp11+(dtau9p(i,k)+dtau9p(k,i))*LlI(k,j)
        temp12=temp12+(dtau99p(i,k)+dtau99p(k,i))*LlI(k,j)
        temp13=temp13+(dtau999p(i,k)+dtau999p(k,i))*LlI(k,j)
        temp14=temp14+(dtau10p(i,k)+dtau10p(k,i))*LlI(k,j)
   
       enddo

      indx=indx+1

      dtau6l(indx)=-temp1-temp2-temp3+temp4
      dtau7l(indx)=-temp5-temp6-temp7+temp8
      dtau8l(indx)=-temp9-temp10
      dtau9l(indx)=-temp11-temp12+temp13
      dtau10l(indx)=-temp14

    enddo
  enddo  
endif

!***************
!
! Kinetic Energy
!
!***************

Tkl=(PI**(THREE*n/TWO))&
*(det_tAkl**(-THREEHALF))&
*(THREE*tau1*nu3+TWO*(tau6+tau7-tau8-tau9+tau10))

!write(*,*)'Tkl',Tkl

!************************
!
! Kinetic Energy Gradient
!
!************************

if (grad_k) then

  do i=1,n*(n+1)/2

    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    temp4=ZERO
    temp5=ZERO
    temp6=ZERO

    do k=1,3*n*(3*n+1)/2

!dnu3
      temp1=temp1+dnu3k(k)*trans(k,i)
!dtau6
      temp2=temp2+dtau6k(k)*trans(k,i)
!dtau7
      temp3=temp3+dtau7k(k)*trans(k,i)
!dtau8
      temp4=temp4+dtau8k(k)*trans(k,i)
!dtau9
      temp5=temp5+dtau9k(k)*trans(k,i)
!dtau10
      temp6=temp6+dtau10k(k)*trans(k,i)

    enddo

     Dk(i)=Dk(i)+(PI**(THREE*n/TWO))&    
                 *(det_tAkl**(-THREEHALF))&
                 *(-THREEHALF*resultAk(i)*(THREE*tau1*nu3+TWO*(tau6+tau7-tau8-tau9+tau10))&
                   +THREE*resultX1k(i)*nu3&
                   -THREE*tau1*temp1&
                   +TWO*(temp2+temp3-temp4-temp5+temp6))

enddo

!write(*,*)'Dk',Dk

endif

if (grad_l) then

  do i=1,n*(n+1)/2

    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    temp4=ZERO
    temp5=ZERO
    temp6=ZERO

    do k=1,3*n*(3*n+1)/2
!dnu3
      temp1=temp1+dnu3l(k)*trans(k,i)
!dtau6
      temp2=temp2+dtau6l(k)*trans(k,i)
!dtau7
      temp3=temp3+dtau7l(k)*trans(k,i)
!dtau8
      temp4=temp4+dtau8l(k)*trans(k,i)
!dtau9=ZERO
      temp5=temp5+dtau9l(k)*trans(k,i)
!dtau10
      temp6=temp6+dtau10l(k)*trans(k,i)

    enddo

      Dl(i)=Dl(i)+(PI**(THREE*n/TWO))&
                 *(det_tAkl**(-THREEHALF))&
                 *(-THREEHALF*resultAl(i)*(THREE*tau1*nu3+TWO*(tau6+tau7-tau8-tau9+tau10))&
                   +THREE*resultX1l(i)*nu3&
                   -THREE*tau1*temp1&
                   +TWO*(temp2+temp3-temp4-temp5+temp6))

 enddo

!write(*,*)'Dl',Dl

endif

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

!write(*,*)'Jm'
!do r=1,n
!  write(*,9292)(Jm(r,t),t=1,n)
!enddo

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

!write(*,*)'betaJm'
!do r=1,n
!  write(*,9292)(betaJm(r,t),t=1,n)
!enddo

! Preforming Kronecker Product of betaJm with I3

betaJmI=ZERO

do r=1,n  

    betaJmI(3*(r-1)+1,3*(i-1)+1)=betaJm(r,i)
    betaJmI(3*(r-1)+2,3*(i-1)+2)=betaJm(r,i)
    betaJmI(3*(r-1)+3,3*(i-1)+3)=betaJm(r,i)
    betaJmI(3*(r-1)+1,3*(j-1)+1)=betaJm(r,j)
    betaJmI(3*(r-1)+2,3*(j-1)+2)=betaJm(r,j)
    betaJmI(3*(r-1)+3,3*(j-1)+3)=betaJm(r,j)

 enddo

!write(*,*)'betaJmI'
!do r=1,3*n
!  write(*,9292)(betaJmI(r,t),t=1,3*n)
!enddo

beta28=ZERO
beta29=ZERO
beta31=ZERO
beta38=ZERO

if (i_ll.eq.j_ll) then

do r=1,n

  beta28(3*(r-1)+1,3*(i-1)+1)=beta3(3*(r-1)+1,nl)*betaJmI(nl,3*(i-1)+1)
  beta28(3*(r-1)+2,3*(i-1)+2)=beta3(3*(r-1)+2,ml)*betaJmI(ml,3*(i-1)+2)
  beta28(3*(r-1)+3,3*(i-1)+3)=beta3(3*(r-1)+3,rl)*betaJmI(rl,3*(i-1)+3)         
  beta28(3*(r-1)+1,3*(j-1)+1)=beta3(3*(r-1)+1,nl)*betaJmI(nl,3*(j-1)+1)
  beta28(3*(r-1)+2,3*(j-1)+2)=beta3(3*(r-1)+2,ml)*betaJmI(ml,3*(j-1)+2)
  beta28(3*(r-1)+3,3*(j-1)+3)=beta3(3*(r-1)+3,rl)*betaJmI(rl,3*(j-1)+3)         

  beta29(3*(r-1)+1,3*(i-1)+1)=beta2(3*(r-1)+1,nl)*betaJmI(nl,3*(i-1)+1)         
  beta29(3*(r-1)+2,3*(i-1)+2)=beta2(3*(r-1)+2,ml)*betaJmI(ml,3*(i-1)+2)
  beta29(3*(r-1)+3,3*(i-1)+3)=beta2(3*(r-1)+3,rl)*betaJmI(rl,3*(i-1)+3)
  beta29(3*(r-1)+1,3*(j-1)+1)=beta2(3*(r-1)+1,nl)*betaJmI(nl,3*(j-1)+1)         
  beta29(3*(r-1)+2,3*(j-1)+2)=beta2(3*(r-1)+2,ml)*betaJmI(ml,3*(j-1)+2)
  beta29(3*(r-1)+3,3*(j-1)+3)=beta2(3*(r-1)+3,rl)*betaJmI(rl,3*(j-1)+3)

  beta31(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta3(3*(i-1)+1,nl)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,nl)
  beta31(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta3(3*(i-1)+2,ml)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,ml)
  beta31(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta3(3*(i-1)+3,rl)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,rl)

  beta38(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl) 
  beta38(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml) 
  beta38(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl) 

  if (i.eq.j) then

     beta31(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta3(3*(i-1)+1,nl)
     beta31(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta3(3*(i-1)+2,ml)
     beta31(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta3(3*(i-1)+3,rl)

     beta38(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)
     beta38(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)
     beta38(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)

  endif
enddo
endif

if (i_ll.ne.j_ll) then

do r=1,n

  beta28(3*(r-1)+1,3*(i-1)+1)=beta3(3*(r-1)+1,nl)*betaJmI(nl,3*(i-1)+1)&
                             +beta3(3*(r-1)+1,sl)*betaJmI(sl,3*(i-1)+1)
  beta28(3*(r-1)+2,3*(i-1)+2)=beta3(3*(r-1)+2,ml)*betaJmI(ml,3*(i-1)+2)&
                             +beta3(3*(r-1)+2,gl)*betaJmI(gl,3*(i-1)+2) 
  beta28(3*(r-1)+3,3*(i-1)+3)=beta3(3*(r-1)+3,rl)*betaJmI(rl,3*(i-1)+3)&
                             +beta3(3*(r-1)+3,hl)*betaJmI(hl,3*(i-1)+3) 
  beta28(3*(r-1)+1,3*(j-1)+1)=beta3(3*(r-1)+1,nl)*betaJmI(nl,3*(j-1)+1)&
                             +beta3(3*(r-1)+1,sl)*betaJmI(sl,3*(j-1)+1)
  beta28(3*(r-1)+2,3*(j-1)+2)=beta3(3*(r-1)+2,ml)*betaJmI(ml,3*(j-1)+2)&
                             +beta3(3*(r-1)+2,gl)*betaJmI(gl,3*(j-1)+2)
  beta28(3*(r-1)+3,3*(j-1)+3)=beta3(3*(r-1)+3,rl)*betaJmI(rl,3*(j-1)+3)&
                             +beta3(3*(r-1)+3,hl)*betaJmI(hl,3*(j-1)+3)

  beta29(3*(r-1)+1,3*(i-1)+1)=beta2(3*(r-1)+1,nl)*betaJmI(nl,3*(i-1)+1)&
                             +beta2(3*(r-1)+1,sl)*betaJmI(sl,3*(i-1)+1) 
  beta29(3*(r-1)+2,3*(i-1)+2)=beta2(3*(r-1)+2,ml)*betaJmI(ml,3*(i-1)+2)&
                             +beta2(3*(r-1)+2,gl)*betaJmI(gl,3*(i-1)+2) 
  beta29(3*(r-1)+3,3*(i-1)+3)=beta2(3*(r-1)+3,rl)*betaJmI(rl,3*(i-1)+3)&
                             +beta2(3*(r-1)+3,hl)*betaJmI(hl,3*(i-1)+3) 
  beta29(3*(r-1)+1,3*(j-1)+1)=beta2(3*(r-1)+1,nl)*betaJmI(nl,3*(j-1)+1)&
                             +beta2(3*(r-1)+1,sl)*betaJmI(sl,3*(j-1)+1) 
  beta29(3*(r-1)+2,3*(j-1)+2)=beta2(3*(r-1)+2,ml)*betaJmI(ml,3*(j-1)+2)&
                             +beta2(3*(r-1)+2,gl)*betaJmI(gl,3*(j-1)+2) 
  beta29(3*(r-1)+3,3*(j-1)+3)=beta2(3*(r-1)+3,rl)*betaJmI(rl,3*(j-1)+3)&
                             +beta2(3*(r-1)+3,hl)*betaJmI(hl,3*(j-1)+3) 

   beta31(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta3(3*(i-1)+1,nl)&
                       +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,nl)
   beta31(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta3(3*(i-1)+2,ml)&
                       +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,ml)
   beta31(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta3(3*(i-1)+3,rl)&
                       +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,rl)
   beta31(3*(r-1)+1,sl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta3(3*(i-1)+1,sl)&
                       +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta3(3*(j-1)+1,sl)
   beta31(3*(r-1)+2,gl)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta3(3*(i-1)+2,gl)&
                       +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta3(3*(j-1)+2,gl)
   beta31(3*(r-1)+3,hl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta3(3*(i-1)+3,hl)&
                       +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta3(3*(j-1)+3,hl)

   beta38(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)&
                       +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl) 
   beta38(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)&
                       +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml) 
   beta38(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)&
                       +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl) 
   beta38(3*(r-1)+1,sl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,sl)&
                       +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,sl) 
   beta38(3*(r-1)+2,gl)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,gl)&
                       +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,gl) 
   beta38(3*(r-1)+3,hl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,hl)&
                       +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,hl) 

  if (i.eq.j) then
 
     beta31(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta3(3*(i-1)+1,nl)
     beta31(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta3(3*(i-1)+2,ml)
     beta31(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta3(3*(i-1)+3,rl)
     beta31(3*(r-1)+1,sl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta3(3*(i-1)+1,sl)
     beta31(3*(r-1)+2,gl)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta3(3*(i-1)+2,gl)
     beta31(3*(r-1)+3,hl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta3(3*(i-1)+3,hl)
  
     beta38(3*(r-1)+1,nl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)
     beta38(3*(r-1)+2,ml)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)
     beta38(3*(r-1)+3,rl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)
     beta38(3*(r-1)+1,sl)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,sl)
     beta38(3*(r-1)+2,gl)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,gl)
     beta38(3*(r-1)+3,hl)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,hl)

  endif
enddo
endif

beta30=ZERO
beta32=ZERO
beta34=ZERO
beta36=ZERO
beta39=ZERO

do r=1,n
  
if (i_kk.eq.j_kk) then

  beta30(3*(r-1)+1,nk)=beta29(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)&
                      +beta29(3*(r-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,nk)
  beta30(3*(r-1)+2,mk)=beta29(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)&
                      +beta29(3*(r-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,mk)
  beta30(3*(r-1)+3,rk)=beta29(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)&
                      +beta29(3*(r-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,rk)

  beta36(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta33(3*(i-1)+1,nk)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,nk)
  beta36(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta33(3*(i-1)+2,mk)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,mk)
  beta36(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta33(3*(i-1)+3,rk)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,rk)

  beta39(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,nk) 
  beta39(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,mk) 
  beta39(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,rk)

  if (i.eq.j) then

    beta30(3*(r-1)+1,nk)=beta29(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)
    beta30(3*(r-1)+2,mk)=beta29(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)
    beta30(3*(r-1)+3,rk)=beta29(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)

    beta36(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta33(3*(i-1)+1,nk)
    beta36(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta33(3*(i-1)+2,mk)
    beta36(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta33(3*(i-1)+3,rk)

    beta39(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)
    beta39(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)
    beta39(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)

  endif

  beta32(3*(r-1)+1,3*(i-1)+1)=beta33(3*(r-1)+1,nk)*betaJmI(nk,3*(i-1)+1)
  beta32(3*(r-1)+2,3*(i-1)+2)=beta33(3*(r-1)+2,mk)*betaJmI(mk,3*(i-1)+2)
  beta32(3*(r-1)+3,3*(i-1)+3)=beta33(3*(r-1)+3,rk)*betaJmI(rk,3*(i-1)+3)
  beta32(3*(r-1)+1,3*(j-1)+1)=beta33(3*(r-1)+1,nk)*betaJmI(nk,3*(j-1)+1)
  beta32(3*(r-1)+2,3*(j-1)+2)=beta33(3*(r-1)+2,mk)*betaJmI(mk,3*(j-1)+2)
  beta32(3*(r-1)+3,3*(j-1)+3)=beta33(3*(r-1)+3,rk)*betaJmI(rk,3*(j-1)+3)

  beta34(3*(r-1)+1,3*(i-1)+1)=beta1(3*(r-1)+1,nk)*betaJmI(nk,3*(i-1)+1)  
  beta34(3*(r-1)+2,3*(i-1)+2)=beta1(3*(r-1)+2,mk)*betaJmI(mk,3*(i-1)+2) 
  beta34(3*(r-1)+3,3*(i-1)+3)=beta1(3*(r-1)+3,rk)*betaJmI(rk,3*(i-1)+3)  
  beta34(3*(r-1)+1,3*(j-1)+1)=beta1(3*(r-1)+1,nk)*betaJmI(nk,3*(j-1)+1) 
  beta34(3*(r-1)+2,3*(j-1)+2)=beta1(3*(r-1)+2,mk)*betaJmI(mk,3*(j-1)+2)  
  beta34(3*(r-1)+3,3*(j-1)+3)=beta1(3*(r-1)+3,rk)*betaJmI(rk,3*(j-1)+3)
 
endif

if (i_kk.ne.j_kk) then

  beta30(3*(r-1)+1,nk)=beta29(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)&
                      +beta29(3*(r-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,nk)
  beta30(3*(r-1)+2,mk)=beta29(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)&
                      +beta29(3*(r-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,mk)
  beta30(3*(r-1)+3,rk)=beta29(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)&
                      +beta29(3*(r-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,rk)
  beta30(3*(r-1)+1,sk)=beta29(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,sk)&
                      +beta29(3*(r-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,sk)
  beta30(3*(r-1)+2,gk)=beta29(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,gk)&
                      +beta29(3*(r-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,gk)
  beta30(3*(r-1)+3,hk)=beta29(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,hk)&
                      +beta29(3*(r-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,hk)

  beta36(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta33(3*(i-1)+1,nk)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,nk)
  beta36(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta33(3*(i-1)+2,mk)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,mk)
  beta36(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta33(3*(i-1)+3,rk)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,rk)
  beta36(3*(r-1)+1,sk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta33(3*(i-1)+1,sk)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta33(3*(j-1)+1,sk)
  beta36(3*(r-1)+2,gk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta33(3*(i-1)+2,gk)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta33(3*(j-1)+2,gk)
  beta36(3*(r-1)+3,hk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta33(3*(i-1)+3,hk)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta33(3*(j-1)+3,hk)

  beta39(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,nk) 
  beta39(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,mk) 
  beta39(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,rk) 
  beta39(3*(r-1)+1,sk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,sk)&
                      +betaJmI(3*(r-1)+1,3*(j-1)+1)*beta1(3*(j-1)+1,sk) 
  beta39(3*(r-1)+2,gk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,gk)&
                      +betaJmI(3*(r-1)+2,3*(j-1)+2)*beta1(3*(j-1)+2,gk) 
  beta39(3*(r-1)+3,hk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,hk)&
                      +betaJmI(3*(r-1)+3,3*(j-1)+3)*beta1(3*(j-1)+3,hk) 

  if (i.eq.j) then

    beta30(3*(r-1)+1,nk)=beta29(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)
    beta30(3*(r-1)+2,mk)=beta29(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)
    beta30(3*(r-1)+3,rk)=beta29(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)
    beta30(3*(r-1)+1,sk)=beta29(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,sk)
    beta30(3*(r-1)+2,gk)=beta29(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,gk)
    beta30(3*(r-1)+3,hk)=beta29(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,hk)

    beta36(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta33(3*(i-1)+1,nk)
    beta36(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta33(3*(i-1)+2,mk)
    beta36(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta33(3*(i-1)+3,rk)
    beta36(3*(r-1)+1,sk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta33(3*(i-1)+1,sk)
    beta36(3*(r-1)+2,gk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta33(3*(i-1)+2,gk)
    beta36(3*(r-1)+3,hk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta33(3*(i-1)+3,hk)

    beta39(3*(r-1)+1,nk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,nk)
    beta39(3*(r-1)+2,mk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,mk)
    beta39(3*(r-1)+3,rk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,rk)
    beta39(3*(r-1)+1,sk)=betaJmI(3*(r-1)+1,3*(i-1)+1)*beta1(3*(i-1)+1,sk)
    beta39(3*(r-1)+2,gk)=betaJmI(3*(r-1)+2,3*(i-1)+2)*beta1(3*(i-1)+2,gk)
    beta39(3*(r-1)+3,hk)=betaJmI(3*(r-1)+3,3*(i-1)+3)*beta1(3*(i-1)+3,hk)

  endif

  beta32(3*(r-1)+1,3*(i-1)+1)=beta33(3*(r-1)+1,nk)*betaJmI(nk,3*(i-1)+1)&
                             +beta33(3*(r-1)+1,sk)*betaJmI(sk,3*(i-1)+1)  
  beta32(3*(r-1)+2,3*(i-1)+2)=beta33(3*(r-1)+2,mk)*betaJmI(mk,3*(i-1)+2)&
                             +beta33(3*(r-1)+2,gk)*betaJmI(gk,3*(i-1)+2)  
  beta32(3*(r-1)+3,3*(i-1)+3)=beta33(3*(r-1)+3,rk)*betaJmI(rk,3*(i-1)+3)&
                             +beta33(3*(r-1)+3,hk)*betaJmI(hk,3*(i-1)+3)  
  beta32(3*(r-1)+1,3*(j-1)+1)=beta33(3*(r-1)+1,nk)*betaJmI(nk,3*(j-1)+1)&
                             +beta33(3*(r-1)+1,sk)*betaJmI(sk,3*(j-1)+1)  
  beta32(3*(r-1)+2,3*(j-1)+2)=beta33(3*(r-1)+2,mk)*betaJmI(mk,3*(j-1)+2)&
                             +beta33(3*(r-1)+2,gk)*betaJmI(gk,3*(j-1)+2)  
  beta32(3*(r-1)+3,3*(j-1)+3)=beta33(3*(r-1)+3,rk)*betaJmI(rk,3*(j-1)+3)&
                             +beta33(3*(r-1)+3,hk)*betaJmI(hk,3*(j-1)+3)  

  beta34(3*(r-1)+1,3*(i-1)+1)=beta1(3*(r-1)+1,nk)*betaJmI(nk,3*(i-1)+1)&
                             +beta1(3*(r-1)+1,sk)*betaJmI(sk,3*(i-1)+1)
  beta34(3*(r-1)+2,3*(i-1)+2)=beta1(3*(r-1)+2,mk)*betaJmI(mk,3*(i-1)+2)&
                             +beta1(3*(r-1)+2,gk)*betaJmI(gk,3*(i-1)+2)   
  beta34(3*(r-1)+3,3*(i-1)+3)=beta1(3*(r-1)+3,rk)*betaJmI(rk,3*(i-1)+3)&  
                             +beta1(3*(r-1)+3,hk)*betaJmI(hk,3*(i-1)+3)   
  beta34(3*(r-1)+1,3*(j-1)+1)=beta1(3*(r-1)+1,nk)*betaJmI(nk,3*(j-1)+1)&  
                             +beta1(3*(r-1)+1,sk)*betaJmI(sk,3*(j-1)+1)   
  beta34(3*(r-1)+2,3*(j-1)+2)=beta1(3*(r-1)+2,mk)*betaJmI(mk,3*(j-1)+2)&  
                             +beta1(3*(r-1)+2,gk)*betaJmI(gk,3*(j-1)+2)
  beta34(3*(r-1)+3,3*(j-1)+3)=beta1(3*(r-1)+3,rk)*betaJmI(rk,3*(j-1)+3)&
                             +beta1(3*(r-1)+3,hk)*betaJmI(hk,3*(j-1)+3)

endif

enddo

beta35=ZERO
beta37=ZERO
beta41=ZERO

do r=1,n

if (i_ll.eq.j_ll) then

  beta35(3*(r-1)+1,nl)=beta34(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)&
                      +beta34(3*(r-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl) 
  beta35(3*(r-1)+2,ml)=beta34(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)&
                      +beta34(3*(r-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml) 
  beta35(3*(r-1)+3,rl)=beta34(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)&
                      +beta34(3*(r-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl) 

  if (i.eq.j) then

    beta35(3*(r-1)+1,nl)=beta34(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)
    beta35(3*(r-1)+2,ml)=beta34(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)
    beta35(3*(r-1)+3,rl)=beta34(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)

  endif
endif

if (i_ll.ne.j_ll) then

  beta35(3*(r-1)+1,nl)=beta34(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)& 
                      +beta34(3*(r-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,nl)
  beta35(3*(r-1)+2,ml)=beta34(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)& 
                      +beta34(3*(r-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,ml)  
  beta35(3*(r-1)+3,rl)=beta34(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)& 
                      +beta34(3*(r-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,rl)    
  beta35(3*(r-1)+1,sl)=beta34(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,sl)& 
                      +beta34(3*(r-1)+1,3*(j-1)+1)*beta2(3*(j-1)+1,sl)
  beta35(3*(r-1)+2,gl)=beta34(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,gl)& 
                      +beta34(3*(r-1)+2,3*(j-1)+2)*beta2(3*(j-1)+2,gl)  
  beta35(3*(r-1)+3,hl)=beta34(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,hl)& 
                      +beta34(3*(r-1)+3,3*(j-1)+3)*beta2(3*(j-1)+3,hl)    

  if (i.eq.j) then

    beta35(3*(r-1)+1,nl)=beta34(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,nl)  
    beta35(3*(r-1)+2,ml)=beta34(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,ml)  
    beta35(3*(r-1)+3,rl)=beta34(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,rl)
    beta35(3*(r-1)+1,sl)=beta34(3*(r-1)+1,3*(i-1)+1)*beta2(3*(i-1)+1,sl)
    beta35(3*(r-1)+2,gl)=beta34(3*(r-1)+2,3*(i-1)+2)*beta2(3*(i-1)+2,gl)
    beta35(3*(r-1)+3,hl)=beta34(3*(r-1)+3,3*(i-1)+3)*beta2(3*(i-1)+3,hl)

  endif
endif

  beta37(3*(r-1)+1,3*(i-1)+1)=beta34(3*(r-1)+1,3*(i-1)+1)*beta29(3*(i-1)+1,3*(i-1)+1)&
                             +beta34(3*(r-1)+1,3*(j-1)+1)*beta29(3*(j-1)+1,3*(i-1)+1)
  beta37(3*(r-1)+2,3*(i-1)+2)=beta34(3*(r-1)+2,3*(i-1)+2)*beta29(3*(i-1)+2,3*(i-1)+2)&
                             +beta34(3*(r-1)+2,3*(j-1)+2)*beta29(3*(j-1)+2,3*(i-1)+2)
  beta37(3*(r-1)+3,3*(i-1)+3)=beta34(3*(r-1)+3,3*(i-1)+3)*beta29(3*(i-1)+3,3*(i-1)+3)&
                             +beta34(3*(r-1)+3,3*(j-1)+3)*beta29(3*(j-1)+3,3*(i-1)+3)
  beta37(3*(r-1)+1,3*(j-1)+1)=beta34(3*(r-1)+1,3*(i-1)+1)*beta29(3*(i-1)+1,3*(j-1)+1)&
                             +beta34(3*(r-1)+1,3*(j-1)+1)*beta29(3*(j-1)+1,3*(j-1)+1)
  beta37(3*(r-1)+2,3*(j-1)+2)=beta34(3*(r-1)+2,3*(i-1)+2)*beta29(3*(i-1)+2,3*(j-1)+2)&
                             +beta34(3*(r-1)+2,3*(j-1)+2)*beta29(3*(j-1)+2,3*(j-1)+2)
  beta37(3*(r-1)+3,3*(j-1)+3)=beta34(3*(r-1)+3,3*(i-1)+3)*beta29(3*(i-1)+3,3*(j-1)+3)&
                             +beta34(3*(r-1)+3,3*(j-1)+3)*beta29(3*(j-1)+3,3*(j-1)+3)

  beta41(3*(r-1)+1,3*(i-1)+1)=beta29(3*(r-1)+1,3*(i-1)+1)*beta34(3*(i-1)+1,3*(i-1)+1)&
                             +beta29(3*(r-1)+1,3*(j-1)+1)*beta34(3*(j-1)+1,3*(i-1)+1) 
  beta41(3*(r-1)+2,3*(i-1)+2)=beta29(3*(r-1)+2,3*(i-1)+2)*beta34(3*(i-1)+2,3*(i-1)+2)&
                             +beta29(3*(r-1)+2,3*(j-1)+2)*beta34(3*(j-1)+2,3*(i-1)+2) 
  beta41(3*(r-1)+3,3*(i-1)+3)=beta29(3*(r-1)+3,3*(i-1)+3)*beta34(3*(i-1)+3,3*(i-1)+3)&
                             +beta29(3*(r-1)+3,3*(j-1)+3)*beta34(3*(j-1)+3,3*(i-1)+3) 
  beta41(3*(r-1)+1,3*(j-1)+1)=beta29(3*(r-1)+1,3*(i-1)+1)*beta34(3*(i-1)+1,3*(j-1)+1)&
                             +beta29(3*(r-1)+1,3*(j-1)+1)*beta34(3*(j-1)+1,3*(j-1)+1) 
  beta41(3*(r-1)+2,3*(j-1)+2)=beta29(3*(r-1)+2,3*(i-1)+2)*beta34(3*(i-1)+2,3*(j-1)+2)&
                             +beta29(3*(r-1)+2,3*(j-1)+2)*beta34(3*(j-1)+2,3*(j-1)+2) 
  beta41(3*(r-1)+3,3*(j-1)+3)=beta29(3*(r-1)+3,3*(i-1)+3)*beta34(3*(i-1)+3,3*(j-1)+3)&
                             +beta29(3*(r-1)+3,3*(j-1)+3)*beta34(3*(j-1)+3,3*(j-1)+3) 

if (i.eq.j) then

  beta37(3*(r-1)+1,3*(i-1)+1)=beta34(3*(r-1)+1,3*(i-1)+1)*beta29(3*(i-1)+1,3*(i-1)+1)
  beta37(3*(r-1)+2,3*(i-1)+2)=beta34(3*(r-1)+2,3*(i-1)+2)*beta29(3*(i-1)+2,3*(i-1)+2)
  beta37(3*(r-1)+3,3*(i-1)+3)=beta34(3*(r-1)+3,3*(i-1)+3)*beta29(3*(i-1)+3,3*(i-1)+3)

  beta41(3*(r-1)+1,3*(i-1)+1)=beta29(3*(r-1)+1,3*(i-1)+1)*beta34(3*(i-1)+1,3*(i-1)+1)
  beta41(3*(r-1)+2,3*(i-1)+2)=beta29(3*(r-1)+2,3*(i-1)+2)*beta34(3*(i-1)+2,3*(i-1)+2)
  beta41(3*(r-1)+3,3*(i-1)+3)=beta29(3*(r-1)+3,3*(i-1)+3)*beta34(3*(i-1)+3,3*(i-1)+3)

endif
enddo

beta40=ZERO
beta42=ZERO

do r=1,n

  if (i_kk.ne.j_kk.and.i_ll.ne.j_ll) then

    beta40(3*(r-1)+1,nk)=beta38(3*(r-1)+1,nl)*beta39(nl,nk)&
                        +beta38(3*(r-1)+1,sl)*beta39(sl,nk)    
    beta40(3*(r-1)+2,mk)=beta38(3*(r-1)+2,ml)*beta39(ml,mk)&
                        +beta38(3*(r-1)+2,gl)*beta39(gl,mk) 
    beta40(3*(r-1)+3,rk)=beta38(3*(r-1)+3,rl)*beta39(rl,rk)&
                        +beta38(3*(r-1)+3,hl)*beta39(hl,rk) 
    beta40(3*(r-1)+1,sk)=beta38(3*(r-1)+1,nl)*beta39(nl,sk)&
                        +beta38(3*(r-1)+1,sl)*beta39(sl,sk) 
    beta40(3*(r-1)+2,gk)=beta38(3*(r-1)+2,ml)*beta39(ml,gk)&
                        +beta38(3*(r-1)+2,gl)*beta39(gl,gk) 
    beta40(3*(r-1)+3,hk)=beta38(3*(r-1)+3,rl)*beta39(rl,hk)&
                        +beta38(3*(r-1)+3,hl)*beta39(hl,hk) 
    
    beta42(3*(r-1)+1,nl)=beta39(3*(r-1)+1,nk)*beta38(nk,nl)&
                        +beta39(3*(r-1)+1,sk)*beta38(sk,nl) 
    beta42(3*(r-1)+2,ml)=beta39(3*(r-1)+2,mk)*beta38(mk,ml)&
                        +beta39(3*(r-1)+2,gk)*beta38(gk,ml) 
    beta42(3*(r-1)+3,rl)=beta39(3*(r-1)+3,rk)*beta38(rk,rl)&
                        +beta39(3*(r-1)+3,hk)*beta38(hk,rl) 
    beta42(3*(r-1)+1,sl)=beta39(3*(r-1)+1,nk)*beta38(nk,sl)&
                        +beta39(3*(r-1)+1,sk)*beta38(sk,sl) 
    beta42(3*(r-1)+2,gl)=beta39(3*(r-1)+2,mk)*beta38(mk,gl)&
                        +beta39(3*(r-1)+2,gk)*beta38(gk,gl) 
    beta42(3*(r-1)+3,hl)=beta39(3*(r-1)+3,rk)*beta38(rk,hl)&
                        +beta39(3*(r-1)+3,hk)*beta38(hk,hl)
  endif

  if (i_kk.ne.j_kk.and.i_ll.eq.j_ll) then

    beta40(3*(r-1)+1,nk)=beta38(3*(r-1)+1,nl)*beta39(nl,nk)
    beta40(3*(r-1)+2,mk)=beta38(3*(r-1)+2,ml)*beta39(ml,mk)
    beta40(3*(r-1)+3,rk)=beta38(3*(r-1)+3,rl)*beta39(rl,rk)
    beta40(3*(r-1)+1,sk)=beta38(3*(r-1)+1,nl)*beta39(nl,sk)
    beta40(3*(r-1)+2,gk)=beta38(3*(r-1)+2,ml)*beta39(ml,gk)
    beta40(3*(r-1)+3,hk)=beta38(3*(r-1)+3,rl)*beta39(rl,hk)

    beta42(3*(r-1)+1,nl)=beta39(3*(r-1)+1,nk)*beta38(nk,nl)&
                        +beta39(3*(r-1)+1,sk)*beta38(sk,nl) 
    beta42(3*(r-1)+2,ml)=beta39(3*(r-1)+2,mk)*beta38(mk,ml)&
                        +beta39(3*(r-1)+2,gk)*beta38(gk,ml) 
    beta42(3*(r-1)+3,rl)=beta39(3*(r-1)+3,rk)*beta38(rk,rl)&
                        +beta39(3*(r-1)+3,hk)*beta38(hk,rl) 

  endif

  if (i_kk.eq.j_kk.and.i_ll.ne.j_ll) then

    beta40(3*(r-1)+1,nk)=beta38(3*(r-1)+1,nl)*beta39(nl,nk)&
                        +beta38(3*(r-1)+1,sl)*beta39(sl,nk) 
    beta40(3*(r-1)+2,mk)=beta38(3*(r-1)+2,ml)*beta39(ml,mk)&
                        +beta38(3*(r-1)+2,gl)*beta39(gl,mk) 
    beta40(3*(r-1)+3,rk)=beta38(3*(r-1)+3,rl)*beta39(rl,rk)&
                        +beta38(3*(r-1)+3,hl)*beta39(hl,rk) 

    beta42(3*(r-1)+1,nl)=beta39(3*(r-1)+1,nk)*beta38(nk,nl)
    beta42(3*(r-1)+2,ml)=beta39(3*(r-1)+2,mk)*beta38(mk,ml)
    beta42(3*(r-1)+3,rl)=beta39(3*(r-1)+3,rk)*beta38(rk,rl)
    beta42(3*(r-1)+1,sl)=beta39(3*(r-1)+1,nk)*beta38(nk,sl)
    beta42(3*(r-1)+2,gl)=beta39(3*(r-1)+2,mk)*beta38(mk,gl)
    beta42(3*(r-1)+3,hl)=beta39(3*(r-1)+3,rk)*beta38(rk,hl)

  endif

  if (i_kk.eq.j_kk.and.i_ll.eq.j_ll) then

    beta40(3*(r-1)+1,nk)=beta38(3*(r-1)+1,nl)*beta39(nl,nk)
    beta40(3*(r-1)+2,mk)=beta38(3*(r-1)+2,ml)*beta39(ml,mk)
    beta40(3*(r-1)+3,rk)=beta38(3*(r-1)+3,rl)*beta39(rl,rk)

    beta42(3*(r-1)+1,nl)=beta39(3*(r-1)+1,nk)*beta38(nk,nl)
    beta42(3*(r-1)+2,ml)=beta39(3*(r-1)+2,mk)*beta38(mk,ml)
    beta42(3*(r-1)+3,rl)=beta39(3*(r-1)+3,rk)*beta38(rk,rl)

  endif

enddo

!**************************
!
! FINDING TRACE OF MATRICES
!
!**************************

!Xi=tr[betaJm]

Xi=betaJm(i,i)+betaJm(j,j)

nu6=beta28(3*(i-1)+1,3*(i-1)+1)&
   +beta28(3*(i-1)+2,3*(i-1)+2)&
   +beta28(3*(i-1)+3,3*(i-1)+3)&
   +beta28(3*(j-1)+1,3*(j-1)+1)&
   +beta28(3*(j-1)+2,3*(j-1)+2)&
   +beta28(3*(j-1)+3,3*(j-1)+3)

nu7=beta32(3*(i-1)+1,3*(i-1)+1)&
   +beta32(3*(i-1)+2,3*(i-1)+2)&
   +beta32(3*(i-1)+3,3*(i-1)+3)&
   +beta32(3*(j-1)+1,3*(j-1)+1)&
   +beta32(3*(j-1)+2,3*(j-1)+2)&
   +beta32(3*(j-1)+3,3*(j-1)+3) 

nu8=beta37(3*(i-1)+1,3*(i-1)+1)&
   +beta37(3*(i-1)+2,3*(i-1)+2)&
   +beta37(3*(i-1)+3,3*(i-1)+3)&
   +beta37(3*(j-1)+1,3*(j-1)+1)&
   +beta37(3*(j-1)+2,3*(j-1)+2)&
   +beta37(3*(j-1)+3,3*(j-1)+3) 

if (i.eq.j) then 

  Xi=betaJm(i,i)

nu6=beta28(3*(i-1)+1,3*(i-1)+1)&
   +beta28(3*(i-1)+2,3*(i-1)+2)&
   +beta28(3*(i-1)+3,3*(i-1)+3)

nu7=beta32(3*(i-1)+1,3*(i-1)+1)&
   +beta32(3*(i-1)+2,3*(i-1)+2)&
   +beta32(3*(i-1)+3,3*(i-1)+3)

nu8=beta37(3*(i-1)+1,3*(i-1)+1)&
   +beta37(3*(i-1)+2,3*(i-1)+2)&
   +beta37(3*(i-1)+3,3*(i-1)+3) 

endif

!write(*,*)'Xi',Xi

!write(*,*)'nu6,nu7,nu8',nu6,nu7,nu8

!**************************
!
! Potential Energy Gradient
!
!**************************

if (grad_k.or.grad_l) then

  alpha25=ZERO

  do r=1,n
    do t=1,n

      alpha25(r,t)=betaJm(r,i)*inv_tAkl(i,t)&
                  +betaJm(r,j)*inv_tAkl(j,t)

      alpha25(r,t)=betaJm(r,i)*inv_tAkl(i,t)&
                  +betaJm(r,j)*inv_tAkl(j,t) 

      if (i.eq.j) then
   
        alpha25(r,t)=betaJm(r,i)*inv_tAkl(i,t)
   
      endif

    enddo
  enddo

!write(*,*)'alpha25'
!do r=1,n
!  write(*,9292)(alpha25(r,t),t=1,n)
!enddo

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

dnu6=ZERO
dnu66=ZERO
dnu666=ZERO
dnu7=ZERO
dnu77=ZERO
dnu777=ZERO
dnu8=ZERO
dnu88=ZERO
dnu888=ZERO
dnu8888=ZERO

do r=1,n

  do t=1,n

    dnu6(3*(r-1)+1,3*(t-1)+1)=beta28(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)&
                             +beta28(3*(r-1)+1,3*(j-1)+1)*inv_tAklI(3*(j-1)+1,3*(t-1)+1)
    dnu6(3*(r-1)+2,3*(t-1)+2)=beta28(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)&
                             +beta28(3*(r-1)+2,3*(j-1)+2)*inv_tAklI(3*(j-1)+2,3*(t-1)+2)
    dnu6(3*(r-1)+3,3*(t-1)+3)=beta28(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)&
                             +beta28(3*(r-1)+3,3*(j-1)+3)*inv_tAklI(3*(j-1)+3,3*(t-1)+3)

    dnu7(3*(r-1)+1,3*(t-1)+1)=beta32(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)&
                             +beta32(3*(r-1)+1,3*(j-1)+1)*inv_tAklI(3*(j-1)+1,3*(t-1)+1) 
    dnu7(3*(r-1)+2,3*(t-1)+2)=beta32(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)&
                             +beta32(3*(r-1)+2,3*(j-1)+2)*inv_tAklI(3*(j-1)+2,3*(t-1)+2) 
    dnu7(3*(r-1)+3,3*(t-1)+3)=beta32(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)&
                             +beta32(3*(r-1)+3,3*(j-1)+3)*inv_tAklI(3*(j-1)+3,3*(t-1)+3) 

    dnu8(3*(r-1)+1,3*(t-1)+1)=beta37(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)&
                             +beta37(3*(r-1)+1,3*(j-1)+1)*inv_tAklI(3*(j-1)+1,3*(t-1)+1)
    dnu8(3*(r-1)+2,3*(t-1)+2)=beta37(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)&
                             +beta37(3*(r-1)+2,3*(j-1)+2)*inv_tAklI(3*(j-1)+2,3*(t-1)+2)
    dnu8(3*(r-1)+3,3*(t-1)+3)=beta37(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)&
                             +beta37(3*(r-1)+3,3*(j-1)+3)*inv_tAklI(3*(j-1)+3,3*(t-1)+3)

    dnu888(3*(r-1)+1,3*(t-1)+1)=beta41(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)&
                               +beta41(3*(r-1)+1,3*(j-1)+1)*inv_tAklI(3*(j-1)+1,3*(t-1)+1)
    dnu888(3*(r-1)+2,3*(t-1)+2)=beta41(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)&
                               +beta41(3*(r-1)+2,3*(j-1)+2)*inv_tAklI(3*(j-1)+2,3*(t-1)+2)
    dnu888(3*(r-1)+3,3*(t-1)+3)=beta41(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)&
                               +beta41(3*(r-1)+3,3*(j-1)+3)*inv_tAklI(3*(j-1)+3,3*(t-1)+3)

if (i.eq.j) then

    dnu6(3*(r-1)+1,3*(t-1)+1)=beta28(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)
    dnu6(3*(r-1)+2,3*(t-1)+2)=beta28(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)
    dnu6(3*(r-1)+3,3*(t-1)+3)=beta28(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)

    dnu7(3*(r-1)+1,3*(t-1)+1)=beta32(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)
    dnu7(3*(r-1)+2,3*(t-1)+2)=beta32(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)
    dnu7(3*(r-1)+3,3*(t-1)+3)=beta32(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)

    dnu8(3*(r-1)+1,3*(t-1)+1)=beta37(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)
    dnu8(3*(r-1)+2,3*(t-1)+2)=beta37(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)
    dnu8(3*(r-1)+3,3*(t-1)+3)=beta37(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)

    dnu888(3*(r-1)+1,3*(t-1)+1)=beta41(3*(r-1)+1,3*(i-1)+1)*inv_tAklI(3*(i-1)+1,3*(t-1)+1)
    dnu888(3*(r-1)+2,3*(t-1)+2)=beta41(3*(r-1)+2,3*(i-1)+2)*inv_tAklI(3*(i-1)+2,3*(t-1)+2)
    dnu888(3*(r-1)+3,3*(t-1)+3)=beta41(3*(r-1)+3,3*(i-1)+3)*inv_tAklI(3*(i-1)+3,3*(t-1)+3)

endif
 
if (i_kk.ne.j_kk) then

   dnu66(3*(r-1)+1,3*(t-1)+1)=beta30(3*(r-1)+1,nk)*inv_tAklI(nk,3*(t-1)+1)&
                             +beta30(3*(r-1)+1,sk)*inv_tAklI(sk,3*(t-1)+1)
   dnu66(3*(r-1)+2,3*(t-1)+2)=beta30(3*(r-1)+2,mk)*inv_tAklI(mk,3*(t-1)+2)&
                             +beta30(3*(r-1)+2,gk)*inv_tAklI(gk,3*(t-1)+2)
   dnu66(3*(r-1)+3,3*(t-1)+3)=beta30(3*(r-1)+3,rk)*inv_tAklI(rk,3*(t-1)+3)&
                             +beta30(3*(r-1)+3,hk)*inv_tAklI(hk,3*(t-1)+3)

   dnu777(3*(r-1)+1,3*(t-1)+1)=beta36(3*(r-1)+1,nk)*inv_tAklI(nk,3*(t-1)+1)&
                              +beta36(3*(r-1)+1,sk)*inv_tAklI(sk,3*(t-1)+1) 
   dnu777(3*(r-1)+2,3*(t-1)+2)=beta36(3*(r-1)+2,mk)*inv_tAklI(mk,3*(t-1)+2)&
                              +beta36(3*(r-1)+2,gk)*inv_tAklI(gk,3*(t-1)+2) 
   dnu777(3*(r-1)+3,3*(t-1)+3)=beta36(3*(r-1)+3,rk)*inv_tAklI(rk,3*(t-1)+3)&
                              +beta36(3*(r-1)+3,hk)*inv_tAklI(hk,3*(t-1)+3)

   dnu88(3*(r-1)+1,3*(t-1)+1)=beta40(3*(r-1)+1,nk)*inv_tAklI(nk,3*(t-1)+1)&
                             +beta40(3*(r-1)+1,sk)*inv_tAklI(sk,3*(t-1)+1)
   dnu88(3*(r-1)+2,3*(t-1)+2)=beta40(3*(r-1)+2,mk)*inv_tAklI(mk,3*(t-1)+2)&
                             +beta40(3*(r-1)+2,gk)*inv_tAklI(gk,3*(t-1)+2)
   dnu88(3*(r-1)+3,3*(t-1)+3)=beta40(3*(r-1)+3,rk)*inv_tAklI(rk,3*(t-1)+3)&
                             +beta40(3*(r-1)+3,hk)*inv_tAklI(hk,3*(t-1)+3)

endif

if (i_kk.eq.j_kk) then

  dnu66(3*(r-1)+1,3*(t-1)+1)=beta30(3*(r-1)+1,nk)*inv_tAklI(nk,3*(t-1)+1)
  dnu66(3*(r-1)+2,3*(t-1)+2)=beta30(3*(r-1)+2,mk)*inv_tAklI(mk,3*(t-1)+2)
  dnu66(3*(r-1)+3,3*(t-1)+3)=beta30(3*(r-1)+3,rk)*inv_tAklI(rk,3*(t-1)+3)

  dnu777(3*(r-1)+1,3*(t-1)+1)=beta36(3*(r-1)+1,nk)*inv_tAklI(nk,3*(t-1)+1)
  dnu777(3*(r-1)+2,3*(t-1)+2)=beta36(3*(r-1)+2,mk)*inv_tAklI(mk,3*(t-1)+2)   
  dnu777(3*(r-1)+3,3*(t-1)+3)=beta36(3*(r-1)+3,rk)*inv_tAklI(rk,3*(t-1)+3)

  dnu88(3*(r-1)+1,3*(t-1)+1)=beta40(3*(r-1)+1,nk)*inv_tAklI(nk,3*(t-1)+1)
  dnu88(3*(r-1)+2,3*(t-1)+2)=beta40(3*(r-1)+2,mk)*inv_tAklI(mk,3*(t-1)+2)
  dnu88(3*(r-1)+3,3*(t-1)+3)=beta40(3*(r-1)+3,rk)*inv_tAklI(rk,3*(t-1)+3)

endif

if (i_ll.ne.j_ll) then

   dnu666(3*(r-1)+1,3*(t-1)+1)=beta31(3*(r-1)+1,nl)*inv_tAklI(nl,3*(t-1)+1)&
                              +beta31(3*(r-1)+1,sl)*inv_tAklI(sl,3*(t-1)+1)
   dnu666(3*(r-1)+2,3*(t-1)+2)=beta31(3*(r-1)+2,ml)*inv_tAklI(ml,3*(t-1)+2)&
                              +beta31(3*(r-1)+2,gl)*inv_tAklI(gl,3*(t-1)+2)
   dnu666(3*(r-1)+3,3*(t-1)+3)=beta31(3*(r-1)+3,rl)*inv_tAklI(rl,3*(t-1)+3)&
                              +beta31(3*(r-1)+3,hl)*inv_tAklI(hl,3*(t-1)+3)

   dnu77(3*(r-1)+1,3*(t-1)+1)=beta35(3*(r-1)+1,nl)*inv_tAklI(nl,3*(t-1)+1)&
                             +beta35(3*(r-1)+1,sl)*inv_tAklI(sl,3*(t-1)+1) 
   dnu77(3*(r-1)+2,3*(t-1)+2)=beta35(3*(r-1)+2,ml)*inv_tAklI(ml,3*(t-1)+2)&
                             +beta35(3*(r-1)+2,gl)*inv_tAklI(gl,3*(t-1)+2) 
   dnu77(3*(r-1)+3,3*(t-1)+3)=beta35(3*(r-1)+3,rl)*inv_tAklI(rl,3*(t-1)+3)&
                             +beta35(3*(r-1)+3,hl)*inv_tAklI(hl,3*(t-1)+3) 

   dnu8888(3*(r-1)+1,3*(t-1)+1)=beta42(3*(r-1)+1,nl)*inv_tAklI(nl,3*(t-1)+1)&
                               +beta42(3*(r-1)+1,sl)*inv_tAklI(sl,3*(t-1)+1) 
   dnu8888(3*(r-1)+2,3*(t-1)+2)=beta42(3*(r-1)+2,ml)*inv_tAklI(ml,3*(t-1)+2)&
                               +beta42(3*(r-1)+2,gl)*inv_tAklI(gl,3*(t-1)+2) 
   dnu8888(3*(r-1)+3,3*(t-1)+3)=beta42(3*(r-1)+3,rl)*inv_tAklI(rl,3*(t-1)+3)&
                               +beta42(3*(r-1)+3,hl)*inv_tAklI(hl,3*(t-1)+3) 

endif

if (i_ll.eq.j_ll) then

  dnu666(3*(r-1)+1,3*(t-1)+1)=beta31(3*(r-1)+1,nl)*inv_tAklI(nl,3*(t-1)+1)
  dnu666(3*(r-1)+2,3*(t-1)+2)=beta31(3*(r-1)+2,ml)*inv_tAklI(ml,3*(t-1)+2)
  dnu666(3*(r-1)+3,3*(t-1)+3)=beta31(3*(r-1)+3,rl)*inv_tAklI(rl,3*(t-1)+3)

  dnu77(3*(r-1)+1,3*(t-1)+1)=beta35(3*(r-1)+1,nl)*inv_tAklI(nl,3*(t-1)+1)
  dnu77(3*(r-1)+2,3*(t-1)+2)=beta35(3*(r-1)+2,ml)*inv_tAklI(ml,3*(t-1)+2)
  dnu77(3*(r-1)+3,3*(t-1)+3)=beta35(3*(r-1)+3,rl)*inv_tAklI(rl,3*(t-1)+3)

  dnu8888(3*(r-1)+1,3*(t-1)+1)=beta42(3*(r-1)+1,nl)*inv_tAklI(nl,3*(t-1)+1)
  dnu8888(3*(r-1)+2,3*(t-1)+2)=beta42(3*(r-1)+2,ml)*inv_tAklI(ml,3*(t-1)+2)
  dnu8888(3*(r-1)+3,3*(t-1)+3)=beta42(3*(r-1)+3,rl)*inv_tAklI(rl,3*(t-1)+3)

endif

  enddo
enddo

endif

if (grad_k) then

  dnu6k=ZERO
  dnu7k=ZERO
  dnu8k=ZERO

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
      temp7=ZERO
      temp8=ZERO
      temp9=ZERO
      temp10=ZERO

      do k=1,3*n

        temp1=temp1+(dnu6(r,k)+dnu6(k,r))*LkI(k,t)
        temp2=temp2+(dnu66(r,k)+dnu66(k,r))*LkI(k,t)
        temp3=temp3+(dnu666(r,k)+dnu666(k,r))*LkI(k,t)
        temp4=temp4+(dnu7(r,k)+dnu7(k,r))*LkI(k,t)
        temp5=temp5+(dnu77(r,k)+dnu77(k,r))*LkI(k,t)
        temp6=temp6+(dnu777(r,k)+dnu777(k,r))*LkI(k,t)
        temp7=temp7+(dnu8(r,k)+dnu8(k,r))*LkI(k,t)
        temp8=temp8+(dnu88(r,k)+dnu88(k,r))*LkI(k,t)
        temp9=temp9+(dnu888(r,k)+dnu888(k,r))*LkI(k,t)
        temp10=temp10+(dnu8888(r,k)+dnu8888(k,r))*LkI(k,t)

     enddo

     dnu6k(indx)=temp1+temp2+temp3
     dnu7k(indx)=temp4+temp5+temp6
     dnu8k(indx)=temp7+temp8+temp9+temp10

  enddo
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
  Z9=ZERO
  Z10=ZERO
  Z11=ZERO
  Z12=ZERO

  do r=1,3*n
    do t=1,3*n

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

      do k=1,3*n

        temp1=temp1+pP(r,k)*dnu6(k,t)
        temp2=temp2+pP(r,k)*dnu66(k,t)
        temp3=temp3+pP(r,k)*dnu666(k,t)
        temp4=temp4+pP(r,k)*dnu7(k,t)
        temp5=temp5+pP(r,k)*dnu77(k,t)
        temp6=temp6+pP(r,k)*dnu777(k,t)
        temp7=temp7+pP(r,k)*dnu8(k,t)
        temp8=temp8+pP(r,k)*dnu88(k,t)  
        temp9=temp9+pP(r,k)*dnu888(k,t) 
        temp10=temp10+pP(r,k)*dnu8888(k,t)

    enddo

      Z3(r,t)=temp1
      Z4(r,t)=temp2
      Z5(r,t)=temp3
      Z6(r,t)=temp4
      Z7(r,t)=temp5
      Z8(r,t)=temp6
      Z9(r,t)=temp7
      Z10(r,t)=temp8
      Z11(r,t)=temp9
      Z12(r,t)=temp10

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
      temp7=ZERO
      temp8=ZERO
      temp9=ZERO
      temp10=ZERO

      do k=1,3*n

        temp1=temp1+Z3(r,k)*Pp(t,k)
        temp2=temp2+Z4(r,k)*Pp(t,k)
        temp3=temp3+Z5(r,k)*Pp(t,k)
        temp4=temp4+Z6(r,k)*Pp(t,k)
        temp5=temp5+Z7(r,k)*Pp(t,k)
        temp6=temp6+Z8(r,k)*Pp(t,k)
        temp7=temp7+Z9(r,k)*Pp(t,k)
        temp8=temp8+Z10(r,k)*Pp(t,k)
        temp9=temp9+Z11(r,k)*Pp(t,k)
        temp10=temp10+Z12(r,k)*Pp(t,k)

      enddo

      dnu6p(r,t)=temp1
      dnu66p(r,t)=temp2
      dnu666p(r,t)=temp3
      dnu7p(r,t)=temp4
      dnu77p(r,t)=temp5
      dnu777p(r,t)=temp6
      dnu8p(r,t)=temp7
      dnu88p(r,t)=temp8  
      dnu888p(r,t)=temp9 
      dnu8888p(r,t)=temp10

    enddo
  enddo

  dnu6l=ZERO
  dnu7l=ZERO
  dnu8l=ZERO

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
      temp7=ZERO
      temp8=ZERO
      temp9=ZERO
      temp10=ZERO

      do k=1,3*n

        temp1=temp1+(dnu6p(r,k)+dnu6p(k,r))*LlI(k,t)
        temp2=temp2+(dnu66p(r,k)+dnu66p(k,r))*LlI(k,t)
        temp3=temp3+(dnu666p(r,k)+dnu666p(k,r))*LlI(k,t)
        temp4=temp4+(dnu7p(r,k)+dnu7p(k,r))*LlI(k,t)
        temp5=temp5+(dnu77p(r,k)+dnu77p(k,r))*LlI(k,t)
        temp6=temp6+(dnu777p(r,k)+dnu777p(k,r))*LlI(k,t)
        temp7=temp7+(dnu8p(r,k)+dnu8p(k,r))*LlI(k,t)
        temp8=temp8+(dnu88p(r,k)+dnu88p(k,r))*LlI(k,t)
        temp9=temp9+(dnu888p(r,k)+dnu888p(k,r))*LlI(k,t)
        temp10=temp10+(dnu8888p(r,k)+dnu8888p(k,r))*LlI(k,t)

     enddo

     dnu6l(indx)=temp1+temp2+temp3
     dnu7l(indx)=temp4+temp5+temp6
     dnu8l(indx)=temp7+temp8+temp9+temp10

  enddo
enddo

endif

!********************
!
! COMPUTING POTENTIAL
! (non-normalized)   
!********************

Rkl=TWO*(PI**((THREE*n-ONE)/TWO))&
*(det_tAkl**(-THREE/TWO))&
*(Xi**(-ONE/TWO))&
*(ONEHALF*nu3&
-(ONE/(SIX))*(Xi**(-ONE))*(nu6+nu7)&
+(ONE/(TEN))*(Xi**(-TWO))*nu8)

if (i.ne.j) Vkl=Vkl+Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*Rkl

if (i.eq.j) Vkl=Vkl+Glob_PseudoCharge0*Glob_PseudoCharge(i)*Rkl

!if (i.ne.j) write(*,*)'i,j,Rkl,Vkl',i,j,Rkl,Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*Rkl
!if (i.eq.j) write(*,*)'i,j,Rkl,Vkl',i,j,Rkl,Glob_PseudoCharge0*Glob_PseudoCharge(i)*Rkl

if (grad_k) then

pk=ZERO

  do r=1,n*(n+1)/2

    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    temp4=ZERO
       
    do k=1,3*n*(3*n+1)/2

!dnu3
      temp1=temp1+dnu3k(k)*trans(k,r)
!dnu6
      temp2=temp2+dnu6k(k)*trans(k,r)
!dnu7
      temp3=temp3+dnu7k(k)*trans(k,r)
!dnu8
      temp4=temp4+dnu8k(k)*trans(k,r)

   enddo

pk(r)=pk(r)+TWO*(PI**((THREE*n-ONE)/TWO))&
    *(det_tAkl**(-THREE/TWO))&
    *(Xi**(-ONE/TWO))&
    *((ONEHALF*(Xi**(-ONE))*resultJk(r)-THREEHALF*resultAk(r))&
      *(ONEHALF*nu3&
        -(ONE/SIX)*(Xi**(-ONE))*(nu6+nu7)&
        +(ONE/TEN)*(Xi**(-TWO))*nu8)&
      -ONEHALF*temp1&
      -(ONE/SIX)*(Xi**(-TWO))*resultJk(r)*(nu6+nu7)&
      +(ONE/SIX)*(Xi**(-ONE))*(temp2+temp3)&
      +(ONE/FIVE)*(Xi**(-THREE))*resultJk(r)*nu8&
      -(ONE/TEN)*(Xi**(-TWO))*temp4)

  enddo

  do r=1,n*(n+1)/2

!    if (i.ne.j) write(*,*)'i,j,pk',i,j,Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*pk(r)
!    if (i.eq.j) write(*,*)'i,j,pk',i,j,Glob_PseudoCharge0*Glob_PseudoCharge(i)*pk(r)

    if (i.ne.j) Dk(r)=Dk(r)+Glob_PseudoCharge(i)*Glob_PseudoCharge(j)*pk(r)

    if (i.eq.j) Dk(r)=Dk(r)+Glob_PseudoCharge0*Glob_PseudoCharge(i)*pk(r)

  enddo

 endif

if (grad_l) then

pl=ZERO

  do r=1,n*(n+1)/2

    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    temp4=ZERO

    do k=1,3*n*(3*n+1)/2

!dnu3
      temp1=temp1+dnu3l(k)*trans(k,r)
!dnu6
      temp2=temp2+dnu6l(k)*trans(k,r)
!dnu7
      temp3=temp3+dnu7l(k)*trans(k,r)  
!dnu8
      temp4=temp4+dnu8l(k)*trans(k,r)  

   enddo

pl(r)=pl(r)+TWO*(PI**((THREE*n-ONE)/TWO))&
    *(det_tAkl**(-THREE/TWO))&
    *(Xi**(-ONE/TWO))&
    *((ONEHALF*(Xi**(-ONE))*resultJl(r)-THREEHALF*resultAl(r))&
      *(ONEHALF*nu3&
        -(ONE/SIX)*(Xi**(-ONE))*(nu6+nu7)&
        +(ONE/TEN)*(Xi**(-TWO))*nu8)&
      -ONEHALF*temp1&        
      -(ONE/SIX)*(Xi**(-TWO))*resultJl(r)*(nu6+nu7)&
      +(ONE/SIX)*(Xi**(-ONE))*(temp2+temp3)&
      +(ONE/FIVE)*(Xi**(-THREE))*resultJl(r)*nu8&
      -(ONE/TEN)*(Xi**(-TWO))*temp4)

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

!write(*,*)'Dk',Dk
!write(*,*)'Dl',Dl

Hkl=Tkl+Vkl

!write(*,*)'Skl,Tkl,Vkl,Hkl',Skl,Tkl,Vkl,Hkl

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
