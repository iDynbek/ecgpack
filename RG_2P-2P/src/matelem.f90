module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=0 and L=1 Gaussians.
  use globvars
  implicit none

contains

  subroutine MatrixElementsL1ForExpcValsP(m_k, m_l, mm_k, mm_l, vechLk, vechLl, Pbra, Pket, &
                                          Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
                                          MVkl, drach_MVkl, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, prvalkl, &
                                          NumCFGridPoints, CFGrid, &
                                          CFkl, NumDensGridPoints, DensGrid, Denskl, AreCorrFuncNeeded, ArePartDensNeeded)

!Arguments
    integer,intent(in)       :: m_k,m_l,mm_k,mm_l
    real(wp),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)   :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
    real(wp),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl
    real(wp),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
    real(wp),intent(out)  :: deltarkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: drach_deltarkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: prvalkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: rmrmkl(Glob_n,Glob_n,Glob_n,Glob_n)
    integer,intent(in)       :: NumCFGridPoints,NumDensGridPoints
    real(wp),intent(in)   :: CFGrid(2,NumCFGridPoints),DensGrid(2,NumDensGridPoints)
    real(wp),intent(out)  :: CFkl(Glob_n*(Glob_n+1)/2,NumCFGridPoints)
    real(wp),intent(out)  :: Denskl(Glob_n+1,NumDensGridPoints)
    logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded

    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
    integer           n,np
    integer           tvk(nn),tvl(nn),tbk(nn),tbl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn),tvl8(nn)
    real(wp)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn),rmkl1(nn,nn),Vkl1,Hkl1
    real(wp)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn), inv_tAkltAl(nn,nn)
    real(wp)       eta2(nn,nn),inv_tAkltAlM(nn,nn),eta22(nn,nn)
    real(wp)       eta(nn,nn)
    real(wp)       W1(nn,nn),W2(nn,nn),W3(nn,nn),W4(nn,nn),W5(nn,nn),W6(nn,nn),W7(nn,nn)
    real(wp)       W44(nn,nn),W55(nn,nn),W77(nn,nn)
    real(wp)       W4b(nn,nn),W5b(nn,nn),W7b(nn,nn)
    real(wp)       W44b(nn,nn),W55b(nn,nn),W77b(nn,nn),temp444,temp4444,temp444b,temp4444b
    real(wp)       inv_tAkltvl(nn),tvkinv_tAkl(nn),tvkinv_tAkltAlM(nn),u1(nn)
    real(wp)       inv_tAkltbl(nn),tbkinv_tAkl(nn),tbkinv_tAkltAlM(nn),u11(nn)
    real(wp)       temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp44
    real(wp)       temp55,temp66,temp55b,temp66b,temp44b,temp4b,temp5b,temp6b,temp1b,temp11b
    real(wp)       temp10,temp11,temp12,temp13,temp14,threshold,tr1, tr2, tr3, tr4,tr4vkvl,tr4bkbl,tr4bkvl,tr4vkbl
    real(wp)       det_Lk, det_Ll, det_tAkl,tau1,tau2,tau3,inv_tau3 ,V2kl,tau22,tau33,m
    integer           i,j,k,t,indx,p,q
    real(wp)       TrAJ(nn,nn),sqrtTrAJ(nn,nn),TrAJAJ(nn,nn,nn,nn)
    real(wp)       jAj(nn,nn,nn,nn),jAtvl(nn,nn),tvkAj(nn,nn),Mass_For_Darwin(0:nn)

    real(wp)   m1,m2,m3,tau331,tau332,tau333,tau334,tau221,tau222,tau223,tau224,temp38,temp338,temp36,temp336,temp37,temp337
    real(wp)   temp31,temp331,temp32,temp332,temp33,temp333,temp34,temp334,temp35,temp335,templast
    real(wp)   eta221(nn,nn),eta222(nn,nn),eta223(nn,nn),eta224(nn,nn),u111(nn)
    real(wp)   temp441,temp442,temp443,temp4440,temp4441,temp4442,h,term1,term2
    real(wp)   inv_tAkltAk(nn,nn),inv_tAkltAkM(nn,nn)
    real(wp)   inv_tAkltbk(nn),tvlinv_tAkl(nn),tvlinv_tAkltAkM(nn),tvlinv_tAkltAlM(nn)
    real(wp)   tbltbk(nn,nn),tvltvk(nn,nn),tbltvk(nn,nn),tvltbk(nn,nn)

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
    m=m1-m3  !look formula 40 in document
!Evaluating overlap
!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2Raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
    temp1=FOUR*det_tAkl*sqrt(det_tAkl)
    Skl=Glob_PiRaised3n2*m/temp1

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
    h=tau3*tau22+tau33*tau2-tau333*tau224-tau334*tau223
!Evaluating the kinetic energy
    Tkl=Skl*(SIX*tau1+FOUR*h/m)
    Vkl=ZERO
    Vkl1=ZERO
    temp5=Skl*TWO
    temp1=temp5/Glob_SqrtPi
    temp8=Skl/(Glob_Pi*Glob_SqrtPi)
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
      eta(i,i)=temp4*temp44-temp443*temp444
      term1=tau3*temp44+tau33*temp4-tau333*temp444-tau334*temp443
      term2=temp4*temp44-temp443*temp444
!rm2kl(i,i)=temp5*(ONE-TWO*ONETHIRD*term1/(m*temp2) + EIGHT*ONEFIFTH*term2/(THREE*m*temp2*temp2))/temp2
      rmkl(i,i)=temp1*(ONE-ONETHIRD*term1/(m*temp2) + ONEFIFTH*term2/(m*temp2*temp2))/temp3
!rmkl(i,i)=ME_over_rij(i,i,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
      Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*rmkl(i,i)
!Vkl1=Vkl1+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge0)*rmkl1(i,i)
      rkl(i,i)= temp1*temp3*(ONE+ONETHIRD*term1/(m*temp2) - ONEFIFTH*term2/(THREE*m*temp2*temp2))
!r2kl(i,i)=Skl*THREEHALF*temp2*(ONE+TWO*ONETHIRD*term1/(m*temp2))
      temp10=temp8/(temp2*temp3)
      deltarkl(i,i)=temp10*(ONE-term1/(m*temp2)+term2/(THREE*m*temp2*temp2))
!prvalkl(i,i)=PI*temp10*( TWO*(Glob_EulerConst+log(temp2))*(ONE-term1/(m*temp2)+term2/(m*temp2*temp2)) &
!+ FOUR*(term1-TWO*term2/temp2)/(THREE*m*temp2)+EIGHT*term2/(15*m*temp2*temp2) )
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
      eta(j,i)=temp4*temp44-temp443*temp444
      eta(i,j)=temp4*temp44-temp443*temp444
      term1=tau3*temp44+tau33*temp4-tau333*temp444-tau334*temp443
      term2=temp4*temp44-temp443*temp444
!rm2kl(j,i)=temp5*(ONE-TWO*ONETHIRD*term1/(m*temp2) + EIGHT*ONEFIFTH*term2/(THREE*m*temp2*temp2))/temp2
!rm2kl(i,j)=rm2kl(j,i)
      rmkl(j,i)=temp1*(ONE-ONETHIRD*term1/(m*temp2)+ONEFIFTH*term2/(m*temp2*temp2))/temp3
!rmkl(j,i)=ME_over_rij(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
      rmkl(i,j)=rmkl(j,i)
!rmkl1(i,j)=rmkl1(j,i)
      Vkl=Vkl+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*rmkl(j,i)
!Vkl1=Vkl1+ScaledChargeProd(Glob_PseudoCharge(i),Glob_PseudoCharge(j))*rmkl1(j,i)
      rkl(j,i)=temp1*temp3*(ONE+ONETHIRD*term1/(m*temp2) - ONEFIFTH*term2/(THREE*m*temp2*temp2))
      rkl(i,j)=rkl(j,i)
!r2kl(j,i)=Skl*THREEHALF*temp2*(ONE+TWO*ONETHIRD*term1/(m*temp2))
!r2kl(i,j)=r2kl(j,i)
      temp10=temp8/(temp2*temp3)
      deltarkl(j,i)=temp10*(ONE-term1/(m*temp2)+term2/(THREE*m*temp2*temp2))
      deltarkl(i,j)=deltarkl(j,i)
!prvalkl(j,i)=PI*temp10*( TWO*(Glob_EulerConst+log(temp2))*(ONE-term1/(m*temp2)+term2/(m*temp2*temp2)) &
!+ FOUR*(term1-TWO*term2/temp2)/(THREE*m*temp2)+EIGHT*term2/(15*m*temp2*temp2) )
!prvalkl(i,j)=prvalkl(j,i)
    enddo
    enddo
    Hkl=Tkl+Vkl

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
    Darwinkl=-Darwinkl*Glob_Pi/2

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

!Evaluation of correlation functions
    if (AreCorrFuncNeeded) then
      temp1=Skl/(Glob_Pi*Glob_SqrtPi)
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
      temp1=Skl/(Glob_Pi*Glob_SqrtPi)
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

  end subroutine MatrixElementsL1ForExpcValsP

  subroutine symmetrize_matrix(W)
!subroutine symmetrize_matrix makes an arbitrary square matrix W
!symmetric by the following procedure:
!W = (1/2)*(W + W')
!Input:
!   W :: n x n real matrix

    integer, parameter :: nn = Glob_MaxAllowedNumOfPseudoParticles
    real(wp)           W(nn, nn), t
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
    real(wp) ScaledChargeProd,q1,q2,x
    x=q1*q2
    if (x<0.0_wp) then
      ScaledChargeProd=x*Glob_AttractionScalingParam
    else
      if ((q1>0.0_wp).and.(q2>0.0_wp)) then
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
    real(wp)   SG_ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn)
    integer       i,j
    real(wp)   t_V,Skl
!Local variables:
    integer       p,q,s,n
    real(wp)   temp1,temp2,temp3,temp4,temp5,temp6
    real(wp)   AX(nn,nn),AY(nn,nn)
    real(wp)   Aj(nn),AjX(nn),AjY(nn),AXAj(nn),AYAj(nn)
    real(wp)   t_J,t_X,t_Y
    real(wp)   t_XJ,t_YJ,t_XY
    real(wp)   t_XYJ,t_YXJ,det_tAkl

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

    temp1=Glob_PiRaised3n2/(Glob_SqrtPi*det_tAkl**(THREEHALF))
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
    real(wp)   SG_ME_rXr_rYr_rZr_over_rij
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),Z(nn,nn),inv_tAkl(nn,nn)
    integer       i,j
!Local variables:
    integer       p,q,s,n
    real(wp)   temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12
    real(wp)   AX(nn,nn),AY(nn,nn),AZ(nn,nn)
    real(wp)   Aj(nn),AjX(nn),AjY(nn),AjZ(nn),AXAj(nn),AYAj(nn),AZAj(nn)
    real(wp)   AZAY(nn,nn),AYAZ(nn,nn),AZAX(nn,nn),AXAZ(nn,nn),AXAY(nn,nn),AYAX(nn,nn)
    real(wp)   AZAYAj(nn),AYAZAj(nn),AZAXAj(nn),AXAZAj(nn),AXAYAj(nn),AYAXAj(nn)
    real(wp)   t_J,t_X,t_Y,t_Z
    real(wp)   t_XJ,t_YJ,t_ZJ
    real(wp)   t_XY,t_ZY,t_ZX,t_YX,t_ZYX,t_YZX
    real(wp)   t_XYJ,t_YXJ,t_ZYJ,t_YZJ,t_XZJ,t_ZXJ
    real(wp)   t_ZYXJ,t_YZXJ,t_YXZJ,t_XYZJ,t_ZXYJ,t_XZYJ
    real(wp)   det_tAkl,term1,term2,term3,term4,term5,term6,term7,term8

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

    temp1=Glob_PiRaised3n2/(Glob_SqrtPi*det_tAkl**(THREEHALF))
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
    real(wp)   ME_rXr_over_rij
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn)
    real(wp)   t_J,t_X,t_XJ
    real(wp)   t_JV1,t_XV1,t_JXV1,t_XJV1
    real(wp)   t_JV2,t_XV2,t_JXV2,t_XJV2
    real(wp)   t_JV5,t_XV5,t_JXV5,t_XJV5
    real(wp)   t_JV6,t_XV6,t_JXV6,t_XJV6
    real(wp)   Ajtvl,Ajtbl
    real(wp)   temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
    real(wp)   temp11,temp22,temp33,mu,mX,mXJ,u,Xmu,muXJ

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
    m=m1-m3  !look formula 40 in document

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

    mu=tau3*t_JV1+tau33*t_JV2-tau333*t_JV5-tau334*t_JV6
    mX=tau3*t_XV1+tau33*t_XV2-tau333*t_XV5-tau334*t_XV6
    u=t_JV1*t_JV2-t_JV5*t_JV6
    mXJ=tau3*(t_XJV1+t_JXV1)+tau33*(t_XJV2+t_JXV2)-tau333*(t_XJV5+t_JXV5)-tau334*(t_XJV6+t_JXV6)
    muXJ=t_JV2*(t_XJV1+t_JXV1)+t_JV1*(t_XJV2+t_JXV2)-t_JV6*(t_XJV5+t_JXV5)-t_JV5*(t_XJV6+t_JXV6)
    Xmu=t_XV2*t_JV1+t_XV1*t_JV2-t_XV6*t_JV5-t_XV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))
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
    real(wp)   ME_over_rij
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)   temp1,temp2,temp3,temp11,temp22,temp33,mu,u

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
    m=m1-m3  !look formula 40 in document

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

    mu=tau3*t_JV1+tau33*t_JV2-tau333*t_JV5-tau334*t_JV6
    u=t_JV1*t_JV2-t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij

  function ME_d_X_over_rij_d(i,j,X,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)

    real(wp)   ME_d_X_over_rij_d
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),tAl(nn,nn),tAk(nn,nn),inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

!Local variables:
    integer       c,s,p,q,n,k
    real(wp)   tvkXtAl(nn),tbkXtAl(nn)
    real(wp)   temp1,temp2,temp3,temp4,temp5
    real(wp)   tAkX(nn,nn),tAkXtAl(nn,nn),tAkXtvl(nn),tAkXtbl(nn),XtAl(nn,nn)

!Doing multiplication Xtvl=X*tvl,tvkX=tvk*X
    n=Glob_n
    do s=1,n
    do c=1,n
      temp1=ZERO
      temp2=ZERO
      do k=1,n
        temp1=temp1+tAk(c,k)*X(k,s)
        temp2=temp2+X(c,k)*tAl(k,s)
      enddo
      tAkX(c,s)=temp1
      XtAl(c,s)=temp2
    enddo
    enddo

!Doing multiplication tAlXtAk=tAlX*tAk,tAkXtAl=tAkX*tAl
    do s=1,n
    do c=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+tAkX(c,k)*tAl(k,s)
      enddo
      tAkXtAl(c,s)=temp1
    enddo
    enddo

    call symmetrize_matrix(tAkXtAl)

    do s=1,n
      temp1=ZERO
      temp2=ZERO
      do c=1,n
        temp1=temp1+tAkX(c,s)*tvl(c)
        temp2=temp2+tAkX(c,s)*tbl(c)
      enddo
      tAkXtvl(s)=temp1
      tAkXtbl(s)=temp2
    enddo

    do p=1,n
      temp1=ZERO
      temp2=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*XtAl(q,p)
        temp2=temp2+tbk(q)*XtAl(q,p)
      enddo
      tvkXtAl(p)=temp1
      tbkXtAl(p)=temp2
    enddo

    temp1=4*ME_rXr_over_rij(i,j,tAkXtAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
    temp2=-2*ME_over_rij_tvl(i,j,inv_tAkl,det_tAkl,tvk,tAkXtvl,tbk,tbl)
    temp3=-2*ME_over_rij_tbl(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tAkXtbl)
    temp4=-2*ME_over_rij_tvk(i,j,inv_tAkl,det_tAkl,tvkXtAl,tvl,tbk,tbl)
    temp5=-2*ME_over_rij_tbk(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbkXtAl,tbl)
    ME_d_X_over_rij_d=temp1+temp2+temp3+temp4+temp5

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
    real(wp)   ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
    integer       p,q,n,s,k
    real(wp)   Aj(nn),AjX(nn),AjY(nn)
    real(wp)   t_J,t_X,t_XJ,t_Y,t_YJ,t_XYJ
    real(wp)   t_JV1,t_XV1,t_JXV1,t_XJV1,t_YV1,t_JYV1,t_YJV1,t_XYV1,t_YXV1,t_XYJV1,t_YXJV1,t_XJYV1,t_YJXV1,t_JXYV1,t_JYXV1
    real(wp)   t_JV2,t_XV2,t_JXV2,t_XJV2,t_YV2,t_JYV2,t_YJV2,t_XYV2,t_YXV2,t_XYJV2,t_YXJV2,t_XJYV2,t_YJXV2,t_JXYV2,t_JYXV2
    real(wp)   t_JV5,t_XV5,t_JXV5,t_XJV5,t_YV5,t_JYV5,t_YJV5,t_XYV5,t_YXV5,t_XYJV5,t_YXJV5,t_XJYV5,t_YJXV5,t_JXYV5,t_JYXV5
    real(wp)   t_JV6,t_XV6,t_JXV6,t_XJV6,t_YV6,t_JYV6,t_YJV6,t_XYV6,t_YXV6,t_XYJV6,t_YXJV6,t_XJYV6,t_YJXV6,t_JXYV6,t_JYXV6
    real(wp)    Ajtvl,Ajtbl,tvkAj,tbkAj,prod
    real(wp)    temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp44,temp77,temp88,temp99
    real(wp)    temp11,temp22,temp33,mu,mX,mXJ,u,Xmu,muXJ,mY,mYJ,Ymu,muYJ,t_XYJV11,t_XYJV22,t_XYJV55,t_XYJV66,YX,mXYJ
    real(wp)    temp2Y,temp22Y,t_XY,temp55,temp66,AXAj(nn),AYAj(nn),muYX,muXYJ,XYJ,YXJ,t_YXJ,XJYJ,mYX
    real(wp)    tvkinv_tAklX(nn),tbkinv_tAklX(nn),tvkinv_tAklY(nn),tbkinv_tAklY(nn),AY(nn,nn),AX(nn,nn)
    real(wp)    AXinv_tAkltvl(nn),AXinv_tAkltbl(nn),AYinv_tAkltvl(nn),AYinv_tAkltbl(nn)
    real(wp)    term1,term2,term3,term4,term5,term6,term7,term8,term9,term10,term11,term12,term13,temp1Y,temp11Y

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
    m=m1-m3  !look formula 40 in document

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
    mu=tau3*t_JV1+tau33*t_JV2-tau333*t_JV5-tau334*t_JV6
    mX=tau3*t_XV1+tau33*t_XV2-tau333*t_XV5-tau334*t_XV6
    mY=tau3*t_YV1+tau33*t_YV2-tau333*t_YV5-tau334*t_YV6
    mYX=tau3*t_YXV1+tau33*t_YXV2-tau333*t_YXV5-tau334*t_YXV6
    u=t_JV1*t_JV2-t_JV5*t_JV6
    mXJ=tau3*(t_XJV1+t_JXV1)+tau33*(t_XJV2+t_JXV2)-tau333*(t_XJV5+t_JXV5)-tau334*(t_XJV6+t_JXV6)
    mYJ=tau3*(t_YJV1+t_JYV1)+tau33*(t_YJV2+t_JYV2)-tau333*(t_YJV5+t_JYV5)-tau334*(t_YJV6+t_JYV6)
    muXJ=t_JV2*(t_XJV1+t_JXV1)+t_JV1*(t_XJV2+t_JXV2)-t_JV6*(t_XJV5+t_JXV5)-t_JV5*(t_XJV6+t_JXV6)
    muYJ=t_JV2*(t_YJV1+t_JYV1)+t_JV1*(t_YJV2+t_JYV2)-t_JV6*(t_YJV5+t_JYV5)-t_JV5*(t_YJV6+t_JYV6)
    Xmu=t_XV2*t_JV1+t_XV1*t_JV2-t_XV6*t_JV5-t_XV5*t_JV6
    Ymu=t_YV2*t_JV1+t_YV1*t_JV2-t_YV6*t_JV5-t_YV5*t_JV6

    mXYJ=tau3*t_XYJV11+tau33*t_XYJV22-tau333*t_XYJV55-tau334*t_XYJV66
    YX=t_YV2*t_XV1+t_YV1*t_XV2-t_YV6*t_XV5-t_YV5*t_XV6
    muYX=t_JV1*t_YXV2+t_JV2*t_YXV1-t_JV6*t_YXV5-t_JV5*t_YXV6
    muXYJ=t_JV1*t_XYJV22+t_JV2*t_XYJV11-t_JV6*t_XYJV55-t_JV5*t_XYJV66
    XYJ=t_XV1*t_YJV2+t_XV2*t_YJV1-t_XV5*t_YJV6-t_XV6*t_YJV5
    YXJ=t_YV1*t_XJV2+t_YV2*t_XJV1-t_YV5*t_XJV6-t_YV6*t_XJV5
XJYJ=(t_XJV1+t_JXV1)*(t_YJV2+t_JYV2)+(t_XJV2+t_JXV2)*(t_YJV1+t_JYV1)-(t_XJV5+t_JXV5)*(t_YJV6+t_JYV6)-(t_XJV6+t_JXV6)*(t_YJV5+t_JYV5)

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))
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
    real(wp)   rPr_rQr
!arguments
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
    integer       tvk(nn),tbk(nn)
    real(wp)   P(nn,nn),Q(nn,nn),inv_tAkl(nn,nn),tau3,tau33,tau333,tau334,det_tAkl
    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),tvkinv_tAkl(nn),tbkinv_tAkl(nn)

!declaration
    integer       i,j,n,k
    real(wp)   temp1,temp2,temp3,temp4,P1,P2,P5,P6,Q1,Q2,Q5,Q6,trP,trQ,trQP
    real(wp)   PQ1,PQ2,PQ5,PQ6,QP1,QP2,QP5,QP6
    real(wp)   inv_tAklP(nn,nn),inv_tAklQ(nn,nn),inv_tAklQP(nn,nn)
    real(wp)   tvkinv_tAklP(nn),tbkinv_tAklP(nn),tvkinv_tAklQ(nn),tbkinv_tAklQ(nn)
    real(wp)   Pgamma,Qgamma,PQgamma,QPgamma

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

    temp1=Glob_PiRaised3n2/(FOUR*det_tAkl**(THREEHALF))
    temp2=THREEHALF*(tau3*tau33-tau333*tau334)*(trQP+THREEHALF*trQ*trP)
    temp3=Q1*P2+Q2*P1-Q5*P6-Q6*P5
    Pgamma=P1*tau3+P2*tau33-P5*tau333-P6*tau334
    Qgamma=Q1*tau3+Q2*tau33-Q5*tau333-Q6*tau334
    PQgamma=PQ1*tau3+PQ2*tau33-PQ5*tau333-PQ6*tau334
    QPgamma=QP1*tau3+QP2*tau33-QP5*tau333-QP6*tau334

    rPr_rQr= temp1*(THREEHALF*trQ*Pgamma+THREEHALF*trP*Qgamma+PQgamma+QPgamma+temp2+temp3)

  end function rPr_rQr

  function dXddYd(X,Y,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)
    real(wp)   dXddYd
!arguments
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
    integer       tvk(nn),tbk(nn),tvl(nn),tbl(nn)
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),tau3,tau33,tau333,tau334,det_tAkl
    real(wp)   tAk(nn,nn),tAl(nn,nn)
    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),tvkinv_tAkl(nn),tbkinv_tAkl(nn)
    real(wp)   inv_tAkltAl(nn,nn),inv_tAkltAk(nn,nn)
!declaration
    integer       i,j,n,k
    real(wp)   P(nn,nn),Q(nn,nn),trAkX,trAlY
    real(wp)   temp1,temp2,temp3,temp4,temp33,temp44,P1,P2,P5,P6,Q1,Q2,Q5,Q6,trP,trQ,trQP
    real(wp)   PQ1,PQ2,PQ5,PQ6,QP1,QP2,QP5,QP6
    real(wp)   inv_tAklP(nn,nn),inv_tAklQ(nn,nn),inv_tAklQP(nn,nn),tvkXtAk(nn),tbkXtAk(nn)
    real(wp)   tvkinv_tAklP(nn),tbkinv_tAklP(nn),tvkinv_tAklQ(nn),tbkinv_tAklQ(nn)
    real(wp)   XtAk(nn,nn),YtAl(nn,nn)
    real(wp)   Pgamma,Qgamma,PQgamma,QPgamma,big_term1,big_term2,big_term3,big_term4
    real(wp)   term1,term2,term3,term4,term5,term6,term7,term8,prod,gamma
    real(wp)   term9,term10,term11,term12,term13,term14,term15,term16
    real(wp)   term17,term18,term19,term20,term21,term22,term23,term24,term25
    real(wp)   tAkX(nn,nn),tAlY(nn,nn),Ytvl(nn),Ytbl(nn),Xtbk(nn),tvkX(nn),tbkX(nn),tvlY(nn)
    real(wp)   term3_1,term3_2,term4_1,term4_2,term8_1,term8_2,term9_1,term9_2,term11_1,term11_2,term16_1,term16_2
    real(wp)   term12_1,term12_2,term17_1,term17_2,term13_1,term14_2,term18_2,term19_1,temp5,temp6

    n=Glob_n
    prod=Glob_PiRaised3n2/(FOUR*det_tAkl**(THREEHALF))
    gamma=tau3*tau33-tau333*tau334

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
    Qgamma=Q1*tau3+Q2*tau33-Q5*tau333-Q6*tau334
    Pgamma=P1*tau3+P2*tau33-P5*tau333-P6*tau334

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
    term3=24*trAkX*prod*(tau33*term3_1-tau333*term3_2)
    term4=24*trAkX*prod*(tau3*term4_1-tau334*term4_2)
    term6=-24*trAlY*prod*(THREEHALF*trP*gamma+Pgamma)
    term7=16*rPr_rQr(P,Q,tvk,tbk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,tvkinv_tAkl,tbkinv_tAkl,inv_tAkltvl,inv_tAkltbl)
    term8=-24*trP*prod*(tau33*term3_1-tau333*term3_2)-16*prod*(tau33*term8_1+P1*term3_1-P6*term3_2-tau333*term8_2)
    term9=-24*trP*prod*(tau3*term4_1-tau334*term4_2)-16*prod*(P2*term4_1+tau3*term9_1-tau334*term9_2-P5*term4_2)
    term11=24*trAlY*prod*(tau33*term11_1-tau334*term11_2)
    term12=-24*trQ*prod*(tau33*term11_1-tau334*term11_2)-16*prod*(tau33*term12_1+Q1*term11_1-tau334*term12_2-Q5*term11_2)
    term13=16*prod*(tau33*term13_1-term11_2*term3_2)
    term14=16*prod*(term11_1*term4_1-tau334*term14_2)
    term16=24*trAlY*prod*(tau3*term16_1-tau333*term16_2)
    term17=-24*trQ*prod*(tau3*term16_1-tau333*term16_2)&
            -16*prod*(Q2*term16_1+tau3*term17_1-Q6*term16_2-tau333*term17_2)
    term18=16*prod*(term3_1*term16_1-tau333*term18_2)
    term19=16*prod*(tau3*term19_1-term4_2*term16_1)

    big_term1=term1+term2+term3+term4+term6
    big_term2=term7+term8+term9+term11
    big_term3=term12+term13+term14
    big_term4=term16+term17+term18+term19
    dXddYd= big_term1+big_term2+big_term3+big_term4

  end function dXddYd
  function ME_over_rij_tbk(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tbk
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvl(nn),tvk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tbk(nn)
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)    temp1,temp2,temp3,temp11,temp22,temp33,mu,u

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
    m=m1-m3  !look formula 40 in document

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

    mu=tau3*t_JV1+tau33*t_JV2-tau333*t_JV5-tau334*t_JV6
    u=t_JV1*t_JV2-t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tbk=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tbk

  function ME_over_rij_tvk(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tvk
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tvk(nn)
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)   temp1,temp2,temp3,temp11,temp22,temp33,mu,u

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
    m=m1-m3  !look formula 40 in document

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

    mu=tau3*t_JV1+tau33*t_JV2-tau333*t_JV5-tau334*t_JV6
    u=t_JV1*t_JV2-t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tvk=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tvk

  function ME_over_rij_tbl(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tbl
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tbl(nn)
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)   temp1,temp2,temp3,temp11,temp22,temp33,mu,u

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
    m=m1-m3  !look formula 40 in document

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

    mu=tau3*t_JV1+tau33*t_JV2-tau333*t_JV5-tau334*t_JV6
    u=t_JV1*t_JV2-t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tbl=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tbl

  function ME_over_rij_tvl(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tvl
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tvl(nn)
!Local variables:
    integer       p,q,n,k
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)   temp1,temp2,temp3,temp11,temp22,temp33,mu,u

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
    m=m1-m3  !look formula 40 in document

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

    mu=tau3*t_JV1+tau33*t_JV2-tau333*t_JV5-tau334*t_JV6
    u=t_JV1*t_JV2-t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tvl=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tvl

  function ME_dXd(X,tvk,tvl,inv_tAkltvl,inv_tAkl,tAk,tAl,inv_tAkltAl,Skl,tau3)
    real(wp)   ME_dXd
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
    real(wp)   X(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn),inv_tAkltAl(nn,nn)
    integer       i,j,n,k,tvk(nn),tvl(nn)
    real(wp)   inv_tAkltAlX(nn,nn),inv_tAkltAlXtAk(nn,nn),tvkinv_tAkltAlX(nn),inv_tAkltvl(nn)
    real(wp)   temp1, Skl,tau,tau1,tau2,tau3
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
!Depending on the kind parameter (wp=4,8,10,16 - double, extended, or quadruple
!precision) for real numbers, a different xmin is used.
!In all cases the accuracy is close to machine precision corresponding to
!that kind parameter (1-2 last significant figures might be inaccurate in
!the worst case).
    real(wp) ftransaux,x
    real(wp),parameter :: xmin_4=0.30_wp !for single precision
    real(wp),parameter :: xmin_8=0.27_wp !for double precision
    real(wp),parameter :: xmin_10=0.2_wp !for extended precision
    real(wp),parameter :: xmin_16=0.065_wp !for quadruple precision
!Local variables
    real(wp) x2,ax,t,xmin

    ax=abs(x)
    selectcase (wp)
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
      t=(524288.0_wp/50702925.0_wp)
      t=(262144.0_wp/22309287.0_wp)+x2*t
      t=(65536.0_wp/4849845.0_wp)+x2*t
      t=(32768.0_wp/2078505.0_wp)+x2*t
      t=(2048.0_wp/109395.0_wp)+x2*t
      t=(1024.0_wp/45045.0_wp)+x2*t
      t=(256.0_wp/9009.0_wp)+x2*t
      t=(128.0_wp/3465.0_wp)+x2*t
      t=(16.0_wp/315.0_wp)+x2*t
      t=(8.0_wp/105.0_wp)+x2*t
      t=(2.0_wp/15.0_wp)+x2*t
      t=(1.0_wp/3.0_wp)+x2*t
      ftransaux=x*t
    else
      t=sqrt(1.0_wp-x*x)
      ftransaux=(ax-t*acos(t))/(ax*x)
    endif

  end function ftransaux

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
    !   m_k, mm_k :: integers that determine which pseudoparticles carry l=1 momentum
    !   vechLk :: Array of length (n(n+1)/2) of
    !     exponential parameters.
    !   P  :: The symmetry permutation matrix of size n x n
    !Output:
    !   Skk         ::        Overlap matrix element (normalized)

    !Arguments
    integer,intent(in)          :: m_k, mm_k
    real(wp),intent(in)      :: vechLk(Glob_np)
    real(wp),intent(in)      :: P(Glob_n,Glob_n)
    real(wp),intent(out)     :: Skk

    !Parameters (These are needed to declare static arrays. Using static
    !arrays makes the function call a little faster in comparison with
    !the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

    !Local variables
    integer           n, np
    real(wp)       vl(nn), bl(nn), vkinv_tAkl(nn), bkinv_tAkl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn)
    real(wp)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)
    real(wp)       W1(nn,nn)
    real(wp)       temp1, temp2
    real(wp)       det_tAkl
    real(wp)       tau3, tau33, tau333, tau334, m, m1, m3
    integer           i,j,k,q,t,indx

    n=Glob_n
    np=Glob_np

    !print*, "m_k = ", m_k
    !print*, "mm_k = ", mm_k
    !print*, "vechLk = ", vechLk
    !print*, "P = ", P

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

    !print*, "det_tAkl = ", det_tAkl

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
    Skk=Glob_PiRaised3n2*m/temp1

    !print*, "temp1 = ", temp1
    !print*, "m = ", m
    !print*, "Skk = ", Skk
    !print*, "Glob_PiRaised3n2 = ", Glob_PiRaised3n2

  end subroutine overlapMatrixElementsLP

  subroutine spinPreCalc(n, nFactorial, parityFactor, SSFmassChargeCoefficient, SSNCmassChargeCoefficient, &
                         SOmassChargeCoefficient, AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, ketMatrix,&
                         spatialYoung0, spatialYoung1,SSNCspinME, SiMinusME, SiPlusME, SziME, spinFreeME, SiSjME)
    use spinStuff
    implicit none

    !input vars:
    character(len = maxLen), intent(in) :: spatialYoung0, spatialYoung1
    integer, intent(in) :: n, nFactorial

    !output vars:
    real(wp), dimension(nFactorial), intent(out) :: parityFactor
    real(wp), dimension(n, n, 4), intent(out) :: SOmassChargeCoefficient, AMMmassChargeCoefficient, &
                                                    AMMFinmassChargeCoefficient
    real(wp), dimension(n, n), intent(out) :: SSNCmassChargeCoefficient, SSFmassChargeCoefficient
    real(wp), dimension(n, n, nFactorial), intent(out) :: ketMatrix
    real(wp), dimension(n, n, nFactorial), intent(out) :: SSNCspinME, SiSjME
    real(wp), dimension(nFactorial, 2), intent(out) :: spinFreeME
    real(kind = wp), dimension(n, nFactorial), intent(out) :: SiMinusME, SiPlusME, SziME

    ! local variables
    integer :: i, j, k, l, m
    character(len = maxLen) :: mySpatialYoung0, mySpatialYoung1
    integer, dimension(nFactorial) :: parities
    integer, dimension(n, n, nFactorial) :: allPermutations
    real(kind = wp), dimension(:), allocatable :: finalSpinFunction0, finalSpinFunction1
    integer, dimension(:, :), allocatable :: primitives0, primitives1
    integer :: numberOfPrimitives0, numberOfPrimitives1

    SOmassChargeCoefficient = ZERO
    SSNCmassChargeCoefficient = ZERO
    SSFmassChargeCoefficient = ZERO

    do i = 1, n
      do j = 1, n
        SSFmassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                         (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * EIGHT * Glob_Pi / THREE
      enddo
    enddo

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

    !Singlet wf
    call getSpinFunction(n, nFactorial, mySpatialYoung1, allPermutations, parities, &
                         finalSpinFunction0, primitives0, numberOfPrimitives0)

    !Triplet wf
    call getSpinFunction(n, nFactorial, mySpatialYoung0, allPermutations, parities, &
                         finalSpinFunction1, primitives1, numberOfPrimitives1)

    call getSpinOpMeanValues(n, nFactorial, allPermutations, finalSpinFunction0, finalSpinFunction1, primitives0, primitives1, &
                             numberOfPrimitives0, numberOfPrimitives1, spinFreeME, SSNCspinME, SiMinusME, SiPlusME, SziME, SiSjME)

    ketMatrix = ZERO
    do i = 1, nFactorial
      do k = 1, n
        do l = 1, n
          !note the transposition here
          ketMatrix(k, l, i) = real(allPermutations(l, k, i), kind=wp)
        enddo
      enddo

    enddo

    do i = 1, nFactorial
      parityFactor(i) = real(parities(i), kind=wp)
    enddo

  end subroutine spinPreCalc

  subroutine spinDependentMatrixElements(m_k, m_l, mm_k, mm_l, vechLk, vechLl, Pket, &
                                         SOspinME, SSNCspinME, SSNCmassChargeCoefficient, SOmassChargeCoefficient, &
                                         AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, &
                                         SSNCkl, SO1kl, SO2kl, AMM1kl, AMM2kl, AMM1finkl, AMM2finkl)
    !This subroutine computes symmetry adapted matrix element
    !with two real L=1 correlated Gaussians. These matrix element
    !is used in calculations of expectation values.

    !Input:
    !   m_k,m_l,mm_k, mm_l :: integers that determine which x or y-components is in the
    !                premultiplier of the Gaussian
    !   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.

    !Output (all matrix elements are computed with normalized functions):

    !   SSNCkl :: Non-contact spin-spin term (without the factor of alpha**2)
    !   SO1kl, SO2kl  :: Spin-Orbit corrections (without the factor of alpha**2)
    !         1 and 2 stay for spin-same orbit and spin-another orbit contributions
    !   AMM1kl, AMM2kl  :: AMM corrections (without the factor of alpha**2)
    !         1 and 2 stay for spin-same orbit and spin-another orbit contributions

    !Arguments
    integer,intent(in)       :: m_k, m_l, mm_k, mm_l
    real(wp),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)   :: Pket(Glob_n,Glob_n)

    real(wp), intent(out)  :: SO1kl, SO2kl, AMM1kl, AMM2kl, AMM1finkl, AMM2finkl, SSNCkl
    !Parameters (These are needed to declare static arrays. Using static
    !arrays makes the function call a little faster in comparison with
    !the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_MaxAllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2
    real(wp),intent(in)   :: SSNCspinME(Glob_n, Glob_n), &
                                SOspinME(Glob_n), &
                                SOmassChargeCoefficient(Glob_n, Glob_n, 4), &
                                AMMmassChargeCoefficient(Glob_n, Glob_n, 4), &
                                AMMFinmassChargeCoefficient(Glob_n, Glob_n, 4), &
                                SSNCmassChargeCoefficient(Glob_n, Glob_n)

    !Local variables
    integer           n, np
    integer           tvk(nn),tvl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
    real(wp)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)

    real(wp)       W1(nn,nn)
    real(wp)       temp1, temp2, temp3, temp1010, temp1001, temp0110, temp0101, t1, t2, det_tAkl
    integer :: i, j, k, indx

    integer :: pm_k, pm_l, pmm_k, pmm_l ! new non-zero components of v_k and v_l
    real(wp) :: commonFactor, gamma, gamma_diag, jiAlAklinvVl, &
                   jjAlAklinvVl, localEps

    ! V-quantities
    real(wp) :: VkAklinvVl, jiAkAklinvVl, jiAlAklinvVk, jiAklinvVk, jiAklinvVl, &
                   jjAkAklinvVl, jjAlAklinvVk, jjAklinvVk, jjAklinvVl

    ! W-quantities
    real(wp) :: WkAklinvWl, jiAklinvWk, &
                   jiAkAklinvWl, jiAlAklinvWk, jiAklinvWl, &
                   jjAkAklinvWl, jjAlAklinvWk, jjAklinvWk, jjAklinvWl

    !mixed quantities
    real(wp) :: VkAklinvWl, WkAklinvVl

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

    !common factor (ONEHALF - for consistent normalization with Skl)

    commonFactor = ONEHALF * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
    SO1kl = ZERO
    SO2kl = ZERO

    AMM1kl = ZERO
    AMM2kl = ZERO
    AMM1finkl = ZERO
    AMM2finkl = ZERO

    SSNCkl = ZERO

    WkAklinvWl = inv_tAkl(pmm_k, pmm_l)
    WkAklinvVl = inv_tAkl(pmm_k, pm_l)
    VkAklinvWl = inv_tAkl(pm_k, pmm_l)
    VkAklinvVl = inv_tAkl(pm_k, pm_l)

    do indexI = 1, n

      ! gamma diagonal coefficient
      gamma_diag = ONE / sqrt(inv_tAkl(indexI, indexI))

      ! calculating all the traces we need
      ! tr(Axy') is computed as (y, Ax) everywhere
      ! variable names: jiAlAklinvVk = (j^i, A_l A_{kl}^(-1) v_k) (names doesn't account for permutations)
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

      temp1 = temp1010 + temp0101 - temp1001 - temp0110

      SO1kl = SO1kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
      AMM1kl = AMM1kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
      AMM1finkl = AMM1finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 1) * temp1

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

        temp1 = temp1010 + temp0101 - temp1001 - temp0110

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
        AMM2finkl = AMM2finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 2) * temp1

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

        temp1 = temp1010 + temp0101 - temp1001 - temp0110

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
        AMM2kl  = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
        AMM2finkl  = AMM2finkl + SOspinME(indexI) * AMMfinmassChargeCoefficient(indexI, indexJ, 3) * temp1

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

        temp1 = temp1010 + temp0101 - temp1001 - temp0110

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
        AMM2kl = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
        AMM2finkl = AMM2finkl + SOspinME(indexI) * AMMfinmassChargeCoefficient(indexI, indexJ, 4) * temp1

        if (Glob_selectTransition == 2 .or. indexJ <= indexI .or. abs(SSNCspinME(indexI, indexJ)) < localEps ) cycle
        !SSNC term
        !sum f1010 + f0101
        t1 = jiAklinvVk * jiAklinvVl + jjAklinvVk * jjAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
        t2 = jiAklinvWk * jiAklinvWl + jjAklinvWk * jjAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
        temp2 = -gamma**5 / FIVE * (t1 * WkAklinvWl + t2 * VkAklinvVl) + &
                TWO * gamma**7 / SEVEN * (t1 * t2)

        !sum f0110 + f1001
        t1 = jiAklinvVk * jiAklinvWl + jjAklinvVk * jjAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
        t2 = jiAklinvWk * jiAklinvVl + jjAklinvWk * jjAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
        temp3 =  -gamma**5 / FIVE * (t1  * WkAklinvVl + t2 * VkAklinvWl) + &
                TWO * gamma**7 / SEVEN * (t1 * t2)

        temp1 = ONETHIRD*(temp2 - temp3) !ONETHIRD - from common factor and spin part (1/sqrt(6))

        SSNCkl = SSNCkl + SSNCspinME(indexI, indexJ) * SSNCmassChargeCoefficient(indexI, indexJ) * temp1

      enddo ! indexJ cycle
    enddo ! indexI cycle

    SSNCkl = SSNCkl * commonFactor
    SO1kl = SO1kl * commonFactor
    SO2kl = SO2kl * commonFactor
    AMM1kl = AMM1kl * commonFactor
    AMM2kl = AMM2kl * commonFactor
    AMM1finkl = AMM1finkl * commonFactor
    AMM2finkl = AMM2finkl * commonFactor

  end subroutine spinDependentMatrixElements
end module matelem

