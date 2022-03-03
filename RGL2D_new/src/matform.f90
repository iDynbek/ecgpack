module matform
!Module matform contains procedures that form Hamiltonian
!and overlap matrices and related routines
use matelem
implicit none
 
contains 

subroutine StoreHS(i,j,Hij,Sij)
!Routine StoreHS stores calculated matrix elements of the
!Hamiltonian and the overlap in proper places of global
!arrays. Upon doing this, the routine nolmalizes 
!matrix elements. 
!Important comment:  i must be greater or equal to j. 
integer,intent(in)      ::  i,j
real(dprec),intent(in)  ::  Hij,Sij
real(dprec)             ::  f

select case (Glob_GSEPSolutionMethod) 
case('I')
  !Only the lower triangle of array Glob_H (including the diagonal) is 
  !used to store H-Glob_ApproxEnergy*S.
  !The entire array Glob_S is used to store S.
  if (i==j) then
	Glob_diagS(i)=Sij
	Glob_S(i,i)=ONE
	Glob_H(i,i)=Hij/Glob_diagS(i)-Glob_ApproxEnergy
  else
    f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
	Glob_S(i,j)=Sij*f
	Glob_S(j,i)=Glob_S(i,j)
	Glob_H(i,j)=(Hij-Glob_ApproxEnergy*Sij)*f
  endif
case('G')
  !In the case when Glob_GSEPSolutionMethod='G' the diagonal
  !of the Hamiltonian matrix is stored in Glob_diagH
  !The diagonal of the overlap is stored in Glob_diagS
  !The lower triangles of arrays Glob_H and 
  !Glob_S are used to store H and S
  if (i==j) then
	Glob_diagS(i)=Sij 
	Glob_diagH(i)=Hij/Glob_diagS(i)
  else
    f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
	Glob_S(i,j)=Sij*f
	Glob_H(i,j)=Hij*f
  endif
endselect
end subroutine StoreHS



subroutine StoreHSD(i,j,Hij,Sij,Di,Dj)
!Routine StoreHSdHdS stores calculated matrix elements of the
!Hamiltonian and the overlap as well as their derivatives in 
!proper places of global arrays. Upon doing this, the routine 
!normalizes matrix elements. There must be i>=j
!when routine is called. Di and Dj are the
!derivatives of Hij and Sij with respect to nonlinear parameters
!of i-th function and j-th functions:
!
!Di(1:Glob_np)              is dHijdvechLi
!Di(Glob_np+1:2*Glob_np)    is dSijdvechLi
!
!Dj(1:Glob_np)              is dHijdvechLj
!Dj(Glob_np+1:2*Glob_np)    is dSijdvechLj

integer,intent(in)     ::  i,j
real(dprec),intent(in) ::  Hij,Sij
real(dprec),intent(in) ::  Di(2*Glob_npt),Dj(2*Glob_npt)
real(dprec)            ::  f

select case (Glob_GSEPSolutionMethod) 
case('I')
  !Only the lower triangle of array Glob_H (including the diagonal) is 
  !used to store H-Glob_ApproxEnergy*S.
  !The entire array Glob_S is used to store S.
  if (i==j) then
	Glob_diagS(i)=Sij
	Glob_S(i,i)=ONE
	Glob_H(i,i)=Hij/Glob_diagS(i)-Glob_ApproxEnergy
    Glob_D(1:2*Glob_npt,i-Glob_nfru,i)=TWO*Di(1:2*Glob_npt)/Glob_diagS(i)
  else
    f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
	Glob_S(i,j)=Sij*f
	Glob_S(j,i)=Glob_S(i,j)
	Glob_H(i,j)=(Hij-Glob_ApproxEnergy*Sij)*f
    Glob_D(1:2*Glob_npt,i-Glob_nfru,j)=Di(1:2*Glob_npt)*f
	if (j>Glob_nfru) Glob_D(1:2*Glob_npt,j-Glob_nfru,i)=Dj(1:2*Glob_npt)*f
  endif
case('G')
  !In the case when Glob_GSEPSolutionMethod='G' the diagonal
  !of the Hamiltonian matrix is stored in Glob_diagH
  !The diagonal of the overlap is stored in Glob_diagS
  !Lower triangles of arrays Glob_H and 
  !Glob_S are used to store H and S
  if (i==j) then
	Glob_diagS(i)=Sij 
	Glob_diagH(i)=Hij/Glob_diagS(i)
    Glob_D(1:2*Glob_npt,i-Glob_nfru,i)=TWO*Di(1:2*Glob_npt)/Glob_diagS(i)
  else
    f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
	Glob_S(i,j)=Sij*f
	Glob_H(i,j)=Hij*f
    Glob_D(1:2*Glob_npt,i-Glob_nfru,j)=Di(1:2*Glob_npt)*f
	if (j>Glob_nfru) Glob_D(1:2*Glob_npt,j-Glob_nfru,i)=Dj(1:2*Glob_npt)*f
  endif
endselect
end subroutine StoreHSD



subroutine ComputeMatElem(Nmin,Nmax)
!Subroutine ComputeMatElem computes matrix elements of the 
!Hamiltonian and the overlap with basis functions whose number
!ranges from Nmin to Nmax. The derivatives of H and S are NOT 
!calculated at all. Routine StoreHS is called to store 
!calculated matrix elements in proper global arrays. It is 
!assumed that matrix elements of the first Nmin-1 functions are 
!already calculated and stored properly when the the routine
!is called. Thus, only those matrix elements are computed that
!are not known yet. If all matrix elements are needed then 
!one should set Nmin=1. 
!  Input parameters :
!   Nmin-1 :: The number of functions whose matrix elements
!             are already known.
!     Nmax :: The number of functions whose matrix elements
!             need to be calculated.
!Arguments :
integer     Nmin,Nmax
!Local variables :
integer     k,l,i,kk,ll,ii,j,q
integer     kstart,lstart,kstop,lstop,n,np,np1,npt,nb
real(dprec) Paramk(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
real(dprec) Paraml(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
integer     mk,ml,mmk,mml
real(dprec) Skl,Hkl,Skl1,Hkl1,Skl2,Hkl2,Skl3,Hkl3,Skl4,Hkl4
real(dprec) Vkl1,Tkl1,Vkl2,Tkl2,Vkl3,Tkl3,Vkl4,Tkl4
real(dprec) Skl5,Hkl5,Skl6,Hkl6,Skl7,Hkl7,Skl8,Hkl8
real(dprec) Vkl5,Tkl5,Vkl6,Tkl6,Vkl7,Tkl7,Vkl8,Tkl8
real(dprec) Ssum,Hsum
!These arrays are not actually used but needed for proper calling
!of subroutine MatrixElements. Thus, one can set some small size
!for them
real(dprec)  Dk(2),Dl(2),Dk1(2),Dl1(2),Dk2(2),Dl2(2),Dk3(2),Dl3(2),Dk4(2),Dl4(2)
real(dprec)  Dk5(2),Dl5(2),Dk6(2),Dl6(2),Dk7(2),Dl7(2),Dk8(2),Dl8(2)

n=Glob_n
np=Glob_np
np1=np+1
npt=Glob_npt
nb=Glob_HSBuffLen
Glob_HklBuff1(1:nb)=ZERO
Glob_SklBuff1(1:nb)=ZERO
i=0
  
do k=Nmin,Nmax
  Paramk(1:npt)=Glob_NonlinParam(1:npt,k)
  mk=Glob_ZIndex(k)
  mmk=1
  !if (mk==1) mmk=2
  do l=k,1,-1
    i=i+1
	if (i==1) then
	  kstart=k
	  lstart=l
	endif	
    Paraml(1:npt)=Glob_NonlinParam(1:npt,l)
    ml=Glob_ZIndex(l)
    mml=1
    !if (ml==1) mml=2  below is the code for S state
	Hsum=ZERO; Ssum=ZERO
	q=(i-1)*Glob_NumYHYTerms-1
	do j=1,Glob_NumYHYTerms
	  if (mod(q+j,Glob_NumOfProcs)==Glob_ProcID) then       
        if (mk==1 .or. ml==1 )    then  
                call MatrixElementsL1(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl1,Skl1,Tkl1,Vkl1,Dk1,Dl1,.false.,.false.)
                       call MatrixElementsL11(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl5,Skl5,Tkl5,Vkl5,Dk5,Dl5,.false.,.false.)
                Hkl=3*Hkl1+6*Hkl5
                Skl=3*Skl1+6*Skl5
                Dk=3*Dk1+6*Dk5
                Dl=3*Dl1+6*Dl5
               
               
        elseif (ml /=1 .or. mk/=1) then !.and. mk/=ml
        call MatrixElementsL1(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl1,Skl1,Tkl1,Vkl1,Dk1,Dl1,.false.,.false.)
!        call MatrixElementsL1(mk,mmk,ml,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!               Hkl2,Skl2,Tkl2,Vkl2,Dk2,Dl2,.false.,.false.)              
!        call MatrixElementsL1(mml,ml,mmk,mk,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!               Hkl3,Skl3,Tkl3,Vkl3,Dk3,Dl3,.false.,.false.)
!        call MatrixElementsL1(mmk,mml,mk,ml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!               Hkl4,Skl4,Tkl4,Vkl4,Dk4,Dl4,.false.,.false.)  
               
        call MatrixElementsL11(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl5,Skl5,Tkl5,Vkl5,Dk5,Dl5,.false.,.false.)
        call MatrixElementsL11(mk,mmk,ml,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl6,Skl6,Tkl6,Vkl6,Dk6,Dl6,.false.,.false.)     

!                Hkl=1.5*Hkl1+1.5*Hkl2+3*Hkl5+3*Hkl6
!                Skl=1.5*Skl1+1.5*Skl2+3*Skl5+3*Skl6
!                Dk=1.5*Dk1+1.5*Dk2+3*Dk5+3*Dk6
!                Dl=1.5*Dl1+1.5*Dl2+3*Dl5+3*Dl6
                Hkl=3*Hkl1+3*Hkl5+3*Hkl6
                Skl=3*Skl1+3*Skl5+3*Skl6
                Dk=3*Dk1+3*Dk5+3*Dk6
                Dl=3*Dl1+3*Dl5+3*Dl6
        endif       
!                Hkl=2*Hkl1+2*Hkl2+2*Hkl3-3*Hkl6-3*Hkl5
!                Skl=2*Skl1+2*Skl2+2*Skl3-3*Skl6-3*Skl5
!                Dk=2*Dk1+2*Dk2+2*Dk3-3*Dk6-3*Dk5
!                Dl=2*Dl1+2*Dl2+2*Dl3-3*Dl6-3*Dl5
            
                
                 !In this comments below I list the  approximately  working codes for S state out of two p electrons   
!                   call MatrixElementsL1(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!                        Hkl1,Skl1,Tkl1,Vkl1,Dk1,Dl1,.false.,.false.)
!                   call MatrixElementsL11(mk,mmk,ml,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!                        Hkl2,Skl2,Tkl2,Vkl2,Dk2,Dl2,.false.,.false.)
!                Hkl=6*Hkl1-6*Hkl2
!                Skl=6*Skl1-6*Skl2
!                Dk=6*Dk1-6*Dk2
!                Dl=6*Dl1-6*Dl2      
                
		Hsum=Hsum+Glob_YHYCoeff(j)*Hkl
		Ssum=Ssum+Glob_YHYCoeff(j)*Skl
	  endif
	enddo
	Glob_HklBuff1(i)=Hsum
	Glob_SklBuff1(i)=Ssum
	if (i==Glob_HSBuffLen) then
       call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
		  MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
	   call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
		  MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
	   ii=0
	   do kk=kstart,k
         if (kk==kstart) then
		   ll=lstart
         else
           ll=kk
		 endif
		 if (kk==k) then
		   lstop=l
         else
		   lstop=1
		 endif
		 do while (ll>=lstop)
		   ii=ii+1
		   call StoreHS(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii))
		   ll=ll-1
		 enddo
	   enddo
	   i=0
	   Glob_HklBuff1(1:nb)=ZERO 
	   Glob_SklBuff1(1:nb)=ZERO  
	endif
  enddo
enddo
if (i>0) then
  call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
	  MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
	  MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  ii=0
  do kk=kstart,Nmax
    if (kk==kstart) then
	  ll=lstart
    else
      ll=kk
	endif
	lstop=1
	do while (ll>=lstop)
	  ii=ii+1
	  call StoreHS(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii))
	  ll=ll-1
	enddo
  enddo    
endif

end subroutine ComputeMatElem 



subroutine ComputeMatElemAndDeriv(Nmin,Nmax)
!Subroutine ComputeMatElem computes matrix elements of the 
!Hamiltonian and the overlap as well as their derivatives with 
!basis functions whose number ranges from Nmin to Nmax. Routine 
!StoreHSD is called to store the calculated values in proper 
!global arrays. It is assumed that matrix elements of the first 
!Nmin-1 functions are already calculated and placed where 
!needed when the the routine is called. Thus, only those matrix 
!elements are computed that are not known yet. If all matrix 
!elements are needed then one should set Nmin=1. Also, one needs
!to make sure that the value of the global variables Glob_nfo and
!Glob_nfru are equal to what they should be.
!  Input parameters :
!   Nmin-1 :: The number of functions whose matrix elements
!             are already known.
!     Nmax :: The number of functions whose matrix elements
!             need to be calculated.
!Arguments :
integer     Nmin,Nmax
!Local variables :
integer     k,l,i,kk,ll,ii,j,q
integer     kstart,lstart,kstop,lstop,n,np,npt,npt2,nb
real(dprec) Paramk(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
real(dprec) Paraml(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
integer     mk,ml,mmk,mml
real(dprec) Skl,Hkl,Skl1,Hkl1,Skl2,Hkl2,Skl3,Hkl3,Skl4,Hkl4
real(dprec) Skl5,Hkl5,Skl6,Hkl6,Skl7,Hkl7,Skl8,Hkl8
real(dprec) Vkl1,Tkl1,Vkl2,Tkl2,Vkl3,Tkl3,Vkl4,Tkl4
real(dprec) Vkl5,Tkl5,Vkl6,Tkl6,Vkl7,Tkl7,Vkl8,Tkl8
real(dprec) Ssum,Hsum
real(dprec) Dk(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dl(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dksum(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dlsum(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dk1(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dl1(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dk2(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dl2(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dk3(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dl3(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dk4(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dl4(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dk5(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dl5(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dk6(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
real(dprec) Dl6(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
!real(dprec) Dk7(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
!real(dprec) Dl7(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
!real(dprec) Dk8(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
!real(dprec) Dl8(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
logical     grad_l

n=Glob_n
np=Glob_np
npt=Glob_npt
npt2=np*2
nb=Glob_HSBuffLen

Glob_HklBuff1(1:nb)=ZERO
Glob_SklBuff1(1:nb)=ZERO
Glob_DkBuff1(1:npt2,1:nb)=ZERO
if (Glob_nfo>1) Glob_DlBuff1(1:npt2,1:nb)=ZERO
i=0  

do k=Nmin,Nmax
  Paramk(1:npt)=Glob_NonlinParam(1:npt,k)
  mk=Glob_ZIndex(k)
  mmk=1
  !if (mk==1) mmk=2
  do l=k,1,-1
    i=i+1
	if (i==1) then
	  kstart=k
	  lstart=l
	endif	
    Paraml(1:npt)=Glob_NonlinParam(1:npt,l)
    ml=Glob_ZIndex(l)
    mml=1
   !if (ml==1) mml=2
	Hsum=ZERO 
	Ssum=ZERO
    Dksum(1:npt2)=ZERO
	if ((l>Glob_nfru).and.(l/=k)) then
      grad_l=.true.
      Dlsum(1:npt2)=ZERO	          
	else 
      grad_l=.false.
	endif    
	q=(i-1)*Glob_NumYHYTerms-1
	do j=1,Glob_NumYHYTerms
	  if (mod(q+j,Glob_NumOfProcs)==Glob_ProcID) then  

               
        if (mk==1 .or. ml==1 )    then  
                call MatrixElementsL1(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl1,Skl1,Tkl1,Vkl1,Dk1,Dl1,.true.,grad_l)
                       call MatrixElementsL11(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl5,Skl5,Tkl5,Vkl5,Dk5,Dl5,.true.,grad_l)
                Hkl=3*Hkl1+6*Hkl5
                Skl=3*Skl1+6*Skl5
                Dk=3*Dk1+6*Dk5
                Dl=3*Dl1+6*Dl5
               
               
        elseif (ml /=1 .or. mk/=1) then !.and. mk/=ml
        call MatrixElementsL1(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl1,Skl1,Tkl1,Vkl1,Dk1,Dl1,.true.,grad_l)
!        call MatrixElementsL1(mk,mmk,ml,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!               Hkl2,Skl2,Tkl2,Vkl2,Dk2,Dl2,.false.,.false.)              
!        call MatrixElementsL1(mml,ml,mmk,mk,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!               Hkl3,Skl3,Tkl3,Vkl3,Dk3,Dl3,.false.,.false.)
!        call MatrixElementsL1(mmk,mml,mk,ml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!               Hkl4,Skl4,Tkl4,Vkl4,Dk4,Dl4,.false.,.false.)  
               
        call MatrixElementsL11(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl5,Skl5,Tkl5,Vkl5,Dk5,Dl5,.true.,grad_l)
        call MatrixElementsL11(mk,mmk,ml,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
               Hkl6,Skl6,Tkl6,Vkl6,Dk6,Dl6,.true.,grad_l)     

!                Hkl=1.5*Hkl1+1.5*Hkl2+3*Hkl5+3*Hkl6
!                Skl=1.5*Skl1+1.5*Skl2+3*Skl5+3*Skl6
!                Dk=1.5*Dk1+1.5*Dk2+3*Dk5+3*Dk6
!                Dl=1.5*Dl1+1.5*Dl2+3*Dl5+3*Dl6
                Hkl=3*Hkl1+3*Hkl5+3*Hkl6
                Skl=3*Skl1+3*Skl5+3*Skl6
                Dk=3*Dk1+3*Dk5+3*Dk6
                Dl=3*Dl1+3*Dl5+3*Dl6
        endif       
!                Hkl=2*Hkl1+2*Hkl2+2*Hkl3-3*Hkl6-3*Hkl5
!                Skl=2*Skl1+2*Skl2+2*Skl3-3*Skl6-3*Skl5
!                Dk=2*Dk1+2*Dk2+2*Dk3-3*Dk6-3*Dk5
!                Dl=2*Dl1+2*Dl2+2*Dl3-3*Dl6-3*Dl5
            
                
                 !In this comments below I list the  approximately  working codes for S state out of two p electrons   
!                   call MatrixElementsL1(mk,ml,mmk,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!                        Hkl1,Skl1,Tkl1,Vkl1,Dk1,Dl1,.false.,.false.)
!                   call MatrixElementsL11(mk,mmk,ml,mml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
!                        Hkl2,Skl2,Tkl2,Vkl2,Dk2,Dl2,.false.,.false.)
!                Hkl=6*Hkl1-6*Hkl2
!                Skl=6*Skl1-6*Skl2
!                Dk=6*Dk1-6*Dk2
!                Dl=6*Dl1-6*Dl2      
                
                 
		Hsum=Hsum+Glob_YHYCoeff(j)*Hkl
		Ssum=Ssum+Glob_YHYCoeff(j)*Skl
        Dksum(1:npt2)=Dksum(1:npt2)+Glob_YHYCoeff(j)*Dk(1:npt2)
        if ((l>Glob_nfru).and.(l/=k)) Dlsum(1:npt2)=Dlsum(1:npt2)+Glob_YHYCoeff(j)*Dl(1:npt2)
	  endif
	enddo
	Glob_HklBuff1(i)=Hsum
	Glob_SklBuff1(i)=Ssum
	Glob_DkBuff1(1:npt2,i)=Dksum(1:npt2)
    if ((l>Glob_nfru).and.(l/=k)) Glob_DlBuff1(1:npt2,i)=Dlsum(1:npt2)
	if (i==Glob_HSBuffLen) then
       call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
                          MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
	   call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
	                      MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
       call MPI_ALLREDUCE(Glob_DkBuff1,Glob_DkBuff2,i*npt2, &
                          MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
       if (Glob_nfo>1) call MPI_ALLREDUCE(Glob_DlBuff1,Glob_DlBuff2,i*npt2, &
                                          MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
	   ii=0
	   do kk=kstart,k
         if (kk==kstart) then
		   ll=lstart
         else
           ll=kk
		 endif
		 if (kk==k) then
		   lstop=l
         else
		   lstop=1
		 endif
		 do while (ll>=lstop)
		   ii=ii+1
		   call StoreHSD(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii), &
		                 Glob_DkBuff2(1:npt2,ii),Glob_DlBuff2(1:npt2,ii))
		   ll=ll-1
		 enddo
	   enddo
	   i=0
       Glob_HklBuff1(1:nb)=ZERO
       Glob_SklBuff1(1:nb)=ZERO
       Glob_DkBuff1(1:npt2,1:nb)=ZERO
       if (Glob_nfo>1) Glob_DlBuff1(1:npt2,1:nb)=ZERO
	endif
  enddo
enddo
if (i>0) then
  call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
                     MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
                     MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  call MPI_ALLREDUCE(Glob_DkBuff1,Glob_DkBuff2,i*npt2, &
                     MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  if (Glob_nfo>1) call MPI_ALLREDUCE(Glob_DlBuff1,Glob_DlBuff2,i*npt2, &
                                     MPI_DPREC,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
  ii=0
  do kk=kstart,Nmax
    if (kk==kstart) then
	  ll=lstart
    else
      ll=kk
	endif
	lstop=1
	do while (ll>=lstop)
	  ii=ii+1
	  call StoreHSD(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii), &
	                Glob_DkBuff2(1:npt2,ii),Glob_DlBuff2(1:npt2,ii))
	  ll=ll-1
	enddo
  enddo    
endif

end subroutine ComputeMatElemAndDeriv



end module matform
