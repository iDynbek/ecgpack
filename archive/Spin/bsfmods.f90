module BSFPrecision
implicit none
! single precision (4 Byte):  selected_real_kind(6, 37)
! double precision (8 Byte):  selected_real_kind(15, 307)
! quad precision (16 Byte): selected_real_kind(33, 4931)
integer, parameter :: dp = selected_real_kind(15, 307)
end module

module BSFGlobalVariables
use BSFPrecision
implicit none

integer :: BSFBinomINFO=0, BSFBinomINFO8=0
end module

module BSFFunctions

use BSFPrecision

implicit none

contains

function BSFBinomial(n,m)

use BSFGlobalVariables, only: BSFBinomINFO

! BSFBinomINFO = -1, if there is incorrect input data (=> BSFBinomial=-1)
! BSFBinomINFO =  1, if the result bigger than 2147483647 (normal integer*4) (=> BSFBinomial=-1)
! BSFBinomINFO =  0, otherwise

implicit none

    integer, intent(in) :: n,m
    integer             :: BSFBinomial
    
    integer  :: i,tmpi
    real(dp) :: tmpr
    real(dp) :: rBSFBinomial
    
    BSFBinomINFO=0
    if((n>m-1).AND.(m>-1)) then
       BSFBinomial=1
       tmpi=Min(m,n-m)
       if((n>1).AND.(tmpi>0)) then
          rBSFBinomial=1._dp
          do i=0,tmpi-1
             tmpr=REAL(n-i,kind=dp)
	     tmpr=tmpr/REAL(tmpi-i,kind=dp)
             rBSFBinomial=rBSFBinomial*tmpr
          end do
	  if(rBSFBinomial>2147483647._dp) then
	     BSFBinomINFO=1
	     BSFBinomial=-1
	  else
	     BSFBinomial=NINT(rBSFBinomial)
	  end if
       end if
    else
       BSFBinomINFO=-1
       BSFBinomial=-1
    end if
    
end function BSFBinomial

subroutine BSFSubUnpairedConfs(n,nup,nconf,conf)
!
! Construct of possible configurations of "up"s, "nup", and "down"s, "n-nup",
! which are equivalent with UNPAIRED configurations.
!
! INPUT
!
! n     :  Sum of "up"s and "down"s. (n > 0)
! nup   :  Number of "up"s. (nup >= 0 and nup <= n)
! nconf :  Number of possible configurations. (nconf >= 1)
!
! OUTPUT
!
! conf  :  Array of configurations with dimensions nconf * MAX(nup,1).
!          If nup=0 then the dimensions are 1*1 and conf(1,1)=0.
!          If nup>0 then the dimensions are nconf*nup and
!          stores the position(s) of "up"s.
!
!          Example: n = 4, nup = 3, nconf = Binomial(4,3) = 4
!
!                   conf(1,1-3) = 1, 2, 3
!                   conf(2,1-3) = 1, 2, 4
!                   conf(3,1-3) = 1, 3, 4
!                   conf(4,1-3) = 2, 3, 4

implicit none

    integer, intent(in)                                   :: n,nup,nconf
    integer, dimension(1:nconf,1:MAX(nup,1)), intent(out) :: conf
    
    logical                          :: logicnum
    integer                          :: num,p,i
    integer, dimension(1:MAX(nup,1)) :: ctmp
    
    if(nup==0) then
       ctmp(1)=0
       conf(1,1)=ctmp(1)
    else
       num=nup+1
       do i=1,nup
          ctmp(i)=num-i
       end do
       conf(1,:)=ctmp
    end if
    
    num=1
    p=1
    if((nup.NE.n).AND.(nup.NE.0)) then
       logicnum=.True.
       do while(logicnum)
          if(ctmp(p)==n-p+1) then
	     if(p<nup) then
	        p=p+1
		do while((ctmp(p)==n-p+1).AND.(p<nup))
		   p=p+1
		end do
	     end if
	     if((p==nup).AND.(ctmp(p)==n-p+1)) then
	        logicnum=.False.
	     else
	        num=num+1
                ctmp(p)=ctmp(p)+1
		do i=1,p-1
		   ctmp(p-i)=ctmp(p)+i
		end do
                conf(num,:)=ctmp
                p=1
	     end if
	  else
	     num=num+1
             ctmp(p)=ctmp(p)+1
             conf(num,:)=ctmp
	  end if
       end do
    end if
    conf(:,1:MAX(nup,1))=conf(:,MAX(nup,1):1:-1)

end subroutine BSFSubUnpairedConfs

subroutine BSFSubPairedConfs(nup,nconf,conf,signConf)
!
! Construct of possible configurations of "up"s, "nup", and "down"s, "= nup",
! which are equivalent with PAIRED configurations and gives the signs of them.
! The possible paired configurations are given by the expression
! Product_i=1^nup [up(2*i-1)down(2*i)-down(2*i-1)up(2*i)].
!
! INPUT
!
! nup      : Number of "up"s. (nup = the half of number of paired particles and
!            in this subroutine nup >= 1 has to be satisfied.)
! nconf    : Number of possible configurations. (nconf >= 1)
!
! OUTPUT
!
! conf     : Array of configurations with dimensions nconf * nup.
!            It stores the position(s) of "up"s.
! signConf : Signs in the above expression.
!
!          Example: nup = 2 (=> 4 paired particle), nconf = 2**nup = 4
!
!                   conf(1,1-2) = 1, 3      signConf(1) =  1
!                   conf(2,1-2) = 1, 4      signConf(2) = -1
!                   conf(3,1-2) = 2, 3      signConf(3) = -1
!                   conf(4,1-2) = 2, 4      signConf(4) =  1

implicit none

    integer, intent(in)                            :: nup,nconf
    integer, dimension(1:nconf,1:nup), intent(out) :: conf
    integer, dimension(1:nconf), intent(out)       :: signConf
    
    integer :: i,j,tmp1,tmp2
    
    signConf=0
    conf=0
    do j=1,nconf
       tmp1=j
       tmp2=nconf/2
       do i=1,nup
          if(tmp1>tmp2) then
	     conf(j,i)=2*i
	     tmp1=tmp1-tmp2
	     signConf(j)=signConf(j)+1
	  else
             conf(j,i)=2*i-1
          end if
          tmp2=tmp2/2          
       end do
       signConf(j)=(-1)**MOD(signConf(j),2)
    end do
    
end subroutine BSFSubPairedConfs

subroutine BSFSubBasicSpinConfs(n,narr,conf,multConf)
!
! Configurations for the Basic Spin Function
!
! Construct the possible configurations of "up"s, "nup", and "down"s, "n-nup", for the
! full system of "n" particles [conf], and the multipliers of them [multConf].
!
! CONDITIONS
!
! We suppose that n > 0, np = n-nunp, and MOD(np,2) = 0 are satisfied.
! If nunp = 0, then nconfunp = 1 also.
! If np = 0, then nconfp = 1 also.
! 
! INPUT
! 
! n        : Number of particles. (Sum of "up"s and "down"s. n > 0)
!
! narr(1:8) = nunp,nupunp,nconfunp,np,nupp,nconfp,nup,nconf
!
! nunp     : Number of unpaired particles. (nunp >= 0)
! nupunp   : Number of unpaired "up"s. (nupunp >= 0 and nupunp <= n)
! nconfunp : Number of possible configurations of unpaired particles. (nconfunp >= 1)
!            If nunp = 0, then nconfunp = 1 also.
! np       : Number of paired particles. (np >= 0)
! nupp     : Number of paired "up"s. (nupp >= 0 and nupp <= n)
! nconfp   : Number of possible configurations of paired particles. (nconfp >= 1)
!            If np = 0, then nconfp = 1 also.
! nup      : Number of "up"s in the n-particle system. (nup >= 0 and nup <= n)
!            [nup = nupp + nupunp]
! nconf    : Number of possible configurations for the n-particle system. (nconf >= 1)
!            [nconf = nconfp * nconfunp]
! 
! OUTPUT
! 
! conf     : If the "up" particle, numbered by "nparticle", nparticle = 0,1,...,n, then conf(i,j) = nparticle,
!            where i = 1,...,nconf, and j = 1,...,nup. Every configuration, i, has nup "up"s, conf(i,j)
!            shows the "nparticle"s in increasing order as j increase. In the case nup = 0, the dimensions of
!            conf are 1*1 and conf(1,1)=0, otherwise the dimensions are nconf * nup.
! multConf : Array of multipliers of the possible configurations in the expression of the basic spin
!            function:
!                     Sum_i=1^nconf multConf(i) * {conf(i,:) -> Product_j=1^n {alpha OR beta}(j)}
! 

implicit none

    integer, intent(in)                                         :: n
    integer, dimension(1:8), intent(in)                         :: narr
    real(kind=dp), dimension(1:narr(8)), intent(out)            :: multConf
    integer, dimension(1:narr(8),1:MAX(narr(7),1)), intent(out) :: conf
    
    integer                              :: i,j,k,tmp
    integer, allocatable, dimension(:,:) :: confunp,confp
    integer, allocatable, dimension(:)   :: signConfp
    integer                              :: nunp,nupunp,nconfunp,np,nupp,nconfp,nup,nconf
    
    nunp=narr(1)
    nupunp=narr(2)
    nconfunp=narr(3)
    np=narr(4)
    nupp=narr(5)
    nconfp=narr(6)
    nup=narr(7)
    nconf=narr(8)
    
    multConf=0._dp
    conf=0
    if(np==0) then
       allocate(confunp(1:nconfunp,1:MAX(nupunp,1)))
       call BSFSubUnpairedConfs(nunp,nupunp,nconfunp,confunp)
       conf=confunp
       deallocate(confunp)
       multConf=nconfunp**(-0.5_dp)
    else
       allocate(confp(1:nconfp,1:nupp),signConfp(1:nconfp))
       if(nunp==0) then
          call BSFSubPairedConfs(nupp,nconfp,confp,signConfp)
          conf=confp
          multConf=2._dp**(-0.5_dp*nupp)
	  do i=1,nconfp
	     multConf(i)=multConf(i)*signConfp(i)
	  end do
       else
          allocate(confunp(1:nconfunp,1:MAX(nupunp,1)))
          call BSFSubUnpairedConfs(nunp,nupunp,nconfunp,confunp)
          call BSFSubPairedConfs(nupp,nconfp,confp,signConfp)
          multConf=2._dp**(-0.5_dp*nupp)*nconfunp**(-0.5_dp)
	  do j=1,nconfunp
	     do i=1,nconfp
	        tmp=(j-1)*nconfp+i
	        multConf(tmp)=multConf(tmp)*signConfp(i)
	        do k=1,nupp
		   conf(tmp,k)=confp(i,k)
		end do
		do k=nupp+1,nup
		   conf(tmp,k)=confunp(j,k-nupp)+2*nupp
		end do
	     end do
	  end do
	  deallocate(confunp)
       end if
       deallocate(confp,signConfp)
    end if
    
end subroutine BSFSubBasicSpinConfs

subroutine BSFSubNumbers(n,Sn2,Mn2,narr)
!
! Calculation of "n"s.
!
! INPUT
! 
! n        : Number of particles. (Sum of "up"s and "down"s. n >= 0)
! Sn2      : Two times total S (Spin)
! Mn2      : Two times total M (S_z)
! 
! OUTPUT
! 
! narr(1) = nunp     : Number of unpaired particles. (nunp >= 0)
! narr(2) = nupunp   : Number of unpaired "up"s. (nupunp >= 0 and nupunp <= n)
! narr(3) = nconfunp : Number of possible configurations of unpaired particles. (nconfunp >= 1)
!                      If nunp = 0, then nconfunp = 1 also.
! narr(4) = np       : Number of paired particles. (np >= 0)
! narr(5) = nupp     : Number of paired "up"s. (nupp >= 0 and nupp <= n)
! narr(6) = nconfp   : Number of possible configurations of paired particles. (nconfp >= 1)
!                      If np = 0, then nconfp = 1 also.
! narr(7) = nup      : Number of "up"s in the n-particle system. (nup >= 0 and nup <= n)
!                      [nup = nupp + nupunp]
! narr(8) = nconf    : Number of possible configurations for the n-particle system. (nconf >= 1)
!                      [nconf = nconfp * nconfunp]

use BSFGlobalVariables, only: BSFBinomINFO

implicit none

    integer, intent(in)                  :: n,Sn2,Mn2
    integer, dimension(1:8), intent(out) :: narr
    
    if(n<1) then
       write(*,*) "BSFSubNumbers: n < 1"
       stop
    end if
    if(Sn2<0) then
       write(*,*) "BSFSubNumbers: Sn < 0"
       stop
    end if
    if(ABS(Mn2)>Sn2) then
       write(*,*) "BSFSubNumbers: ABS(Mn) > Sn"
       stop
    end if
    narr(1)=Sn2
    narr(2)=(Sn2+Mn2)/2
    narr(3)=BSFBinomial(narr(1),narr(2))
    select case(BSFBinomINFO)
    case(-1)
       write(*,*) "BSFSubNumbers: Incorrect input data in the function BSFBinomial. Binomial(",narr(1),",",narr(2),")"
       stop
    case(1)
       write(*,*) "BSFSubNumbers: Integer, the result of BSFBinomial, too big for its kind!"
       stop
    case(0)
       narr(4)=n-narr(1)
       if(MOD(narr(4),2).NE.0) then
          write(*,*) "BSFSubNumbers: nunp=",narr(1),", np should be an even number! (n=nunp+np)"
          stop
       end if
       narr(5)=narr(4)/2
       narr(6)=2**narr(5)
       narr(7)=narr(5)+narr(2)
       narr(8)=narr(3)*narr(6)
    end select
    
end subroutine BSFSubNumbers

end module
