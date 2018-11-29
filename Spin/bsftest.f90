program BasicSpinFuncMain

! Basic Spin Function

use BSFPrecision
use BSFGlobalVariables
use BSFFunctions

implicit none

integer                                  :: Sn2,Mn2
integer                                  :: n,nunp,nupunp,nconfunp,np,nupp,nconfp,&
                                            nup,nconf
integer, dimension(1:8)                  :: narr
integer, allocatable, dimension(:,:)     :: conf
integer, allocatable, dimension(:)       :: signConf
real(kind=dp), allocatable, dimension(:) :: multConf

integer                                  :: i,j,m,res,maxSn2,startSn2


! *** Binomial Test ***

write(*,*) "Binomial Test"
n=51
m=25
res=BSFBinomial(n,m)
if(BSFBinomINFO==-1) then
   write(*,*) "Incorrect input data in the function BSFBinomial. Binomial(",n,",",m,")"
else if(BSFBinomINFO==1) then
   write(*,*) "Integer too big for its kind!"
else
   write(*,*) "BSFBinomINFO=",BSFBinomINFO
   write(*,*) "BSFBinomial(Big)=",res
end if

! *** Unpaired Configurations ***
! Condition: n > 0

write(*,*) "Unpaired Configurations"
nunp=5
nupunp=3
nconfunp=BSFBinomial(nunp,nupunp)
allocate(conf(1:nconfunp,1:MAX(nupunp,1)))
call BSFSubUnpairedConfs(nunp,nupunp,nconfunp,conf)
write(*,*) "BSFBinomINFO=",BSFBinomINFO
do i=1,nconfunp
write(*,*) conf(i,:)
end do
deallocate(conf)

! *** Paired Configurations ***
! Condition: n > 0

write(*,*) "Paired Configurations"
nupp=4
nconfp=2**nupp
allocate(conf(1:nconfp,1:nupp),signConf(1:nconfp))
call BSFSubPairedConfs(nupp,nconfp,conf,signConf)
do i=1,nconfp
write(*,*) signConf(i),conf(i,:)
end do
deallocate(conf,signConf)

! *** Basic Spin Configurations ***

write(*,*) "Basic Spin Configurations"

n=3
write(*,*) "Basic Spin Configurations for n = ",n
if(MOD(n,2)==0) then
   maxSn2=n/2
   startSn2=0
else
   maxSn2=(n-1)/2
   startSn2=1
end if
do j=0,maxSn2
   Sn2=startSn2+2*j
   do Mn2=-Sn2,Sn2,2
      write(*,*) " "
      write(*,*) "n    = ",n
      write(*,*) "2*Sn = ",Sn2
      write(*,*) "2*Mn = ",Mn2
      write(*,*) " "
      call BSFSubNumbers(n,Sn2,Mn2,narr)
      allocate(conf(1:narr(8),1:MAX(narr(7),1)),multConf(1:narr(8)))
      call BSFSubBasicSpinConfs(n,narr,conf,multConf)
      do i=1,narr(8)
      write(*,*) multConf(i),conf(i,:)
      end do
      deallocate(conf,multConf)
   end do
end do

end program
