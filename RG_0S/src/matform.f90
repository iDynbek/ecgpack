module matform
!Module matform contains procedures that form Hamiltonian
!and overlap matrices and related routines
  use matelem
#ifdef USE_CUDA
  use matelem_gpu
#endif
#ifdef USE_OPENACC
  use matelem_acc
#endif
  implicit none

contains

  subroutine StoreHS(i,j,Hij,Sij)
!Routine StoreHS stores calculated matrix elements of the
!Hamiltonian and the overlap in proper places of global
!arrays. Upon doing this, the routine normalizes
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
!already calculated and stored properly when the routine
!is called. Thus, only those matrix elements are computed that
!are not known yet. If all matrix elements are needed then
!one should set Nmin=1.
!  Input parameters :
!   Nmin-1 :: The number of functions whose matrix elements
!             are already known.
!     Nmax :: The number of functions whose matrix elements
!             need to be calculated.
!Arguments :
    integer         Nmin,Nmax
!Local variables :
    integer     k,l,i,kk,ll,ii,j,q
    integer     kstart,lstart,kstop,lstop,n,np,np1,npt,nb
    real(dprec) Paramk(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
    real(dprec) Paraml(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
    real(dprec) Skl,Hkl
    real(dprec) Ssum,Hsum
!These arrays are not actually used but needed for proper calling
!of subroutine MatrixElements. Thus, one can set some small size
!for them
    real(dprec)  Dk(2),Dl(2)
#ifdef USE_CUDA
!GPU batch temporaries (allocated only when Glob_UseGPU)
    integer              :: npairs_gpu, ipair_gpu
    integer, allocatable :: k_list_gpu(:), l_list_gpu(:)
    real(dprec), allocatable :: Hout_gpu(:), Sout_gpu(:)
    real(dprec), allocatable :: Hout_gpu_r(:), Sout_gpu_r(:)
#endif
#ifdef USE_OPENACC
!OpenACC batch temporaries
    integer              :: npairs_acc, ip_acc, qa, ja, ka, la
    integer, allocatable :: k_list_acc(:), l_list_acc(:)
    real(dprec), allocatable :: Hout_acc(:), Sout_acc(:)
    real(dprec), allocatable :: Hout_acc_r(:), Sout_acc_r(:)
    real(dprec) :: Hkl_a, Skl_a, Hs_a, Ss_a
    real(dprec) :: massM(Glob_n,Glob_n), chg(Glob_n)
    real(dprec) :: chg0, pir, ats, rps, rpp, rpm
#endif

    n=Glob_n
    np=Glob_np
    np1=np+1
    npt=Glob_npt
    nb=Glob_HSBuffLen

#ifdef USE_CUDA
    if (Glob_UseGPU) then
      !One GPU call covers the entire [Nmin,Nmax] range.
      npairs_gpu = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
      allocate(k_list_gpu(npairs_gpu), l_list_gpu(npairs_gpu))
      allocate(Hout_gpu(npairs_gpu), Sout_gpu(npairs_gpu))
      allocate(Hout_gpu_r(npairs_gpu), Sout_gpu_r(npairs_gpu))
      ipair_gpu = 0
      do k=Nmin,Nmax
        do l=k,1,-1
          ipair_gpu = ipair_gpu + 1
          k_list_gpu(ipair_gpu) = k
          l_list_gpu(ipair_gpu) = l
        enddo
      enddo
      call gpu_compute_matelem_batch( &
          Glob_NonlinParam(1:np,1:Nmax), Nmax, &
          k_list_gpu, l_list_gpu, npairs_gpu, &
          Glob_YHYMatr(1:n,1:n,1), Glob_YHYCoeff, Glob_NumYHYTerms, &
          n, np, Glob_NumOfProcs, Glob_ProcID, &
          Glob_Piraised3n2, Glob_MassMatrix(1:n,1:n), &
          Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
          Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
          Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
          Hout_gpu, Sout_gpu)
      call MPI_ALLREDUCE(Hout_gpu, Hout_gpu_r, npairs_gpu, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      call MPI_ALLREDUCE(Sout_gpu, Sout_gpu_r, npairs_gpu, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      ipair_gpu = 0
      do k=Nmin,Nmax
        do l=k,1,-1
          ipair_gpu = ipair_gpu + 1
          call StoreHS(k, l, Hout_gpu_r(ipair_gpu), Sout_gpu_r(ipair_gpu))
        enddo
      enddo
      deallocate(k_list_gpu, l_list_gpu, Hout_gpu, Sout_gpu, Hout_gpu_r, Sout_gpu_r)
    else
#endif
#ifdef USE_OPENACC
    if (Glob_UseGPU) then
      !OpenACC single-source GPU path: the same Fortran MatrixElements_energy_acc
      !runs on the device, one gang/vector per (k,l) pair.
      npairs_acc = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
      allocate(k_list_acc(npairs_acc), l_list_acc(npairs_acc))
      allocate(Hout_acc(npairs_acc), Sout_acc(npairs_acc))
      allocate(Hout_acc_r(npairs_acc), Sout_acc_r(npairs_acc))
      ip_acc = 0
      do ka=Nmin,Nmax
        do la=ka,1,-1
          ip_acc = ip_acc + 1
          k_list_acc(ip_acc) = ka
          l_list_acc(ip_acc) = la
        enddo
      enddo
      massM = Glob_MassMatrix(1:n,1:n)
      chg   = Glob_PseudoCharge(1:n)
      chg0  = Glob_PseudoCharge0
      pir   = Glob_Piraised3n2
      ats   = Glob_AttractionScalingParam
      rps   = Glob_RepulsionScalingParam
      rpp   = Glob_RepulsionScalingParamPlus
      rpm   = Glob_RepulsionScalingParamMinus
      !One-time: park the constant permutation data on the device so it is not
      !re-copied on every call (it never changes during a run).
      if (.not. Glob_AccConstLoaded) then
        !$acc enter data copyin(Glob_YHYMatr, Glob_YHYCoeff)
        Glob_AccConstLoaded = .true.
      endif
      !Parallelize BOTH axes: gang over (k,l) pairs, vector over permutation terms
      !(reduction), mirroring the hand-CUDA block-per-pair / thread-per-term layout.
      !$acc parallel loop gang &
      !$acc   copyin(Glob_NonlinParam(1:np,1:Nmax), k_list_acc, l_list_acc, massM, chg) &
      !$acc   copyout(Hout_acc, Sout_acc) &
      !$acc   present(Glob_YHYMatr, Glob_YHYCoeff) &
      !$acc   private(ka, la, Hs_a, Ss_a, qa)
      do ip_acc = 1, npairs_acc
        ka = k_list_acc(ip_acc)
        la = l_list_acc(ip_acc)
        Hs_a = ZERO; Ss_a = ZERO
        qa = (ip_acc-1)*Glob_NumYHYTerms - 1
        !$acc loop vector reduction(+:Hs_a,Ss_a) private(Hkl_a, Skl_a)
        do ja = 1, Glob_NumYHYTerms
          if (mod(qa+ja, Glob_NumOfProcs) == Glob_ProcID) then
            call MatrixElements_energy_acc( &
                Glob_NonlinParam(1:np,ka), Glob_NonlinParam(1:np,la), &
                Glob_YHYMatr(1:n,1:n,ja), n, np, &
                pir, massM, chg, chg0, ats, rps, rpp, rpm, Hkl_a, Skl_a)
            Hs_a = Hs_a + Glob_YHYCoeff(ja)*Hkl_a
            Ss_a = Ss_a + Glob_YHYCoeff(ja)*Skl_a
          endif
        enddo
        Hout_acc(ip_acc) = Hs_a
        Sout_acc(ip_acc) = Ss_a
      enddo
      !$acc end parallel loop
      call MPI_ALLREDUCE(Hout_acc, Hout_acc_r, npairs_acc, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      call MPI_ALLREDUCE(Sout_acc, Sout_acc_r, npairs_acc, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      ip_acc = 0
      do ka=Nmin,Nmax
        do la=ka,1,-1
          ip_acc = ip_acc + 1
          call StoreHS(ka, la, Hout_acc_r(ip_acc), Sout_acc_r(ip_acc))
        enddo
      enddo
      deallocate(k_list_acc, l_list_acc, Hout_acc, Sout_acc, Hout_acc_r, Sout_acc_r)
    else
#endif

    Glob_HklBuff1(1:nb)=ZERO
    Glob_SklBuff1(1:nb)=ZERO
    i=0

    do k=Nmin,Nmax
      Paramk(1:npt)=Glob_NonlinParam(1:npt,k)
      do l=k,1,-1
        i=i+1
        if (i==1) then
          kstart=k
          lstart=l
        endif
        Paraml(1:npt)=Glob_NonlinParam(1:npt,l)
        Hsum=ZERO; Ssum=ZERO
        q=(i-1)*Glob_NumYHYTerms-1
        do j=1,Glob_NumYHYTerms
          if (mod(q+j,Glob_NumOfProcs)==Glob_ProcID) then
            call MatrixElements(Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j),Hkl,Skl,Dk,Dl,.false.,.false.)
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

#ifdef USE_OPENACC
    endif  !Glob_UseGPU
#endif
#ifdef USE_CUDA
    endif  !Glob_UseGPU
#endif

  end subroutine ComputeMatElem

  subroutine ComputeMatElemAndDeriv(Nmin,Nmax)
!Subroutine ComputeMatElem computes matrix elements of the
!Hamiltonian and the overlap as well as their derivatives with
!basis functions whose number ranges from Nmin to Nmax. Routine
!StoreHSD is called to store the calculated values in proper
!global arrays. It is assumed that matrix elements of the first
!Nmin-1 functions are already calculated and placed where
!needed when the routine is called. Thus, only those matrix
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
    integer         Nmin,Nmax
!Local variables :
    integer     k,l,i,kk,ll,ii,j,q
    integer     kstart,lstart,kstop,lstop,n,np,npt,npt2,nb
    real(dprec) Paramk(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
    real(dprec) Paraml(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1)/2)
    real(dprec) Skl,Hkl
    real(dprec) Ssum,Hsum
    real(dprec) Dk(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
    real(dprec) Dl(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
    real(dprec) Dksum(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
    real(dprec) Dlsum(Glob_MaxAllowedNumOfPseudoParticles*(Glob_MaxAllowedNumOfPseudoParticles+1))
    logical     grad_l
#ifdef USE_CUDA
    integer              :: npairs_gpu, ipair_gpu
    integer, allocatable :: k_list_gpu(:), l_list_gpu(:), grad_l_gpu(:)
    real(dprec), allocatable :: Hout_gpu(:), Sout_gpu(:)
    real(dprec), allocatable :: Hout_gpu_r(:), Sout_gpu_r(:)
    real(dprec), allocatable :: Dkout_gpu(:,:), Dlout_gpu(:,:)
    real(dprec), allocatable :: Dkout_gpu_r(:,:), Dlout_gpu_r(:,:)
#endif

    n=Glob_n
    np=Glob_np
    npt=Glob_npt
    npt2=np*2
    nb=Glob_HSBuffLen

#ifdef USE_CUDA
    if (Glob_UseGPU) then
      npairs_gpu = Nmax*(Nmax+1)/2 - (Nmin-1)*Nmin/2
      allocate(k_list_gpu(npairs_gpu), l_list_gpu(npairs_gpu), grad_l_gpu(npairs_gpu))
      allocate(Hout_gpu(npairs_gpu), Sout_gpu(npairs_gpu))
      allocate(Hout_gpu_r(npairs_gpu), Sout_gpu_r(npairs_gpu))
      allocate(Dkout_gpu(npt2,npairs_gpu), Dlout_gpu(npt2,npairs_gpu))
      allocate(Dkout_gpu_r(npt2,npairs_gpu), Dlout_gpu_r(npt2,npairs_gpu))
      ipair_gpu = 0
      do k=Nmin,Nmax
        do l=k,1,-1
          ipair_gpu = ipair_gpu + 1
          k_list_gpu(ipair_gpu) = k
          l_list_gpu(ipair_gpu) = l
          if ((l>Glob_nfru).and.(l/=k)) then
            grad_l_gpu(ipair_gpu) = 1
          else
            grad_l_gpu(ipair_gpu) = 0
          endif
        enddo
      enddo
      Dkout_gpu = ZERO; Dlout_gpu = ZERO
      call gpu_compute_matelem_and_deriv_batch( &
          Glob_NonlinParam(1:np,1:Nmax), Nmax, &
          k_list_gpu, l_list_gpu, grad_l_gpu, npairs_gpu, &
          Glob_YHYMatr(1:n,1:n,1), Glob_YHYCoeff, Glob_NumYHYTerms, &
          n, np, Glob_NumOfProcs, Glob_ProcID, &
          Glob_Piraised3n2, Glob_MassMatrix(1:n,1:n), &
          Glob_PseudoCharge(1:n), Glob_PseudoCharge0, &
          Glob_AttractionScalingParam, Glob_RepulsionScalingParam, &
          Glob_RepulsionScalingParamPlus, Glob_RepulsionScalingParamMinus, &
          Hout_gpu, Sout_gpu, Dkout_gpu, Dlout_gpu)
      call MPI_ALLREDUCE(Hout_gpu,  Hout_gpu_r,  npairs_gpu, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      call MPI_ALLREDUCE(Sout_gpu,  Sout_gpu_r,  npairs_gpu, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      call MPI_ALLREDUCE(Dkout_gpu, Dkout_gpu_r, npairs_gpu*npt2, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      if (Glob_nfo>1) &
      call MPI_ALLREDUCE(Dlout_gpu, Dlout_gpu_r, npairs_gpu*npt2, &
          MPI_DPREC, MPI_SUM, MPI_COMM_WORLD, Glob_MPIErrCode)
      ipair_gpu = 0
      do k=Nmin,Nmax
        do l=k,1,-1
          ipair_gpu = ipair_gpu + 1
          call StoreHSD(k, l, Hout_gpu_r(ipair_gpu), Sout_gpu_r(ipair_gpu), &
              Dkout_gpu_r(1:npt2,ipair_gpu), Dlout_gpu_r(1:npt2,ipair_gpu))
        enddo
      enddo
      deallocate(k_list_gpu, l_list_gpu, grad_l_gpu)
      deallocate(Hout_gpu, Sout_gpu, Hout_gpu_r, Sout_gpu_r)
      deallocate(Dkout_gpu, Dlout_gpu, Dkout_gpu_r, Dlout_gpu_r)
    else
#endif

    Glob_HklBuff1(1:nb)=ZERO
    Glob_SklBuff1(1:nb)=ZERO
    Glob_DkBuff1(1:npt2,1:nb)=ZERO
    if (Glob_nfo>1) Glob_DlBuff1(1:npt2,1:nb)=ZERO
    i=0

    do k=Nmin,Nmax
      Paramk(1:npt)=Glob_NonlinParam(1:npt,k)
      do l=k,1,-1
        i=i+1
        if (i==1) then
          kstart=k
          lstart=l
        endif
        Paraml(1:npt)=Glob_NonlinParam(1:npt,l)
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
            call MatrixElements(Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j),Hkl,Skl,Dk,Dl,.true.,grad_l)
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

#ifdef USE_CUDA
    endif  !Glob_UseGPU
#endif

  end subroutine ComputeMatElemAndDeriv

end module matform