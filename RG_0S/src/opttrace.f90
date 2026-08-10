module opttrace
!Optimizer trajectory tracer -- OFF unless ECG_OPTTRACE=1 is set in the
!environment, in which case rank 0 writes one CSV row per objective evaluation
!performed by the nonlinear optimizer (DRMNG). The point is visualization: the
!resulting file is a frame-by-frame record of a basis function being optimized,
!ready to be animated.
!
!This is an INSTRUMENTATION module. It is deliberately not on the production
!branches; keep it on the visualization branch.
!
!Environment:
!  ECG_OPTTRACE=1        enable tracing
!  ECG_OPTTRACE_FILE     output path (default: opttrace.csv)
!
!Output: two files.
!  <file>       CSV, one row per (objective evaluation x optimized function)
!  <file>.meta  key=value metadata needed to interpret the CSV (see below)
!
!CSV columns:
!  step    global counter, increments once per objective evaluation
!  phase   which optimizer routine produced the row (e.g. OPTCYCLE_I)
!  cbs     current basis size when the row was written
!  slot    index of this function within the set being optimized (1..nfo)
!  func    the basis function's own number in the basis
!  kind    1 = the optimizer asked for an energy, 2 = it asked for a gradient
!  status  0 = the evaluation succeeded, nonzero = it failed (energy unusable)
!  energy  the variational energy at this point (only meaningful if status=0)
!  p1..pN  the function's N=npt nonlinear parameters at this point
!
!Reconstructing the Gaussian from p1..pN (documented in the .meta file too):
!the parameters are the lower triangle of L, stored column by column, i.e.
!  indx=0; do i=1,n; do j=i,n; indx=indx+1; L(j,i)=p(indx); enddo; enddo
!and the Gaussian's exponent matrix is A = L*L' (symmetric, positive definite),
!so the basis function is exp(-r'*A*r) in pseudoparticle coordinates r.
  use globvars
  implicit none
  private

  public :: opttrace_init      !open the file and write headers (collective-safe)
  public :: opttrace_finalize  !close the file
  public :: opttrace_active    !.true. if tracing is on (guards the call sites)
  public :: opttrace_record    !write the rows for one objective evaluation

  logical, save            :: on    = .false.
  integer, save            :: iu    = 0
  integer(kind=8), save    :: step  = 0
  character(256), save     :: fname = 'opttrace.csv'

contains

  logical function opttrace_active()
    opttrace_active = on
  end function opttrace_active

  subroutine opttrace_init()
  !Read the env knobs and, if enabled, open the CSV on rank 0 and write both
  !the column header and the companion .meta file. Only rank 0 ever writes, so
  !no broadcast is needed -- the trace is pure output and never affects results.
    character(256) :: val
    integer        :: ios, i, k, indx, mu
    if (Glob_ProcID/=0) return
    call get_environment_variable('ECG_OPTTRACE', val, status=ios)
    on = (ios==0) .and. (val(1:1)=='1')
    if (.not.on) return
    call get_environment_variable('ECG_OPTTRACE_FILE', val, status=ios)
    if ((ios==0).and.(len_trim(val)>0)) fname = trim(val)

    open(newunit=iu, file=trim(fname), status='replace', action='write', iostat=ios)
    if (ios/=0) then
      write(*,'(1x,a,a)') 'opttrace: cannot open ',trim(fname)
      on = .false.
      return
    endif
    write(iu,'(a)',advance='no') 'step,phase,cbs,slot,func,kind,status,energy'
    do i=1,Glob_npt
      write(iu,'(a,i0)',advance='no') ',p',i
    enddo
    write(iu,'(a)') ''

    !Companion metadata: everything needed to turn p1..pN back into a Gaussian.
    open(newunit=mu, file=trim(fname)//'.meta', status='replace', action='write', iostat=ios)
    if (ios==0) then
      write(mu,'(a,i0)')  'nparticles=', Glob_n+1
      write(mu,'(a,i0)')  'npseudo=',    Glob_n
      write(mu,'(a,i0)')  'nparams=',    Glob_npt
      write(mu,'(a,a)')   'basis_type=', trim(Glob_BasisType)
      write(mu,'(a)')     'param_packing=lower triangle of L, column by column:'
      write(mu,'(a)')     '  indx=0; do i=1,n; do j=i,n; indx=indx+1; L(j,i)=p(indx); enddo; enddo'
      write(mu,'(a)')     'exponent_matrix=A = L*transpose(L)  (basis function is exp(-r^T A r))'
      write(mu,'(a)')     'coordinates=pseudoparticle (relative) coordinates r_i = x_(i+1) - x_1'
      if (allocated(Glob_Mass)) then
        write(mu,'(a)',advance='no') 'masses='
        do i=1,Glob_n+1
          write(mu,'(1x,es23.16)',advance='no') Glob_Mass(i)
        enddo
        write(mu,'(a)') ''
      endif
      if (allocated(Glob_PseudoCharge)) then
        write(mu,'(a,es23.16)') 'pseudocharge0=', Glob_PseudoCharge0
        write(mu,'(a)',advance='no') 'pseudocharges='
        do i=1,Glob_n
          write(mu,'(1x,es23.16)',advance='no') Glob_PseudoCharge(i)
        enddo
        write(mu,'(a)') ''
      endif
      !index map: which (row,col) of L each parameter fills
      indx=0
      do i=1,Glob_n
        do k=i,Glob_n
          indx=indx+1
          write(mu,'(a,i0,a,i0,a,i0,a)') 'p',indx,'=L(',k,',',i,')'
        enddo
      enddo
      close(mu)
    endif
    write(*,'(1x,a)')   '--------------- optimizer trace --------------'
    write(*,'(1x,a,a)') '  writing per-evaluation trajectory to ',trim(fname)
    write(*,'(1x,a)')   '----------------------------------------------'
  end subroutine opttrace_init

  subroutine opttrace_record(phase, cbs, nfo, nfru, npt, kind, status, energy, x)
  !Write one row per optimized function for a single objective evaluation.
  !x is the optimizer's packed parameter vector (nfo blocks of npt).
    character(*), intent(in) :: phase
    integer,      intent(in) :: cbs, nfo, nfru, npt, kind, status
    real(wp),     intent(in) :: energy
    real(wp),     intent(in) :: x(nfo*npt)
    integer :: i, j, fnum
    if (.not.on) return
    if (Glob_ProcID/=0) return
    step = step + 1
    do i=1,nfo
      fnum = 0
      if (allocated(Glob_FuncNum)) fnum = Glob_FuncNum(nfru+i)
      write(iu,'(i0,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,es23.16)',advance='no') &
        step,',',trim(phase),',',cbs,',',i,',',fnum,',',kind,',',status,',',energy
      do j=1,npt
        write(iu,'(a,es23.16)',advance='no') ',',x((i-1)*npt+j)
      enddo
      write(iu,'(a)') ''
    enddo
  end subroutine opttrace_record

  subroutine opttrace_finalize()
    if (.not.on) return
    if (Glob_ProcID/=0) return
    close(iu)
    on = .false.
  end subroutine opttrace_finalize

end module opttrace
