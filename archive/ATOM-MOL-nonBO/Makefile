F90 = mpiifort
F77 = mpiifort

OPT=-O3 -ip -ftz -align all -fp-model precise -I$(MKLROOT)/include

FILES = globvars.f90 misc.f90 dmng.f linalg.f90 matelem.f90 matform.f90 workproc.f90 main.f90
OBJS = globvars.o misc.o dmng.o linalg.o matelem.o matform.o workproc.o main.o
LIBS =  -L$(MKLROOT) -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

main	 : $(OBJS)
	$(F90) -o main $(OPTL) $(OBJS) $(LIBS) 

main.o	 : $(FILES)
	$(F90) -c $(OPT) $?

globvars.o : globvars.f90
	$(F90) -c $(OPT) $(LIBS) $?

misc.o : misc.f90
	$(F90) -c $(OPT) $?

dmng.o : dmng.f
	$(F77) -c $(OPT) $?

linalg.o : linalg.f90 
	$(F90) -c $(OPT) $?

matelem.o : matelem.f90 globvars.f90
	$(F90) -c $(OPT) $?

matform.o : matform.f90 matelem.f90 
	$(F90) -c $(OPT) $?

workproc.o : workproc.f90 
	$(F90) -c $(OPT) $?
