#####################################################################
#####################################################################
##                                                                 ##
##               Makefile for 2D MPI/OpenMP UPIC-EMMA              ##
##                                                                 ##
#####################################################################
#####################################################################

#################
# GNU compilers #
#################

#OpenMPI
MPIFC = mpif90 -fopenmp
MPICC = mpicc -fopenmp

CC = gcc

#OPTSF90 = -O3
OPTSF90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTSF90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8
#OPTSF90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95
#OPTSF90 = -O3 -fcheck=bounds
#OPTSF90 = -O3 -g -fbacktrace -fcheck=all -Wall

#OPTSF77 = -O3
OPTSF77 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTSF77 = -O3 -fcheck=bounds -fdefault-real-8 -Wall
#OPTSF77 = -O1 -g -fbacktrace #-fdefault-real-8 -Wall

OPTSCC = -O3 -std=c99
#OPTSCC = -O3 -Wall -std=c99

###################
# Intel compilers #
###################

## MPI
#MPIFC = mpiifort -openmp
#MPICC = mpiicc -openmp
#
#CC = icc
#
#OPTSF90 = -O3
#OPTSF90 = -O3 -r8
#OPTSF90 = -O3 -CB -r8 -warn all -std90
#
#OPTSF77 = -O3
#OPTSF77 = -O3 -r8
#OPTSF77 = -O3 -CB -r8 -warn all -std77
#
#OPTSCC = -O3 -std=c99
#OPTSCC = -O3 -no-vec -Wall -std=c99
#LEGACY = -nofor_main

###############
# Directories #
###############

SRC_PATH_F90 = source/f90/
SRC_PATH_PY  = source/py/
SRC_PATH_F77 = source/f77/
SRC_PATH_C   = source/c/

OBJ_PATH =

MOD_PATH = 

TRG_PATH = 

#############
# shortcuts #
#############
		
ESOBJS = $(addprefix $(OBJ_PATH), libmpinit2.o \
        libmpsort2.o libmpgard2.o \
        libmpfft2.o  libmpfield2.o \
		libmpfieldpml2.o )

EMOBJS = $(addprefix $(OBJ_PATH), libmpbpush2.o \
        libmpcurd2.o )

ESHOBJS = $(addprefix $(OBJ_PATH), libmpinit2_h.o \
        libmpsort2_h.o libmpgard2_h.o \
        libmpfft2_h.o  libmpfield2_h.o \
		libmpfieldpml2_h.o )

EMHOBJS = $(addprefix $(OBJ_PATH), libmpbpush2_h.o libmpcurd2_h.o )

ESMODS = $(addprefix $(OBJ_PATH), rvm.o modprofile2.o modmpinit2.o \
        modmpsort2.o modmpgard2.o  \
        modmpfft2.o  modmpfield2.o \
        modmpfieldpml2.o antenna.o\
        diag.o  input.o )

EMMODS = $(addprefix $(OBJ_PATH), modmpbpush2.o modmpcurd2.o )

#####################
# Compilation rules #
#####################

all : upic-emma clean_mod clean_obj

upic-emma : $(OBJ_PATH)upic-emma.o $(ESOBJS) $(EMOBJS) dtimer.o 
	$(MPIFC) $(OPTSF90) -o $(TRG_PATH)upic-emma $(OBJ_PATH)upic-emma.o \
	$(ESOBJS) $(EMOBJS) $(ESMODS) $(EMMODS) \
	$(OBJ_PATH)mpplib2.o $(OBJ_PATH)mppmod2.o \
	$(OBJ_PATH)omplib.o $(OBJ_PATH)ompplib2.o \
	$(ESHOBJS) $(EMHOBJS) dtimer.o

dtimer.o : $(SRC_PATH_C)dtimer.c
	$(CC) $(OPTSCC) -c $(SRC_PATH_C)dtimer.c

$(OBJ_PATH)omplib.o : $(SRC_PATH_F90)omplib.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)omplib.o -c $(SRC_PATH_F90)omplib.f90

$(OBJ_PATH)mpplib2.o : $(SRC_PATH_F90)mpplib2.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)mpplib2.o -c $(SRC_PATH_F90)mpplib2.f90

$(OBJ_PATH)mppmod2.o : $(SRC_PATH_F90)mppmod2.f90 $(OBJ_PATH)mpplib2.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)mppmod2.o -c $(SRC_PATH_F90)mppmod2.f90

$(OBJ_PATH)libmpinit2.o : $(SRC_PATH_F77)libmpinit2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpinit2.o -c $(SRC_PATH_F77)libmpinit2.f

# $(OBJ_PATH)libmppush2.o : $(SRC_PATH_F77)libmppush2.f
# 	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmppush2.o -c $(SRC_PATH_F77)libmppush2.f

$(OBJ_PATH)libmpbpush2.o : $(SRC_PATH_F77)libmpbpush2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpbpush2.o -c $(SRC_PATH_F77)libmpbpush2.f

$(OBJ_PATH)libmpcurd2.o : $(SRC_PATH_F77)libmpcurd2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpcurd2.o -c $(SRC_PATH_F77)libmpcurd2.f

$(OBJ_PATH)libmpsort2.o : $(SRC_PATH_F77)libmpsort2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpsort2.o -c $(SRC_PATH_F77)libmpsort2.f

$(OBJ_PATH)libmpgard2.o : $(SRC_PATH_F77)libmpgard2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpgard2.o -c $(SRC_PATH_F77)libmpgard2.f

$(OBJ_PATH)libmpfft2.o : $(SRC_PATH_F77)libmpfft2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpfft2.o -c $(SRC_PATH_F77)libmpfft2.f

$(OBJ_PATH)libmpfield2.o : $(SRC_PATH_F77)libmpfield2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpfield2.o -c $(SRC_PATH_F77)libmpfield2.f
	
$(OBJ_PATH)libmpfieldpml2.o : $(SRC_PATH_F77)libmpfieldpml2.f
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpfieldpml2.o -c $(SRC_PATH_F77)libmpfieldpml2.f

$(OBJ_PATH)libmpinit2_h.o : $(SRC_PATH_F90)libmpinit2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpinit2_h.o -c $(SRC_PATH_F90)libmpinit2_h.f90

# $(OBJ_PATH)libmppush2_h.o : $(SRC_PATH_F90)libmppush2_h.f90
# 	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmppush2_h.o -c $(SRC_PATH_F90)libmppush2_h.f90

$(OBJ_PATH)libmpbpush2_h.o : $(SRC_PATH_F90)libmpbpush2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpbpush2_h.o -c $(SRC_PATH_F90)libmpbpush2_h.f90

$(OBJ_PATH)libmpcurd2_h.o : $(SRC_PATH_F90)libmpcurd2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpcurd2_h.o -c $(SRC_PATH_F90)libmpcurd2_h.f90

$(OBJ_PATH)libmpsort2_h.o : $(SRC_PATH_F90)libmpsort2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpsort2_h.o -c $(SRC_PATH_F90)libmpsort2_h.f90

$(OBJ_PATH)libmpgard2_h.o : $(SRC_PATH_F90)libmpgard2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpgard2_h.o -c $(SRC_PATH_F90)libmpgard2_h.f90

$(OBJ_PATH)libmpfft2_h.o : $(SRC_PATH_F90)libmpfft2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpfft2_h.o -c $(SRC_PATH_F90)libmpfft2_h.f90

$(OBJ_PATH)libmpfield2_h.o : $(SRC_PATH_F90)libmpfield2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpfield2_h.o -c $(SRC_PATH_F90)libmpfield2_h.f90
	
$(OBJ_PATH)libmpfieldpml2_h.o : $(SRC_PATH_F90)libmpfieldpml2_h.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)libmpfieldpml2_h.o -c $(SRC_PATH_F90)libmpfieldpml2_h.f90

$(OBJ_PATH)rvm.o : $(SRC_PATH_F90)parser.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)rvm.o -c $(SRC_PATH_F90)parser.f90

$(OBJ_PATH)modprofile2.o : $(SRC_PATH_F90)modprofile2.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modprofile2.o -c $(SRC_PATH_F90)modprofile2.f90

$(OBJ_PATH)modmpinit2.o : $(SRC_PATH_F90)modmpinit2.f90 $(OBJ_PATH)libmpinit2_h.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpinit2.o -c $(SRC_PATH_F90)modmpinit2.f90

$(OBJ_PATH)modmpbpush2.o : $(SRC_PATH_F90)modmpbpush2.f90 $(OBJ_PATH)libmpbpush2_h.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpbpush2.o -c $(SRC_PATH_F90)modmpbpush2.f90

$(OBJ_PATH)modmpcurd2.o : $(SRC_PATH_F90)modmpcurd2.f90 $(OBJ_PATH)libmpcurd2_h.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpcurd2.o -c $(SRC_PATH_F90)modmpcurd2.f90

$(OBJ_PATH)modmpsort2.o : $(SRC_PATH_F90)modmpsort2.f90 $(OBJ_PATH)libmpsort2_h.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpsort2.o -c $(SRC_PATH_F90)modmpsort2.f90

$(OBJ_PATH)modmpgard2.o : $(SRC_PATH_F90)modmpgard2.f90 $(OBJ_PATH)libmpgard2_h.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpgard2.o -c $(SRC_PATH_F90)modmpgard2.f90

$(OBJ_PATH)modmpfft2.o : $(SRC_PATH_F90)modmpfft2.f90 $(OBJ_PATH)libmpfft2_h.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpfft2.o -c $(SRC_PATH_F90)modmpfft2.f90

$(OBJ_PATH)modmpfield2.o : $(SRC_PATH_F90)modmpfield2.f90 $(OBJ_PATH)libmpfield2_h.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpfield2.o -c $(SRC_PATH_F90)modmpfield2.f90

$(OBJ_PATH)ompplib2.o : $(SRC_PATH_F90)ompplib2.f90 $(OBJ_PATH)modmpsort2.o $(OBJ_PATH)modmpfft2.o \
	$(OBJ_PATH)modmpgard2.o $(OBJ_PATH)mppmod2.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)ompplib2.o -c $(SRC_PATH_F90)ompplib2.f90

$(OBJ_PATH)modmpfieldpml2.o : $(SRC_PATH_F90)modmpfieldpml2.f90 $(OBJ_PATH)libmpfieldpml2_h.o $(OBJ_PATH)ompplib2.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)modmpfieldpml2.o -c $(SRC_PATH_F90)modmpfieldpml2.f90

$(OBJ_PATH)diag.o : $(SRC_PATH_F90)diag.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)diag.o -c $(SRC_PATH_F90)diag.f90

$(OBJ_PATH)antenna.o : $(SRC_PATH_F90)antenna.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)antenna.o -c $(SRC_PATH_F90)antenna.f90
	
$(OBJ_PATH)input.o : $(SRC_PATH_F90)input.f90
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)input.o -c $(SRC_PATH_F90)input.f90

$(OBJ_PATH)upic-emma.o : $(SRC_PATH_F90)upic-emma.f90 $(ESMODS) $(EMMODS) $(OBJ_PATH)mppmod2.o $(OBJ_PATH)omplib.o $(OBJ_PATH)ompplib2.o
	$(MPIFC) $(OPTSF90) -o $(OBJ_PATH)upic-emma.o -c $(SRC_PATH_F90)upic-emma.f90

extract : $(SRC_PATH_PY)extract.py
	python $(SRC_PATH_PY)extract.py

clean_mod :
	rm -f $(MOD_PATH)*.mod

clean_obj :
	rm -f $(OBJ_PATH)*.o

clean :
	rm -f $(TRG_PATH)upic-emma $(OBJ_PATH)*.o $(MOD_PATH)*.mod

clean_results :
	rm -rf results/
	
clean_plots :
	rm -rf plots/

clean_all : clean clean_results clean_plots
