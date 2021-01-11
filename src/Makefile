#######################################################################
##                          DYNAMON MAKEFILE                         ##
#######################################################################

# name of the final executable, destination and fDynamo path
EXE     ?= dynamon
EXEPATH ?= ../bin
DYNPATH ?= ../fDynamo

# compiler and flags
FC       = gfortran
Olevel   = -O3
FFLAGS   = $(Olevel) -finit-local-zero -frecord-marker=4 \
           -J. -I$(DYNPATH) -I$(PWD)

#----------------------------------------------------------------------

OBJS := utils.f90 \
        initialization.f90 \
        nofix_qm.f90 \
        common.f90 \
        panadero.f90 \
        calculation_modes.f90 \
        dynamon.f90 \
        $(DYNPATH)/dynamo.a

#----------------------------------------------------------------------

.DEFAULT_GOAL := dynamon
.PHONY        := dynamon fdynamo clean clean_all

dynamon: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(EXEPATH)/$(EXE)

fdynamo: $(DYNPATH)/dynamo.a

$(DYNPATH)/dynamo.a:
	@printf "Compiling fDynamo... "
	@make -s -C $(DYNPATH) && printf "DONE \n"

clean:
	rm -f *.o *.mod

clean_all: clean
	@make clean -C $(DYNPATH)
