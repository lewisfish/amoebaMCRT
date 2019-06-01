ifeq ($(FCOMP),mpifort)
FCFLAGS = -freal-4-real-8 -O2
FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
FCBUILD = -Wall -Wextra -pedantic -std=f2008
else ifeq ($(FCOMP),mpiifort)
FCFLAGS = -O2 -xHOST -r8 
FCDEBUG = -g -debug -traceback -check all -CB -CU -CA -fpe0
#-stand f08 -assume realloc_lhs  -check all  -traceback  -warn all  -fstack-protector  -assume protect_parens  -implicitnone
FCBUILD = -warn all -stand f08
endif
PROGRAM =  mcgrid

SRCS =      ran2.f \
		utils.f90 \
            constants.f90 \
            photon_vars.f90 \
            iarray.f90 \
            opt_prop.f90 \
            subs.f90 \
            fluorophores.f90 \
            absorbers_mod.f90 \
            skin_layers_mod.f90 \
            ch_opt.f90 \
            gridset.f90 \
            reader.f90 \
            vector_class.f90 \
            inttau2.f90 \
            stokes.f90 \
            sourceph.f90 \
            writer.f90 \
            mcpolar.f90 \
            2Dsimplex_mod.f90 \
            simpRun.f90

#                         trackPacket_mod.f90 \

#             getFluence_mod.f90 \

#             main2.f90 

                      
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)
debug:	FCFLAGS += $(FCDEBUG)
debug:	$(PROGRAM)
build:	FCFLAGS += $(FCBUILD)
build:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 
	
%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgridf mcgrid
