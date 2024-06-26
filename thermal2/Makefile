# Makefile for thermal2

include ../../make.inc

STRIP = true #strip

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../UtilXlib \
          $(MOD_FLAG)../../FFTXlib/src $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules $(MOD_FLAG)../../PHonon/PH \
	  $(MOD_FLAG)../src \
	  $(MOD_FLAG).

IFLAGS+=-I../../FFTXlib -I../../include -I../../PHonon/PH 

# Remove the following line if you want to use intel compiler but not MKL 
# vectorized vdCos and vdSin
DFLAGS+= -D__HASVTRIG

#F90FLAGS += -heap-arrays 10
#-check all -traceback -g
#F90FLAGS += -check all -traceback -g

#F90FLAGS += -qopenmp 
#LDFLAGS  += -qopenmp #-pg -g

# For profiling
#F90FLAGS += -p
#LDFLAGS += -p

# For PGI
#F90FLAGS += -g77libs
#F90FLAGS += -Mbounds -Mchkptr -Mchkstk -traceback


PHOBJS = ../../PHonon/PH/libph.a ../../PHonon/PH/libphaux.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/Davidson/libdavid.a ../../KS_Solvers/CG/libcg.a \
         ../../FFTXlib/src/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../upflib/libupf.a


ifdef SCALAPACK_LIBS
MINPACK = lapackified-libs distributed-libs
EXTOBJS = ../src/libd3q.a ../minpack/distributed/libminpack.a ../minpack/lapackified/libminpack.a \
	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
else
MINPACK = lapackified-libs
EXTOBJS = ../src/libd3q.a ../minpack/lapackified/libminpack.a \
	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
endif



TLDEPS= bindir mods libs pw-lib lrmods ph-lib

LIBOBJS = ../../FFTXlib/src/libqefft.a ../../LAXlib/libqela.a

SRC = $(shell ls *.f90 *.c|grep -v PROGRAM_)
OBJ = $(shell ls *.f90 *.c|grep -v PROGRAM_|sed -e 's/\.f90$$/.o/' -e 's/\.c$$/.o/')
print-%  : ; @echo $* = $($*)

THERMALK2_OBJS = \
mpi_thermal.o \
asr2.o \
casimir.o \
code_input.o \
decompose_d2.o \
decompose_zstar.o \
dynbubble.o \
eos.o \
f3_bwfft.o \
fc2_interp.o \
fc3_interp.o \
final_state.o \
functions.o \
gruneisen.0 \
harmonic.o \
input_fc.o \
isotopes.o \
linewidth.o \
lmdif_p.o \
mc_grids.o \
merge_degenerate.o \
more_constants.o \
nanoclock.o \
nanosec.o \
neutrons.o \
nist_isotopes_db.o \
overlap.o \
ph_dos.o \
ph_system.o \
ph_velocity.o \
posix_signal.o \
q_grids.o \
quter.o \
read_md.o \
rigid_d3.o \
timers.o \
variational_tk.o \
variational_tk_symq.o

#EXECUTABLES = d3_lw.x d3_tk.x d3_asr3.x d3_sqom.x d3_sparse.x d3_q2r.x d3_r2q.x d3_qq2rr.x d3_import_shengbte.x d3_recenter.x

MYNAME       = $(shell basename $(CURDIR))
FATHERNAME   = $(shell cd ..; basename `pwd -P`)
THELINKPLACE = $(FATHERNAME)/$(MYNAME)

EXEC = db tk lw qha asr3 sqom sparse q2r r2q qq2rr import_shengbte recenter interpolate2 rmzeu tdph mtx impurity
#make_wedge decompose_d2
#gruneisen
#interpolate2

all : tldeps d3q-libs thermal2-libs $(EXEC)

.PHONY : tldeps d3q-libs thermal2-libs

# Note to self: $^ expands to all dependencies, $< expands to the first dependency
$(EXEC) : % : PROGRAM_%.o tldeps d3q-libs thermal2-libs $(MINPACK)
	$(LD) $(LDFLAGS) -o d3_$@.x $< libthermal2.a \
	    $(EXTOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	$(STRIP)  d3_$@.x
	- ( cd ../../bin && ln -fs ../D3Q/thermal2/d3_$@.x d3_$@.x )

d3q-libs :
	( cd ../src/ && make libd3q.a)

lapackified-libs :
	( cd ../minpack/lapackified/ && make)

distributed-libs :
	( cd ../minpack/distributed/ && make)

thermal2-libs : libthermal2.a

libthermal2.a : $(OBJ)
#$(THERMALK2_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@ 

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean : 
	- /bin/rm -f *.x *.o *~ *.F90 *.mod *.d *.i *.L *.a

include make.depend

