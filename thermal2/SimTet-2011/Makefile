include ../../../make.inc

%.o: %.f
	$(MPIF90) $(F90FLAGS) -c $< -o $(@)

OBJ = $(shell ls *.f | sed -e 's/\.f$$/.o/')


libsimtet.a : $(OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.mod *.d *.i *.L *.a