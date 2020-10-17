FC=gfortran 
OPT=-O1 

PROG_OBJ=otter_cylinders_beck1.o
MODS=otter_input.o otter_cyl_rotate.o otter_cyl_contact.o 
OBJ=$(MODS) $(PROG_OBJ)
PROG=otter_cyl

default : $(PROG)

$(OBJ) : %.o : %.f90
	$(FC) $(OPT) -c $<

$(PROG) : $(OBJ)
	$(FC) $(OPT) -o $@ $^

$(PROG_OBJ) : $(MODS)

clean:
	-rm -f *.o *~
