### OTTER - makefile
#
# OTTER: Complex Structure Generation Toolkit
#	Maintained by the Computational Materials Science Research Group at the University of Kentucky, Dr. Matthew J. Beck, PI.
#	https://www.beckdt.engr.uky.edu
#
# OTTER master: Matthew Beck, m.beck@uky.edu
# OTTER developers: OTTER-dev team at GitHub.org
#
# https://github.com/Computational-Materials-Science-UK/OTTER
#
###

#####
# makfile - base makefile for GNU make
#	Vb.1 - Working version

# External Dependencies:
#	GNU Make
# 	Fortran compiler (default: gfortran)
# Internal Dependencies:
#	Full OTTER codebase

FC=gfortran 
OPT=-O2 -funroll-all-loops 

MODS= 

OBJ=otter_spheres.o otter_fibers.o otter_ligaments.o otter_main.o

otter.x : $(OBJ)
	$(FC) -o $@ $^

%.mod: %.f90
	$(FC) -c $^

%.o: %.f90 $(MODS)
	$(FC) -c -o $@ $^

clean:
	-rm -f *.o *~
