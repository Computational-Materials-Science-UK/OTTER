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
.SUFFIXES:

FC=gfortran 
OPT=-O2 -funroll-all-loops 
EXE := otter.x
SRC_DIR := src
OBJ_DIR := obj

SRC := $(wildcard $(SRC_DIR)/*.f90)
OBJ := $(SRC:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
MAIN := $(OBJ_DIR)/otter_main.o
OBJ := $(filter-out $(MAIN),$(OBJ)) $(MAIN) 

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) -o $@ $^
	mv *.mod $(OBJ_DIR)/.

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(OPT) -c -o $@ $^ 

$(OBJ_DIR):
	mkdir $@ 

clean:
	$(RM) $(OBJ) $(OBJ_DIR)/*.mod *~