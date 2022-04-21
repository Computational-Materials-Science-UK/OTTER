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
#	Vb.2 - Added incorporated ligaments code 02/15/2021

# External Dependencies:
#	GNU Make
# 	Fortran compiler (default: gfortran)
# Internal Dependencies:
#	Full OTTER codebase
.SUFFIXES:

FC=gfortran 
#OPT=-c -fbacktrace -Wall -Wextra -Wno-maybe-uninitialized -Wno-unused-function -pedantic -std=f2008ts -O2
OPT= -O2 -funroll-all-loops -fdefault-real-8 -fdefault-double-8 
EXE := otter.x
SRC_DIR := src
OBJ_DIR := obj

SRC := $(wildcard $(SRC_DIR)/*.f90)
#OBJ := $(SRC:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
MAIN := $(OBJ_DIR)/otter_main.o
#OBJ := $(filter-out $(MAIN),$(OBJ))
TEST := $(OBJ_DIR)/otter_test.o
#OBJ := $(filter-out $(TEST),$(OBJ))
OBJ := $(OBJ_DIR)/bspline_kinds_module.o $(OBJ_DIR)/bspline_sub_module.o \
    $(OBJ_DIR)/bspline_oo_module.o $(OBJ_DIR)/bspline_module.o \
    $(OBJ_DIR)/otter_01_globals.o \
	$(OBJ_DIR)/otter_02_math.o \
	$(OBJ_DIR)/otter_02_ligaments_globals.o\
    $(OBJ_DIR)/otter_03_spheres_globals.o \
	$(OBJ_DIR)/otter_03_fibers_globals.o \
	$(OBJ_DIR)/otter_03_input.o \
    $(OBJ_DIR)/otter_fibers.o $(OBJ_DIR)/otter_ligaments.o \
    $(OBJ_DIR)/otter_spheres.o \
    
    	
 
.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) $(MAIN)
	$(FC) -o $@ $^	
	mv *.mod $(OBJ_DIR)/.
#Smv $(SRC_DIR)/*.mod $(OBJ_DIR)/.

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(OPT) -c -o $@ $^ 

$(OBJ_DIR):
	mkdir $@ 

test: otter_test.x

otter_test.x: $(OBJ) $(TEST)
	$(FC) -o $@ $^
	mv *.mod $(OBJ_DIR)/.

clean:
	$(RM) $(OBJ) $(MAIN) $(TEST) *.x $(OBJ_DIR)/*.mod *~
