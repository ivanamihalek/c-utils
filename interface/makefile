CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= gcc   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = vec

OPTM = -g


VPATH = ./

PLATFORM = $(shell uname)



LIBS =  -lm -L/home/protean5/imihalek/LAPACK  -llapack -lblas

ifeq (${PLATFORM},Linux)
	PD =  -DLINUX
endif

DEFINES =  $(PD)   -DVERBOSE
SOURCES.F = 
SOURCES.f = 
SOURCES.c =  bind_pot.c ifc_vec.c  main.c momin.c  process_pdb.c  neighbors.c read_pdb.c  utils.c
IFC_SRC =  ifc_main.c ifc.c  read_pdb.c 
AUX =  makefile ifc.h utils.h

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)
IFC_OBJ = $(IFC_SRC:.c=.o)

$(NAME): $(OBJECTS)
	$(FCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)

ifc: $(IFC_OBJ)
	$(FCOMPL) $(FFLAGS) -o $@  $(IFC_OBJ) $(OPTM) \
	$(LIBS)    $(INCLUDE)


ls:
	 ls $(SOURCES) $(AUX) 


