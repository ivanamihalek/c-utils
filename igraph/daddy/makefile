CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++ 

HERE = $(shell pwd)

NAME = daddy

OPTM =   -O # -g -O0 # recommended combo for valgrind

VPATH = ./

PLATFORM = $(shell uname)


LIBS =  -lm -lc -lpthread  -ligraph
INCLUDE =  -I/usr/local/include/igraph

ifeq (${PLATFORM},Linux)
	PD =  -DLINUX 
endif

DEFINES =  $(PD)  -DVERBOSE
SOURCES.F = 
SOURCES.f = 
SOURCES.c =  igraph_daddy.c utils.c threading.c solver.c

AUX =  makefile igraph_daddy.h

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)
SOURCES_OC.c = output_contmat_main.c utils.c read_msf.c find_distmat.c 
OBJECTS_OC = $(SOURCES_OC.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)


.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


