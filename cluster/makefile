CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = genclust

OPTM =  -O


VPATH = ./

PLATFORM = $(shell uname)



LIBS =  -lm -lc

ifeq (${PLATFORM},Linux)
	PD =  -DLINUX 
endif

DEFINES =  $(PD)  -DVERBOSE
SOURCES.F = 
SOURCES.f = 
SOURCES.c =  general_clustering.c  clustercounter.c utils.c



AUX =  makefile utils.h 

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
SOURCES_ANN = $(SOURCES_ANN.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)
OBJECTS_ANN =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES_ANN.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

anneal:$(OBJECTS_ANN)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS_ANN) $(OPTM) \
	$(LIBS)    $(INCLUDE)


.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


