CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = sqc

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
SOURCES.c =  main.c utils.c read_msf.c find_distmat.c clustercounter.c



AUX =  makefile utils.h sqc.h

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)
SOURCES_OC.c = output_contmat_main.c utils.c read_msf.c find_distmat.c 
OBJECTS_OC = $(SOURCES_OC.c:.c=.o)

SOURCES_SN.c = seqnet_main.c utils.c read_msf.c find_distmat.c clustercounter.c net_stats.c
OBJECTS_SN = $(SOURCES_SN.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

outcont:$(OBJECTS_OC)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS_OC) $(OPTM) \
	$(LIBS)    $(INCLUDE)

seqnet:$(OBJECTS_SN)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS_SN) $(OPTM) \
	$(LIBS)    $(INCLUDE)


.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


