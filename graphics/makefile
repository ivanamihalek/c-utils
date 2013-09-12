CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f77   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = cpg

OPTM =  -O3

VPATH = ./

PLATFORM = $(shell uname)

INCLUDE = -I/home/protean5/imihalek/pgplot/

CLIBPATH = -L/usr/lib/gcc-lib/i386-redhat-linux/3.2 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2/../../..
LIBS =   -L/home/protean5/imihalek/pgplot -lcpgplot -lpgplot  -L/usr/X11R6/lib -lX11 $(CLIBPATH) -lgcc -lm -lc
 
ifeq (${PLATFORM},Linux)
	PD =  -DLINUX
endif

DEFINES =  $(PD)  -DVERBOSE 
SOURCES.F = palett.F
SOURCES.f = 
SOURCES.c = cpgmatrix.c



AUX =  

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(FCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

.F.o:
	$(FCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


