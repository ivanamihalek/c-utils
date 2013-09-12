
CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= g77   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = lsf3

OPTM =  -O

VPATH = ./

PLATFORM = $(shell uname)



LIBS =  -lm  -llapack -lblas

ifeq (${PLATFORM},Linux)
	PD =  -DLINUX
endif

DEFINES =  $(PD)  
SOURCES.F = 
SOURCES.f = 
SOURCES.c =   leastsq3d.c


AUX = 
SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)



$(NAME): $(OBJECTS)
	$(FCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)


.F.o:
	$(FCOMPL)   $(FFLAGS) -c   $(DEFINES)  $<  $(OPTM) $(INCLUDE)
.f.o:
	$(FCOMPL)  $(FFLAGS) -c   $(DEFINES)  $<  $(OPTM) $(INCLUDE)
.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


