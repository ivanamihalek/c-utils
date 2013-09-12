CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = ../smallig
OPTM =  -g


VPATH = ../src


PLATFORM = $(shell uname)


LIBS =  -lm -lc  -L/usr/lib  -llapack
INCLUDE = -I ../src


DEFINES =   # -DDMALLOC # -DDEBUG
SOURCES.F = 
SOURCES.f = 


SOURCES.c =  smallig.c smallig_read_pdb.c smallig_utils.c


AUX =  makefile  smallig.h  smallig_pdb.h  README 


SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

$(PREP_NAME): $(PREP_OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(PREP_OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

$(TEST_NAME): $(TEST_OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(TEST_OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)



.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(PREP_OBJECTS) $(TEST_OBJECTS) $(OBJECTS) 

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


