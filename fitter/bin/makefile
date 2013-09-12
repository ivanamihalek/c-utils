CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = ../fitter
OPTM =  -g


VPATH = ../src


PLATFORM = $(shell uname)


LIBS =  -lm -lc  -L/usr/lib  -llapack
INCLUDE = -I ../src


DEFINES =   # -DDMALLOC # -DDEBUG
SOURCES.F = 
SOURCES.f = 


SOURCES.c =  fitter.c fitter_read_list.c fitter_read_pdb.c fitter_map2rotation.c fitter_output.c fitter_utils.c


AUX =  makefile  fitter.h  fitter_pdb.h  README 


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


