CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = ../postp
OPTM =  -O


VPATH =                       ../02_parsers:../04_residue_scoring:../06_main:\
        ../01_data_structures:../03_output:../05_cluster_scoring: ../07_utils


PLATFORM = $(shell uname)



LIBS =  -lm -lc
INCLUDE = -I ../00_include


DEFINES =  
SOURCES.F = 
SOURCES.f = 
SOURCES.c = postp.c  postp_aa_code.c postp_adj_matrix.c  postp_alignment.c  postp_area.c \
	postp_clustering_out.c  postp_cluster_weight.c  \
	postp_clustercounter.c  postp_entropy.c postp_hybrid.c  postp_int_trace.c postp_logger.c postp_mapping.c postp_node.c\
	postp_read_cmd_file.c postp_read_cw.c postp_read_pdb.c postp_read_score.c  postp_score2coverage.c \
        postp_score_out.c  postp_scoring.c postp_score2rank.c  postp_seq_pw_dist.c postp_spearman.c  postp_tree.c  postp_utils.c



AUX =  makefile  al.h postp_alignment.h postp_geometry.h postp_utils.h

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


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


