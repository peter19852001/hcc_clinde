.SUFFIXES: .o .c
CC = gcc 
#DEBUG = -g -Wall
DEBUG = -Wall
FLAGS = $(DEBUG) -O2 -static
LIBS = -lgsl -lgslcblas -lm

CORE = hcc_clinde.o parse_option.o tsv.o

CORE2 = grn_cmp_hcc.o parse_option.o

EXEC = hcc_clinde

EXEC2 = grn_cmp_hcc

all: $(EXEC) $(EXEC2)

.c.o: 
	$(CC) $(FLAGS) -c $< -o $@

core : $(CORE)

$(EXEC) : $(CORE)
	$(CC) $(FLAGS) $(CORE) -o $(EXEC) $(LIBS)
	@echo 'Made '

$(EXEC2) : $(CORE2)
	$(CC) $(FLAGS) $(CORE2) -o $(EXEC2) -lm
	@echo 'Made '

clean:
	rm -f *.o *.bak *.*~
	rm -f $(EXEC)
	@echo 'Made'
