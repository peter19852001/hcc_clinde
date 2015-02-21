.SUFFIXES: .o .c
CC = gcc 
#DEBUG = -g -Wall
DEBUG = -Wall
FLAGS = $(DEBUG) -O2 -static
LIBS = -lgsl -lgslcblas -lm

CORE = clinde.o parse_option.o tsv.o

EXEC = clinde

all: $(EXEC)

.c.o: 
	$(CC) $(FLAGS) -c $< -o $@

core : $(CORE)

$(EXEC) : $(CORE) $(OLSLIB)
	$(CC) $(FLAGS) $(CORE) -o $(EXEC) $(LIBS)
	@echo 'Made '

clean:
	rm -f *.o *.bak *.*~
	rm -f $(EXEC)
	@echo 'Made'
