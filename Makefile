#
# COMET - SDSC 
#
# Load Intel compiler prior to make any target, e.g., module load intel
# For other compilers please check the module list
#
CC = icpc
MPCC = mpicxx
OPENMP = -openmp
CFLAGS = -O3
LIBS =

TARGETS = graph_compute

all:	$(TARGETS)

#openmp: openmp.o common.o
#	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o

graph_compute: graph_compute.o utility.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) graph_compute.o utility.o 

graph_compute.o: graph_compute.cpp graph_compute.h utility.h
	$(MPCC) -c $(CFLAGS) graph_compute.cpp

utility.o: utility.cpp utility.h
	$(CC) -c $(CFLAGS) utility.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout
