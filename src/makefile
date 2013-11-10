CXX	= g++
CFLAGS	= -g -Wall -O3 -std=c++0x
OMP	= -fopenmp
INCLUDE	= -I eigen/

all: dse

dse: dse.h dse.cpp DisjointSet.h DisjointSet.cpp
	$(CXX) $(CFLAGS) $(OMP) $(INCLUDE) dse.cpp DisjointSet.cpp -o dse

clean:
	rm -f dse
