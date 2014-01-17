CXX=g++
CC=gcc

# for GLPK
CXXFLAGS_IP=-DWITH_GLPK
LDFLAGS_IP=-lglpk

# for CPLEX
#CPLEX_BASE=/opt/ibm/ILOG/CPLEX_Studio126
#ARCH=x86-64_linux
#CXXFLAGS_IP=-DWITH_CPLEX -DIL_STD -I$(CPLEX_BASE)/concert/include -I$(CPLEX_BASE)/cplex/include 
#LDFLAGS_IP=-L$(CPLEX_BASE)/concert/lib/$(ARCH)/static_pic -L$(CPLEX_BASE)/cplex/lib/$(ARCH)/static_pic -lconcert -lilocplex -lcplex -lpthread

#CXXFLAGS=-g -std=c++11 -Wall -D_GLIBCXX_DEBUG $(CXXFLAGS_IP)
CXXFLAGS=-O2 -std=c++11 -Wall -DNDEBUG $(CXXFLAGS_IP)
LDFLAGS=$(LDFLAGS_IP)

all: practip

clean:
	rm -f practip.o ip.o cmdline.o

practip: practip.o ip.o cmdline.o
	$(CXX) -g -o practip practip.o ip.o cmdline.o $(LDFLAGS)

practip.o: practip.cpp ip.h cmdline.h
ip.o: ip.cpp
cmdline.o: cmdline.c
