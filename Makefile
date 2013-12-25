CXX=g++
CC=gcc
CXXFLAGS=-g -Wall -DWITH_GLPK

all: practip

clean:
	rm -f practip.o ip.o cmdline.o

practip: practip.o ip.o cmdline.o
	$(CXX) -o practip practip.o ip.o cmdline.o -lglpk

practip.o: practip.cpp ip.h cmdline.h
ip.o: ip.cpp
cmdline.o: cmdline.c
