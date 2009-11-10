# Makefile for a FEA finite-strains project

CC = gcc
EXPAT_ROOT = /Users/alexey/Sources/expat-2.0.1
CFLAGS = -O3 -g -ansi -Wall -pedantic -DUSE_EXPAT
INCLUDES = -I $(EXPAT_ROOT)/lib
LINKFLAGS =  $(EXPAT_ROOT)/.libs/libexpat.a
objects = fea_solver.o
output = feasolver

fea_solver : $(objects)
	$(CC) $(objects) $(LINKFLAGS) -o $(output) 

fea_solver.o: fea_solver.c
	$(CC) -c $(CFLAGS) $(INCLUDES) fea_solver.c

.PHONY : clean
clean :
	rm $(objects) $(output)