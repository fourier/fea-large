# Makefile for a FEA finite-strains project

CC = gcc
CFLAGS = -O3 -g -ansi -Wall 
objects = fea_solver.o
output = feasolver

fea_solver : $(objects)
	$(CC) $(objects) -o $(output) 

fea_solver.o: fea_solver.c
	$(CC) -c $(CFLAGS) fea_solver.c

.PHONY : clean
clean :
	rm $(objects) $(output)