# Makefile for a FEA finite-strains project

CC = gcc
EXPAT_ROOT = ../../../../expat-2.0.1
CFLAGS = -ggdb -pg -ansi -pedantic -Wall -Wextra -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wdeclaration-after-statement -Wmissing-declarations -DUSE_EXPAT -DCURRENT_SHAPE_GRADIENTS
#-DDUMP_DATA

INCLUDES = -I $(EXPAT_ROOT)/lib -I $(DMALLOC_ROOT)
LINKFLAGS =  $(EXPAT_ROOT)/.libs/libexpat.a -pg
objects = fea_solver.o
output = feasolver

fea_solver : $(objects)
	$(CC) $(objects) $(LINKFLAGS) -o $(output) 

fea_solver.o: fea_solver.c
	$(CC) -c $(CFLAGS) $(INCLUDES) fea_solver.c

.PHONY : clean
clean :
	rm $(objects) $(output)

check-syntax: 
	gcc -o nul -S ${CHK_SOURCES} 
