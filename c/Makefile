# Makefile for a FEA finite-strains project

CC = gcc
EXPAT_ROOT = ../../../../expat-2.0.1
DMALLOC_ROOT = ../../../../dmalloc-5.5.2
CFLAGS =  -g -O3 -ansi -Wall -pedantic -DUSE_EXPAT -DDMALLOC -DDMALLOC_FUNC_CHECK -DDUMP_DATA
#-DCURRENT_SHAPE_GRADIENTS


INCLUDES = -I $(EXPAT_ROOT)/lib -I $(DMALLOC_ROOT)
LINKFLAGS =  $(EXPAT_ROOT)/.libs/libexpat.a $(DMALLOC_ROOT)/libdmalloc.a
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