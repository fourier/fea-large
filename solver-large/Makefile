# -*- Mode: makefile; -*-
# Makefile for a FEA finite-strains project

CC = gcc

EXPAT_ROOT = ../../../../expat-2.0.1

CFLAGS = -ggdb -pg -ansi -pedantic -Wall -Wextra -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wdeclaration-after-statement -Wmissing-declarations -fopenmp
DEFINES = -DUSE_EXPAT -DCURRENT_SHAPE_GRADIENTS
#-DDUMP_DATA

INCLUDES = -I $(EXPAT_ROOT)/lib
LINKFLAGS =  $(EXPAT_ROOT)/.libs/libexpat.a -pg

SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
OUTPUT = feasolver

%.o : %.c %.h defines.h
	$(CC) -c $(CFLAGS) $(DEFINES) $(INCLUDES) $< -o $@

$(OUTPUT): $(OBJECTS)
	$(CC) -o $(OUTPUT) $(OBJECTS) $(LINKFLAGS)

all: $(OUTPUT)

lint:
	splint +booltype=BOOL -boolfalse=FALSE -booltrue=TRUE *.c

.PHONY : clean
clean :
	rm $(OBJECTS) $(OUTPUT)

check-syntax: 
	gcc -o nul -S ${CHK_SOURCES} 
