# -*- Mode: makefile; -*-
# Makefile for a FEA finite-strains project

CC = gcc

#EXPAT_ROOT = ../../../../expat-2.0.1
LIBSEXP_PATH = ../../../libsexp
LIBSPM_PATH = ../../../libspmatrix
LOGGER_PATH = ../../../liblogger

CFLAGS = -ggdb  --std=c99 -pedantic -Wall -Wextra -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wdeclaration-after-statement -Wmissing-declarations 
DEFINES = -DCURRENT_SHAPE_GRADIENTS 
#-DDUMP_DATA -DUSE_EXPAT 

INCLUDES = -I $(EXPAT_ROOT)/lib -I $(LIBSEXP_PATH) -I $(LIBSPM_PATH)/inc -I $(LOGGER_PATH)
LINKFLAGS =  -L $(LIBSEXP_PATH) -lsexp -L $(LIBSPM_PATH)/lib -lspmatrix -L $(LOGGER_PATH) -llogger  -lm -rdynamic -lrt 
# $(EXPAT_ROOT)/.libs/libexpat.a 

SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
OUTPUT = feasolver

%.o : %.c %.h defines.h
	$(CC) -c $(CFLAGS) $(DEFINES) $(INCLUDES) $< -o $@

$(OUTPUT): $(OBJECTS)
	$(CC) $(OBJECTS) $(LINKFLAGS) -o $(OUTPUT) 

all: $(OUTPUT)

lint:
	splint +booltype=BOOL -boolfalse=FALSE -booltrue=TRUE *.c

.PHONY : clean
clean :
	rm $(OBJECTS) $(OUTPUT)

check-syntax: 
	gcc -o nul -S ${CHK_SOURCES} 
