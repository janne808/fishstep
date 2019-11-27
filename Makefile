# compile flags
SDL=1
TIFF_ENABLE=0
OPTIMIZATION_LEVEL=3

# object files
OBJ=fishstep.o

# compilers
CC=gcc

# SDL options
ifeq ($(SDL),1)
	SDL_OPTS=-DSDL=1 `sdl-config --libs` `sdl-config --cflags`
else
	SDL_OPTS=
endif

# TIFFLib options
ifeq ($(TIFF_ENABLE),1)
	TIFFLIB_OPTS=-DTIFF_ENABLE=1 -ltiff
else
	TIFFLIB_OPTS=
endif

# compiler options
OPTS=-Wall -pthread 
#CFLAGS=-O$(OPTIMIZATION_LEVEL) -lm -lrt -lfftw3
CFLAGS=-g -lm -lrt -lfftw3

fishstep: 	$(OBJ)
	$(CC) -o $@ $+ $(OPTS) $(CFLAGS) $(SDL_OPTS) $(TIFFLIB_OPTS)

fishstep.o:	fishstep.c
	$(CC) $(OPTS) $(CFLAGS) $(SDL_OPTS) $(TIFFLIB_OPTS) -c $<

.PHONY: clean
clean:
	rm *.o fishstep

