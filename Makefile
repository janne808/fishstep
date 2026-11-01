# compile flags
SDL=1
TIFF_ENABLE=0
OPTIMIZATION_LEVEL=3
FFTW3_THREADS=1
# DEBUG=1: debug symbols, no optimization; DEBUG=2: also address/UB sanitizers
DEBUG=0

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

# FFTW3 options
ifeq ($(FFTW3_THREADS),1)
	FFTW3_OPTS=-DFFTW3_THREADS=1 -lfftw3_threads
else
	FFTW3_OPTS=
endif

# debug options
ifeq ($(DEBUG),1)
	OPTIMIZATION_LEVEL=0
	DEBUG_OPTS=-g3
else ifeq ($(DEBUG),2)
	OPTIMIZATION_LEVEL=0
	DEBUG_OPTS=-g3 -fsanitize=address,undefined -fno-omit-frame-pointer
else
	DEBUG_OPTS=
endif

# compiler options
OPTS=-Wall -pthread $(DEBUG_OPTS)
CFLAGS=-O$(OPTIMIZATION_LEVEL) -lm -lrt -lfftw3

fishstep: 	$(OBJ)
	$(CC) -o $@ $+ $(OPTS) $(CFLAGS) $(SDL_OPTS) $(TIFFLIB_OPTS) $(FFTW3_OPTS)

fishstep.o:	fishstep.c
	$(CC) $(OPTS) $(CFLAGS) $(SDL_OPTS) $(TIFFLIB_OPTS) $(FFTW3_OPTS) -c $<

.PHONY: clean shocktube gascloud stardisk
clean:
	rm -f *.o fishstep fishstep_shocktube fishstep_gascloud fishstep_stardisk

# Sod shock tube test: gas only, no gravity, no window, writes shocktube.dat
fishstep_shocktube:	fishstep.c fishstep.h
	$(CC) -o $@ fishstep.c $(OPTS) -DSHOCK_TUBE=1 -DSDL=0 $(CFLAGS) $(FFTW3_OPTS)

shocktube:	fishstep_shocktube
	./fishstep_shocktube
	python3 tests/shocktube.py shocktube.dat


# self-gravitating polytropic gas cylinder: gravity and gas, no window,
# writes gascloud.dat
fishstep_gascloud:	fishstep.c fishstep.h
	$(CC) -o $@ fishstep.c $(OPTS) -DGAS_CLOUD=1 -DSDL=0 $(CFLAGS) $(FFTW3_OPTS)

gascloud:	fishstep_gascloud
	./fishstep_gascloud
	python3 tests/gascloud.py gascloud.dat


# cooling protostellar disk, the default setup without a window for about
# 2.5 outer orbits: with beta = 2 cooling it breaks up into spiral arms and
# gas clumps from t ~ 40, writes stardisk.dat and stardisk_rho.bin (maps)
fishstep_stardisk:	fishstep.c fishstep.h
	$(CC) -o $@ fishstep.c $(OPTS) -DSDL=0 -DT_END=64.0 $(CFLAGS) $(FFTW3_OPTS)

stardisk:	fishstep_stardisk
	./fishstep_stardisk
	python3 tests/stardisk.py stardisk.dat
