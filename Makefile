CC = g++
CFLAGS = -fopenmp -I.
OPTFLAG = -O3
DEPS = fd_der.h interpolate.h tools.h functions.h BICGSTAB.h HelmholtzKernel.h frame.h PDE.h constants.h
OBJ = fd_der.o interpolate.o tools.o functions.o BICGSTAB.o HelmholtzKernel.o largeexample2.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(OPTFLAG) -c $< -o $@

largeexample2: $(OBJ)
	$(CC) $(CFLAGS) $(OPTFLAG) -o largeexample2 $(OBJ)

.PHONY:	clean

clean:
	rm *.o *~
