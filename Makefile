default: all

# compiling files

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: model.cpp model.h
	mpicxx -std=c++11 -Wall -o model.o -c model.cpp

burgers.o: burgers.cpp burgers.h model.h
	mpicxx -std=c++11 -Wall -o burgers.o -c burgers.cpp

compile: main.o  model.o burgers.o
	mpicxx -o my_prog main.o  model.o burgers.o -O3 -ffast-math -funroll-loops  -march=native -ftree-vectorize


# Format for running
# mpiexec -np <number of processes> my_prog <T> <L> <Nx> <Ny> <Nt> <ax> <ay> <b> <c> <Px> <Py>

diff: compile cleanO # serial diffusion
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 0 0 0 1 1 1

advx: compile cleanO # serial x advection
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 1 0 0 0 1 1

advy: compile cleanO # serial y advection
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 0 1 0 0 1 1

burg: compile cleanO # serial burgers
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 1 1

diffp: compile cleanO # px=2 py=1 parallel diffusion
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 0 0 0 1 2 1

advxp: compile cleanO # px=2 py=1 parallel x advection
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 1 0 0 0 2 1

advyp: compile cleanO # px=2 py=1 parallel y advection
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 0 1 0 0 2 1

burgp: compile cleanO # px=2 py=1 parallel burgers
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 2 1
	
burgP: compile cleanO # px=4 py=4 parallel burgers
	mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 4 4

.PHONY: cleanO # Specify that ’clean’ is not a real file
	target

.PHONY: cleanProg # Specify that ’clean’ is not a real file
	target

cleanO:
	rm -f *.o # Clean up object files (and ignore any errors)

cleanProg:
	rm -f *.o # Clean up program (and ignore any errors)

all: diff advx advy burg diffp advxp advyp burgp burgP cleanO cleanProg
