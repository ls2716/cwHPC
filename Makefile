default: all

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: model.cpp model.h
	mpicxx -std=c++11 -Wall -o model.o -c model.cpp

burgers.o: burgers.cpp burgers.h model.h
	mpicxx -std=c++11 -Wall -o burgers.o -c burgers.cpp

compile: main.o  model.o burgers.o
	mpicxx -o my_prog main.o  model.o burgers.o

.PHONY: clean # Specify that ’clean’ is not a real file
	target

.PHONY: cleaner
	target

diff: compile
	mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 0 0 0 1 4 4

advx: compile
	mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 1 0 0 0 4 4

advy: compile
	mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 0 1 0 0 4 4

burg: compile
	mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 4 4

clean:
	rm -f *.o my_prog   # Clean up (and ignore any errors)

cleaner:
	rm -r *.er

all: diff advx advy burg clean

padvx: compile
	collect -o initial.er mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 1 0 0 0 4 4

anadvx: padvx 
	analyzer initial.er

ANadvx: anadvx cleaner
	