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

diff: compile
	mpiexec -np 3 my_prog 1.0 10.0 10 10 400 0 0 0 1 3 1

advx: compile
	mpiexec -np 6 my_prog 1.0 10.0 200 200 400 1 0 0 0 3 2
	
advy: compile
	mpiexec -np 6 my_prog 1.0 10.0 200 200 400 0 1 0 0 3 2

burg: compile
	mpiexec -np 1 my_prog 1.0 10.0 200 200 400 1 0.5 1 0.02 3 1

clean:
	-rm -f *.o my_prog# Clean up (and ignore any errors)

all: diff advx advy burg clean