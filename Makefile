default: all

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: model.cpp model.h
	mpicxx -std=c++11 -Wall -o model.o -c model.cpp

burgers.o: burgers.cpp burgers.h model.h
	mpicxx -std=c++11 -Wall -o burgers.o -c burgers.cpp

myprog: main.o  model.o burgers.o
	mpicxx -o myprog main.o  model.o burgers.o

.PHONY: clean # Specify that ’clean’ is not a real file
	target
	
simple_diff: myprog
	mpiexec -np 6 myprog 1.0 10.0 200 200 400 0 0 0 1 3 2 

clean:
	-rm -f *.o # Clean up (and ignore any errors)

all: simple_diff clean