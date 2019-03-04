default: all

main.o: main.cpp
	g++ -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: model.cpp model.h
	g++ -std=c++11 -Wall -o model.o -c model.cpp

burgers.o: burgers.cpp burgers.h model.h
	g++ -std=c++11 -Wall -o burgers.o -c burgers.cpp

myprog: main.o burgers.o model.o
	g++ -o myprog main.o burgers.o model.o

.PHONY: clean # Specify that ’clean’ is not a real file
	target

clean:
	-rm -f *.o # Clean up (and ignore any errors)

all: myprog clean