default: all

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: model.cpp model.h
	mpicxx -std=c++11 -Wall -o model.o -c model.cpp

burgers.o: burgers.cpp burgers.h model.h avxFun.h
	mpicxx -std=c++11 -Wall -o burgers.o -c burgers.cpp

compile: main.o  model.o burgers.o avxFun.o
	mpicxx -o my_prog main.o  model.o burgers.o avxFun.o -O3 -ffast-math -funroll-loops  -march=native -ftree-vectorize
	

diff: compile
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 0 0 0 1 1 1

advx: compile
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 1 0 0 0 1 1

advy: compile
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 0 1 0 0 1 1

burg: compile
	mpiexec -np 1 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 1 1

diffp: compile
	mpiexec -np 2 my_prog 1.0 10.0 21 21 4000 0 0 0 1 2 1

advxp: compile
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 1 0 0 0 2 1

advyp: compile
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 0 1 0 0 2 1

burgp: compile
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 2 1
	
burgP: compile
	mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 4 4

.PHONY: clean # Specify that ’clean’ is not a real file
	target

clean:
	rm -f *.o my_prog   # Clean up (and ignore any errors)

all: burgP clean

diffP: compile
	mpiexec -np 36 my_prog 1.0 10.0 2001 2001 4000 0 0 0 1 6 6

burgcheck: compile
	mpiexec -np 2 my_prog 1.0 10.0 2001 2001 4000 1 0.5 1 0.02 1 2

diffcheck: compile
	mpiexec -np 2 my_prog 1.0 10.0 21 21 4000 0 0 0 1 2 1


#Below are not used so much

.PHONY: cleaner
	target
	
cleaner:
	rm -r *.er

padvx: compile
	collect -o initial3.er mpiexec -np 16 my_prog 1.0 10.0 2001 2001 4000 1 0 0 0 4 4

anadvx: padvx 
	analyzer initial3.er

ANadvx: anadvx cleaner

avxFun.o: avxFun.cpp avxFun.h
	mpicxx -mavx -mavx2 -mfma -std=c++11 -Wall -o avxFun.o -c avxFun.cpp
	
avxMain.o: avxMain.cpp avxFun.h
	mpicxx -mavx -mavx2 -mfma -std=c++11 -Wall -o avxMain.o -c avxMain.cpp


avxComp: avxMain.o avxFun.o
	mpicxx -o avxProg avxMain.o avxFun.o

avx: avxComp
	mpiexec -np 1 avxProg


