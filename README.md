# cwHPC
HPC coursework repository

The generated MPI programm is used to solve the Burgers equation in 2D.

The Makefile contains targets as specified in the assignment.

Program contains two classes:
1. Model class to read parametrs of the simulation.
2. Burgers class which divides the domain among the processes, does time integration, prints out the energy and save the field to a file.
3. In the main.cpp the simulation is performed using both classes and the integration time is measured and printed out.

Both the makefile and the code is thouroughly commented.