#include <chrono>
#include <iostream>
#include "model.h"
#include "burgers.h"

using namespace std;
// This is a testing script

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	
	//Creating Model
    Model m(argc, argv);
	
	//Creating simulation
	Burgers b(m);
	
	//Initialising the simulation
	b.Initialize();
	
	//Starting the timer
	typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
	
	//Intergrating with time
	b.Integrate();

	//Stopping the timer and printing the duration
	hrc::time_point end = hrc::now();
	chrono::duration<double> elapsed_seconds = end-start;
	
	cout << "Time elapsed: "<<elapsed_seconds.count()<<"s"<<endl;
	
	//Wrapping up - Getting the energy and printing it
	b.WrapUp();
	
//	cout << "Finished"<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;
}