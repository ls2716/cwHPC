#include <chrono>
#include <iostream>
#include "model.h"
#include "burgers.h"

using namespace std;

int main(int argc, char* argv[])
{
    //Initialising MPI
	MPI_Init(&argc, &argv);

	//Creating Model to parse parameters
    Model m(argc, argv);

	//Creating simulation class
	Burgers b(m);

	//Initialising the simulation
	b.Initialize();

	//Starting the timer
	typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

	//Integrating with time
	b.Integrate();

	//Stopping the timer and printing the duration
	hrc::time_point end = hrc::now();
	chrono::duration<double> elapsed_seconds = end-start;

	cout << "Time elapsed: "<<elapsed_seconds.count()<<"s"<<endl;

	//Wrapping up - Getting the energy and printing it to file "grid.txt"
	b.WrapUp();

    //Finalizing MPI
	MPI_Finalize();
    return 0;
}
