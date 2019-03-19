#include <chrono>
#include <iostream>
#include "model.h"
#include "burgers.h"

using namespace std;

int main(int argc, char* argv[])
{
    //Initialising MPI
    int err = MPI_Init(&argc, &argv);
    if(err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI." << endl;
        return -1;
    }

	//Creating Model to parse parameters
    Model m(argc, argv);

//	Creating simulation class
	Burgers b(m);

//	Initialising the simulation
	b.Initialize();

//	Starting the timer
	typedef std::chrono::high_resolution_clock hrc;
    hrc::time_point start = hrc::now();

//	Integrating with time
	b.Integrate();

//	Stopping the timer and printing the duration
	hrc::time_point end = hrc::now();
	chrono::duration<double> elapsed_seconds = end-start;
//	Printing time
	cout << "Time elapsed: "<<elapsed_seconds.count()<<"s"<<endl;
//	Wrapping up
	b.WrapUp();
	
	MPI_Finalize();
    return 0;
}
