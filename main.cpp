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

//	Creating simulation class
	Burgers b(m);

//	Initialising the simulation
	b.Initialize();

//	Starting the timer
	typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

//	Integrating with time
	b.Integrate();

//	Stopping the timer and printing the duration
	hrc::time_point end = hrc::now();
	chrono::duration<double> elapsed_seconds = end-start;

	cout << "Time elapsed: "<<elapsed_seconds.count()<<"s for normal boundary"<<endl;
	b.WrapUp();
//	b.~Burgers();
	cout<<"Got here"<<endl;
	Burgers b2(m);
	b2.Bon=true;
	//Initialising the simulation
	b2.Initialize();

	//Starting the timer
	
    typedef std::chrono::milliseconds ms;
    start = hrc::now();

	//Integrating with time
	b2.Integrate();

	//Stopping the timer and printing the duration
	end = hrc::now();
	chrono::duration<double> elapsed_seconds2 = end-start;

	cout << "Time elapsed: "<<elapsed_seconds2.count()<<"s for streamlined boundary"<<endl;
//	
//	cout<<"Got here"<<endl;
	//Wrapping up - Getting the energy and printing it to file "grid.txt"
	b2.WrapUp();
//	b2.~Burgers();
	MPI_Barrier(MPI_COMM_WORLD);
	cout << "Got here2" <<endl;
    //Finalizing MPI
	MPI_Finalize();
    return 0;
}
