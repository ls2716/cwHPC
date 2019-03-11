#include <chrono>
#include <iostream>
#include "model.h"
#include "burgers.h"

using namespace std;
// This is a testing script

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	
    Model m(argc, argv);
	
	//cout << "Still okay" << endl;
//	for (int i=0; i<argc; i++)
//		cout<<argv[i]<<endl;
    //m.PrintParameters();
	Burgers b(m);
	
//	cout<< "done"<<endl;
	b.Run();
//	b.PrintSubDomain(1);
//	b.PrintBound(1,'r');
	
//	double* ures = b.GetResU();
//	//cout << ures[995+1100*m.GetNx()] << endl;
//	b.PrintGrid();
//	//b.WriteToFile("out_t_1.txt");
//	b.~Burgers();
	cout << "Finished"<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;
}
