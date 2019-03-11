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
	
//	double* ures = b.GetResU();
//	//cout << ures[995+1100*m.GetNx()] << endl;
//	b.PrintGrid();
//	//b.WriteToFile("out_t_1.txt");
	cout << "Almost end"<<endl;
	MPI_Finalize();
    return 0;
}
