#include <chrono>
#include <iostream>
#include "model.h"
#include "burgers.h"

using namespace std;
// This is a testing script

int main(int argc, char* argv[])
{

    Model m(argc, argv);
    m.PrintParameters();
    Burgers b(m);
	b.Run();
	double* ures = b.GetResU();
	//cout << ures[995+1100*m.GetNx()] << endl;
	b.PrintGrid();
	//b.WriteToFile("out_t_1.txt");
    return 0;
}
