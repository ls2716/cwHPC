#include <chrono>
#include <iostream>
#include "../model.h"

using namespace std;
// This is a testing script

int main(int argc, char* argv[])
{

    cout << argc << endl;
    cout << argv[0] << endl;
    Model m(argc, argv);
    double x=5.0;
    m.SetX0(x);
    cout << "x0:  " << m.GetX0()<< endl;
    return 0;
}
