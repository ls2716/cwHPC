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

    return 0;
}
