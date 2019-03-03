#include <chrono>
#include <iostream>
#include "../model.h"

using namespace std;
// This is a testing script

int main(int argc, char* argv[])
{


    Model m(argc, argv);
    m.PrintParameters();

    return 0;
}
