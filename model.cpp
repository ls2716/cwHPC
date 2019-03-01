#include "model.h"
#include <iostream>

using namespace std;


void Model::SetX0(double& x)
{
    x0=x;

}

Model::Model(int argc, char* argv[]) {};
Model::~Model() {};

//int main()
//{
//    char* b;
//    Model modelT(0,0);
//    double x=5.0;
//    modelT.SetX0(x);
//    cout << modelT.GetX0() << endl;
//    return 0;
//}
