#include "model.h"
#include "burgers.h"

using namespace std;

Burgers::Burgers(Model& m)
{
    ax=m.GetAx();
    cout << ax << endl;
};

Burgers::~Burgers(){};
