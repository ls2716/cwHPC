#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include<cmath>
#include<iostream>

class Burgers {

public:
    Burgers(Model& m);
    ~Burgers();

    void Run();
    double* GetRes();


private:

    //Pointer to the grid (results)
    double* grid;


    //Same parameters as model has
    //Numerics
    double L;
    double T;
    int    Nx;
    int    Ny;
    int    Nt;
    double dx;
    double dy;
    double dt;

    //Physics
    double ax;
    double ay;
    double b;
    double c;



};


#endif
