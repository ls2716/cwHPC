#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include<cmath>
#include<iostream>
#include<iomanip>

class Burgers {

public:
    Burgers(Model& m);
    ~Burgers();

    void Run();
    double* GetResU();
	double* GetResV();
	void PrintGrid();


private:

    //Pointer to the grid (results)
    double* ugrid;
	double* vgrid;
	void Initialize();
	void Integrate();

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

	//bools
	bool small;


};


#endif
