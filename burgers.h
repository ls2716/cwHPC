#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include<cmath>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>

class Burgers {

public:
    Burgers(Model& m);
    ~Burgers();

    void Run();
    double* GetResU();
	double* GetResV();
	void PrintGrid();
	void WriteToFile(std::string outname);


private:

	//outputname
	std::string outname;
    //Pointers to the grid (results)
    double* ugrid;
	double* vgrid;
	//Pointers to update value
	double* ugriddt;
	double* vgriddt;
	void Initialize();
	void Integrate();
	void ddt(int is, int js, int ie, int je);
	void Update(int is, int js, int ie, int je);
	
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
	double t=0;

    //Physics
    double ax;
    double ay;
    double b;
    double c;

	//bools
	bool small;
	bool quick;


};


#endif
