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
		int Px;
		int Py;
		int my_rank;
		int P;
		int my_Nx;
		int my_Ny;
		int my_grid_pos_x;
		int my_grid_pos_y;
		

        //Physics
        double ax;
        double ay;
        double b;
        double c;
		double axf;
		double ayf;
		double bxf;
		double byf;
		double cxf;
		double cyf;
		double cxf2;
		double cyf2;
		double t = 0;


};


#endif
