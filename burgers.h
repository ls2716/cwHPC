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
	void WriteToFile();


private:

    //Pointers to the grids odd and even
    double* ugrid_o;
	double* vgrid_o;
	double* ugrid_e;
	double* vgrid_e;
	double* ugrid_l_vertB;
	double* ugrid_r_vertB;
	double* vgrid_l_vertB;
	double* vgrid_r_vertB;
	double* ugrid_myl_vertB;
	double* ugrid_myr_vertB;
	double* vgrid_myr_vertB;
	double* vgrid_myl_vertB;
	double* full_ugrid;
	double* full_vgrid;
	
	void Initialize();
	void Integrate();
	void ddt(int is, int js, int ie, int je);
	void Update(int is, int js, int ie, int je);

	void Assemble();
	
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
		int my_grid_i0;
		int my_grid_j0;
		int my_grid_ie;
		int my_grid_je;

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
