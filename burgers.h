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
	//void PrintGrid();
	void WriteToFile(); // Should be put as private


private:
    //Pointers to the grids odd and even
    /* The integration uses two grids to decrease the number of instructions.
    / Using calculation of dudt the integration would have to Go through the whole array twice.
    / Now the update is done in the same step as calculation of dudt.
    */
    double* ugrid_o;
	double* vgrid_o;
	double* ugrid_e;
	double* vgrid_e;
	double* ugrid_t_horB;
    double* vgrid_t_horB;
    double* ugrid_b_horB;
    double* vgrid_b_horB;
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
	//Pointer for calculations
	double* ugrid_in;
	double* vgrid_in;
	double* ugrid_out;
	double* vgrid_out;
	double* urob_pointer;
	double* vrob_pointer;
	double* uout_pointer;
	double* vout_pointer;

    //Methods for performing simulation
	void Initialize();
	void Integrate();
	void NextStep();
    void BoundaryUpdate();
    void CalculateMyBoundaries();
    void CalculateCorners();
    void CalculateCenter();


//	void ddt(int is, int js, int ie, int je);
//	void Update(int is, int js, int ie, int je);

    // Assembly of full grid
	void Assemble();

    //Same parameters as model has
    //Numerics global
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
    //MPI numerics characterising the domain of the single process
    int my_Nx; // Size of the process domain in x
    int my_Ny; // Size of the process domain in y
    int smy_Nx;
    int smy_Ny;
    int my_size;
    int my_grid_pos_x; // x Position on the proceses gird i.e. <Px;
    int my_grid_pos_y; // y Position on the proceses gird i.e. <Py;
    int my_grid_i0; // Global i coordinate of the left bottom corner of the sub domain
    int my_grid_j0; // Global j coordinate of the left bottom corner of the sub domain
    int my_grid_ie; // Global i coordinate of the right top corner of the sub domain
    int my_grid_je; // Global i coordinate of the right top corner of the sub domain
    int j_offset;
    int cor;

    //Boundary type booleans
    bool upB=false;
    bool downB=false;
    bool rightB=false;
    bool leftB=false;


    MPI_Status   status[8];
    MPI_Request send_request[8],recv_request[8];
    int          ierr;
    int num_B=0;;
    int counter;

    //Physics
    double ax;
    double ay;
    double b;
    double c;
    //Faster calculations coefficients
    double Cij;
    double Cinj;
    double Cipj;
    double Cijn;
    double Cijp;
    double Cbx;
    double Cby;

    int nt; //current time step

};


#endif
