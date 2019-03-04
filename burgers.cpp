#include "model.h"
#include "burgers.h"

using namespace std;

Burgers::Burgers(Model& m)
{
	T = m.GetT();
	L = m.GetL();
	Nx = m.GetNx();
	Ny = m.GetNy();
	Nt = m.GetNt();
	dx = m.GetDx();
	dy = m.GetDy();
	dt = m.GetDt();
    ax = m.GetAx();
	ay = m.GetAy();
	b = m.GetB();
	c = m.GetC();
	small = m.GetSmall();
};

Burgers::~Burgers(){};



void Burgers::Initialize()
{
	ugrid = new double[Nx*Ny];
	vgrid = new double[Nx*Ny];
	double x,y,r;
	double Lo2 = L/2;
	for (int j=0; j<Ny; j++)
		for (int i=0; i<Nx; i++)
		{	
			x=(double)i*dx-Lo2;
			y=(double)j*dy-Lo2;
			r=x*x+y*y;
			if (r>1)
			{
				ugrid[i+j*Nx]=0.0;
				vgrid[i+j*Nx]=0.0;
			}
			else
			{
				r = 2*pow(1-r,4)*(4*r+1);
				ugrid[i+j*Nx]=r;
				vgrid[i+j*Nx]=r;
			}
			
		}
	if (Lo2<1)
	{
		for (int j=0; j<Ny; j++)
		{
			ugrid[0+j*Nx]=0.0;
			vgrid[0+j*Nx]=0.0;
			ugrid[Nx-1+j*Nx]=0.0;
			vgrid[Nx-1+j*Nx]=0.0;
		}
		for (int i=0; i<Nx; i++)
		{
			ugrid[i]=0.0;
			vgrid[i]=0.0;
			ugrid[i+(Ny-1)*Nx]=0.0;
			vgrid[i+(Ny-1)*Nx]=0.0;
		}
	}
	
	cout << "Grid successfully initialized" <<endl;
}

void Burgers::Integrate()
{
	cout << "This doesn't do anything" <<endl;
}

void Burgers::Run()
{
	cout<< endl << "Running simulation" << endl;
	Initialize();
	Integrate();
}

double* Burgers::GetResU()
{
	return ugrid;
}

double* Burgers::GetResV()
{
	return vgrid;
}

void Burgers::PrintGrid()
{
	if (small)
	{
		cout << "Printing u" <<endl << endl;
		for (int j=Ny-1; j>=0; j--)
		{
			for (int i=0; i<Nx; i++)
				cout << setw(10) << ugrid[i+j*Nx];
			cout << endl;
		}
		cout << "Printing v" <<endl << endl;
		for (int j=Ny-1; j>=0; j--)
		{
			for (int i=0; i<Nx; i++)
				cout << setw(10) << vgrid[i+j*Nx];
			cout << endl;
		}
	}
	else
		cout << "U Mad? Grid is too big."<<endl;
}