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
	quick = m.GetQuick();
};

Burgers::~Burgers(){};



void Burgers::Initialize()
{
	ugrid = new double[Nx*Ny];
	vgrid = new double[Nx*Ny];
	ugriddt = new double[Nx*Ny];
	vgriddt = new double[Nx*Ny];
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

void Burgers::WriteToFile(string outname)
{
	ofstream file;
	file.open(outname, ios::out | ios::trunc);
	cout << "Printing u" <<endl << endl;
	for (int j=Ny-1; j>=0; j--)
	{
		for (int i=0; i<(Nx-1); i++)
			file << ugrid[i+j*Nx]<<"	";
		file << ugrid[Nx-1+j*Nx]<<endl;
	}
	cout << "Printing v" <<endl << endl;
	for (int j=Ny-1; j>=0; j--)
	{
		for (int i=0; i<(Nx-1); i++)
			file << vgrid[i+j*Nx]<<"	";
		file << vgrid[Nx-1+j*Nx]<<endl;
	}
	file.close();
}


void Burgers::ddt(int is, int js, int ie, int je)
{
	//cout << "Calculating differences" <<setw(5)<<is<<setw(5)<<js<<setw(5)<<ie<<setw(5)<<je << endl;
	for (int i=is; i<ie; i++)
		for (int j=js; j<je; j++)
		{
			ugriddt[i+j*Nx]=dt*(c*((ugrid[i+1+j*Nx]-2*ugrid[i+j*Nx]+ugrid[i-1+j*Nx])/dx/dx+(ugrid[i+(j+1)*Nx]-2*ugrid[i+j*Nx]+ugrid[i+(j-1)*Nx])/dy/dy)
			-(ax+b*ugrid[i+j*Nx])*(ugrid[i+j*Nx]-ugrid[i-1+j*Nx])/dx-(ay+b*vgrid[i+j*Nx])*(ugrid[i+j*Nx]-ugrid[i+(j-1)*Nx])/dy);
			vgriddt[i+j*Nx]=dt*(c*((vgrid[i+1+j*Nx]-2*vgrid[i+j*Nx]+vgrid[i-1+j*Nx])/dx/dx+(vgrid[i+(j+1)*Nx]-2*vgrid[i+j*Nx]+vgrid[i+(j-1)*Nx])/dy/dy)
			-(ax+b*ugrid[i+j*Nx])*(vgrid[i+j*Nx]-vgrid[i-1+j*Nx])/dx-(ay+b*vgrid[i+j*Nx])*(vgrid[i+j*Nx]-vgrid[i+(j-1)*Nx])/dy);
		}
}

void Burgers::Update(int is, int js, int ie, int je)
{
	t+=dt;
	cout << "Updating t=" << t<<endl;
	for (int i=is; i<ie; i++)
		for (int j=js; j<je; j++)
		{
			ugrid[i+j*Nx]+=ugriddt[i+j*Nx];
			vgrid[i+j*Nx]+=vgriddt[i+j*Nx];
			
		if (small)
			PrintGrid();
		}
}

void Burgers::Integrate()
{
	cout << "Intergrating with time..." <<endl;
	while (t<T)
	{
		if ((T-t)<dt)
			dt=T-t;
		ddt(1,1,Nx-1,Ny-1);
		Update(1,1,Nx-1,Ny-1);
	}
	cout << "Finished." << endl;
}

void Burgers::Run()
{
	cout<< endl << "Running simulation - output filename: " << outname << endl;
	Initialize();
	WriteToFile("out_t_0.txt");
	Integrate();
	WriteToFile("out_t_1.txt");
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