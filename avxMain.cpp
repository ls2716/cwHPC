//#include <mpi.h>
#include <chrono>
#include <iostream>
#include "avxFun.h"
#include <iomanip>



using namespace std;



int main()
{

	int Nx=10;
	int Ny=10;
	double Cij = 0.3;
	double Cinj = 0.2;
	double Cipj = 0.25;
	double Cijp = 0.25;
	double Cijn = 0.2;
	double* u = new double[Nx*Ny];
	for (int i=0; i<Nx; i++)
		for (int j=0; j<Ny; j++)
			u[j*Nx+i] = (double)i+ j*0.05;
	
	cout.precision(5);
	for (int j=Ny-1;j>=0;j--)
	{
		for (int i=0;i<Nx;i++)
			cout <<setw(8)<<fixed<<u[j*Nx+i];
		cout<<endl;
	}
	
	avxFun a;
	
	
	a.calculateSing(&u[2*Nx+4],&u[2*Nx+5],&u[2*Nx+3],&u[3*Nx+4],&u[1*Nx+4]);
	
	
	
	return 0;
}