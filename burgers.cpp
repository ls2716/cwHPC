#include "model.h"
#include "burgers.h"
#include <algorithm>

using namespace std;

Burgers::Burgers(Model& m)
{
	cout << "Created Burgers" << endl;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	T = m.GetT();
	L = m.GetL();
	Nx = m.GetNx()-2;
	Ny = m.GetNy()-2;
	Nt = m.GetNt();
	dx = m.GetDx();
	dy = m.GetDy();
	dt = m.GetDt();
    ax = m.GetAx();
	ay = m.GetAy();
	b = m.GetB();
	c = m.GetC();
	Px = m.GetPx();
	Py = m.GetPy();
	P=Px*Py;
	
	axf=ax*dt/dx;
	ayf=ay*dt/dy;
	bxf=b*dt/dx;
	byf=b*dt/dy;
	cxf=c*dt/dx/dx;
	cyf=c*dt/dy/dy;
	
	//This could be slower
	cxf2=c*dt*2/dx/dx;
	cyf2=c*dt*2/dy/dy;
	
	int modX = Nx%Px;
	int modY = Ny%Py;
	my_grid_pos_x = my_rank%Px;
	my_grid_pos_y = my_rank/Px;
	my_Ny=Ny/Py;
	my_Nx=Nx/Px;
	my_grid_i0 = my_grid_pos_x*my_Nx+max(modX+my_grid_pos_x-Px,0);
	my_grid_j0 = my_grid_pos_y*my_Ny+max(modY+my_grid_pos_y-Py,0);
	if (my_grid_pos_y+modY>=Py)
		my_Ny++;
	if (my_grid_pos_x+modX>=Px)
		my_Nx++;
	my_grid_ie=my_grid_i0+my_Nx;
	my_grid_je=my_grid_j0+my_Ny;
	

	cout << "My rank: " << my_rank << " My pos x: " << my_grid_pos_x << " My pos y: " << my_grid_pos_y<< " My Nx: " << my_Nx<< " My Ny: " << my_Ny<< endl; 
	cout << "My rank: " << my_rank << " My pos i0: " << my_grid_i0 << " My pos j0: " << my_grid_j0<< " My pos ie: " << my_grid_ie << " My pos je: " << my_grid_je<<endl;
};

Burgers::~Burgers(){};



void Burgers::Initialize()
{
	
	int size = (my_Ny+2)*my_Nx;
	ugrid_o = new double[size];
	vgrid_o = new double[size];
	ugrid_e = new double[size];
	vgrid_e = new double[size];
	ugrid_l_vertB = new double[my_Ny];
	ugrid_r_vertB = new double[my_Ny];
	vgrid_l_vertB = new double[my_Ny];
	vgrid_r_vertB = new double[my_Ny];
	ugrid_myl_vertB = new double[my_Ny];
	ugrid_myr_vertB = new double[my_Ny];
	vgrid_myr_vertB = new double[my_Ny];
	vgrid_myl_vertB = new double[my_Ny];
	double x,y,r;
	double Lo2 = L/2;
	int j_offset;
	for (int j=0; j<my_Ny; j++)
	{
		r=(double)(my_grid_j0+j-1)*dy-Lo2;
		r*=r;
		j_offset=(j+1)*my_Nx;
		for (int i=0; i<my_Nx; i++)
		{	
			x=(double)(my_grid_i0+i+1)*dx-Lo2;
			r+=x*x;
			if (r>1)
			{
				ugrid_e[i+j_offset]=0.0;
				vgrid_e[i+j_offset]=0.0;
			}
			else
			{
				r = 2*pow(1-r,4)*(4*r+1);
				ugrid_e[i+j_offset]=r;
				vgrid_e[i+j_offset]=r;
			}			
		}
		ugrid_myr_vertB[j] = ugrid_e[my_Nx-1+j_offset];
		ugrid_myl_vertB[j] = ugrid_e[j_offset];
		vgrid_myr_vertB[j] = vgrid_e[my_Nx-1+j_offset];
		vgrid_myl_vertB[j] = vgrid_e[j_offset];
	}

//		for (int j=0; j<size_vertB; j++)
//		{
//			
//			ugrid_myr_vertB[j] = 0;
//			vgrid_myl_vertB[j] = 0;
//		}
//		j_offset=size-my_Nx;
//		for (int i=0; i<my_Nx; i++)
//		{
//			ugrid_e[i]=0.0;
//			vgrid_e[i]=0.0;
//			ugrid_e[i+size]=0.0;
//			vgrid_e[i+size]=0.0;
//		}
//	
	cout << "Grid successfully initialized" <<endl;
}

void Burgers::WriteToFile()
{
	cout << "Writing to a file" <<endl;
	int full_Nx=Nx+2;
	int full_Ny=Ny+2;
	ofstream file;
	file.open("grid.txt", ios::out | ios::trunc);
	cout << "Printing u" <<endl << endl;
	for (int j=full_Ny-1; j>=0; j--)
	{
		for (int i=0; i<(full_Nx-1); i++)
			file << full_ugrid[i+j*full_Nx]<<"	";
		file << full_ugrid[full_Nx-1+j*full_Nx]<<endl;
	}
	cout <<endl<<endl<< "Printing v" <<endl << endl;
	for (int j=full_Ny-1; j>=0; j--)
	{
		for (int i=0; i<(full_Nx-1); i++)
			file << full_vgrid[i+j*full_Nx]<<"	";
		file << full_vgrid[full_Nx-1+j*full_Nx]<<endl;
	}
	file.close();
}

void Burgers::Assemble()
{
	if (my_rank==0)
	{
		
		int full_Nx=Nx+2;
		int full_Ny=Ny+2;
		full_ugrid = new double[full_Nx*full_Ny];
		full_vgrid = new double[full_Nx*full_Ny];
	
		for (int j=0; j<full_Ny; j++)
		{
			full_ugrid[j*full_Nx] = 0;
			full_vgrid[j*full_Nx+full_Nx] = 0;
		}
		for (int i=0; i<full_Nx; i++)
		{
			full_ugrid[i] = 0;
			full_vgrid[i+(full_Ny-1)*full_Nx] = 0;
		}
		int m_Nx;
		int m_Ny;
		int m_i0;
		int m_j0;
		int m_size;
		double* m_ugrid = new double[(my_Nx+1)*(my_Ny+1)];
		double* m_vgrid = new double[(my_Nx+1)*(my_Ny+1)];
		for (int m=1; m<P; m++)
		{
			cout << "Receiving from: "<< m << endl;
			// Latter to single messages - down to two - maybe will be faster -it will be
			MPI_Recv(&m_Nx, 1, MPI_INT, m, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&m_Ny, 1, MPI_INT, m, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&m_i0, 1, MPI_INT, m, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&m_j0, 1, MPI_INT, m, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			m_size=m_Nx*m_Ny;
			cout << "Still okay" << endl;
			MPI_Recv(m_ugrid, m_size, MPI_DOUBLE, m, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(m_vgrid, m_size, MPI_DOUBLE, m, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			cout << "Received everything from: "<<m << " i0: "<<m_i0<< " j0 "<<m_j0<<endl;
			for (int j=0; j<m_Ny; j++)
			{
				for (int i=0; i<m_Nx; i++)
				{
					full_ugrid[(m_i0+i)+(m_j0)*full_Ny+j*m_Ny]=m_ugrid[i+j*m_Ny];
					full_vgrid[(m_i0+i)+(m_j0)*full_Ny+j*m_Ny]=m_vgrid[i+j*m_Ny];
					//cout << i+j*m_Ny << endl;
				}
			}
			cout << "Assembled from: "<<m <<endl;
		}
		for (int j=0; j<m_Ny; j++)
		{
			for (int i=0; i<m_Nx; i++)
			{
				full_ugrid[(1+i)+full_Ny+j*my_Ny]=ugrid_e[i+j*my_Ny];
				full_vgrid[(1+i)+full_Ny+j*my_Ny]=vgrid_e[i+j*my_Ny];
			}
		}
	}
	else //Sending
	{
		cout << "My rank: " << my_rank << " Sending. " << endl;
		MPI_Send(&my_Nx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&my_Ny, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&my_grid_i0, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		MPI_Send(&my_grid_j0, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		MPI_Send(ugrid_e, (my_Nx*my_Ny), MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
		MPI_Send(vgrid_e, (my_Nx*my_Ny), MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
		cout << "My rank: " << my_rank << " Sent. " << endl;
	}
}

//void Burgers::ddt(int is, int js, int ie, int je)
//{
//	//cout << "Calculating differences" <<setw(5)<<is<<setw(5)<<js<<setw(5)<<ie<<setw(5)<<je << endl;
//	for (int i=is; i<ie; i++)
//		for (int j=js; j<je; j++)
//		{
//			ugriddt[i+j*Nx]=dt*(c*((ugrid[i+1+j*Nx]-2*ugrid[i+j*Nx]+ugrid[i-1+j*Nx])/dx/dx+(ugrid[i+(j+1)*Nx]-2*ugrid[i+j*Nx]+ugrid[i+(j-1)*Nx])/dy/dy)
//			-(ax+b*ugrid[i+j*Nx])*(ugrid[i+j*Nx]-ugrid[i-1+j*Nx])/dx-(ay+b*vgrid[i+j*Nx])*(ugrid[i+j*Nx]-ugrid[i+(j-1)*Nx])/dy);
//			vgriddt[i+j*Nx]=dt*(c*((vgrid[i+1+j*Nx]-2*vgrid[i+j*Nx]+vgrid[i-1+j*Nx])/dx/dx+(vgrid[i+(j+1)*Nx]-2*vgrid[i+j*Nx]+vgrid[i+(j-1)*Nx])/dy/dy)
//			-(ax+b*ugrid[i+j*Nx])*(vgrid[i+j*Nx]-vgrid[i-1+j*Nx])/dx-(ay+b*vgrid[i+j*Nx])*(vgrid[i+j*Nx]-vgrid[i+(j-1)*Nx])/dy);
//		}
//}
//
//void Burgers::Update(int is, int js, int ie, int je)
//{
//	t+=dt;
//	cout << "Updating t=" << t<<endl;
//	for (int i=is; i<ie; i++)
//		for (int j=js; j<je; j++)
//		{
//			ugrid[i+j*Nx]+=ugriddt[i+j*Nx];
//			vgrid[i+j*Nx]+=vgriddt[i+j*Nx];
//			
//		}
//}

void Burgers::Integrate()
{
	cout << "Intergrating with time..." <<endl;
	while (t<T)
	{
		if ((T-t)<dt)
			dt=T-t;
//		ddt(1,1,Nx-1,Ny-1);
//		Update(1,1,Nx-1,Ny-1);
	}
	cout << "Finished." << endl;
}

void Burgers::Run()
{
	cout<< endl << "Running simulation. " << endl;
	Initialize();
//	WriteToFile("out_t_0.txt");
//	Integrate();
//	WriteToFile("out_t_1.txt");
	Assemble();
//	if (my_rank==0)
//		WriteToFile();
}

