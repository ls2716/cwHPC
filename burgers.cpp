#include "model.h"
#include "burgers.h"
#include <algorithm>

using namespace std;


Burgers::Burgers(Model& m)
{
	cout << "Created Burgers" << endl;
	//Rank info
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //Getting parameters from the model
	T = m.GetT();
	L = m.GetL();
	Nx = m.GetNx()-2; // The subdomains do not span the boundary
	Ny = m.GetNy()-2; // The subdomains do not span the boundary
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
	nt=0; //Setting time to zero

    // Calculating coefficients for integration
    Cij = 1 - dt*(ax/dx + ay/dy +2*c*(1/dx/dx+1/dy/dy));
    Cinj = dt/dx*(ax+c/dx);
    Cijn = dt/dy*(ay+c/dy);
    Cipj = c*dt/dx/dx;
    Cijp = c*dt/dy/dy;
    Cbx = b*dt/dx;
    Cby = b*dt/dy;
	cout << " " << Cij << " " << Cinj << " " << Cipj << " " << Cijn << " " << Cijp << " " << Cbx << " " << Cby << endl;


    /* Inferring data about the subdomain grid
    / Partitioning so that if the grid is not divisible by Px or Py the subdomains span whole domain
    / and so that the grids get bigger from the top right corner of the domain.
    / i.e if Nx mod Px == 1 then all the top row subdomain will be bigger in x by one
    / if Nx mod Px == 2 then first two top subdomains will be bigger by one.
    */
	int modX = Nx%Px;
	int modY = Ny%Py;
	my_grid_pos_x = my_rank%Px;
	my_grid_pos_y = my_rank/Px;
	my_Ny=Ny/Py;
	my_Nx=Nx/Px;
	
	//Two lines below could be faster using if statements
	my_grid_i0 = my_grid_pos_x*my_Nx+max(modX+my_grid_pos_x-Px,0)+1;
	my_grid_j0 = my_grid_pos_y*my_Ny+max(modY+my_grid_pos_y-Py,0)+1;

	if (my_grid_pos_y+modY>=Py)
		my_Ny++;
	if (my_grid_pos_x+modX>=Px)
		my_Nx++;
	my_grid_ie=my_grid_i0+my_Nx;
	my_grid_je=my_grid_j0+my_Ny;
	
	smy_Nx=my_Nx-1;
	smy_Ny=my_Ny-1;
	//Setting boundary type
	if (my_grid_pos_x!=0)
	{
        leftB=true;
        num_B+=2;
    }
    if (my_grid_pos_x!=(Px-1))
    {
        rightB=true;
        num_B+=2;
    }
    if (my_grid_pos_y!=0)
    {
        downB=true;
        num_B+=2;
    }
    if (my_grid_pos_y!=(Py-1))
    {
        upB=true;
        num_B+=2;
    }
//    status = new MPI_Status[num_B];
//    send_request = new MPI_Request[num_B];
//    recv_request = new MPI_Request[num_B];


    //Printing the subdomain info
	cout << "My rank: " << my_rank << " My pos x: " << my_grid_pos_x << " My pos y: " << my_grid_pos_y<< " My Nx: " << my_Nx<< " My Ny: " << my_Ny<< endl
	<< "My rank: " << my_rank << " My pos i0: " << my_grid_i0 << " My pos j0: " << my_grid_j0<< " My pos ie: " << my_grid_ie << " My pos je: " << my_grid_je<<endl;
};

//Need to release the memory
Burgers::~Burgers()
{
  delete[] ugrid_o;
	delete[] vgrid_o;
	delete[] ugrid_e;
	delete[] vgrid_e;
	delete[] ugrid_t_horB;
	delete[] ugrid_b_horB;
	delete[] vgrid_t_horB;
	delete[] vgrid_b_horB;
	delete[] ugrid_l_vertB;
	delete[] ugrid_r_vertB;
	delete[] vgrid_l_vertB;
	delete[] vgrid_r_vertB;
	delete[] ugrid_myl_vertB;
	delete[] ugrid_myr_vertB;
	delete[] vgrid_myr_vertB;
	delete[] vgrid_myl_vertB;
//	delete[] full_ugrid;
//	delete[] full_vgrid;
//	delete[] ugrid_in;
//	delete[] ugrid_out;
//	delete[] vgrid_in;
//	delete[] vgrid_out;
//	delete[] urob_pointer;
//	delete[] vrob_pointer;
//	delete[] uout_pointer;
//	delete[] vout_pointer;

};


/* Grid initialisation
/  the grid for each process consists of the following for both u and v (same):
/  1. actual subdomain in the next my_Nx*my_Ny entries using row major
/  2a. left boundary in array *grid_l_vertB and right boundary in array *grid_r_vertB (my_Ny entries)
/  2b. top boundary in array *grid_t_vertB and bottom boundary in array *grid_b_vertB (my_Ny entries)
/  3. own right and left boundaries to send to neightbouring subdomains (my_Ny entries) - these will be copied from actual sub domain
/  All of them are initialised using the formula given in the assignment
*/
void Burgers::Initialize()
{

	my_size = my_Ny*my_Nx;
	ugrid_o = new double[my_size];
	vgrid_o = new double[my_size];
	ugrid_e = new double[my_size];
	vgrid_e = new double[my_size];
	ugrid_t_horB = new double[my_Nx];
    vgrid_t_horB = new double[my_Nx];
    ugrid_b_horB = new double[my_Nx];
    vgrid_b_horB = new double[my_Nx];
	ugrid_l_vertB = new double[my_Ny];
	ugrid_r_vertB = new double[my_Ny];
	vgrid_l_vertB = new double[my_Ny];
	vgrid_r_vertB = new double[my_Ny];
	ugrid_myl_vertB = new double[my_Ny];
	ugrid_myr_vertB = new double[my_Ny];
	vgrid_myr_vertB = new double[my_Ny];
	vgrid_myl_vertB = new double[my_Ny];
	vgrid_out = vgrid_e;
	ugrid_out = ugrid_e;

    // Initial value calculation
	double x,y,r; // numerics
	double Lo2 = L/2;
	for (int j=0; j<my_Ny; j++)
	{
		y=(double)(my_grid_j0+j)*dy-Lo2; // y position is position of the left bottom corner + j of the loop + 1 from boundary offset
		y*=y;
		j_offset=(j)*my_Nx; // Offset j to target the actual subdomain in the domain matrix (rows 1 to my_Ny inculsive)
		for (int i=0; i<my_Nx; i++)
		{
			x=(double)(my_grid_i0+i)*dx-Lo2;
			r=pow(y+x*x,0.5);
			if (r>1)
			{
				
				ugrid_e[i+j_offset]=0.0;
				vgrid_e[i+j_offset]=0.0;
			}
			else
			{
//				cout << "Condition "<< (i+my_grid_i0) <<" " <<(j+my_grid_j0)<<"  "<<((double)(my_grid_i0+i)*dx-Lo2) <<" "<<((double)(my_grid_j0+j)*dy-Lo2)<<endl;
				r = 2*pow(1-r,4)*(4*r+1);
				ugrid_e[i+j_offset]=r;
				vgrid_e[i+j_offset]=r;
//				cout<< "  "<<r <<endl;
			}
		}
		//Updating own boundaries of process 0 - have to be sent at the beginning of each time step
		ugrid_myr_vertB[j] = ugrid_e[smy_Nx+j_offset];
		ugrid_myl_vertB[j] = ugrid_e[j_offset];
		vgrid_myr_vertB[j] = vgrid_e[smy_Nx+j_offset];
		vgrid_myl_vertB[j] = vgrid_e[j_offset];
	}

	cout << "Grid successfully initialized" <<endl;
}

void Burgers::PrintSubDomain(int mrank)
{
	if (my_rank==mrank)
	{
		cout<<"Printing domain of "<<mrank<<endl;
		for (int j=smy_Ny; j>=0;j--)
		{
			for (int i=0;i<my_Nx;i++)
				cout<<" "<<ugrid_out[i+j*my_Nx];
			cout<<endl;
		}
		
	}
}


void Burgers::PrintBound(int mrank, char w)
{
	if (my_rank==mrank)
	{
		cout.precision(4);
		if (w=='r')
			for (int j=smy_Ny; j>=0; j--)
				cout<<" r "<<fixed<<ugrid_myr_vertB[j];
		cout<<endl;
	}
}

void Burgers::BoundaryUpdate()
{
    counter=0;
    if (downB)
    {
        ierr=MPI_Isend(&ugrid_in[0],my_Nx,MPI_DOUBLE, my_rank-Px,2,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(ugrid_b_horB,my_Nx,MPI_DOUBLE,my_rank-Px,0,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
        ierr=MPI_Isend(&vgrid_in[0],my_Nx,MPI_DOUBLE, my_rank-Px,6,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(vgrid_b_horB,my_Nx,MPI_DOUBLE,my_rank-Px,4,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
		
    }
//	MPI_Barrier(MPI_COMM_WORLD);
    if (upB)
    {
        ierr=MPI_Isend(&ugrid_in[smy_Ny*my_Nx],my_Nx,MPI_DOUBLE, my_rank+Px,0,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(ugrid_t_horB,my_Nx,MPI_DOUBLE,my_rank+Px,2,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
        ierr=MPI_Isend(&vgrid_in[smy_Ny*my_Nx],my_Nx,MPI_DOUBLE, my_rank+Px,4,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(vgrid_t_horB,my_Nx,MPI_DOUBLE,my_rank+Px,6,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
		
    }
	MPI_Barrier(MPI_COMM_WORLD);
    if (leftB)
    {
        ierr=MPI_Isend(ugrid_myl_vertB,my_Ny,MPI_DOUBLE, my_rank-1,1,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(ugrid_l_vertB,my_Ny,MPI_DOUBLE,my_rank-1,3,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
        ierr=MPI_Isend(vgrid_myl_vertB,my_Ny,MPI_DOUBLE, my_rank-1,5,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(vgrid_l_vertB,my_Ny,MPI_DOUBLE,my_rank-1,7,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
    }
    if (rightB)
    {
        ierr=MPI_Isend(ugrid_myr_vertB,my_Ny,MPI_DOUBLE, my_rank+1,3,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(ugrid_r_vertB,my_Ny,MPI_DOUBLE,my_rank+1,1,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
        ierr=MPI_Isend(vgrid_myr_vertB,my_Ny,MPI_DOUBLE, my_rank+1,7,MPI_COMM_WORLD,&send_request[counter]);
        ierr=MPI_Irecv(vgrid_r_vertB,my_Ny,MPI_DOUBLE,my_rank+1,5,MPI_COMM_WORLD,&recv_request[counter]);
        counter++;
    }
    ierr=MPI_Waitall(counter,send_request,status);
    ierr=MPI_Waitall(counter,recv_request,status);
//    cout << "Updated boundaries" << endl;
	//Send here
	MPI_Barrier(MPI_COMM_WORLD);
//	if (my_rank==2)
//	{
//		for (int j=smy_Ny; j>=0; j--)
//		{
//			for (int i=0;i<my_Nx; i++)
//				cout<<" "<<fixed<<vgrid_e[i+j*my_Nx];
//			cout<<endl;
//			
//		}
//	if (downB)
//	{
//		cout<<"My rank: " << my_rank << " Sent down" << endl;
//		printB(vgrid_in,my_Nx,'h');
//		cout<<"My rank: " << my_rank << " Received down" << endl;
//		printB(vgrid_b_horB,my_Nx,'h');
//	}
//	if (upB)
//	{
//		cout<<"My rank: " << my_rank << " Sent up" << endl;
//		printB(&vgrid_in[smy_Ny*my_Nx],my_Nx,'h');
//		cout<<"My rank: " << my_rank << " Received up" << endl;
//		printB(vgrid_t_horB,my_Nx,'h');
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	if (leftB)
//	{
//		cout<<"My rank: " << my_rank << " Sent left" << endl;
//		printB(vgrid_myl_vertB,my_Ny,'v');
//		cout<<"My rank: " << my_rank << " Received left" << endl;
//		printB(vgrid_l_vertB,my_Ny,'v');
//	}
//	if (rightB)
//	{
//		cout<<"My rank: " << my_rank << " Sent right" << endl;
//		printB(vgrid_myr_vertB,my_Ny,'v');
//		cout<<"My rank: " << my_rank << " Received right" << endl;
//		printB(vgrid_r_vertB,my_Ny,'v');
//	}
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
}

void Burgers::CalculateMyBoundaries()
{
    if (downB)
    {
        for (int i=1;i<(smy_Nx);i++)
        {
            ugrid_out[i] = Cij*ugrid_in[i]+Cinj*ugrid_in[i-1]+Cipj*ugrid_in[i+1]+Cijp*ugrid_in[i+my_Nx]+Cijn*ugrid_b_horB[i]+Cbx*ugrid_in[i]*(ugrid_in[i-1]-ugrid_in[i])+Cby*vgrid_in[i]*(ugrid_b_horB[i]-ugrid_in[i]);
            vgrid_out[i] = Cij*vgrid_in[i]+Cinj*vgrid_in[i-1]+Cipj*vgrid_in[i+1]+Cijp*vgrid_in[i+my_Nx]+Cijn*vgrid_b_horB[i]+Cbx*ugrid_in[i]*(vgrid_in[i-1]-vgrid_in[i])+Cby*vgrid_in[i]*(vgrid_b_horB[i]-vgrid_in[i]);
        }
    }
    else
    {
        for (int i=1;i<(smy_Nx);i++)
        {
            ugrid_out[i] = Cij*ugrid_in[i]+Cinj*ugrid_in[i-1]+Cipj*ugrid_in[i+1]+Cijp*ugrid_in[i+my_Nx]+Cbx*ugrid_in[i]*(ugrid_in[i-1]-ugrid_in[i])+Cby*vgrid_in[i]*(-ugrid_in[i]);
            vgrid_out[i] = Cij*vgrid_in[i]+Cinj*vgrid_in[i-1]+Cipj*vgrid_in[i+1]+Cijp*vgrid_in[i+my_Nx]+Cbx*ugrid_in[i]*(vgrid_in[i-1]-vgrid_in[i])+Cby*vgrid_in[i]*(-vgrid_in[i]);
        }
    }
    if (upB)
    {
        urob_pointer = &ugrid_in[my_Nx*(smy_Ny)];
        vrob_pointer = &vgrid_in[my_Nx*(smy_Ny)];
        uout_pointer = &ugrid_out[my_Nx*(smy_Ny)];
        vout_pointer = &vgrid_out[my_Nx*(smy_Ny)];
        for (int i=1;i<(smy_Nx);i++)
        {
            uout_pointer[i] = Cij*urob_pointer[i]+Cinj*urob_pointer[i-1]+Cipj*urob_pointer[i+1]+Cijp*ugrid_t_horB[i]+Cijn*urob_pointer[i-my_Nx]+Cbx*urob_pointer[i]*(urob_pointer[i-1]-urob_pointer[i])+Cby*vrob_pointer[i]*(urob_pointer[i-my_Nx]-urob_pointer[i]);
            vout_pointer[i] = Cij*vrob_pointer[i]+Cinj*vrob_pointer[i-1]+Cipj*vrob_pointer[i+1]+Cijp*vgrid_t_horB[i]+Cijn*vrob_pointer[i-my_Nx]+Cbx*urob_pointer[i]*(vrob_pointer[i-1]-vrob_pointer[i])+Cby*vrob_pointer[i]*(vrob_pointer[i-my_Nx]-vrob_pointer[i]);
        }
    }
    else
    {
        urob_pointer = &ugrid_in[my_Nx*(smy_Ny)];
        vrob_pointer = &vgrid_in[my_Nx*(smy_Ny)];
        uout_pointer = &ugrid_out[my_Nx*(smy_Ny)];
        vout_pointer = &vgrid_out[my_Nx*(smy_Ny)];
        for (int i=1;i<(smy_Nx);i++)
        {
            uout_pointer[i] = Cij*urob_pointer[i]+Cinj*urob_pointer[i-1]+Cipj*urob_pointer[i+1]+Cijn*urob_pointer[i-my_Nx]+Cbx*urob_pointer[i]*(urob_pointer[i-1]-urob_pointer[i])+Cby*vrob_pointer[i]*(urob_pointer[i-my_Nx]-urob_pointer[i]);
            vout_pointer[i] = Cij*vrob_pointer[i]+Cinj*vrob_pointer[i-1]+Cipj*vrob_pointer[i+1]+Cijn*vrob_pointer[i-my_Nx]+Cbx*urob_pointer[i]*(vrob_pointer[i-1]-vrob_pointer[i])+Cby*vrob_pointer[i]*(vrob_pointer[i-my_Nx]-vrob_pointer[i]);
        }
    }
    if (leftB)
    {
        for (int j=1;j<(smy_Ny);j++)
        {
            j_offset = j*my_Nx;
            ugrid_out[j_offset] = ugrid_myl_vertB[j] = Cij*ugrid_in[j_offset]+Cinj*ugrid_l_vertB[j]+Cipj*ugrid_in[j_offset+1]+Cijp*ugrid_in[j_offset+my_Nx]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_l_vertB[j]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
            vgrid_out[j_offset] = vgrid_myl_vertB[j] = Cij*vgrid_in[j_offset]+Cinj*vgrid_l_vertB[j]+Cipj*vgrid_in[j_offset+1]+Cijp*vgrid_in[j_offset+my_Nx]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_l_vertB[j]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
        }
    }
    else
    {
        for (int j=1;j<(smy_Ny);j++)
        {
            j_offset = j*my_Nx;
            ugrid_out[j_offset] = ugrid_myl_vertB[j] = Cij*ugrid_in[j_offset]+Cipj*ugrid_in[j_offset+1]+Cijp*ugrid_in[j_offset+my_Nx]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
            vgrid_out[j_offset] = vgrid_myl_vertB[j] = Cij*vgrid_in[j_offset]+Cipj*vgrid_in[j_offset+1]+Cijp*vgrid_in[j_offset+my_Nx]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
        }
    }
    if (rightB)
    {
        for (int j=1;j<smy_Ny;j++)
        {
            j_offset = j*my_Nx+smy_Nx;
            ugrid_out[j_offset] = ugrid_myr_vertB[j] = Cij*ugrid_in[j_offset]+Cinj*ugrid_in[j_offset-1]+Cipj*ugrid_r_vertB[j]+Cijp*ugrid_in[j_offset+my_Nx]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_in[j_offset-1]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
            vgrid_out[j_offset] = vgrid_myr_vertB[j] = Cij*vgrid_in[j_offset]+Cinj*vgrid_in[j_offset-1]+Cipj*vgrid_r_vertB[j]+Cijp*vgrid_in[j_offset+my_Nx]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_in[j_offset-1]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
        }
    }
    else
    {
        for (int j=1;j<smy_Ny;j++)
        {
            j_offset = j*my_Nx+smy_Nx;
            ugrid_out[j_offset] = ugrid_myr_vertB[j] = Cij*ugrid_in[j_offset]+Cinj*ugrid_in[j_offset-1]+Cijp*ugrid_in[j_offset+my_Nx]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_in[j_offset-1]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
            vgrid_out[j_offset] = vgrid_myr_vertB[j] = Cij*vgrid_in[j_offset]+Cinj*vgrid_in[j_offset-1]+Cijp*vgrid_in[j_offset+my_Nx]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_in[j_offset-1]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
        }
    }
}
//TODO - put this into calculation of boundaries
void Burgers::CalculateCorners()
{
//	cout<<"My rank: "<<my_rank<<" Calculating corners"<<endl;
    //Bottom left
    if (leftB&&downB)
    {
        ugrid_out[0]= ugrid_myl_vertB[0] = Cij*ugrid_in[0]+Cinj*ugrid_l_vertB[0]+Cipj*ugrid_in[1]+Cijp*ugrid_in[my_Nx]+Cijn*ugrid_b_horB[0]+Cbx*ugrid_in[0]*(ugrid_l_vertB[0]-ugrid_in[0])+Cby*vgrid_in[0]*(ugrid_b_horB[0]-ugrid_in[0]);
        vgrid_out[0]= vgrid_myl_vertB[0] = Cij*vgrid_in[0]+Cinj*vgrid_l_vertB[0]+Cipj*vgrid_in[1]+Cijp*vgrid_in[my_Nx]+Cijn*vgrid_b_horB[0]+Cbx*ugrid_in[0]*(vgrid_l_vertB[0]-vgrid_in[0])+Cby*vgrid_in[0]*(vgrid_b_horB[0]-vgrid_in[0]);
//		cout<<"My rank: "<<my_rank<<" Complex left bottom corner"<<endl;
	}
    else if (leftB)
    {
        ugrid_out[0]= ugrid_myl_vertB[0] = Cij*ugrid_in[0]+Cinj*ugrid_l_vertB[0]+Cipj*ugrid_in[1]+Cijp*ugrid_in[my_Nx]+Cbx*ugrid_in[0]*(ugrid_l_vertB[0]-ugrid_in[0])+Cby*vgrid_in[0]*(-ugrid_in[0]);
        vgrid_out[0]= vgrid_myl_vertB[0] = Cij*vgrid_in[0]+Cinj*vgrid_l_vertB[0]+Cipj*vgrid_in[1]+Cijp*vgrid_in[my_Nx]+Cbx*ugrid_in[0]*(vgrid_l_vertB[0]-vgrid_in[0])+Cby*vgrid_in[0]*(-vgrid_in[0]);
//		cout<<"My rank: "<<my_rank<<" Complex left (bottom)corner"<<endl;
	}
    else if (downB)
    {
        ugrid_out[0]= ugrid_myl_vertB[0] = Cij*ugrid_in[0]+Cipj*ugrid_in[1]+Cijp*ugrid_in[my_Nx]+Cijn*ugrid_b_horB[0]+Cbx*ugrid_in[0]*(-ugrid_in[0])+Cby*vgrid_in[0]*(ugrid_b_horB[0]-ugrid_in[0]);
        vgrid_out[0]= vgrid_myl_vertB[0] = Cij*vgrid_in[0]+Cipj*vgrid_in[1]+Cijp*vgrid_in[my_Nx]+Cijn*vgrid_b_horB[0]+Cbx*ugrid_in[0]*(-vgrid_in[0])+Cby*vgrid_in[0]*(vgrid_b_horB[0]-vgrid_in[0]);
//		cout<<"My rank: "<<my_rank<<" Complex bottom (left)corner"<<endl;
	}
    else
    {
        ugrid_out[0]= ugrid_myl_vertB[0] = Cij*ugrid_in[0]+Cipj*ugrid_in[1]+Cijp*ugrid_in[my_Nx]+Cbx*ugrid_in[0]*(-ugrid_in[0])+Cby*vgrid_in[0]*(-ugrid_in[0]);
        vgrid_out[0]= vgrid_myl_vertB[0] = Cij*vgrid_in[0]+Cipj*vgrid_in[1]+Cijp*vgrid_in[my_Nx]+Cbx*ugrid_in[0]*(-vgrid_in[0])+Cby*vgrid_in[0]*(-vgrid_in[0]);
//		cout<<"My rank: "<<my_rank<<" Grid bottom left corner"<<endl;
	}
    //Top left
    if (leftB&&upB)
    {
        j_offset=my_Nx*smy_Ny;
        ugrid_out[j_offset]= ugrid_myl_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cinj*ugrid_l_vertB[smy_Ny]+Cipj*ugrid_in[j_offset+1]+Cijp*ugrid_t_horB[0]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_l_vertB[smy_Ny]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myl_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cinj*vgrid_l_vertB[smy_Ny]+Cipj*vgrid_in[j_offset+1]+Cijp*vgrid_t_horB[0]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_l_vertB[smy_Ny]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Complex left top corner"<<endl;
	}
    else if (leftB)
    {
        j_offset=my_Nx*smy_Ny;
        ugrid_out[j_offset]= ugrid_myl_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cinj*ugrid_l_vertB[smy_Ny]+Cipj*ugrid_in[j_offset+1]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_l_vertB[smy_Ny]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myl_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cinj*vgrid_l_vertB[smy_Ny]+Cipj*vgrid_in[j_offset+1]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_l_vertB[smy_Ny]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Complex left (top) corner"<<endl;
	}
    else if (upB)
    {
        j_offset=my_Nx*smy_Ny;
        ugrid_out[j_offset]= ugrid_myl_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cipj*ugrid_in[j_offset+1]+Cijp*ugrid_t_horB[0]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myl_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cipj*vgrid_in[j_offset+1]+Cijp*vgrid_t_horB[0]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Complex (left) top corner"<<endl;
	}
    else
    {
        j_offset=my_Nx*smy_Ny;
        ugrid_out[j_offset]= ugrid_myl_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cipj*ugrid_in[j_offset+1]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myl_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cipj*vgrid_in[j_offset+1]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Grid left top corner"<<endl;
	}
    //Top right
    if (rightB&&upB)
    {
        j_offset = my_Nx*my_Ny-1;
        ugrid_out[j_offset]= ugrid_myr_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cinj*ugrid_in[j_offset-1]+Cipj*ugrid_r_vertB[smy_Ny]+Cijp*ugrid_t_horB[smy_Nx]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_in[j_offset-1]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myr_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cinj*vgrid_in[j_offset-1]+Cipj*vgrid_r_vertB[smy_Ny]+Cijp*vgrid_t_horB[smy_Nx]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_in[j_offset-1]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Complex right top corner"<<endl;
	}
    else if (rightB)
    {
        j_offset = my_Nx*my_Ny-1;
        ugrid_out[j_offset]= ugrid_myr_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cinj*ugrid_in[j_offset-1]+Cipj*ugrid_r_vertB[smy_Ny]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_in[j_offset-1]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myr_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cinj*vgrid_in[j_offset-1]+Cipj*vgrid_r_vertB[smy_Ny]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_in[j_offset-1]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Complex right (top) corner"<<endl;
	}
    else if (upB)
    {
        j_offset = my_Nx*my_Ny-1;
        ugrid_out[j_offset]= ugrid_myr_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cinj*ugrid_in[j_offset-1]+Cijp*ugrid_t_horB[smy_Nx]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_in[j_offset-1]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myr_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cinj*vgrid_in[j_offset-1]+Cijp*vgrid_t_horB[smy_Nx]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_in[j_offset-1]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Complex (right) top corner"<<endl;
	}
    else
    {
        j_offset = my_Nx*my_Ny-1;
        ugrid_out[j_offset]= ugrid_myr_vertB[smy_Ny] = Cij*ugrid_in[j_offset]+Cinj*ugrid_in[j_offset-1]+Cijn*ugrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(ugrid_in[j_offset-1]-ugrid_in[j_offset])+Cby*vgrid_in[j_offset]*(ugrid_in[j_offset-my_Nx]-ugrid_in[j_offset]);
        vgrid_out[j_offset]= vgrid_myr_vertB[smy_Ny] = Cij*vgrid_in[j_offset]+Cinj*vgrid_in[j_offset-1]+Cijn*vgrid_in[j_offset-my_Nx]+Cbx*ugrid_in[j_offset]*(vgrid_in[j_offset-1]-vgrid_in[j_offset])+Cby*vgrid_in[j_offset]*(vgrid_in[j_offset-my_Nx]-vgrid_in[j_offset]);
//		cout<<"My rank: "<<my_rank<<" Grid right top corner"<<endl;
	}
    //Bottom right
    if (downB&&rightB)
    {
        ugrid_out[smy_Nx]= ugrid_myr_vertB[0] = Cij*ugrid_in[smy_Nx]+Cinj*ugrid_in[smy_Nx-1]+Cipj*ugrid_r_vertB[0]+Cijp*ugrid_in[smy_Nx+my_Nx]+Cijn*ugrid_b_horB[smy_Nx]+Cbx*ugrid_in[smy_Nx]*(ugrid_in[smy_Nx-1]-ugrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(ugrid_b_horB[smy_Nx]-ugrid_in[smy_Nx]);
        vgrid_out[smy_Nx]= vgrid_myr_vertB[0] = Cij*vgrid_in[smy_Nx]+Cinj*vgrid_in[smy_Nx-1]+Cipj*vgrid_r_vertB[0]+Cijp*vgrid_in[smy_Nx+my_Nx]+Cijn*vgrid_b_horB[smy_Nx]+Cbx*ugrid_in[smy_Nx]*(vgrid_in[smy_Nx-1]-vgrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(vgrid_b_horB[smy_Nx]-vgrid_in[smy_Nx]);
//		cout<<"My rank: "<<my_rank<<" Complex right bot corner"<<endl;
	}
    else if (downB)
    {
        ugrid_out[smy_Nx]= ugrid_myr_vertB[0] = Cij*ugrid_in[smy_Nx]+Cinj*ugrid_in[smy_Nx-1]+Cijp*ugrid_in[smy_Nx+my_Nx]+Cijn*ugrid_b_horB[smy_Nx]+Cbx*ugrid_in[smy_Nx]*(ugrid_in[smy_Nx-1]-ugrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(ugrid_b_horB[smy_Nx]-ugrid_in[smy_Nx]);
        vgrid_out[smy_Nx]= vgrid_myr_vertB[0] = Cij*vgrid_in[smy_Nx]+Cinj*vgrid_in[smy_Nx-1]+Cijp*vgrid_in[smy_Nx+my_Nx]+Cijn*vgrid_b_horB[smy_Nx]+Cbx*ugrid_in[smy_Nx]*(vgrid_in[smy_Nx-1]-vgrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(vgrid_b_horB[smy_Nx]-vgrid_in[smy_Nx]);
//		cout<<"My rank: "<<my_rank<<" Complex (right) bot corner"<<endl;
	}
    else if (rightB)
    {
        ugrid_out[smy_Nx]= ugrid_myr_vertB[0] = Cij*ugrid_in[smy_Nx]+Cinj*ugrid_in[smy_Nx-1]+Cipj*ugrid_r_vertB[0]+Cijp*ugrid_in[smy_Nx+my_Nx]+Cbx*ugrid_in[smy_Nx]*(ugrid_in[smy_Nx-1]-ugrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(-ugrid_in[smy_Nx]);
        vgrid_out[smy_Nx]= vgrid_myr_vertB[0] = Cij*vgrid_in[smy_Nx]+Cinj*vgrid_in[smy_Nx-1]+Cipj*vgrid_r_vertB[0]+Cijp*vgrid_in[smy_Nx+my_Nx]+Cbx*ugrid_in[smy_Nx]*(vgrid_in[smy_Nx-1]-vgrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(-vgrid_in[smy_Nx]);
//		cout<<"My rank: "<<my_rank<<" Complex right (bot) corner"<<endl;
	}
    else
    {
        ugrid_out[smy_Nx]= ugrid_myr_vertB[0] = Cij*ugrid_in[smy_Nx]+Cinj*ugrid_in[smy_Nx-1]+Cijp*ugrid_in[smy_Nx+my_Nx]+Cbx*ugrid_in[smy_Nx]*(ugrid_in[smy_Nx-1]-ugrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(-ugrid_in[smy_Nx]);
        vgrid_out[smy_Nx]= vgrid_myr_vertB[0] = Cij*vgrid_in[smy_Nx]+Cinj*vgrid_in[smy_Nx-1]+Cijp*vgrid_in[smy_Nx+my_Nx]+Cbx*ugrid_in[smy_Nx]*(vgrid_in[smy_Nx-1]-vgrid_in[smy_Nx])+Cby*vgrid_in[smy_Nx]*(-vgrid_in[smy_Nx]);
//		cout<<"My rank: "<<my_rank<<" Grid right bot corner"<<endl;
	}
}

void Burgers::CalculateCenter()
{
    for (int j=1;j<smy_Ny;j++)
    {
        j_offset = j*my_Nx;
        for (int i=1;i<smy_Nx;i++)
        {
            cor = i+j_offset;
            ugrid_out[cor] = Cij*ugrid_in[cor]+Cinj*ugrid_in[cor-1]+Cipj*ugrid_in[cor+1]+Cijp*ugrid_in[cor+my_Nx]+Cijn*ugrid_in[cor-my_Nx]+Cbx*ugrid_in[cor]*(ugrid_in[cor-1]-ugrid_in[cor])+Cby*vgrid_in[cor]*(ugrid_in[cor-my_Nx]-ugrid_in[cor]);
            vgrid_out[cor] = Cij*vgrid_in[cor]+Cinj*vgrid_in[cor-1]+Cipj*vgrid_in[cor+1]+Cijp*vgrid_in[cor+my_Nx]+Cijn*vgrid_in[cor-my_Nx]+Cbx*ugrid_in[cor]*(vgrid_in[cor-1]-vgrid_in[cor])+Cby*vgrid_in[cor]*(vgrid_in[cor-my_Nx]-vgrid_in[cor]);
        }
    }
}

void Burgers::NextStep()
{
    nt++;
	if ((nt%10==0)&&(my_rank==0))
		cout<<"Step: "<<nt<<"/"<<Nt<<endl;
    //Changing pointers
    if ((nt%2)==0)
    {
        ugrid_in=ugrid_o;
        vgrid_in=vgrid_o;
        ugrid_out=ugrid_e;
        vgrid_out=vgrid_e;
    }
    else
    {
        ugrid_in=ugrid_e;
        vgrid_in=vgrid_e;
        ugrid_out=ugrid_o;
        vgrid_out=vgrid_o;
    }
    //Setting boundaries
    BoundaryUpdate();
    CalculateCorners();
    CalculateMyBoundaries();
    CalculateCenter();
	if ((nt%10==0))
		Energy();
}

void Burgers::Energy()
{
	MPI_Barrier(MPI_COMM_WORLD);
//	cout << "Calculating energy" << endl;
	energy=0;
	for (int j=0;j<my_Ny;j++)
	{
		j_offset=j*my_Nx;
		for (int i=0; i<my_Nx; i++)
			energy+=dx*dy*(pow(ugrid_out[i+j_offset],2)+pow(vgrid_out[i+j_offset],2));
	}
//	cout << "My rank : " <<my_rank<<" My energy "<<energy<<endl;
	if (my_rank==0)
	{
		double m_energy;
		for (int m=1; m<P; m++)
		{
//			cout<<"Receiving energy from " <<m<<endl;
			MPI_Recv(&m_energy, 1, MPI_DOUBLE, m, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//			cout <<"Got energy from "<<m<<endl;
			energy+=m_energy;
		}
		cout << "Total Energy:  "<<(energy/2.0)<<endl;
	}
	else
	{
//		cout << "My rank: "<<my_rank<<" I am sending energy"<<endl; 
		MPI_Send(&energy, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
		
	}
	MPI_Barrier(MPI_COMM_WORLD);
}


void Burgers::Integrate()
{
//	cout << "Intergrating with time..." <<endl;
//	while (nt<Nt)
	for (int k=0; k<2000;k++)
	{
        NextStep();
//		ddt(1,1,Nx-1,Ny-1);
//		Update(1,1,Nx-1,Ny-1);
		
	}
//	cout << "Finished." << endl;
}

void Burgers::Run()
{
	cout<< endl << "Running simulation. " << endl;
	Initialize();
//	WriteToFile("out_t_0.txt");
	Integrate();
//	WriteToFile("out_t_1.txt");
	
	Assemble();
	MPI_Barrier(MPI_COMM_WORLD);
	Energy();
	if (my_rank==0)
		WriteToFile();
	MPI_Barrier(MPI_COMM_WORLD);
}






// Below are ready








void Burgers::WriteToFile()
{
	cout << "Writing to a file" <<endl;
	int full_Nx=Nx+2;
	int full_Ny=Ny+2;
	ofstream file0;
	file0.open("grid.txt", ios::out | ios::trunc);
	cout.precision(3);
//	cout << "Printing u" <<endl << endl;
	file0<<full_Nx<<"	"<<full_Ny<<endl;
	for (int j=full_Ny-1; j>=0; j--)
	{
		for (int i=0; i<(full_Nx-1); i++)
			file0 << fixed<<full_ugrid[i+j*full_Nx]<<"	";
		file0 << full_ugrid[full_Nx-1+j*full_Nx]<<endl;
	}
//	cout << "Done printing u" << endl;
	file0<<full_Nx<<"	"<<full_Ny<<endl;
//	cout <<endl<<endl<< "Printing v" <<endl << endl;
	for (int j=full_Ny-1; j>=0; j--)
	{
		for (int i=0; i<(full_Nx-1); i++)
			file0 <<fixed<< full_vgrid[i+j*full_Nx]<<"	";
		file0 << full_vgrid[full_Nx-1+j*full_Nx]<<endl;
	}
	file0.close();
	cout << "Done printing v" << endl;
}

//Assembling full grid for final Energy calculation and writing to a file
void Burgers::Assemble()
{
	if (my_rank==0)
	{
		cout << "Starting assembly"<<endl;
        //Defining new grid and parameters
		int full_Nx=Nx+2;
		int full_Ny=Ny+2;
		full_ugrid = new double[full_Nx*full_Ny+2];
		full_vgrid = new double[full_Nx*full_Ny+2];

        //Setting boundary data - coudl be faster if do it during printing - without mentioning them here
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
		
		//Parameters of receiving processes
		int m_Nx;
		int m_Ny;
		int m_i0;
		int m_j0;
		int m_size;
		double* m_ugrid = new double[(my_Nx+1)*(my_Ny+1)+1];
		double* m_vgrid = new double[(my_Nx+1)*(my_Ny+1)+1];
//		double* m_ugrid;
//		double* m_vgrid;
//		cout<<"Okay"<<endl;
		//Receive data from each process - coulbd be faster by streamlining the ints
		for (int m=1; m<P; m++)
		{
//			cout << "Receiving from: "<< m << endl;
			// Latter to single messages - down to two - maybe will be faster -it will be
			MPI_Recv(&m_Nx, 1, MPI_INT, m, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&m_Ny, 1, MPI_INT, m, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&m_i0, 1, MPI_INT, m, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&m_j0, 1, MPI_INT, m, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			m_size=m_Nx*m_Ny;
//			cout << "Still okay" << endl;
			MPI_Recv(m_ugrid, m_size, MPI_DOUBLE, m, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(m_vgrid, m_size, MPI_DOUBLE, m, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//			cout << "Received everything from: "<<m << " i0: "<<m_i0<< " j0 "<<m_j0<<endl;
			//Assembling on the full grid
			for (int j=0; j<m_Ny; j++)
			{
				for (int i=0; i<m_Nx; i++)
				{
					full_ugrid[(m_i0+i)+(m_j0)*full_Ny+j*full_Ny]=m_ugrid[i+j*m_Nx];
					full_vgrid[(m_i0+i)+(m_j0)*full_Ny+j*full_Ny]=m_vgrid[i+j*m_Nx];
//					if (full_ugrid[(m_i0+i)+(m_j0)*full_Ny+j*full_Ny]>0)
//					cout << "Condition asm "<< (i+m_i0) <<" " <<(j+m_j0)<<"  "<<((double)(m_i0+i)*dx-L/2) <<" "<<((double)(m_j0+j)*dy-L/2)<<endl;
				}
			}
//			cout << "Assembled from: "<<m <<endl;
		}
//		cout<<"Still okay"<<endl;
		delete[] m_ugrid;
		delete[] m_vgrid;
		for (int j=0; j<my_Ny; j++)
		{
			for (int i=0; i<my_Nx; i++)
			{
				full_ugrid[(1+i)+full_Ny+j*full_Ny]=ugrid_out[i+j*my_Nx];
				full_vgrid[(1+i)+full_Ny+j*full_Ny]=vgrid_out[i+j*my_Nx];
//				if (ugrid_out[i+j*my_Ny]>0)
//					cout<<"Still ok "<<i<<" "<<j<<endl;
			}
		}
//		cout<<"Assembled mine"<<endl;
	}
	else //Sending
	{
//		cout << "My rank: " << my_rank << " Sending. " << endl;
		MPI_Send(&my_Nx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&my_Ny, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&my_grid_i0, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		MPI_Send(&my_grid_j0, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		MPI_Send(ugrid_out, (my_Nx*my_Ny), MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
		MPI_Send(vgrid_out, (my_Nx*my_Ny), MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
//		cout << "My rank: " << my_rank << " Sent. " << endl;
	}
}
