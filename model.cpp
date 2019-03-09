#include "model.h"

using namespace std;

//Function to check validity
void Model::IsValid()
{
    int def=0;
    //Time cannot be negative or too big
    if (T<0)
    {
        cout << "Error. Final time cannot be negative." << endl;
        def = 1;
    }
    if (dt>0.02)
    {
        cout << "Warning. Final time is too big for accuracy." << endl;
        def = 2;
    }
    if (L<0)
    {
        cout << "Error. Domain length cannot be negative." << endl;
        def = 1;
    }
    if (dx>0.05)
    {
        cout << "Warning. Domain is too big for accuracy." << endl;
        def = 2;
    }
	if (dy>0.05)
    {
        cout << "Warning. Domain is too big for accuracy." << endl;
        def = 2;
    }
    if ((b<0)||(b>1.0))
    {
        cout << "Warning. Value of b is strange." << endl;
        def = 2;
    }
    if (c<0)
    {
        cout << "Warning. c in negative - negative diffusion - instabilities!." << endl;
        def = 2;
    }
	if (P!=(Px*Py))
	{
		cout << "The grid is wrongly subdivided" <<endl;
		def = 1;
	}

	switch (def)
    {
        case 0:     cout << "Model is valid." << endl;
                    break;
        case 1:     cout << "Model is invalid." << endl;
                    PrintParameters();
                    exit(EXIT_FAILURE);
                    break;
        case 2:     cout << "Warning. Parameters could be invalid." << endl;
                    PrintParameters();
                    cout << "Do you want to continue [y/n] ?" << endl;
                    char ans='n';
                    cin >> ans;
                    if (ans!='y')
                        exit(EXIT_FAILURE);
                    break;
    }
}

//Filling missing model parameters
void Model::ParameterFill()
{
	dt=T/Nt;
	dx=L/Nx;
    dy=L/Ny;
	cout << "My rank: "<< my_rank << " I have filled my parameters." << endl;
}

//Printing parameters function
void Model::PrintParameters()
{
    cout << "Printing parameters:" << endl;
    cout << "T = " << T << endl;
    cout << "L = " << L << endl;
    cout << "ax = " << ax << endl;
    cout << "ay = " << ay << endl;
    cout << "b = " << b << endl;
    cout << "c = " << c << endl;
    cout << "Nx = " << Nx << endl;
    cout << "Ny = " << Ny << endl;
    cout << "Nt = " << Nt << endl;
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "dt = " << dt << endl;
	cout << "P = " << P << endl;
	cout << "Px = " << Px << endl;
	cout << "Py = " << Py << endl;
}


void Model::ParseParameters(int argc, char* argv[])
{
	// Reading the data
	T = strtod(argv[1],NULL);
	L = strtod(argv[2],NULL);
	Nx = strtod(argv[3],NULL);
	Ny = strtod(argv[4],NULL);
	Nt = strtod(argv[5],NULL);
	ax = strtod(argv[6],NULL);
	ay = strtod(argv[7],NULL);
	b = strtod(argv[8],NULL);
	c = strtod(argv[9],NULL);
	Px = strtod(argv[10],NULL);
	Py = strtod(argv[11],NULL);


//    cout << "My rank: "<< my_rank << " I have all inputs." << " My pos x: " << my_grid_pos_x << " My pos y: " << my_grid_pos_y<< " My Nx: " << my_Nx<< " My Ny: " << my_Ny<< endl;

}


/* Constructor has to get parameters both of the simulation
/ as well as numerics for handling the sizes of each grid.
/ Each of the processes has to construct their own model and read parameters.
/ It is faster than broadcast.
*/
Model::Model(int argc, char* argv[])
{
    //Getting rank and number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);

    cout << endl << "Starting execution. - My rank = "<<my_rank <<endl <<endl;

	//Reading parameters
	ParseParameters(argc, argv);

	//Filling the rest of the parameters
	cout << "My rank: "<<my_rank << " Filling rest of the parameters" << endl << endl;
	ParameterFill();

	//Checking if valid and asking if problems
	if (my_rank==0)
	{
		IsValid();
		PrintParameters();
	}
	MPI_Barrier(MPI_COMM_WORLD);
};

//Nothing to be done on the deletion
Model::~Model() {};
