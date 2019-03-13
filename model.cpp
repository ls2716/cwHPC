#include "model.h"

using namespace std;

//Function to check validity
void Model::IsValidVal()
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

//Check if valid input - a float or an integer where applicable
void Model::IsValidInp(int argc, char* argv[])
{
	if (argc!=12)
	{
		cout<<"Not enough or too many parameters!"<<endl;
		exit(EXIT_FAILURE);
	}
	//Regular expressions
	regex integer("(\\+|-)?[[:digit:]]+");
	regex floating("((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?))?");
	
	if (!regex_match(argv[1],floating))
	{
		cout<<"Argument 1 (T) not a float."<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[2],floating))
	{
		cout<<"Argument 2 (L) not a float."<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[3],integer))
	{
		cout<<"Argument 3 (Nx) not an integer"<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[4],integer))
	{
		cout<<"Argument 4 (Ny) not an integer"<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[5],integer))
	{
		cout<<"Argument 5 (Nt) not an integer"<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[6],floating))
	{
		cout<<"Argument 6 (ax) not a float."<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[7],floating))
	{
		cout<<"Argument 7 (ay) not a float."<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[8],floating))
	{
		cout<<"Argument 8 (b) not a float."<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[9],floating))
	{
		cout<<"Argument 9 (c) not a float."<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[10],integer))
	{
		cout<<"Argument 10 (Px) not an integer"<<endl;
		exit(EXIT_FAILURE);
	}
	if (!regex_match(argv[5],integer))
	{
		cout<<"Argument 11 (Py) not an integer"<<endl;
		exit(EXIT_FAILURE);
	}
}


//Filling missing model parameters
void Model::ParameterFill()
{
	dt=T/Nt;
	dx=L/(Nx-1);
    dy=L/(Ny-1);
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
	Nx = strtol(argv[3],NULL,10);
	Ny = strtol(argv[4],NULL,10);
	Nt = strtol(argv[5],NULL,10);
	ax = strtod(argv[6],NULL);
	ay = strtod(argv[7],NULL);
	b = strtod(argv[8],NULL);
	c = strtod(argv[9],NULL);
	Px = strtol(argv[10],NULL,10);
	Py = strtol(argv[11],NULL,10);


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

	//Check if input is valid (actual integer or floating point number)
	if (my_rank==0)
		IsValidInp(argc, argv);
	MPI_Barrier(MPI_COMM_WORLD);
    

    cout << endl << "Starting execution. - My rank = "<<my_rank <<endl <<endl;

	//Reading parameters
	ParseParameters(argc, argv);

	//Filling the rest of the parameters
	cout << "My rank: "<<my_rank << " Filling rest of the parameters" << endl << endl;
	ParameterFill();

	//Checking if valid and asking if problems
	MPI_Barrier(MPI_COMM_WORLD);
	if (my_rank==0)
	{
		IsValidVal();
		PrintParameters();
	}
	MPI_Barrier(MPI_COMM_WORLD);
};

//Nothing to be done on the deletion
Model::~Model() {};
