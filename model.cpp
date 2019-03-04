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
    if (dx>0.02)
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
    Nx=2001;
    Ny=2001;
    Nt=4000;
    dx=L/Nx;
    dy=L/Ny;
    dt=T/Nt;
	if (small)
	{
		cout << "Small filling." <<endl;
		Nx=11;
		Ny=11;
		Nt=4000;
		dx=L/Nx;
		dy=L/Ny;
		dt=T/Nt;
	}
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
}

// Method for reading input parameters from a file
bool Model::readInputFile(string filename)
{
    cout << "Reading parameters from: " <<filename<< endl;
    ifstream inputfile;
    inputfile.open(filename);
    if (inputfile.fail())
    {
        cout << "Couldn't read the input file." << endl;
        cout << "Check if '" << filename << "' exists." << endl;
        exit(EXIT_FAILURE);
    }

    //Checking if appropriate number of lines
    string line;
    int lTotal=0;
    while(!inputfile.eof())
    {
        getline(inputfile, line);
        lTotal ++;
    }
    if (lTotal<6)
    {
        cout << "Not enough data. " <<endl;
        cout << "Input parameters should be written in separate lines in following order:" << endl;
        cout << "T, L, ax, ay, b, c" << endl;
        exit(EXIT_FAILURE);
    }

    //Going to the beginning of the file
    inputfile.clear();
    inputfile.seekg(0, ios::beg);

    //Reading parameters
    inputfile >> T;
    inputfile >> L;
    inputfile >> ax;
    inputfile >> ay;
    inputfile >> b;
    inputfile >> c;

    //Reading file
    return true;
}

//Function for reading parameters from the command line
void Model::readInputCmd()
{
    cout << "Reading parameters from the command line." << endl;
    cout << "Please input final time T." << endl;
    cin >> T;
    cout << "Please input domain size L." << endl;
    cin >> L;
    cout << "Please input parameter ax." << endl;
    cin >> ax;
    cout << "Please input parameter ay." << endl;
    cin >> ay;
    cout << "Please input parameter b." << endl;
    cin >> b;
    cout << "Please input parameter c." << endl;
    cin >> c;
}

void Model::readTest(char testname)
{
	switch (testname)
	{
		case '1':	T=1;
						L=10;
						ax=0;
						ay=0;
						b=0;
						c=1;
						cout << "Test 1 loaded" << endl;
						break;
		default:		cout << "Invalid test." << endl;
						exit(EXIT_FAILURE);
		
	}
}


void Model::ParseParameters(int argc, char* argv[])
{
	// Looking for argument defining input parameters data
    int inputType = 0;
	string filename;
	char test;
    //Creating default filename and string to hold the arguments
    string option;
    //Iterating through the arguments to get options
    for (int i=1; i<argc; i++)
    {
        option = argv[i];
        if ((option[0]!='-'))
        {
            cout << "Invalid flag: "<< option << endl;
            exit(EXIT_FAILURE);
        }
        switch (option[1])
        {
            case 'p':   filename=option.substr(2);
                        cout << "Reading parameters from a file." << endl;
                        cout << "Input filename set to: "<<filename<<endl;
                        inputType = 1;
                        break;
			case 't':   test = option[2];
                        cout << "Reading parameters from a test." << endl;
                        cout << "Test set to: test"<<test<<endl;
                        inputType = 2;
                        break;
			case 's':   cout << "Small grid testing." << endl;
                        small = true;
                        break;
            //Place for other options
        }
    }
	
    //Reading parameters
    bool readSuccess;
    switch (inputType)
	{
		case 0:	readInputCmd();
				readSuccess = true;
				break;
		case 1:	readSuccess = readInputFile(filename);
				break;
		case 2:	readSuccess = true;
				readTest(test);
				break;
	}
	
	
    if (readSuccess)
        cout << "Parameters were read successfully." << endl << endl;

}

Model::Model(int argc, char* argv[])
{
    cout << endl << "Starting execution." <<endl <<endl;
    ParseParameters(argc, argv);
	
    //Filling the rest of the parameters
    cout << "Filling rest of the parameters" << endl << endl;
    ParameterFill();

    //Checking if valid and asking if problems
    IsValid();
    
};

Model::~Model() {};
