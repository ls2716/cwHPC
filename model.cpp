#include "model.h"


using namespace std;

//Function to check validity
int Model::IsValid()
{
    int def=0;
    //Time cannot be negative or too big
    if (T<0)
    {
        cout << "Final time cannot be negative." << endl;
        return 1;
    }
    if (dt>0.02)
    {
        cout << "Final time is too big for accuracy." << endl;
        def = 2;
    }
    if (L<0)
    {
        cout << "Domain length cannot be negative." << endl;
        return 1;
    }
    if (dx>0.02)
    {
        cout << "Domain is too big for accuracy." << endl;
        def = 2;
    }
    if ((b<0)||(b>1.0))
    {
        cout << "Value of b is strange." << endl;
        def = 2;
    }
    if (c<0)
    {
        cout << "c in negative - negative diffusion - instabilities!." << endl;
        def = 2;
    }
    return def;
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
        cout << "Check if " << filename << "exists." << endl;
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


Model::Model(int argc, char* argv[])
{
    // Looking for argument defining input parameters data
    bool comm = true;
    //Creating default filename and string to hold the arguments
    string filename="input.txt";
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
                        comm = false;
                        break;
            //Place for other options
        }
    }
    //Reading parameters
    bool readSuccess;
    //From a file if overriden with flag
    if (!comm)
        readSuccess = readInputFile(filename);
    else
    {
        //From the command line
        readInputCmd();
        readSuccess = true;
    }
    if (readSuccess)
        cout << "Parameters were read successfully." << endl << endl;
    //Filling the rest of the parameters
    cout << "Filling rest of the parameters" << endl << endl;
    ParameterFill();

    //Checking if valid and asking if problems
    int info = IsValid();
    switch (info)
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

};

Model::~Model() {};
