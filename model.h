#ifndef CLASS_MODEL
#define CLASS_MODEL
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>
#include <mpi.h>

class Model {
     public:
		Model(int argc, char* argv[]);
        ~Model();

        

        //Getters
        double GetL()      const { return L; }
        double GetT()      const { return T; }
        int    GetNx()     const { return Nx; }
        int    GetNy()     const { return Ny; }
        int    GetNt()     const { return Nt; }
        double GetDx()     const { return dx; }
        double GetDy()     const { return dy; }
        double GetDt()     const { return dt; }
        double GetAx()     const { return ax; }
        double GetAy()     const { return ay; }
        double GetB()      const { return b; }
        double GetC()      const { return c; }
		int    GetPx()     const { return Px; }
        int    GetPy()     const { return Py; }



	private:

		//Parsing the parameters from program arguments
		void ParseParameters(int argc, char* argv[]);

//		Check if parameters are valid
//        void ValidateParameters();
		
		
		//Parameter print
        void PrintParameters();
        //Validity checks
        void IsValidVal();
		void IsValidInp(int argc, char* argv[]);
		
        //Filing rest of the parameters
        void ParameterFill();
		
		

        //Numerics
        double L;
        double T;
        int    Nx;
        int    Ny;
        int    Nt;
        double dx;
        double dy;
        double dt;

        //Numerics for MPI
		int Px;
		int Py;
		int my_rank;
		int P;



        //Physics
        double ax;
        double ay;
        double b;
        double c;
		double axf;
		double ayf;
		double bxf;
		double byf;
		double cxf;
		double cyf;
		double cxf2;
		double cyf2;


};

#endif
