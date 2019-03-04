#ifndef CLASS_MODEL
#define CLASS_MODEL
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include mpi.h

class Model {
     public:
		Model(int argc, char* argv[]);
        ~Model();

        void PrintParameters();


        void IsValid();

        //Getters
        bool   IsVerbose() const { return verbose; }
        bool   IsHelp()    const { return help; }
//        double GetX0()     const { return x0; }
//        double GetY0()     const { return y0; }
//        double GetLx()     const { return Lx; }
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
		double GetSmall()  const { return small; }
		double GetQuick()  const { return small; }
		std::string GetOutname() const { return outname; }
        //Add any other getters here...



	private:
	
		//outputname
		std::string outname="out.txt";
        //Input reader
        bool readInputFile(std::string filename);
        void readInputCmd();
		void readTest(char testname);
		
		void ParseParameters(int argc, char* argv[]);
		
		//Check if parameters are valid
        void ValidateParameters();
        void ParameterFill();

        bool verbose=false;
        bool help=false;
		bool small=false;
		bool quick=false;

        //Numerics
        double L;
        double T;
        int    Nx;
        int    Ny;
        int    Nt;
        double dx;
        double dy;
        double dt;
		

        //Physics
        double ax;
        double ay;
        double b;
        double c;

//        //Bools to check if initialised
//        bool isL=false;
//        bool isT=false;
//        bool isax=false;
//        bool isay=false;
//        bool isb=false;
//        bool isc=false;



        //Add any additional parameters here...
};

#endif
