#ifndef CLASS_AVX
#define CLASS_AVX
#include <immintrin.h>

class avxFun{

	
public:

	avxFun(double& nCij, double& nCinj, double& nCipj, double& nCijn, double& nCijp, double& nCbx, double& nCby){
	 Cij = nCij;
	 Cinj = nCinj;
	 Cipj = nCipj;
	 Cijp = nCijp;
	 Cijn = nCijn;
	 Cbx = nCbx;
	 Cby = nCby;};
	~avxFun(){};
	
	void calculateSing(double* uij, double *uipj, double* uinj, double* uijp, double* uijn);
	void calculateMat(double* uij, double* uout,double* vij, double* vout, int& mNx, int mNy);
	
	
private:
	
	
	
	double Cij;
	double Cinj;
	double Cipj;
	double Cijp;
	double Cijn;
	double Cbx;
	double Cby;
	double* f;
	

	
	

};

























#endif