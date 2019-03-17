#include <immintrin.h>
#include <stdio.h>
//#include <mpi.h>
#include <chrono>
#include <iostream>
#include "avxFun.h"

using namespace std;

void avxFun::calculateSing(double* uij, double *uipj, double* uinj, double* uijp, double* uijn)
{
	__m256d vec1, vec2, vec3, vec4, vec5,vecW, veCij,veCipj,veCinj,veCijp,veCijn,veCbx,veCby,vecS,vecS1,vecS2;
	vec1 = _mm256_loadu_pd(uij);
	vec2 = _mm256_loadu_pd(uipj);
	vec3 = _mm256_loadu_pd(uinj);
	vec4 = _mm256_loadu_pd(uijp);
	vec5 = _mm256_loadu_pd(uijn);
	veCij = _mm256_set1_pd(Cij);
	veCipj = _mm256_set1_pd(Cipj);
	veCinj = _mm256_set1_pd(Cinj);
	veCijp = _mm256_set1_pd(Cijp);
	veCijn = _mm256_set1_pd(Cijn);
	veCbx = _mm256_set1_pd(Cbx);
	veCby = _mm256_set1_pd(Cby);
//	ugrid_out[cor] = Cij*ugrid_in[cor]+Cinj*ugrid_in[cor-1]+Cipj*ugrid_in[cor+1]+Cijp*ugrid_in[cor+my_Nx]
//	+Cijn*ugrid_in[cor-my_Nx]+Cbx*ugrid_in[cor]*(ugrid_in[cor-1]-ugrid_in[cor])+Cby*vgrid_in[cor]*(ugrid_in[cor-my_Nx]-ugrid_in[cor]);
            
	vecS1 = _mm256_mul_pd(veCbx,_mm256_mul_pd(vec1,_mm256_sub_pd(vec3,vec1)));
	vecS2 = _mm256_mul_pd(veCby,_mm256_mul_pd(vec1,_mm256_sub_pd(vec5,vec1)));
	vecS = _mm256_add_pd(_mm256_mul_pd(veCbx,_mm256_mul_pd(vec1,_mm256_sub_pd(vec3,vec1))),_mm256_mul_pd(veCby,_mm256_mul_pd(vec1,_mm256_sub_pd(vec5,vec1))));
	//Below is the formula that works 
	vecW = _mm256_fmadd_pd(veCij,vec1,_mm256_fmadd_pd(veCinj,vec3,_mm256_fmadd_pd(veCipj,vec2,_mm256_fmadd_pd(veCijp,vec4,_mm256_fmadd_pd(veCijn,vec5,_mm256_add_pd(_mm256_mul_pd(veCbx,_mm256_mul_pd(vec1,_mm256_sub_pd(vec3,vec1))),_mm256_mul_pd(veCby,_mm256_mul_pd(vec1,_mm256_sub_pd(vec5,vec1)))))))));
	
	double* f = (double*)&vecW;
	cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<" "<<endl;
//	return vec1;
}
