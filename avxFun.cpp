#include <immintrin.h>
#include <stdio.h>
//#include <mpi.h>
#include <chrono>
#include <iostream>
#include "avxFun.h"

using namespace std;

void avxFun::calculateSing(double* uij, double *uipj, double* uinj, double* uijp, double* uijn)
{
	__m256d vec1, vec2, vec3, vec4, vec5,vecW, veCij,veCipj,veCinj,veCijp,veCijn,veCbx,veCby;
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
            
//	vecS1 = _mm256_mul_pd(veCbx,_mm256_mul_pd(vec1,_mm256_sub_pd(vec3,vec1)));
//	vecS2 = _mm256_mul_pd(veCby,_mm256_mul_pd(vec1,_mm256_sub_pd(vec5,vec1)));
//	vecS = _mm256_add_pd(_mm256_mul_pd(veCbx,_mm256_mul_pd(vec1,_mm256_sub_pd(vec3,vec1))),_mm256_mul_pd(veCby,_mm256_mul_pd(vec1,_mm256_sub_pd(vec5,vec1))));
	//Below is the formula that works 
	vecW = _mm256_fmadd_pd(veCij,vec1,_mm256_fmadd_pd(veCinj,vec3,_mm256_fmadd_pd(veCipj,vec2,_mm256_fmadd_pd(veCijp,vec4,_mm256_fmadd_pd(veCijn,vec5,_mm256_add_pd(_mm256_mul_pd(veCbx,_mm256_mul_pd(vec1,_mm256_sub_pd(vec3,vec1))),_mm256_mul_pd(veCby,_mm256_mul_pd(vec1,_mm256_sub_pd(vec5,vec1)))))))));
	
	double* f = (double*)&vecW;
	cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<" "<<endl;
//	return vec1;
}


void avxFun::calculateMat(double* uij, double* uout, int& mNx, int& DoTimes)
{
	//Initialisation
	__m256d vec1u, vec2u, vec3u, vec4u, vec5u,vecWu, veCij,veCipj,veCinj,veCijp,veCijn,veCbx,veCby;
	
	veCij = _mm256_set1_pd(Cij);
	veCipj = _mm256_set1_pd(Cipj);
	veCinj = _mm256_set1_pd(Cinj);
	veCijp = _mm256_set1_pd(Cijp);
	veCijn = _mm256_set1_pd(Cijn);
	veCbx = _mm256_set1_pd(Cbx);
	veCby = _mm256_set1_pd(Cby);
	vec1u = _mm256_loadu_pd(uij);
	
//	vec2u = _mm256_loadu_pd(uipj);
//	vec3u = _mm256_loadu_pd(uinj);
	vec4u = _mm256_loadu_pd(uij-mNx);
//	vec5u = _mm256_loadu_pd(uijn);

	int offset;
	for (int j=0;j<DoTimes;j++)
	{
		offset=j*mNx*3;
		//First iteration
//		vec1u = _mm256_loadu_pd(uij);
		vec2u = _mm256_loadu_pd(uij+1+offset);
		vec3u = _mm256_loadu_pd(uij-1+offset);
//		vec4u = _mm256_loadu_pd(uijp);
		vec5u = _mm256_loadu_pd(uij+offset+mNx);
		vecWu = _mm256_fmadd_pd(veCij,vec1u,_mm256_fmadd_pd(veCinj,vec3u,_mm256_fmadd_pd(veCipj,vec2u,_mm256_fmadd_pd(veCijp,vec5u,_mm256_fmadd_pd(veCijn,vec4u,_mm256_add_pd(_mm256_mul_pd(veCbx,_mm256_mul_pd(vec1u,_mm256_sub_pd(vec3u,vec1u))),_mm256_mul_pd(veCby,_mm256_mul_pd(vec1u,_mm256_sub_pd(vec4u,vec1u)))))))));
		f=(double*)&vecWu;
		cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<" "<<endl;
		
		
		//Second iteration
//		vec1u = _mm256_loadu_pd(uij);
		vec2u = _mm256_loadu_pd(uij+1+offset+mNx);
		vec3u = _mm256_loadu_pd(uij-1+offset+mNx);
		vec4u = _mm256_loadu_pd(uij+offset+mNx+mNx);
//		vec5u = _mm256_loadu_pd(uijn+offset);
		vecWu = _mm256_fmadd_pd(veCij,vec5u,_mm256_fmadd_pd(veCinj,vec3u,_mm256_fmadd_pd(veCipj,vec2u,_mm256_fmadd_pd(veCijp,vec4u,_mm256_fmadd_pd(veCijn,vec1u,_mm256_add_pd(_mm256_mul_pd(veCbx,_mm256_mul_pd(vec5u,_mm256_sub_pd(vec3u,vec5u))),_mm256_mul_pd(veCby,_mm256_mul_pd(vec5u,_mm256_sub_pd(vec1u,vec5u)))))))));
		f=(double*)&vecWu;
		cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<" "<<endl;
		
		//Third iteration
		vec1u = _mm256_loadu_pd(uij+offset+mNx+mNx+mNx);
		vec2u = _mm256_loadu_pd(uij+1+offset+mNx+mNx);
		vec3u = _mm256_loadu_pd(uij-1+offset+mNx+mNx);
//		vec4u = _mm256_loadu_pd(uijp+offset+mNx);
//		vec5u = _mm256_loadu_pd(uijn+offset);
		vecWu = _mm256_fmadd_pd(veCij,vec4u,_mm256_fmadd_pd(veCinj,vec3u,_mm256_fmadd_pd(veCipj,vec2u,_mm256_fmadd_pd(veCijp,vec1u,_mm256_fmadd_pd(veCijn,vec5u,_mm256_add_pd(_mm256_mul_pd(veCbx,_mm256_mul_pd(vec4u,_mm256_sub_pd(vec3u,vec4u))),_mm256_mul_pd(veCby,_mm256_mul_pd(vec4u,_mm256_sub_pd(vec5u,vec4u)))))))));
		f=(double*)&vecWu;
		cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<" "<<endl;
	}
}