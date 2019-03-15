#include <immintrin.h>
#include <stdio.h>
//#include <mpi.h>
#include <chrono>
#include <iostream>
#include "avxFun.h"

using namespace std;

void avxFun::calculateSing(double* uij, double *uipj, double* uinj, double* uijp, double* uijn)
{
	__m256d vec1, vec2, vec3, vec4, vec5,vecW;
	vec1 = _mm256_loadu_pd(uij);
	
	double* f = (double*)&vec1;
	cout<<f[0]<<endl;
//	return vec1;
}
