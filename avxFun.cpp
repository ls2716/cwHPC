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
//	cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<" "<<endl;
//	return vec1;
}


void avxFun::calculateMat(double* uij, double* uout,double* vij, double* vout, int& mNx, int mNy)
{
	//Initialisation
	__m256d vec1u, vec2u, vec3u, vec4u, vec5u,vec1v, vec2v, vec3v, vec4v, vec5v,vecWu, veCij,veCipj,veCinj,veCijp,veCijn,veCbx,veCby;
	
	veCij = _mm256_set1_pd(Cij);
	veCipj = _mm256_set1_pd(Cipj);
	veCinj = _mm256_set1_pd(Cinj);
	veCijp = _mm256_set1_pd(Cijp);
	veCijn = _mm256_set1_pd(Cijn);
	veCbx = _mm256_set1_pd(Cbx);
	veCby = _mm256_set1_pd(Cby);
	
	
	//Doing the rest 
	int cor;
	int j_offset;
	uij--;
	vij--;
	uout--;
	vout--;
			for (int j=0;j<mNy;j++)
			{
				j_offset = j*mNx;
				for (int i=mNx-1-((mNx-2)%4);i<mNx-1;i++)
				{
					cor = i+j_offset;
					uout[cor] = Cij*uij[cor]+Cinj*uij[cor-1]+Cipj*uij[cor+1]+Cijp*uij[cor+mNx]+Cijn*uij[cor-mNx]+Cbx*uij[cor]*(uij[cor-1]-uij[cor])+Cby*vij[cor]*(uij[cor-mNx]-uij[cor]);
					vout[cor] = Cij*vij[cor]+Cinj*vij[cor-1]+Cipj*vij[cor+1]+Cijp*vij[cor+mNx]+Cijn*vij[cor-mNx]+Cbx*uij[cor]*(vij[cor-1]-vij[cor])+Cby*vij[cor]*(vij[cor-mNx]-vij[cor]);
			}
		}
		for (int j=mNy-(mNy%3);j<mNy;j++)
			{
				j_offset = j*mNx;
				for (int i=1;i<(mNx-1);i++)
				{
					cor = i+j_offset;
					uout[cor] = Cij*uij[cor]+Cinj*uij[cor-1]+Cipj*uij[cor+1]+Cijp*uij[cor+mNx]+Cijn*uij[cor-mNx]+Cbx*uij[cor]*(uij[cor-1]-uij[cor])+Cby*vij[cor]*(uij[cor-mNx]-uij[cor]);
					vout[cor] = Cij*vij[cor]+Cinj*vij[cor-1]+Cipj*vij[cor+1]+Cijp*vij[cor+mNx]+Cijn*vij[cor-mNx]+Cbx*uij[cor]*(vij[cor-1]-vij[cor])+Cby*vij[cor]*(vij[cor-mNx]-vij[cor]);
			}
		}
	
	uij=uij-3;
	uout=uout-3;
	vij=vij-3;
	vout=vout-3;
	int DoTimes = mNy/3;
	
	

        for (int ioff = 0; ioff < (mNx - 5); ioff += 4) {
	    uij = uij + 4;
		vij = vij + 4;
		uout = uout+ 4;
		vout = vout+ 4;
	    vec1u = _mm256_loadu_pd(uij);
		vec1v = _mm256_loadu_pd(vij);
	    //	vec2u = _mm256_loadu_pd(uipj);
	    //	vec3u = _mm256_loadu_pd(uinj);
	    vec4u = _mm256_loadu_pd(uij - mNx);
		vec4v = _mm256_loadu_pd(vij - mNx);
	    //	vec5u = _mm256_loadu_pd(uijn);

	    int offset;
	    for(int j = 0; j < DoTimes; j++) {
	        offset = j * mNx * 3;
	        // First iteration u
	        //		vec1u = _mm256_loadu_pd(uij);
	        vec2u = _mm256_loadu_pd(uij + 1 + offset);
	        vec3u = _mm256_loadu_pd(uij - 1 + offset);
	        //		vec4u = _mm256_loadu_pd(uijp);
	        vec5u = _mm256_loadu_pd(uij + offset + mNx);
	        vecWu = _mm256_fmadd_pd(veCij, vec1u,
	            _mm256_fmadd_pd(veCinj, vec3u,
	                _mm256_fmadd_pd(veCipj, vec2u,
	                    _mm256_fmadd_pd(veCijp, vec5u,
	                        _mm256_fmadd_pd(veCijn, vec4u,
	                            _mm256_add_pd(
	                                _mm256_mul_pd(veCbx, _mm256_mul_pd(vec1u, _mm256_sub_pd(vec3u, vec1u))),
	                                _mm256_mul_pd(veCby, _mm256_mul_pd(vec1v, _mm256_sub_pd(vec4u, vec1u)))))))));
	        f = (double*)&vecWu;
	        copy(f, f + 4, uout + offset);
//	        cout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << endl;
			// First iteration v
	        //		vec1u = _mm256_loadu_pd(uij);
	        vec2v = _mm256_loadu_pd(vij + 1 + offset);
	        vec3v = _mm256_loadu_pd(vij - 1 + offset);
	        //		vec4u = _mm256_loadu_pd(uijp);
	        vec5v = _mm256_loadu_pd(vij + offset + mNx);
	        vecWu = _mm256_fmadd_pd(veCij, vec1v,
	            _mm256_fmadd_pd(veCinj, vec3v,
	                _mm256_fmadd_pd(veCipj, vec2v,
	                    _mm256_fmadd_pd(veCijp, vec5v,
	                        _mm256_fmadd_pd(veCijn, vec4v,
	                            _mm256_add_pd(
	                                _mm256_mul_pd(veCbx, _mm256_mul_pd(vec1u, _mm256_sub_pd(vec3v, vec1v))),
	                                _mm256_mul_pd(veCby, _mm256_mul_pd(vec1v, _mm256_sub_pd(vec4v, vec1v)))))))));
	        f = (double*)&vecWu;
	        copy(f, f + 4, vout + offset);
//	        cout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << endl;

	        // Second iteration u 
	        //		vec1u = _mm256_loadu_pd(uij);
	        vec2u = _mm256_loadu_pd(uij + 1 + offset + mNx);
	        vec3u = _mm256_loadu_pd(uij - 1 + offset + mNx);
	        vec4u = _mm256_loadu_pd(uij + offset + mNx + mNx);
	        //		vec5u = _mm256_loadu_pd(uijn+offset);
	        vecWu = _mm256_fmadd_pd(veCij, vec5u,
	            _mm256_fmadd_pd(veCinj, vec3u,
	                _mm256_fmadd_pd(veCipj, vec2u,
	                    _mm256_fmadd_pd(veCijp, vec4u,
	                        _mm256_fmadd_pd(veCijn, vec1u,
	                            _mm256_add_pd(
	                                _mm256_mul_pd(veCbx, _mm256_mul_pd(vec5u, _mm256_sub_pd(vec3u, vec5u))),
	                                _mm256_mul_pd(veCby, _mm256_mul_pd(vec5v, _mm256_sub_pd(vec1u, vec5u)))))))));
	        f = (double*)&vecWu;
	        copy(f, f + 4, uout + offset + mNx);
//	        cout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << endl;
			// Second iteration v 
	        //		vec1u = _mm256_loadu_pd(uij);
	        vec2v = _mm256_loadu_pd(vij + 1 + offset + mNx);
	        vec3v = _mm256_loadu_pd(vij - 1 + offset + mNx);
	        vec4v = _mm256_loadu_pd(vij + offset + mNx + mNx);
	        //		vec5u = _mm256_loadu_pd(uijn+offset);
	        vecWu = _mm256_fmadd_pd(veCij, vec5v,
	            _mm256_fmadd_pd(veCinj, vec3v,
	                _mm256_fmadd_pd(veCipj, vec2v,
	                    _mm256_fmadd_pd(veCijp, vec4v,
	                        _mm256_fmadd_pd(veCijn, vec1v,
	                            _mm256_add_pd(
	                                _mm256_mul_pd(veCbx, _mm256_mul_pd(vec5u, _mm256_sub_pd(vec3u, vec5u))),
	                                _mm256_mul_pd(veCby, _mm256_mul_pd(vec5v, _mm256_sub_pd(vec1u, vec5u)))))))));
	        f = (double*)&vecWu;
	        copy(f, f + 4, vout + offset + mNx);
//	        cout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << endl;

	        // Third iteration u
	        vec1u = _mm256_loadu_pd(uij + offset + mNx + mNx + mNx);
	        vec2u = _mm256_loadu_pd(uij + 1 + offset + mNx + mNx);
	        vec3u = _mm256_loadu_pd(uij - 1 + offset + mNx + mNx);
	        //		vec4u = _mm256_loadu_pd(uijp+offset+mNx);
	        //		vec5u = _mm256_loadu_pd(uijn+offset);
	        vecWu = _mm256_fmadd_pd(veCij, vec4u,
	            _mm256_fmadd_pd(veCinj, vec3u,
	                _mm256_fmadd_pd(veCipj, vec2u,
	                    _mm256_fmadd_pd(veCijp, vec1u,
	                        _mm256_fmadd_pd(veCijn, vec5u,
	                            _mm256_add_pd(
	                                _mm256_mul_pd(veCbx, _mm256_mul_pd(vec4u, _mm256_sub_pd(vec3u, vec4u))),
	                                _mm256_mul_pd(veCby, _mm256_mul_pd(vec4v, _mm256_sub_pd(vec5u, vec4u)))))))));
	        f = (double*)&vecWu;
	        copy(f, f + 4, uout + offset + mNx + mNx);
//	        cout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << endl;
			// Third iteration v
	        vec1v = _mm256_loadu_pd(vij + offset + mNx + mNx + mNx);
	        vec2v = _mm256_loadu_pd(vij + 1 + offset + mNx + mNx);
	        vec3v = _mm256_loadu_pd(vij - 1 + offset + mNx + mNx);
	        //		vec4u = _mm256_loadu_pd(uijp+offset+mNx);
	        //		vec5u = _mm256_loadu_pd(uijn+offset);
	        vecWu = _mm256_fmadd_pd(veCij, vec4v,
	            _mm256_fmadd_pd(veCinj, vec3v,
	                _mm256_fmadd_pd(veCipj, vec2v,
	                    _mm256_fmadd_pd(veCijp, vec1v,
	                        _mm256_fmadd_pd(veCijn, vec5v,
	                            _mm256_add_pd(
	                                _mm256_mul_pd(veCbx, _mm256_mul_pd(vec4u, _mm256_sub_pd(vec3u, vec4u))),
	                                _mm256_mul_pd(veCby, _mm256_mul_pd(vec4v, _mm256_sub_pd(vec5u, vec4u)))))))));
	        f = (double*)&vecWu;
	        copy(f, f + 4, vout + offset + mNx + mNx);
//	        cout << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << endl;
			
			
	    }
        }
}