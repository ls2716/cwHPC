#include <immintrin.h>
#include <stdio.h>
#include <mpi.h>
#include <chrono>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
	
	cout << "Doing AVX" << endl;
	typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
	__m256 result;
	
	for (int i=0; i<100000000; i++)
	{
    /* Initialize the two argument vectors */
    __m256 evens = _mm256_set_ps(2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0);
    __m256 odds = _mm256_set_ps(1.0, 4.0, 5.0, 7.0, 9.0, 11.5, 13.0, 15.0);

    /* Compute the difference between the two vectors */
    result = _mm256_sub_ps(evens, odds);
	}
	hrc::time_point end = hrc::now();
	chrono::duration<double> elapsed_seconds = end-start;

	cout << "Time elapsed for AVX: "<<elapsed_seconds.count()<<"s for normal boundary"<<endl;
    /* Display the elements of the result vector */
    float* f = (float*)&result;
    printf("AVX result %f %f %f %f %f %f %f %f\n", f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);

	start = hrc::now();
	double result2[8];
	for (int i=0; i<100000000; i++)
	{
    /* Initialize the two argument vectors */
    double evens[8] = {16.0,14.0,12.0,10.0,8.0,6.0,4.0,2.0};
	double odds[8] = {15.0,13.0,11.5,9.0,7.0,5.0,4.0,1.0};


    /* Compute the difference between the two vectors */
    
	for (int i=0;i<8;i++)
		result2[i]=evens[i]-odds[i];
	}
	
	end = hrc::now();
	elapsed_seconds = end-start;
	
	cout << "Time elapsed for normal: "<<elapsed_seconds.count()<<"s for normal boundary"<<endl;
    /* Display the elements of the result vector */
    f = (float*)&result2;
    printf("AVX result %f %f %f %f %f %f %f %f\n", f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
	
    MPI_Finalize();
    return 0;
}