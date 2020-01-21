#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include <chrono>
#include <cuda_profiler_api.h>
#include <fstream>
using namespace std;
// Fractional order model calculated on GPU with nvcc

#define ThreadPerBlock 1024
__global__
void add(int m, int n, double *Y, double *dx1, double *dx2, double *temp2)
{  
	dx1[m-1] = temp2[0];
	dx2[m-1] = temp2[1];
	__shared__ double cache[ThreadPerBlock]; 
	double temp = 0;
	int sindex = threadIdx.x;
	int stride = blockDim.x;
	for(int s = sindex; s <min(m,n); s+=stride)
	{
		if(s==0 || s==1)
		{
		
		}
		else
		{
			if(blockIdx.x == 0) temp += (Y[s] * dx1[m - s]);
			if(blockIdx.x == 1) temp += (Y[s] * dx2[m - s]);
		}
	}
	cache[sindex] = temp;
	__syncthreads();	
	int i = blockDim.x / 2; 
	while(i!=0)
	{
		if(sindex < i)
			cache[sindex] += cache[sindex + i];
		__syncthreads();
		i /= 2;
	}
	if(sindex == 0)
	{
		if(blockIdx.x == 0) temp2[0] = cache[0] * (-1); //dx1[m] = cache[0] * (-1);
		if(blockIdx.x == 1) temp2[1] = cache[0] * (-1); //dx2[m] = cache[0] * (-1);
	}
}
int main()
{
    //10-1024, 15-32768
    int N = 1 << 18; //19 //T
    int n = 1 << 15;      //L
    double a = 0.8;
    int ngpus = 4;
    const int NGPUS = 4;
    int const S = 8;
    int SS, P, s;
    double *d_Y[NGPUS], *d_dx[NGPUS * 2];
    double *h_Y, *h_dx[NGPUS * 2], *u, *y;
    double *A[S], *B, *C;
    double *d_temp[4], *h_temp[4];
    cudaStream_t stream[NGPUS];
    // memory allocation
    for (int i = 0; i < ngpus; i++)
    {
        cudaSetDevice(i);
        cudaStreamCreate(&stream[i]);
        cudaMalloc((void **)&d_Y[i], n * sizeof(double));
        cudaMalloc((void **)&d_temp[i], 2 * sizeof(double));
        cudaMallocHost((void **)&h_temp[i], 2 * sizeof(double));
        // two d_dx and h_dx for each calculation block on each GPU
        cudaMalloc((void **)&d_dx[i * 2], N * sizeof(double));
        cudaMalloc((void **)&d_dx[i * 2 + 1], N * sizeof(double));
        cudaMallocHost((void **)&h_dx[i * 2], N * sizeof(double));
        cudaMallocHost((void **)&h_dx[i * 2 + 1], N * sizeof(double));
    }
    for (int st = 0; st < S; st++)
    {
        cudaMallocHost((void **)&A[st], S * sizeof(double));
    }
    cudaMallocHost((void **)&B, S * sizeof(double));
    cudaMallocHost((void **)&C, S * sizeof(double));

    // set matrice data
    switch (S)
    {
    case 2: //2 order system parameters!!!!!!!!!!!!!!!!alfa
        A[0][0] = 0.7;
        A[0][1] = 0;
        A[1][0] = 1;
        A[1][1] = 0.4;
        B[0] = 1;
        B[1] = 0;
        C[0] = 0;
        C[1] = 1;
        break; //end of system/*/
    case 4:    //4 order system parameters
        A[0][0] = -2.1624 + 0.8, A[0][1] = -1.9225, A[0][2] = -0.7318, A[0][3] = -0.0859;
        A[1][0] = 1.0, A[1][1] = 0.8, A[1][2] = 0.0, A[1][3] = 0.0;
        A[2][0] = 0.0, A[2][1] = 1.0, A[2][2] = 0.8, A[2][3] = 0.0;
        A[3][0] = 0.0, A[3][1] = 0.0, A[3][2] = 1.0, A[3][3] = 0.8;
        B[0] = 1, B[1] = 0, B[2] = 0, B[3] = 0;
        C[0] = 0, C[1] = -1, C[2] = 0, C[3] = 0;
        break; //end of system/*/
    case 8:    //8order system parameters in progresss
        A[0][0] = -4.4557 + 0.8, A[0][1] = -8.5928, A[0][2] = -8.9663, A[0][3] = -5.2783, A[0][4] = -1.6870, A[0][5] = -0.2600, A[0][6] = -0.0171, A[0][7] = -0.0012;
        A[1][0] = 1.0, A[1][1] = 0.8, A[1][2] = 0.0, A[1][3] = 0.0, A[1][4] = 0.0, A[1][5] = 0.0, A[1][6] = 0.0, A[1][7] = 0.0;
        A[2][0] = 0.0, A[2][1] = 1.0, A[2][2] = 0.8, A[2][3] = 0.0, A[2][4] = 0.0, A[2][5] = 0.0, A[2][6] = 0.0, A[2][7] = 0.0;
        A[3][0] = 0.0, A[3][1] = 0.0, A[3][2] = 1.0, A[3][3] = 0.8, A[3][4] = 0.0, A[3][5] = 0.0, A[3][6] = 0.0, A[3][7] = 0.0;
        A[4][0] = 0.0, A[4][1] = 0.0, A[4][2] = 0.0, A[4][3] = 1.0, A[4][4] = 0.8, A[4][5] = 0.0, A[4][6] = 0.0, A[4][7] = 0.0;
        A[5][0] = 0.0, A[5][1] = 0.0, A[5][2] = 0.0, A[5][3] = 0.0, A[5][4] = 1.0, A[5][5] = 0.8, A[5][6] = 0.0, A[5][7] = 0.0;
        A[6][0] = 0.0, A[6][1] = 0.0, A[6][2] = 0.0, A[6][3] = 0.0, A[6][4] = 0.0, A[6][5] = 1.0, A[6][6] = 0.8, A[6][7] = 0.0;
        A[7][0] = 0.0, A[7][1] = 0.0, A[7][2] = 0.0, A[7][3] = 0.0, A[7][4] = 0.0, A[7][5] = 0.0, A[7][6] = 1.0, A[7][7] = 0.8;
        B[0] = 0, B[1] = 0, B[2] = 0, B[3] = 1, B[4] = 0, B[5] = 1, B[6] = 0, B[7] = 0;
        C[0] = 0, C[1] = -1, C[2] = 0, C[3] = 0, C[4] = 1, C[5] = 0, C[6] = 0, C[7] = 0;
        break; //end of system/*/
    }
    cudaMallocHost((void **)&h_Y, n * sizeof(double));
    cudaMallocHost((void **)&u, N * sizeof(double));
    cudaMallocHost((void **)&y, N * sizeof(double));

    // set u, and Y reset y
    for (int j = 0; j < N; j++)
    {
        u[j] = 1;
        y[j] = 0;
        h_dx[0][j] = 0;
        if (j < n)
        {
            if (j == 0)
            {
                h_Y[0] = 1;
            }
            else
            {
                if (j == 1)
                {
                    h_Y[1] = a;
                }
                else
                {
                    h_Y[j] = h_Y[j - 1] * ((a - j + 1) / j);
                    h_Y[j - 1] = h_Y[j - 1] * pow(-1, j + 1);
                }
            }
        }
    }
    h_Y[n - 1] = h_Y[n - 1] * (-1);
    u[0] = 0;

    //values for t=1
    for (SS = 0; SS < S; SS++)
    {
        h_dx[SS][1] = B[SS] * u[1];
        y[1] += C[SS] * h_dx[SS][1];
    }

    // data file
    ofstream myfile;
    myfile.open("FIND.txt", ios::app);

    // load Y to devices
    for (int i = 0; i < ngpus; i++)
    {
        cudaSetDevice(i);
        if (i == 0)
        {
            cudaMemcpyAsync(d_Y[0], h_Y, n * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
        }
        else
        {
            cudaMemcpyAsync(d_Y[i], d_Y[0], n * sizeof(double), cudaMemcpyDeviceToDevice, stream[i]);
        }
    }
    cudaDeviceSynchronize(); //synchornize davices
    for (int i = 0; i < ngpus; i++)
    {
        h_temp[i][0] = h_dx[i * 2][1];
        h_temp[i][1] = h_dx[i * 2 + 1][1];
    }
    for (int i = 0; i < ngpus; i++)
    {
        cudaSetDevice(i);
        cudaMemcpyAsync(d_dx[i * 2], h_dx[i * 2], sizeof(double), cudaMemcpyHostToDevice, stream[i]);
        cudaMemcpyAsync(d_dx[i * 2 + 1], h_dx[i * 2 + 1], sizeof(double), cudaMemcpyHostToDevice, stream[i]);
        cudaMemcpyAsync(d_temp[i], h_temp[i], 2 * sizeof(double), cudaMemcpyHostToDevice, stream[i]);
    }
    cudaDeviceSynchronize();
    double t1, t2;
    // start calculating
    for (int m = 2; m < N; m++)
    {
        if (m >= N-10)
        { 
            t1 = omp_get_wtime();
        }
        for (int i = 0; i < ngpus; i++)
        {
            cudaSetDevice(i);
            add<<<2, ThreadPerBlock, 0, stream[i]>>>(m, n, d_Y[i], d_dx[i * 2], d_dx[i * 2 + 1], d_temp[i]);
        }
        cudaDeviceSynchronize();
        for (int i = 0; i < ngpus; i++)
        {
            cudaSetDevice(i);
            cudaMemcpyAsync(h_temp[i], d_temp[i], 2 * sizeof(double), cudaMemcpyDeviceToHost, stream[i]);
        }
        cudaDeviceSynchronize();
        for (int i = 0; i < ngpus; i++)
        {
            h_dx[i * 2][m] = h_temp[i][0];
            h_dx[i * 2 + 1][m] = h_temp[i][1];
        }
        for (P = 0; P < S; P++)
        {
            for (s = 0; s < S; s++)
            {
                h_dx[P][m] += (A[P][s] * h_dx[s][m - 1]);
            }
            h_dx[P][m] += B[P] * u[m];
        }
        for (P = 0; P < S; P++)
        {
            y[m] += (C[P] * h_dx[P][m]);
        }
        for (int i = 0; i < ngpus; i++)
        {
            h_temp[i][0] = h_dx[i * 2][m];
            h_temp[i][1] = h_dx[i * 2 + 1][m];
        }
        for (int i = 0; i < ngpus; i++)
        {
            cudaSetDevice(i);
            cudaMemcpyAsync(d_temp[i], h_temp[i], 2 * sizeof(double), cudaMemcpyHostToDevice, stream[i]);
        }
        cudaDeviceSynchronize();
        if (m >= N-10)
        { 
            t2 = omp_get_wtime();
            cout << t2-t1 << endl; 
            myfile << t2-t1 << endl;
        }
    }

     //"debug" purpose
     for (int k = 0; k < 10; k++)
     {
         cout << y[k] << ' ' << k << endl;
     }
     cout << endl;
     for (int k = N / 2; k < N / 2 + 10; k++)
     {
         cout << y[k] << ' ' << k << endl;
     }
     cout << endl;
     for (int k = N - 10; k < N; k++)
     {
         cout << y[k] << ' ' << k << endl;
     }

    cudaFree(h_Y);
    cudaFree(d_Y);
    cudaFree(h_dx);
    cudaFree(d_dx);
    cudaFree(u);
    cudaFree(A);
    cudaFree(B);
    cudaFree(C);
    cudaFree(y);

    // data save
     cout << "N=" << N << " n=" << n << " TpB=" << ThreadPerBlock << " NGPUS=" << NGPUS << endl;
    myfile.close();
    return 0;
}
