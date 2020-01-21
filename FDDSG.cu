#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include <cuda_profiler_api.h>
using namespace std;
//prototype single gpu divided dot product
//N-t
//n-L
#define ThreadPerBlock 256
__global__
void add(int m, int n, double *Y, double *u, double *p_dx)
{
	__shared__ double cache[ThreadPerBlock];
	double temp = 0;
	
  	int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
  	int cacheIndex = threadIdx.x;
  	for (int i = index; i < min(m,n); i += stride)
  	{
  		temp += Y[i] * u[m-i];
	}
	cache[cacheIndex] = temp;
	__syncthreads();	
	int i = blockDim.x / 2; 
	while(i!=0)
	{
		if(cacheIndex < i)
			cache[cacheIndex] += cache[cacheIndex + i];
		__syncthreads();
		i /= 2;
	}
	if(cacheIndex == 0)
		p_dx[blockIdx.x] = cache[0];
}

__global__
void red(double *p_dx, int num)
{
	int idx = threadIdx.x;
	int i = num/ 2; 
	int ii = num % 2;
	while(i!=0)
	{
		if(ii!=0 && idx == 0)
		{
			p_dx[0] += p_dx[num-1];
		}
		if(idx < i)
			p_dx[idx] += p_dx[idx + i];
		__syncthreads();
		ii = i % 2;
		i /= 2;
	}
	
}
int main()
{
	double t1,t2,t3,t4;
	t1 = omp_get_wtime(); 
	
	int N = 1 << 19;//19
	int n = 1 << 17;//17
	double a = 0.5;
	
	int blockSize = 256;
	int numBlocks = (N + blockSize - 1) / blockSize;
	
	double *h_Y, *h_u, *ph_dx, *dx, *tep, *h_t;
	double *d_Y, *d_u, *pd_dx, *d_t;
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	cudaMalloc((void **) &d_Y, n*sizeof(double));
	cudaMalloc((void **) &d_u, N*sizeof(double));
	cudaMalloc((void **) &pd_dx, numBlocks*sizeof(double));
	cudaMalloc((void **) &d_t, 5*sizeof(double));
	
	
	cudaMallocHost((void **) &h_t, 5*sizeof(double));
	cudaMallocHost((void **) &tep, 2*sizeof(double));
	cudaMallocHost((void **) &ph_dx, numBlocks*sizeof(double));
	cudaMallocHost((void **) &h_Y, n*sizeof(double));
	cudaMallocHost((void **) &h_u, N*sizeof(double));
	cudaMallocHost((void **) &dx, N*sizeof(double));

	for(int j=0;j<N;j++)
	{
		h_u[j] = 1;
		dx[j] = 0;
		if(j<n)
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
	h_Y[n-1] = h_Y[n-1] * (-1);
	h_u[0] = 0;
    t3 = omp_get_wtime();
    cudaMemcpyAsync(d_Y, h_Y, n*sizeof(double), cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(d_u, h_u, N*sizeof(double), cudaMemcpyHostToDevice, stream);
    for(int m=0;m<N;m++)
    {
    	numBlocks  = (m + blockSize) / blockSize;
   		add<<<numBlocks, blockSize>>>(m,n,d_Y,d_u,pd_dx);
		cudaDeviceSynchronize();
		cudaMemcpyAsync(ph_dx, pd_dx, numBlocks*sizeof(double), cudaMemcpyDeviceToHost, stream);
		cudaDeviceSynchronize();
		if(numBlocks == 1)
	   	{
	   		dx[m] = ph_dx[0];
	    }
	    else
	    {
			for(int k=0;k<numBlocks;k++)
			{
				dx[m]+=ph_dx[k];		
			}
		}
	}
    t4 = omp_get_wtime();
    for(int k=0;k<10;k++)
    {
		cout<<dx[k]<<' '<<k<<endl;//h_Y[k]<<' '<<h_u[k]<<endl;
	}
	cout<<endl;
	for(int k=N/2;k<N/2+10;k++)
    {
		cout<<dx[k]<<' '<<k<<endl;
	}
	cout<<endl;
	for(int k=N-10;k<N;k++)
    {
		cout<<dx[k]<<' '<<k<<endl;
	}
	cout<<endl;
    cudaStreamDestroy(stream);
    cudaFree(h_Y);
    cudaFree(h_u);
    cudaFree(dx);
    cudaFree(d_Y);
    cudaFree(d_u);
    cudaFree(pd_dx);
    cudaFree(ph_dx);
	t2 = omp_get_wtime();
	cout<<"N="<<N<<"\t"<<"n="<<n<<endl;
	cout<<t4-t3<<' '<<t2-t1<<endl;
	return 0;
}


