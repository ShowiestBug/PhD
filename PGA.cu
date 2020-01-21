#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include <chrono>
#include <cuda_profiler_api.h>
#include <fstream>
#include <ctime>
using namespace std;

int minimum(double *in, int n)
{	
	int k = 1;
 	int a = in[0];
	for(int i=1;i<n;i++)
	{
		if(a<in[i])
		{
			a = in[i];
			k = i;
		}
	}
	return k;
}
//__global__
//void add(int N, int n, double *Y, double *dx, double *u)
//{
//
//
//}
double fun(double *vec)
{
	double temp=0;
	for(int i=0;i<5;i++)
	{
		temp += vec[i];
	}
	return temp;
}
int main()
{
	int const nvars = 5;
	int numPop = 20;
	int numBlocks = 10;
	int IteCount = 0;
	double *Pop[nvars], *Score;
	double BestScore, *BestGenes;
	cudaMallocManaged(&Score, numPop*numBlocks*sizeof(double));
	for(int n=0;n<nvars;n++)
	{
		cudaMallocManaged(&Pop[n], numPop*numBlocks*sizeof(double));
	}
	int EliteKids = 3;
	int MutaKids = 6;
	int CrossParents = (2 * (numPop - EliteKids - MutaKids)) + MutaKids;
	int Pass = numBlocks - 1;
	if(Pass > 0.2 * numPop)
	{
        Pass = 0.2 * numPop;
	}
	srand((int)time(0));
	//random population
    for(int i=0;i<nvars;i++)
    {
    	for(int j=0;j<numPop*numBlocks;j++)
    	{
    		Pop[i][j] = (rand() % 1000 + 1);
    		Pop[i][j] /= 1000;
    	}
	}
    //score population
    double *vec;
    cudaMallocManaged(&vec, nvars*sizeof(double));
    for(int i=0;i<numPop*numBlocks;i++)
    {
    	for(int j=0;j<nvars;j++)
    	{
    		vec[j] = Pop[j][i];
    	}
    	Score[i] = fun(vec);
	}
	//main loop
    int exitFlag = 0;
//    do
//    {
//    IteCount++;
	double a[] = {3,7,2,5,6,4,9,10};
	cout<<minimum(a,8);
//
//    exitFlag++;
//	}while(exitFlag==0);
	return 0;
}
