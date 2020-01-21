#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include <math.h>
 
using namespace std; 

int main()
{
	double *Y, *u, *dx;
	double a = 0.5;
	int N = 1 << 17;//19
	int n = 1 << 16;//17
	long l1 = n;//100k
	long l2 = N;//400k
	double Time = 0, Time2 = 0,t1,t2;
	ofstream myfile;
	Y = (double*)malloc(l1 * sizeof(double));
	u = (double*)malloc(l2 * sizeof(double));
	dx = (double*)malloc(l2 * sizeof(double));
	int m, s, f;
	for (int i = 0; i < l1; i++)
	{
		if (i == 0)
		{
			Y[0] = 1;
		}
		else
		{
			if (i == 1)
			{
				Y[1] = a;
			}
			else
			{
				Y[i] = Y[i - 1] * ((a - i + 1) / i);
				Y[i - 1] = Y[i - 1] * pow(-1, i + 1);
			}
		}
	}
	Y[n-1] = Y[n-1] *(-1);
	u[0] = 0;
	dx[0] = 0;
	for (int k = 1; k < l2; k++)
	{
		u[k] = 1;
		dx[k] = 0;
	}
	for(int main = 0; main < 1; main++)//main loop
	{
		for (int k = 0; k < l2; k++)//reset input, output and states
		{
			if(k == 0)
			{
				u[k] = 0;	
			}
			else
			{
				u[k] = 1;
			}
			dx[k] = 0;
		}
		t1 = omp_get_wtime();
		#pragma omp parallel for shared(Y,dx,u)  private(m,s,f)
		for (m = 0; m < l2; m++)
		{
			if (m < l1)
			{
				for (s = 0; s < m; s++)
				{
					dx[m] = dx[m] + (Y[s] * u[m - s]);
				}
			}
			else
			{			
				for (f = 0; f < l1; f++)
				{
					dx[m] = dx[m] + (Y[f] * u[m - f]);
				} 
			}
		}
		t2 = omp_get_wtime();	
		Time = Time + (t2-t1);
	}
	//myfile.open("FDD.txt",ios::app);
	//Time=Time/10.0;	
	//myfile<<Time<<" ";
	//myfile.close();
	for(int k=0;k<10;k++)
    {
		cout<<dx[k]<<' '<<k<<endl;
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
	free (Y);
	free (u);
	free (dx);
	return 0;
}
