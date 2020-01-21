#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <cstdlib>
#include <fstream>

using namespace std; 
//computing fracsystem divaided per states
//array[x][y] x - rows, y - columns!!!!!!!
//2n,4n and 8n order system
/*
//DNS po stanach
test l1=l2=50000 a =0.8 S=2 serial= 16.82000 parallel= 16.660000 speedup=1.01 
test l1=l2=100000 a =0.8 S=2 serial= 78.722000 parallel= 78.444000 speedup=1.004062
test l1=l2=100000 a =0.8 S=4 serial= 133.733000 parallel= 90.713000 speedup=1.474
test l1=l2=100000 a =0.8 S=4 serial= 135.101000 parallel= 102.774000 speedup=1.314
test l1=100000 l2=100010 Last10iterations time serial 0.03... parallel 0.02... speedup 1.32
changelog 07.11.2017 22.00------------------------------------------------------------------------

Ver FD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!10.12.2017
*/ 
int main(int argc, char** argv) {
	double *Y, *u, **dx, *y;
	double a = 0.8;
	long l1 = 1<<17;//L
	long l2 = 1<<17;//T
	int N = l2;
	int n = l1;
	int S = 8,SS,k;// states number
	double Time = 0, Time2 = 0,t1,t2,t3;
	ofstream myfile;
	Y = (double*)malloc(l1 * sizeof(double));//coeficients
	y = (double*)malloc(l2 * sizeof(double));//outpt vector *****l1->l2********
	u = (double*)malloc(l2 * sizeof(double));
	dx = (double**)malloc(S * sizeof(double*));// first declare rows then columns \/
	for(int XX=0;XX<S;XX++)
	{
		dx[XX] = (double*)malloc(l2*sizeof(double)); // <- each column
	}
	double **A,*B,*C,D;
	A = (double**)malloc(S*sizeof(double*));	
	for(int aa=0;aa<S;aa++)
	{
		A[aa] = (double*)malloc(S*sizeof(double));
	}  
	B = (double*)malloc(S*sizeof(double));
	C = (double*)malloc(S*sizeof(double));
	switch (S)
	{
		case 2:		//2 order system parameters!!!!!!!!!!!!!!!!alfa
			A[0][0] = 0.7;
			A[0][1] = 0;
			A[1][0] = 1;
			A[1][1] = 0.4;
			B[0] = 1;
			B[1] = 0;
			C[0] = 0;
			C[1] = 1;
			break;//end of system/*/
		case 4:		//4 order system parameters
			A[0][0] = -2.1624+0.8, A[0][1] = -1.9225, A[0][2] = -0.7318, A[0][3] = -0.0859;
			A[1][0] = 1.0,	   A[1][1] = 0.8,     A[1][2] = 0.0,     A[1][3] = 0.0;
			A[2][0] = 0.0,     A[2][1] = 1.0,     A[2][2] = 0.8,     A[2][3] = 0.0;
			A[3][0] = 0.0,     A[3][1] = 0.0,     A[3][2] = 1.0,     A[3][3] = 0.8;
			B[0] = 1, B[1] = 0,  B[2] = 0, B[3] = 0;
			C[0] = 0, C[1] = -1, C[2] = 0, C[3] = 0;
			break; 	//end of system/*/
		case 8:		//8order system parameters in progresss
			A[0][0] = -4.4557+0.8, A[0][1] = -8.5928, A[0][2] = -8.9663, A[0][3] = -5.2783, A[0][4] = -1.6870, A[0][5] = -0.2600, A[0][6] = -0.0171, A[0][7] = -0.0012;
			A[1][0] = 1.0,	   A[1][1] = 0.8,    A[1][2] = 0.0,     A[1][3] = 0.0,     A[1][4] = 0.0,     A[1][5] = 0.0,     A[1][6] = 0.0,     A[1][7] = 0.0;
			A[2][0] = 0.0,     A[2][1] = 1.0,     A[2][2] = 0.8,    A[2][3] = 0.0,     A[2][4] = 0.0,     A[2][5] = 0.0,     A[2][6] = 0.0,     A[2][7] = 0.0;
			A[3][0] = 0.0,     A[3][1] = 0.0,     A[3][2] = 1.0,     A[3][3] = 0.8,    A[3][4] = 0.0,     A[3][5] = 0.0,     A[3][6] = 0.0,     A[3][7] = 0.0;
			A[4][0] = 0.0,     A[4][1] = 0.0,     A[4][2] = 0.0,     A[4][3] = 1.0,     A[4][4] = 0.8,    A[4][5] = 0.0,     A[4][6] = 0.0,     A[4][7] = 0.0;
			A[5][0] = 0.0,     A[5][1] = 0.0,     A[5][2] = 0.0,     A[5][3] = 0.0,     A[5][4] = 1.0,     A[5][5] = 0.8,    A[5][6] = 0.0,     A[5][7] = 0.0;
			A[6][0] = 0.0,     A[6][1] = 0.0,     A[6][2] = 0.0,     A[6][3] = 0.0,     A[6][4] = 0.0,     A[6][5] = 1.0,     A[6][6] = 0.8,    A[6][7] = 0.0;
			A[7][0] = 0.0,     A[7][1] = 0.0,     A[7][2] = 0.0,     A[7][3] = 0.0,     A[7][4] = 0.0,     A[7][5] = 0.0,     A[7][6] = 1.0,     A[7][7] = 0.8;
			B[0] = 0, B[1] = 0,  B[2] = 0, B[3] = 1, B[4] = 0, B[5] = 1,  B[6] = 0, B[7] = 0;
			C[0] = 0, C[1] = -1, C[2] = 0, C[3] = 0, C[4] = 1, C[5] = 0, C[6] = 0, C[7] = 0;
			break; 	//end of system/*/
	}	
	int m,s,f;
	for (int i = 0; i < l1; i++)//coefficients
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
	Y[n-1] = Y[n-1] * (-1);//!!!!!!!
	int P,PP;
	for(int main = 0; main < 1; main++)//main loop
	{
		for (k = 0; k < l2; k++)//reset input, output and states
		{	
			if(k == 0)
			{
				u[k] = 0;	
			}
			else
			{
				u[k] = 1;
			}
			for(SS = 0;SS<S;SS++) 	dx[SS][k] = 0;
			y[k] = 0;
		}
		for(SS = 0; SS<S; SS++)//t=1  	
		{
			dx[SS][1] = B[SS] * u[1];
			y[1] = y[1] + C[SS] * dx[SS][1];
    	}
		t1 = omp_get_wtime();
		for (m = 2; m < l2; m++)
		{
			if(m == (l2-1)) 
			{		
					t3 = omp_get_wtime(); 
			}					
			if (m < l1)
			{
				#pragma omp parallel for shared(Y,dx,u,m,A,B,S) private(P,s)
				for(P = 0;P<S;P++)
				{
					for (s = 2; s < m; s++)
					{
						dx[P][m] = dx[P][m] + (Y[s] * dx[P][m - s]);
					}
					dx[P][m] = dx[P][m] * (-1);
					for(s = 0;s < S; s++)
					{
						dx[P][m] = dx[P][m] + (A[P][s] * dx[s][m-1]);
					}
					dx[P][m] = dx[P][m] + B[P]*u[m];
				}
				for(P = 0; P<S; P++)
				{
					y[m] = y[m]	+ (C[P]*dx[P][m]);
				}
			}
			else
			{
				#pragma omp parallel for shared(Y,dx,u,m,A,B)  private(PP,f)
				for(PP = 0;PP<S;PP++)
				{
					for (f = 2; f < l1; f++)
					{
						dx[PP][m] = dx[PP][m] + (Y[f] * dx[PP][m - f]);
					}
					dx[PP][m] = dx[PP][m] * (-1);
					for(f = 0;f < S; f++)
					{
						dx[PP][m] = dx[PP][m] + (A[PP][f]*dx[f][m-1]);
					}
					dx[PP][m] = dx[PP][m] + B[PP]*u[m];
				}
				for(PP = 0; PP<S; PP++)
				{
					y[m] = y[m]	+ (C[PP]*dx[PP][m]);
				}
			}
		}
		t2 = omp_get_wtime();
		Time2 = Time2 + (t2-t3);
		Time = Time + (t2-t1);
	}
	// for(int i=0;i<10;i++)
	// {
	// 	for(int j=0;j<8;j++)
	// 	{
	// 		cout<<dx[j][i]<<' ';
	// 	}
	// 	cout<<endl;
	// }
	for(int k=0;k<10;k++)
    {
		cout<<y[k]<<' '<<k<<endl;
	}
	cout<<endl;
	for(int k=N/2;k<N/2+5;k++)
    {
		cout<<y[k]<<' '<<k<<endl;
	}
	cout<<endl;
	for(int k=N-5;k<N;k++)
    {
		cout<<y[k]<<' '<<k<<endl;
	}
//	myfile.open("DNS.txt",ios::app);
	Time=Time/1.0;
	Time2=Time2/1.0;
	cout<<l1<<' '<<l2<<endl;
	cout<<Time<<' '<<Time2<<endl;
//	myfile<<Time<<" "<<Time2<<" ";
//	myfile.close();
	return 0;
}
