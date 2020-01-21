#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <windows.h>

//sub-vector fractional-order model

using namespace std;
#define NumTH 4
void set_thread_affinity()
{
    DWORD_PTR mask[8];
    int z;
    for (int i = 0; i < 8; i++)
    {
        z = i + 48;
        mask[i] = 1 << z;
    }
#pragma omp parallel default(shared)
    {
        SetThreadAffinityMask(GetCurrentThread(), mask[omp_get_thread_num()]);
    }
}
int main(int argc, char **argv)
{
    double *Y, *u, **dx, *y;
    double a = 0.8;
    long l1 = 1 << 14; //L 17-132k, 14-16k
    long l2 = 1 << 16; //T
    const int S = 8; 
	int SS, k;  // states number
    const int c = NumTH;
    int NN = l1 / c;
    int *Ki, *Pi;
    Pi = (int *)malloc(c * sizeof(int));
    Ki = (int *)malloc(c * sizeof(int));
    double **SUM;
    SUM = (double **)malloc(S * sizeof(double *));
    for (SS = 0; SS < S; SS++)
    {
        SUM[SS] = (double *)malloc(c * sizeof(double));
    }
    for (SS = 0; SS < S; SS++)
    {
        for (k = 0; k < c; k++)
        {
            SUM[SS][k] = 0;
        }
    }
    for (int x = 0; x < c; x++)
    {
        if (x == 0)
        {
            Pi[x] = 2;
        }
        else
        {
            Pi[x] = x * NN;
        }
        if (x == (c - 1))
        {
            Ki[x] = l1 - 1;
        }
        else
        {
            Ki[x] = ((x + 1) * NN) - 1;
        }
    }

    Y = (double *)malloc(l1 * sizeof(double)); //coeficients
    y = (double *)malloc(l2 * sizeof(double)); //outpt vector
    u = (double *)malloc(l2 * sizeof(double));
    dx = (double **)malloc(S * sizeof(double *)); // first declare rows then columns \/
    for (int XX = 0; XX < S; XX++)
    {
        dx[XX] = (double *)malloc(l2 * sizeof(double)); // <- each column
    }
    double **A, *B, *C; // D;
    A = (double **)malloc(S * sizeof(double *));
    for (int aa = 0; aa < S; aa++)
    {
        A[aa] = (double *)malloc(S * sizeof(double));
    }
    B = (double *)malloc(S * sizeof(double));
    C = (double *)malloc(S * sizeof(double));
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
    case 8:    //8order system parameters
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
        break; //end of system/?*/
    }
    for (int i = 0; i < l1; i++) //coefficients
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
    int m, g, h;
    int s, p, d;
    int P, U, J, test, pp, kk;
    double t1, t2, t3;
    ofstream myfile;
    myfile.open("SBFDMtest.txt", ios::app);
    myfile << endl;
    for (int l = 0; l < c; l++)
    {
        cout << Pi[l] << ' ' << Ki[l] << endl;
    }
    for (k = 0; k < l2; k++) //reset input, output and states
    {
        if (k == 0)
        {
            u[k] = 0;
        }
        else
        {
            u[k] = 1;
        }
        for (SS = 0; SS < S; SS++)
        {
            dx[SS][k] = 0;
        }
        y[k] = 0;
    }
    for (SS = 0; SS < S; SS++) //first and second states
    {
        dx[SS][0] = 0;
        dx[SS][1] = B[SS] * u[1];
        y[1] = y[1] + C[SS] * dx[SS][1];
    }
    omp_set_num_threads(NumTH);
    set_thread_affinity();
    double temp[S];
    for (m = 2; m < l2; m++)
    {
        if (m >= l2 - 10)
        {
            t1 = omp_get_wtime();
        }
        if (m < l1) //sum for T<L
        {
            test = m / NN + 1; //how many parts are we having now
			#pragma omp parallel for shared(Y, dx, u, m, S, NN, test, Ki, Pi, SUM) private(h, d, U, pp, kk,temp) default(none)
            for (h = 0; h < test; h++) //start of sum
            {
                for (int t = 0; t < c; t++) temp[t] = 0;
                pp = Pi[h];
                kk = min(m, Ki[h]);
                for (U = 0; U < S; U++)
                {
                    for (d = pp; d < kk; d++)
                    {
                        temp[U] += (Y[d] * dx[U][m - d]);
                        //SUM[U][h] += (Y[d] * dx[U][m - d]);
                    }
                    SUM[U][h] = temp[U];
                }
            }
        } //end of sum for T<L
        else
        { //->sum for T>=L
#pragma omp parallel for shared(Y, dx, u, m, S, NN, c, Ki, Pi, SUM) private(g, s, J, pp, kk,temp) default(none)
            for (g = 0; g < c; g++) //-> divide sum (parallel)
            {
                for (int t = 0; t < c; t++)
				{
					temp[t] = 0;
            	}
				pp = Pi[g];
                kk = Ki[g];
                for (J = 0; J < S; J++)
                {
                    for (s = pp; s < kk; s++)
                    {
                        temp[J] += (Y[s] * dx[J][m - s]);
                        //SUM[J][g] += (Y[s] * dx[J][m - s]);
                    }
                    SUM[J][g] = temp[J];
                }
            }                   //end of whole sum
        }                       //end of sum for T>=L
        if (m >= l2 - 10)
        {
            t3 = omp_get_wtime();
        }
        for (P = 0; P < S; P++) //->turn sum into negative, start adding (Af+Ia) and B*Input
        {
            for (s = 0; s < c; s++)
            {
                dx[P][m] += SUM[P][s];
                SUM[P][s] = 0;
            }
            dx[P][m] = dx[P][m] * (-1);
            for (p = 0; p < S; p++) //->(Af+Ia)
            {
                dx[P][m] += (A[P][p] * dx[p][m - 1]);
            }
            dx[P][m] += B[P] * u[m];
        }
        for (P = 0; P < S; P++) //->output
        {
            y[m] += (C[P] * dx[P][m]);
        }
        if (m >= l2 - 10)
        {
            t2 = omp_get_wtime();
            cout << t2 - t1 << t3-t1 << t2-t3 << endl;
            myfile << t2 - t1 << endl;
        }
    }
    for (int k = 0; k < 10; k++)
    {
        cout << y[k] << ' ' << k << endl;
    }
    cout << endl;
    for (int k = l2 / 2; k < l2 / 2 + 10; k++)
    {
        cout << y[k] << ' ' << k << endl;
    }
    cout << endl;
    for (int k = l2 - 10; k < l2; k++)
    {
        cout << y[k] << ' ' << k << endl;
    }
    cout << "N=" << l2 << " n=" << l1 << " Threads=" << NumTH << endl;
    myfile.close();
    free(Y);
    free(SUM);
    free(dx);
    free(y);
    free(u);
    return 0;
}