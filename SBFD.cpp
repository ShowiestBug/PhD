#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <windows.h>//temp

//sub-vector fractional-order difference
#define NumTH 8
using namespace std;
void set_thread_affinity() {
	DWORD_PTR mask[8];
	int z;
	for(int i=0;i<8;i++)
	{
		z = i + 48;
		mask[i] = 1<< z;	
	}
	
	#pragma omp parallel default(shared)
	{
		//DWORD_PTR mask = (1 << (omp_get_thread_num()+40));
		SetThreadAffinityMask(GetCurrentThread(), mask[omp_get_thread_num()]);
	}
}
int main(int argc, char **argv) {
	double *Y, *u, *dx;
	double a = 0.5;
	long l1 = 1<<19; //L 16-64k, 19-512k
	long l2 = 1<<19; //T
	Y = (double *)malloc(l1 * sizeof(double));
	u = (double *)malloc(l2 * sizeof(double));
	dx = (double *)malloc(l2 * sizeof(double));
	int m, s, k, d, h, g, test, pp, kk;
	int c = NumTH; // number of cores/threads
	int NN = l1 / c;
	int *Ki, *Pi;
	Pi = (int *)malloc(c * sizeof(int));
	Ki = (int *)malloc(c * sizeof(int));
	for (int x = 0; x < c; x++) {
		if (x == 0) {
			Pi[x] = 0;
		} else {
			Pi[x] = x * NN;
		}
		if (x == (c - 1)) {
			Ki[x] = l1 - 1;
		} else {
			Ki[x] = ((x + 1) * NN) - 1;
		}
	}
	double t1, t2;
	for (int i = 0; i < l1; i++) { //coefficients
		if (i == 0) {
			Y[0] = 1;
		} else {
			if (i == 1) {
				Y[1] = a;
			} else {
				Y[i] = Y[i - 1] * ((a - i + 1) / i);
				Y[i - 1] = Y[i - 1] * pow(-1, i + 1);
			}
		}
	}
	ofstream myfile;
	myfile.open("SBFD19.txt", ios::app);
	myfile << endl;
	for (k = 0; k < l2; k++) { //reset input, output and states
		if (k == 0) {
			u[k] = 0;
		} else {
			u[k] = 1;
		}
		dx[k] = 0;
	}
	double xx = 0.0, yy;

	omp_set_num_threads(NumTH);
	set_thread_affinity();


	for(int l=0;l<c;l++)
	{
		cout << Pi[l] << ' ' << Ki[l] << endl;
	}
	for (m = 0; m < l2; m++) { // calculation of FFD
		if (m >= l2-10) {
			t1 = omp_get_wtime();
		}
		if (m < l1) { //sum for T<L
			test = m / NN + 1; //how many parts are we having now
			#pragma omp parallel for shared(Y, dx, u, m, NN, test, Ki, Pi) private(h, d, yy, kk, pp) default(none) reduction(+:xx)
			for (h = 0; h < test; h++) { //start of sum
				yy = 0.0;
				pp = Pi[h];
				kk = min(m, Ki[h]);
				for (d = pp; d <= kk; d++) {
					yy += (Y[d] * u[m - d]);
				}
				xx += yy;
			}
		} //end of sum for T<L
		else {
			//->sum for T>=L
			#pragma omp parallel for shared(Y, dx, u, m, c, NN, Ki, Pi) private(g, s, yy, kk, pp) default(none) reduction(+:xx)
			for (g = 0; g < c; g++) { //-> divide sum (parallel)
				yy = 0.0;
				pp = Pi[g];
				kk = Ki[g];
				for (s = pp; s <= kk; s++) {
					yy += (Y[s] * u[m - s]);
				} //->end of first part of sum*/
				xx += yy;
			}
		}
		dx[m] = xx;
		xx = 0.0;
		if (m >= l2-10) {
			t2 = omp_get_wtime();
			cout << t2-t1 << endl;
			myfile << t2-t1 << endl;
		}
	}
	for (int k = 0; k < 10; k++)
     {
         cout << dx[k] << ' ' << k << endl;
     }
     cout << endl;
     for (int k = (1<<13)-5;k< (1<<13)+5; k++)
     {
         cout << dx[k] << ' ' << k << endl;
     }
     cout<<endl;
     for (int k = l2-10;k< l2; k++)
     {
         cout << dx[k] << ' ' << k << endl;
     }
     cout << endl;
    // data save
     cout << "N=" << l2 << " n=" << l1 << " Threads=" << NumTH << endl;
	myfile.close();
	free(Y);
	free(u);
	free(dx);
}