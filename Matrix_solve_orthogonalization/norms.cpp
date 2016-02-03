#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#include <pthread.h>

double nev(int n, double *A, double *x){
	double m,l = 0.;
	int w,j;
	for (w = 0; w < n; w++){
		m = 0.;
		for (j = 0; j < n; j++){
			m += A[w*(n+1) + j]*x[j];
		}
		m += A[w*(n+1)+n];
		l += m*m;
	}
	return sqrt(l);
}

double error_rate(int n, double* x){
	double l = 0.;
	for (int i = 0; i < n; i += 2){
		l += (x[i]-1)*(x[i]-1);
	}
	for (int i = 1; i < n; i += 2){
		l += x[i]*x[i];
	}
	return sqrt(l);
}

