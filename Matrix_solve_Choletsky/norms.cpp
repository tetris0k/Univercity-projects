#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

double residual_norm1(int n, double *A, double *b, double *x){
	double m,l = 0.;
	for (int i = 0; i < n; i++){
		m = 0.;
		for (int j = 0; j < n; j++){
			m += A[i*n+j]*x[j];
		}
		m -= b[i];
		l += m*m;
	}
	return sqrt(l);
}
double residual_norm(int n, double *b, double *x){
	double m,l = 0.;
	for (int i = 0; i < n; i++){
		m = 0.;
		for (int j = 0; j < n; j++){
			m += Aij(i,j)*x[j];
		}
		m -= b[i];
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
	
