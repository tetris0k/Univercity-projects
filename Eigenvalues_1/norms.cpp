#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 10e-25

double matrix_norm(int n, double* A){
	double sum, max = 0.;
	for (int i = 0; i < n; i++){
		sum = 0.;
		for (int j = 0; j < n; j++){
			sum += fabs(A[i*n + j]);
		}
		if (max < sum)
			max = sum;
	}
	return max;
}	

double vector_norm(double* v, int n){
	double res = 0.;
	for (int i = 0; i < n; i++){
		res += v[i]*v[i];
	}
	return sqrt(res);	
}
	
int issym(double* A, int n){
	for (int i = 0; i < n; i++){
		for (int j = i; j < n; j++){
			if ( fabs (A[i*n+j] - A[j*n+i]) > EPS){
				return 2;
			}
		}
	}
	return 1;
}

double inv_trace(double* A, int n){
	double res = 0.;
	for (int i = 0; i < n; i++){
		res += A[i*n+i];
	}
	return res;
}

double inv_length(double* A, int n){
	double res = 0.;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			res += A[i*n + j]*A[i*n + j];
		}
	}
	return sqrt(res);
}
