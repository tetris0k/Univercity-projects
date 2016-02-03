#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

void Solve(double* A , double* b, int n, double* Y){
	double *X = new double [n];
	double *d = new double [n];
	double l,tk;
	int i,j;
	d[0] = A[0] / fabs(A[0]);
	A[0] = sqrt(fabs(A[0]));
	for(i = 1; i < n; i++){
		A[i] /= A[0]*d[0];
	}
	for(i = 1; i < n; i++){
		

		l = A[i*(n+1) - i*(i+1)/2];
		for (int k = 0; k < i; k++){
			tk = k*(k+1)/2;
			l -= A[k*n + i - tk]*A[k*n + i - tk]*d[k];
		}	
		d[i] = l/fabs(l);
		A[i*(n+1) - i*(i+1)/2] = sqrt(fabs(l));


		for(int t = 1; t < n-i; t++){
			j=i+t;			
			l = 0;
			for (int k = 0; k < i; k++){
				tk = k*(k+1)/2;
				l -= A[k*n + i - tk]*A[k*n + j - tk]*d[k];
			}
			tk = i*(i+1)/2;
			l += A[i*n + j - tk];
			l /= A[i*(n+1) - tk]*d[i];
			A[i*(n+1) + t - tk] = l;
		}
	}

	Y[0] = b[0]/A[0];
	for (int k = 1; k < n; k++){
		Y[k] = b[k];
		for (int i = 0; i <	k; i++){
			tk = i*(i+1)/2;
			Y[k] -= Y[i]*A[k + i*n - tk];
		}
		Y[k] /= A[k*n + k - k*(k+1)/2];
	}

	for (int i = 0; i < n; i++){
		X[i] = Y[i]*d[i];
	}

	Y[n-1] = X[n-1] / A[n*(n+1)/2 - 1];
	for (int k = 2; k < n+1; k++){
		Y[n-k] = X[n-k];
		for (int i = 1; i < k; i++){
			Y[n-k] -= Y[n-i]*A[n*(n+1)/2 - i - k*(k-1)/2];
		}
		Y[n-k] /=  A[n*(n+1)/2 - k  - k*(k-1)/2];
	}
//	R_solve(A,b,n,Y);
//	Dsolve(d,Y,n,X);
//	Rsolve(A,X,n,Y);
	delete [] X;
	delete [] d;
}
/*
double dii(double *A, double* d, int i, int n){
	double l = 0;
	for (int k=0; k<i; k++){
		l -= A[k*n + i - k*(k+1)/2]*A[k*n + i - k*(k+1)/2]*d[k];
	}
	
	l+=A[i*(n+1) - i*(i+1)/2];
	return l/fabs(l);

}

double rii(double *A, double *d, int i, int n){
	double l = 0;
	for (int k = 0; k<i; k++){
		l -= A[k*n + i - k*(k+1)/2]*A[k*n + i - k*(k+1)/2]*d[k];
	}
	l += A[i*(n+1) - i*(i+1)/2];
	return sqrt(fabs(l));
}

double rij(double *A, double *d, int i, int j, int n){
	double l = 0;
	for (int k = 0; k < i; k++){
		l -= A[k*n + i - k*(k+1)/2]*A[k*n + j - k*(k+1)/2]*d[k];
	}
	l += A[i*n + j - i*(i+1)/2];
	l /= A[i*(n+1) - i*(i+1)/2]*d[i];
	return l;
}

void Rsolve(double *A, double *b, int n, double *x){
	x[n-1] = b[n-1] / A[n*(n+1)/2 - 1];
	for (int k = 2; k < n+1; k++){
		x[n-k] = b[n-k];
		for (int i = 1; i < k; i++){
			x[n-k] -= x[n-i]*A[n*(n+1)/2 - i - k*(k-1)/2];
		}
		x[n-k] /=  A[n*(n+1)/2 - k  - k*(k-1)/2];
	}
}

void Dsolve(double *d, double *b, int n, double *x){
	for (int i = 0; i < n; i++){
		x[i] = b[i]*d[i];
	}
}

void R_solve(double *A, double *b, int n, double *x){
	x[0] = b[0]/A[0];
	for (int k = 1; k < n; k++){
		x[k] = b[k];
		for (int i = 0; i <	k; i++){
			x[k] -= x[i]*A[k + i*n - i*(i+1)/2];
		}
		x[k] /= A[k*n + k - k*(k+1)/2];
	}
}
*/
