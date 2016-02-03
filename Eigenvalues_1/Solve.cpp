#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 10e-100

int values(double *A, int n, double *Q1, double *Q2, double *v, double eps){
	double norm = matrix_norm(n,A);
	double qe = fabs(eps*norm);
//	printf("need %e\n", qe);
	if (n > 2){
		TDiag(A,n);
	}
	double* A1;
	double sk;
	double t = fabs(A[(n-1)*n + (n-2)]);
	int r = n;
	int c = 0;
	while (r > 2)
	{
		sk = A[(r-1)*n + (r-1)];
		for (int i = 0; i < r; i++){
			A[i*n + i] -= sk;
		}
		A1 = new double [r*r];
		for (int i = 0; i < r; i++){
			for (int j = 0; j < r; j++){
				A1[i*r+j] = A[i*n + j];
			}
		}
		if (QR(A1,r,Q1,Q2)==-1){
			delete [] A1;
			return -1;
		}
		for (int i = 0; i < r; i++){
			for (int j = 0; j < r; j++){
				A[i*n+j] = A1[i*r+j];
			}
		}
		delete [] A1;
		for (int i = 0; i < r; i++){
			A[i*n + i] += sk;
		}
		if ((fabs(A[(r-1)*n + (r-2)])-t)<EPS){
			c++;
		}
		t = fabs(A[(r-1)*n + (r-2)]);
		if (t < qe){
			v[r-1] = A[(r-1)*n + (r-1)];
			r--;
		}
		if (c > 10){
		   return 3;
		}
	}
	double D = (A[0]+A[n+1])*(A[0]+A[n+1]) - 4.*(A[0]*A[n+1] - A[1]*A[n]);
	v[1] = (A[0]+A[n+1] - sqrt(D))/2.;
	v[0] = (A[0]+A[n+1] + sqrt(D))/2.;
	return 1;
}

void TDiag(double* A, int n){
	for (int i = 1; i < n-1; i++){
		for (int j = i+1; j < n; j++){
			double x = A[i*n + i-1];
			double y = A[j*n + i-1];
			double R = sqrt(x*x + y*y);
			if (fabs(y) > EPS){
				if (R > EPS){				
					double cos = x/R;
					double sin = -y/R;
					double ii = A[i*n+i]*cos - A[j*n+i]*sin;
					double ij = A[i*n+j]*cos - A[j*n+j]*sin;
					double ji = A[i*n+i]*sin + A[j*n+i]*cos;
					double jj = A[i*n+j]*sin + A[j*n+j]*cos;
					A[(i-1)*n+i] = R;
					A[(i-1)*n+j] = 0.;			
					A[i*n+i-1] = R;
					A[j*n+i-1] = 0.;
					A[i*n+i] = ii*cos - ij*sin;
					A[i*n+j] = ii*sin + ij*cos;
					A[j*n+i] = A[i*n + j];
					A[j*n+j] = ji*sin + jj*cos;					
					for (int k = i+1; k < n; k++){
						if (k != j){
							x = A[i*n+k];
							y = A[j*n+k];
							A[k*n+i] = x*cos - y*sin;
							A[i*n+k] = A[k*n + i]; 
							A[k*n+j] = x*sin + y*cos;
							A[j*n+k] = A[k*n + j];
						}
					}
				}
			}
		}
	}
//	PrintMatrix(n, A, 5);
}

int QR(double* A , int n, double* Q1, double* Q2){
	for (int k = 1; k < n; k++){
		Xotr(A,n,k,Q1,Q2);
		UA(A,Q1,Q2,n,k);
	}
	for (int k = 1; k < n; k++){
		if (AU(A,Q1,Q2,n,k)==-1){
			return -1;
		}
	}
	return 1;
}

void Xotr(double* A, int n, int k, double* Q1, double* Q2){
	double s = A[k*n + k-1]*A[k*n + k-1];
	double na = sqrt(A[(k-1)*n + k-1]*A[(k-1)*n + k-1] + s);
	if (na < EPS){
		Q1[k-1] = 0.;
		Q2[k-1] = 0.;
		A[k*n + k-1] = 0.;
		A[(k-1)*n + k-1] = 0.;
	} else {
		Q1[k-1] = A[(k-1)*n + k-1] - na;
		Q2[k-1] = A[k*n + k-1];
		double nx = sqrt(Q1[k-1]*Q1[k-1] + s);
		if (nx > EPS){
			Q1[k-1] /= nx;
			Q2[k-1] /= nx;
		} else {
			Q1[k-1] = 0.;
			Q2[k-1] = 0.;
			A[k*n + k-1] = 0.;
		}
	}
}

void UA(double* A, double* Q1, double* Q2, int n, int k){
	int q = n-k+1;
	double sk;	
	double y1,y2;
	for (int i = 0; i < q; i++){
		y1 = A[(k-1)*n + k-1+i];
		y2 = A[k*n + k-1+i];
		sk = 2.*(y1*Q1[k-1] + y2*Q2[k-1]);
		A[(k-1)*n + k-1+i] -= sk*Q1[k-1];
		A[k*n + k-1+i] -= sk*Q2[k-1];		
	}
}


int AU(double* A, double* Q1, double* Q2, int n, int k){
	double sk;
	double y1,y2;
	for (int i=0; i<n; i++){
		y1 = A[i*n + k-1];
		y2 = A[i*n + k];
		if (fabs((y1/4.)*Q1[k-1]) > 9223372036854775807/2.){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		if (fabs((y2/4.)*Q2[k-1]) > 9223372036854775807/2.){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		if (fabs(y1*(Q1[k-1]/2.) + y2*(Q2[k-1]/2.)) > 9223372036854775807/2.){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		if (fabs((y1*Q1[k-1] + y2*Q2[k-1])) > 9223372036854775807/2.){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		sk = 2.*(y1*Q1[k-1] + y2*Q2[k-1]);
		if (fabs(sk*(Q1[k-1]/2.)) > 9223372036854775807/2.){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		if (fabs(sk*(Q2[k-1]/2.)) > 9223372036854775807/2.){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		if (A[i*n + k-1] < (sk*Q1[k-1] - 9223372036854775807)){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		if (A[i*n + k] < (sk*Q2[k-1] - 9223372036854775807)){
			printf("can't work with such a big numbers.\n");
			return -1;
		}
		A[i*n + k-1] -= sk*Q1[k-1];
		A[i*n + k] -= sk*Q2[k-1];
	}
	return 1;
}


	

