#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 1e-50

void TrigOtra(int n, double* a){
	int i,j,t;
	double Sk, normA, normX,l;
	double* x;
	double* y;
	double* z;
	x = new double [n];
	y = new double [n];
	z = new double [n];
	for (i = 0; i < n-2; i++)
	{
		normA = 0.;
		Sk = 0.;
		normX = 0.;
		l = 0.;
		for (j = 0; j < i+1; j++){
			x[j] = 0.;
		}
		x[i+1] = a[(i+1)*n + i]; 
		for (j = i + 2; j < n; j++) {
			x[j] = a[j*n + i]; 
			Sk = Sk + x[j]*x[j]; 
		}

		normA = sqrt(x[i+1]*x[i+1] + Sk);
		if (normA < 1e-100){
			a[(i+1)*n + i] = 0.;
			a[(i+2)*n + i] = 0.;
			continue;
		}
		if (Sk < 1e-100){
			a[(i+2)*n + i] = 0.;
			continue;
		}
		

		x[i+1] =x[i+1] - normA;

		normX = sqrt(x[i+1]*x[i+1] + Sk);
		
		for (j = i + 1; j < n; ++j){
			x[j] = x[j] /  normX;
		}	                      //получили вектор отраж
		for (j = 0; j < n; j++){
			y[j] = 0.;
		}
		for (j = 0; j < n; j++){
			for (t = i+1; t < n; t++){
				y[j] = y[j] + x[t]*a[j*n + t];
			}
		}

		for (j = i+1; j < n; j++){
			l = l + x[j]*y[j];
		}
		for (j = 0; j < n; j++){
			z[j] = 2.*y[j] - 2.*l*x[j];
		}
		for (j = 0; j < n; j++){
			for (t = i+1; t < n; t++){
				a[j*n+t] -= z[j]*x[t];
			}
		}
		for (j = i+1; j < n; j++){
			for (t = 0; t < n; t++){
				a[j*n+t] -= x[j]*z[t];
			}
		}

		for (j = i+2; j < n; j++){
			a[j*n + i] = 0.;
		}
		for (j = i+2; j < n; j++){
			a[i*n + j] = 0.;
		}
	}
	delete [] x;
	delete [] y;
	delete [] z;
}


void QR(int n, double* a, int k, double* cos, double* sin){
	int i,j;
	double x,y,r;
	for (i = 0; i < k-1 ; i++){ 
		x = a[i*n + i];
		y = a[(i+1)*n + i];
		r = sqrt(x*x + y*y);                                 
		if (r < EPS){
			if (a[i * n + i] > 0.) { 
				cos[i] = 1.;
			} else {
				cos[i] = (-1.);
			}
			sin[i] = 0.;
		} else {
			cos[i] = x / r;
			sin[i] = -y / r;
		}
		a[i*n + i] = r;
		a[(i+1)*n + i] = 0.;		
		for (j = i + 1; j < k; j++){
			x = a[i*n + j];
			y = a[(i+1)*n + j];
			a[i*n + j] = x*cos[i] - y*sin[i];
			a[(i+1)*n + j] = x*sin[i] + y*cos[i];
		}
	}
}

void RQ(int n, double* a, int k, double* cos, double* sin){
	int i,j;
	double x,y;
	for (i = 0; i < k-1; i++){
		for (j = 0; j < i+2; j++){
			x = a[j*n + i];
			y = a[j*n + i+1];
			a[j*n + i] = x*cos[i] - y*sin[i];
			a[j*n + i+1] = x*sin[i] + y*cos[i];
		}
	}
}

void FindValues(int n, double* a, double* values, double eps, double* sin, double* cos){
	double e = eps*Norm(n,a);
	if (n > 2){
		TrigOtra(n, a);
	}
	double sk;
	int k = n;
	while (k > 2) {
		sk = a[(k-1)*n + (k-1)];
		for (int i = 0; i < k; i++){
			a[i*n + i] -= sk;
		}
		QR(n, a, k, cos, sin);
		RQ(n,a,k,cos,sin);	
		for (int i = 0; i < k; i++){
			a[i*n + i] += sk;
		}
		if (fabs(a[(k-1)*n + (k-2)]) < e){
			values[k-1] = a[(k-1)*n + (k-1)];
			k--;
		}		
	}
	double d = (a[0]+a[n+1])*(a[0]+a[n+1]) - 4.*(a[0]*a[n+1] - a[1]*a[n]);
	values[0] = (a[0]+a[n+1] + sqrt(d))/2.;
	values[1] = (a[0]+a[n+1] - sqrt(d))/2.;
}

double Norm(int n, double* a){
	int i,j;
	double t;
	double res = 0.;
	for (i = 0; i < n; i++){
		t = 0.;
		for (j = 0; j < n; j++){
			t = t + fabs(a[i * n + j]);
		}
		if (res < t){
			res = t;
		}
	}
	return res;
}

void PrintVector(int n, double* x, int MAX_OUTPUT_SIZE){
	int i;
	int nPrint;
	nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;
	for (i = 0; i < nPrint; ++i)
		printf("%5.8g ", x[i]);
	printf("\n");
}

void PrintMatrix(int n, double * A, int MAX_OUTPUT_SIZE){
    int i;
    int j;
    int nPrint;

    nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;

    for (i = 0; i < nPrint; ++i) 	{
        printf("| ");
        for (j = 0; j < nPrint; ++j)
            printf("%5.3g ", A[i * n + j]);
        if (nPrint < n)
            printf("\t...\t");
        printf("|");
		printf ("\n");
	}
}
