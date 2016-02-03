#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include <pthread.h>

void Solve(int n, double* A, double* x){	
	double *E = new double [(n+1)*(n+1)];
	double *v = new double [n+1];
	for (int i = 0; i < n+1; i++){
		for (int j = 0; j < n+1; j++){
			if (i==j) E[i*(n+1)+j]=1; else E[i*(n+1)+j]=0;
		}
	}
	for (int k=1;k<n+1;k++){
		nowv(n,E,v,k);
	    for(int i = k; i < n+1; i++){
			master(n,v,i,k,A,E);
		}
	}
	for (int z=0;z<n;z++){
		x[z] = E[(n+1)*z+n];
	}
	delete []E;
	delete []v;
}


void master(int n, double *v, int i , int k, double *A, double *E){
	double x,y;
	x = smult(n,i,k,A,E);
	y = smult2(n,v,k,A);
	double a;
	a = x/y;
	for (int j = 0; j < n+1; j++){
		E[j*(n+1)+i] -= a*v[j];
	}
}

double smult2(int n, double *v , int k, double *A){
	double res=0.;
	for (int j = 0; j < n+1; j++){
		res = res + A[(k-1)*(n+1)+j]*v[j];
	}	    
	return res;
}

double smult(int n, int i , int k, double *A, double *E){
	double res=0.;
	for (int j = 0; j < k; j++){
		res = res + A[(k-1)*(n+1)+j]*E[j*(n+1)+i];
	}	
	res = res + A[(k-1)*(n+1)+i]*E[i*(n+1)+i];    
	return res;
}

void nowv(int n, double *E, double *v, int k){
	for(int i=0;i<n+1;i++){
		v[i]=E[i*(n+1)+k-1];
	}
}

void PrintMatrix(int n, double * A, double * b, int MAX_OUTPUT_SIZE)
{
    int i;
    int j;
    int nPrint;

    nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;

    for (i = 0; i < nPrint; ++i) 	{
        printf("| ");
        for (j = 0; j < nPrint; ++j)
            printf("%5.3g ", A[i * (n+1) + j]);
        if (nPrint < n)
            printf("\t...\t");
        printf("|");
        printf("%5.3g\n", b[i]);
    }
}
