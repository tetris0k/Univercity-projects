#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 10E-10

double Aij(int i, int j){
    return fabs(i-j)+1.;
}
	
int func_input(double* A, int n){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			A[i*(n+1) + j] = Aij(i,j);
		}
        if ( fabs(A[i*n + i]) < EPS ){
			printf("на диагонали нулевой элемент\n");
			return -1;
		}
		A[i*(n+1)+n]=0.;
		for (int k = 0; k < n; k+=2){
			A[i*(n+1) + n]-=Aij(i,k);
		}
	}
	return 0;
}
