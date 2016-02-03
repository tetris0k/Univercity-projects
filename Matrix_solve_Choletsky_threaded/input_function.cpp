#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 10E-15

double Aij(int i, int j){
    return fabs(i-j)+1.;
}
	
int func_input(double* A, double* b, int n){
	for (int i = 0; i < n; i++){
		for (int j = i; j < n; j++){
			A[i*n + j - i*(i+1)/2] = Aij(i,j);
		}
		if ( fabs(A[i*n + i - i*(i+1)/2]) < EPS ){
			printf("There are zero elements on diagonal.\n");
			return -1;
		}
		b[i] = 0.;
		for (int k = 0; k < n; k+=2){
			b[i]+=Aij(i,k);
		} 
	}
	return 0;
}
