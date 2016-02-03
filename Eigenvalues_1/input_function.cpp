#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 10E-35

double Aij(int i, int j){
	return 1./(i+j+1.);
}
	
int func_input(double* A, int n){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			A[i*n + j] = Aij(i,j);
			if (fabs(Aij(i,j) - Aij(j,i)) > EPS){
				printf("It is not diagonal matrix.\n");
				return -1;
			}
		}
	}
	return 0;
}
