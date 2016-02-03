#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 10e-15
double Aij(int i, int j){
	return fabs(i-j)+1.;
}
	
int func_input(double* A, int n){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			A[i*n+ j] = Aij(i,j);
		}
	}
	return 0;
}
