#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

int file_input(FILE* f, int n, double* A){
	for (int o = 0; o < n; o++){
		for (int p = 0; p < n; p++){
			if(fscanf(f, "%lf", &A[o*n+p]) == 0){
				printf("что-то не так с одним из элементов матрицы\n");
				return -1;
			}
		}		
	}
	return 0;
}
	


void file_output(FILE* g, int n, double* x){
	for (int i = 0; i < n; i++){
		fprintf(g, "%f ", x[i]);
	}
}
