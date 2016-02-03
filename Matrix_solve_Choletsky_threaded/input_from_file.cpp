#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define EPS 10E-10

int file_input(FILE* f, int n, double* A, double* A1, double* b){
	for (int i = 0; i < n; i++){
		for (int j = 0 ; j < i ; j++){
			if (fscanf(f, "%lf", &A1[i*n+j]) == 0){
				printf("Problem with elements in matrix.\n");
				return -1;
			}
		}
		for (int j = i; j < n; j++){
			if(fscanf(f, "%lf", &A[i*n+j-i*(i+1)/2]) == 0){
				printf("Problem with elements in matrix.\n");
				return -1;
			}
			A1[i*n+j]=A[i*n+j-i*(i+1)/2];
		}
		if ( fabs(A[i*n + i - i*(i+1)/2]) < EPS ){
			printf("There are zero elements on diagonal.\n");
			return -1;
		}
		if(fscanf(f, "%lf", &b[i]) == 0){
			printf("Problem with elements in matrix.\n");
			return -1;
		}
	}
	return 0;
}


void file_output(FILE* g, int n, int max, double* Y){
	printf("Solution:\n");
	for (int i = 0; i < n; i++){
		fprintf(g, "%f ", Y[i]);
		if (i < max){
			printf("%f ", Y[i]);
		}
	}
	printf("\n");
}

void PrintMatrix(int n, double * A, double * b, int max_size){
    int i;
    int j;
    int nPrint;
	if (n > max_size){
		nPrint = max_size;
	} else {
		nPrint = n;
	}

    for (i = 0; i < nPrint; ++i) 	{
        printf("| ");
        for (j = 0; j < nPrint; ++j)
			if (j<i){
            	printf("%5.3g ", A[j * n + i-j*(j+1)/2]);
			} else {
				printf("%5.3g ", A[i * n + j-i*(i+1)/2]);
			}
        if (nPrint < n)
            printf("\t...\t");
        printf("|");
        printf("%5.3g\n", b[i]);
    }
}
