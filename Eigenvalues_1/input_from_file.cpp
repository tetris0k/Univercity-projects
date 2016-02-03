#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"


int file_input(FILE* f, int n, double* A){
	for (int i = 0; i < n; i++){
		for (int j = 0 ; j < n ; j++){
			if (fscanf(f, "%lf", &A[i*n+j]) == 0){
				printf("Problem with elements in matrix.\n");
				return -1;
			}
		}
	}
	return 0;
}


void file_output(FILE* g, int n, int max, double* Y){
	printf("Solution:\n");
	for (int i = 0; i < n; i++){
		fprintf(g, "%3.8g ", Y[i]);
		if (i < max){
			printf("%3.8g ", Y[i]);
		}
	}
	printf("\n");
}
void PrintVector(int n, double* x, int MAX_OUTPUT_SIZE){
	int i;
	int nPrint;
	nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;
	for (i = 0; i < nPrint; ++i)
		printf("%10.3g ", x[i]);
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
