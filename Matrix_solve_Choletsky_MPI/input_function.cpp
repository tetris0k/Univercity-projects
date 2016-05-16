#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define EPS 10E-15

double Aij(int i, int j){
    return 1./(i+j+1.);
}
	
int func_input(double* A, double* b, int n, int p, int my_rank){
	int t = (int)(n/p) + 1;
    for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
            if (j % p == my_rank){
                A[i * t + (int)(j/p)] = Aij(i,j);    
            }
		}
		if (my_rank == i % p){
			if ( fabs(A[i * t + (int)(i/p)]) < EPS ){
				printf("There are zero elements on diagonal.\n");
				MPI_Abort (MPI_COMM_WORLD, 1);
				return -1;
			}
		}	
		b[i] = 0.;
		for (int k = 0; k < n; k+=2){
			b[i] += Aij(i,k);
		} 
	}
	return 0;
}
