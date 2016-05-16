#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "functions.h"
#define EPS 10E-10

int file_input(FILE* f, int n, double* A, double* A1, double* b, int p, int my_rank, double* c){
    int t = (int)(n/p) + 1;
    MPI_Status status;
	for (int i = 0; i < n; i++){
        int ip = (int)(i/p);
        if (my_rank == 0){
            for (int j = 0 ; j < n ; j++){
                if (fscanf(f, "%lf", &c[j]) == 0){
				    printf("Problem with elements in matrix.\n");
				    return -1;
                }
            }
            if ( fabs(c[i]) < EPS ){
                printf("There are zero elements on diagonal.\n");
                return -1;
            }
        }
        if (i % p == 0){
            if (my_rank == 0){
                for (int j = 0; j < n; j++){ 
                    A[j*t + ip] = c[j];
                    A1[j*t + ip] = c[j];
                }
            }
        } else {
			MPI_Bcast(c, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Sendrecv_replace(c, n, MPI_DOUBLE, i % p, i, 0, i, MPI_COMM_WORLD, &status);
            if (my_rank == i % p){
                for (int j = 0; j < n; j++){
                    A[j * t + ip] = c[j];
                    A1[j * t + ip] = c[j];
                }
            }
        }
		if (my_rank == 0){
            if(fscanf(f, "%lf", &b[i]) == 0){
                printf("Problem with elements in matrix.\n");
                return -1;
            }
        }
        MPI_Bcast(&b[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	return 0;
}


void file_output(FILE* g, int n, int max, double* Y ){
    printf("Solution:\n");
    for (int i = 0; i < n; i++){
        fprintf(g, "%f ", Y[i]);
        if (i < max){
            printf("%f ", Y[i]);
        }
    }
    printf("\n");
}

void PrintMatrix(int n, double * A, double * b, int max_size, int p, int my_rank, double* c){
    int i;
    int j;
    int nPrint;
    if (n > max_size){
	nPrint = max_size;
    } else {
	nPrint = n;
    }
    int t = (int)(n/p) + 1;
   // MPI_Status status; 
    for (i = 0; i < nPrint; ++i) {
        if (my_rank == 0)
            printf("| ");
        if (my_rank == i % p){
            for (j = 0; j < n; j++)
                c[j] = A[j * t + (int)(i/p)];
        }
        MPI_Bcast(c, n, MPI_DOUBLE, i%p, MPI_COMM_WORLD);
      //  MPI_Sendrecv_replace(c, n, MPI_DOUBLE, 0 , 3, i % p, 3, MPI_COMM_WORLD, &status);
        if (my_rank == 0){
            for (j = 0; j < nPrint; ++j){
                printf("%5.3g ", c[j]);
            }
	    if (nPrint < n)
		printf("\t...\t");
	    printf("|");
	    printf("%5.3g\n", b[i]);
	}
	MPI_Barrier (MPI_COMM_WORLD);    
    }
}
