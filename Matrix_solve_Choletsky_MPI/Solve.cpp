#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std;
#include "functions.h"

void Solve(double* A , double* b, int n, double* x, double* c, double* d, int p, int my_rank){
    double l,tmp;
    int np = (int)(n/p) + 1;
    int i,j,k,ip,jp;
    if (my_rank == 0){
        d[0] = A[0] / fabs(A[0]);
        A[0] = sqrt(fabs(A[0]));
        tmp = A[0];
    }
    MPI_Bcast (&d[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&tmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(i = 1; i < n; i++){
        if (my_rank == i % p)
            A[(int)(i/p)] /= tmp * d[0];
	}

	for(i = 1; i < n; i++){
		ip = (int)(i/p);
	
		if (my_rank == i % p){
			l = A[i * np + ip];              //dii:
			for (k = 0; k < i; k++){
				l -= A[k*np + ip] * A[k*np + ip] * d[k];
			}
			d[i] = l/fabs(l);
			A[i*np + ip] = sqrt(fabs(l));        //done
			for (j = 0; j < i+1; j++){
				c[j] = A[j*np + ip];
			}
		}
		
		MPI_Bcast(&d[i], 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD);
		MPI_Bcast(c, n, MPI_DOUBLE, i%p, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		for(j = i+1; j < n; j++){
			if (my_rank == j % p){
				jp = (int)(j/p);
				l = 0.;                            //rij:
				for (k = 0; k < i; k++){
					l -= c[k] * A[k*np + jp] * d[k];
				}
				l += A[i*np + jp];
				l /= c[i] * d[i];
				A[i*np + jp] = l;                  //done
			}
		}
		MPI_Barrier (MPI_COMM_WORLD);
	}
	MPI_Barrier (MPI_COMM_WORLD);

	if (my_rank == 0){              //R_Solve
		x[0] = b[0]/A[0];
	}
	MPI_Bcast(&x[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	for (k = 1; k < n; k++){
		jp = (int)(k/p);
		x[k] = b[k];
		if (my_rank == k % p){
			for (i = 0; i <	k; i++){
				x[k] -= x[i] * A[i*np + jp];
			}
			x[k] /= A[k*np + jp];
		}
		MPI_Bcast(&x[k], 1, MPI_DOUBLE, k%p, MPI_COMM_WORLD);			
	}
	MPI_Barrier (MPI_COMM_WORLD);

	for (i = 0; i < n; i++){   //D_Solve
		x[i] *= d[i];
	}
	MPI_Barrier (MPI_COMM_WORLD);

	if (my_rank == (n - 1) % p){
		x[n-1] /= A[(n - 1) * np + (int)((n - 1)/p)];              //RSolve
		for (i = n - 2; i >= 0; i--)
			x[i] -= A[i * np + (int)((n - 1)/p)] * x[n-1];
	}
	MPI_Bcast (x, n, MPI_DOUBLE, (n - 1) % p, MPI_COMM_WORLD);
	MPI_Barrier (MPI_COMM_WORLD);
	for (i = n - 2;	i >= 0; i--){
		if (my_rank == i % p){
			x[i] /= A[i * np + (int)(i/p)];
			for (j = i-1; j >= 0; j--)
				x[j] -= A[j * np + (int)(i/p)] * x[i]; 
		}
		MPI_Bcast (x, n, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
		MPI_Barrier (MPI_COMM_WORLD);
	}
}
