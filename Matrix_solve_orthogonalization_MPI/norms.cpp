#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

double nev(int n, double *A, double *x, int count, int my_rank){
	double m,l = 0.;
	double sum = 0.;
	int w,j;
	for (w = 0; w < n; w++){
		m = 0.;
		if (my_rank == w % count){                     					//каждый процесс считает свой кусочек невязки
			for (j = 0; j < n; j++){
				m += A[(n+1)*(int)(w/count) + j]*x[j];
			}
			m += A[(n+1)*(int)(w/count) + n];
			l += m*m;
		}	
	}
	MPI_Reduce(&l, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);    //суммируем все кусочки из всех процессов
	return sqrt(sum);                                                   //и пихаем в 0 процесс
}

double error_rate(int n, double* x){
	double l = 0.;
	for (int i = 0; i < n; i += 2){
		l += (x[i]-1)*(x[i]-1);
	}
	for (int i = 1; i < n; i += 2){
		l += x[i]*x[i];
	}
	return sqrt(l);
}

