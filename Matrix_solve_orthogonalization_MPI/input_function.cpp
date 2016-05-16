#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

double Aij(int i, int j){
    return fabs(i-j) + 1.;
}
	
int func_input(double* A, int n, int count, int my_rank){
	for (int i = 0; i < n; i++){
		if (my_rank == i % count){                     		//храним А в процессах по строкам
			for (int j = 0; j < n; j++) 
				A[(int)(i/count)*(n+1) + j] = Aij(i,j);
			for (int k = 0; k < n; k+=2){
				A[(int)(i/count)*(n+1) + n]-=Aij(i,k);  	//храним b[i] последними элементами в строках А,
			}												//со знаком "-"
		}
	}
	return 0;
}
