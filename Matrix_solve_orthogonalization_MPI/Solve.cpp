#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include "math.h"
#define EPS 10E-10

void Solve(int n, double* A, double* x, int my_rank, int count, double* w){	
	int i,j,k;
	double smult_1, smult_2;						//матрица Е храница по столбцам в процессах
	int nc = (int)((n+1)/count) + 1;               	//количество столбцов в E внутри каждого процесса
	double *E = new double [(n+1)*nc];				//равно = цел_часть((n+1)/кол-во_стоблцов) + 1
	double *v = new double [n+1];
	for (i = 0; i < n+1; i++){
		for (j = 0; j < n+1; j++){
			if (my_rank == j % count){
				if (i==j) E[i*nc + (int)(j/count)]=1; else E[i*nc + (int)(j/count)]=0;	//Е изначально единичная
			}		
		}
	}
	for (k=1;k<n+1;k++){
		smult_2 = 0.;
		if (my_rank == (k-1) % count){
			for (j = 0; j < n+1; j++){							//вычисляем сначала все нужное в определенном процессе
				v[j] = E[j * nc + (int)((k-1)/count)];       	//это nowv
				w[j] = A[(n+1) * (int)((k-1)/count) + j];	 	//это вектор из A, который нужен в smult
				smult_2 += v[j] * w[j];			             	//сразу считаем smult2 внутри процесса, где все хранится
			}
			if (fabs(smult_2) < EPS){
				printf("Невозможно решить систему. появилось деление на 0\n");
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		MPI_Bcast(v, n+1, MPI_DOUBLE, (k-1) % count, MPI_COMM_WORLD);         	//рассылаем всем процессам nowv,
		MPI_Bcast(w, n+1, MPI_DOUBLE, (k-1) % count, MPI_COMM_WORLD);         	//Вектор из A, который нужен в smult,
		MPI_Bcast(&smult_2 , 1, MPI_DOUBLE, (k-1) % count, MPI_COMM_WORLD); 	//и smult_2
		MPI_Barrier(MPI_COMM_WORLD);											//синхронизируем				
		
	    for(i = k; i < n+1; i++){
			smult_1 = 0.;
			if (my_rank == i % count){
				smult_1 = smult(i,k,w,E,count,nc) / smult_2;          			//вычисляем коэфициент a,
				for (j = 0; j < n+1; j++){										//равный отношению скалярных произведений
					E[j * nc + (int)(i/count)] -= v[j] * smult_1;     			//изменяем E
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);											//синхронизируем
	}
	if (my_rank == n % count){											//получили решение системы в последнем столбце E
		for (int z=0;z<n;z++){                                         	//то есть в процессе номер ост_от_дел(n на кол-во_процессов)
			x[z] = E[z*nc+(int)(n/count)];
		}
	}
	MPI_Bcast(x, n, MPI_DOUBLE, n % count, MPI_COMM_WORLD);            	//отправляем решение во все процессы
	delete []E;
	delete []v;
}
/*			//Старые функции (для наглядности сходства)
void master(int n, double *v, int i , int k, double *A, double *E){
	double x,y;
	x = smult(n,i,k,A,E);
	y = smult2(n,v,k,A);
	double a;
	a = x/y;
	for (int j = 0; j < n+1; j++){
		E[j*(n+1)+i] -= a*v[j];
	}
}

double smult2(int n, double *v , int k, double *A){
	double res=0.;
	for (int j = 0; j < n+1; j++){
		res = res + A[(k-1)*(n+1)+j]*v[j];
	}	    
	return res;
}
*/
double smult( int i , int k, double *w, double *E, int count, int nc){
	double res = 0.;
	for (int j = 0; j < k; j++){
		res = res + w[j]*E[j*nc + (int)(i/count)];
	}	
	res = res + w[i]*E[i*nc + (int)(i/count)];    
	return res;
}
/*
void nowv(int n, double *E, double *v, int k){
	for(int i=0;i<n+1;i++){
		v[i]=E[i*(n+1)+k-1];
	}
}
*/

