#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

int file_input(FILE* f, int n, double* A, int my_rank, int count, double* c){
	for (int o = 0; o < n; o++){                                            //Храним А по строкам в процессах
		if (my_rank == 0){
			for (int p = 0; p < n+1; p++){                                  //Считываем информацию из файла внутри 0 процесса
				if(fscanf(f, "%lf", &c[p]) == 0){							//в промежуточный массив c
					printf("что-то не так с одним из элементов матрицы\n");
					MPI_Abort (MPI_COMM_WORLD, 1);
					return -1;
				}
			}		
			c[n] *= -1;			                                            //Умножаем последний элемент на -1, так как это b[i], 
		}																    //и по алгоритму нам надо его хранить в А со знаком "-"
		if (o % count == 0){                                                //Определяем, в каком процессе должна храниться текущая строка
			if (my_rank == 0){											    //Если она должна остаться в 0, 
				for (int j = 0; j < n+1; j++){ 							    //то просто переписываем элементы из c в А внутри 0 процесса
				    A[(n+1)*(int)(o/count) + j] = c[j];
				}
			}
		} else {                                                            //Если она должна быть в другом процессе, 
			MPI_Bcast(c, n + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);			    //пересылаем ее туда (пересылается всем процессам,
			if (my_rank == o % count){									    //потому что так проще, а на скорость работы это не влияет)
			    for (int j = 0; j < n+1; j++){							    //И внутри нужного процесса переписываем из c в A
					A[(n+1)*(int)(o/count) + j] = c[j];
			    }
			}
		}
		MPI_Barrier (MPI_COMM_WORLD);                                       //Синхронизация
	}
	return 0;
}
	


void file_output(FILE* g, int n, int max, double* x){
	printf("решение:\n");
	for (int i = 0; i < n; i++){
		fprintf(g, "%f ", x[i]);
		if (i < max){
			printf("%f ", x[i]);
		}
	}
	printf("\n");
}

void PrintMatrix(int n, double * A, int max_size, int my_rank, int count, double* c){ 
    int i;                       
    int j;                                                         //все печатает только процесс 0
    int nPrint;													   //поэтому будем пересылать все ему
	if (n > max_size){
		nPrint = max_size;
	} else {
		nPrint = n;
	} 
    for (i = 0; i < nPrint; ++i) {
        if (my_rank == 0)
            printf("| ");
        if (my_rank == i % count){                                
            for (j = 0; j < n+1; j++)                               //считываем внутри k процесса нужную для вывода 
                c[j] = A[(n+1)*(int)(i/count) + j];                 //строку в промежуточный массив С                
        }
		MPI_Bcast (c, n+1, MPI_DOUBLE, i % count, MPI_COMM_WORLD);  //отправляем массив С в процесс 0 (опять
        if (my_rank == 0){											//отправляем всем, потому что так проще,  
            for (j = 0; j < nPrint; ++j){							//а на скорость не влияет)
                printf("%5.3g ", c[j]);
            }
        	if (nPrint < n)
           		printf("\t...\t");									//0 процесс все печатает
        	printf("|");
        	printf("%5.3g\n", -c[n]);                               //последний элемент в строке печатает со знаком "-",
		}															//потому что это b[i], и в матрице он хранился с "-"
		MPI_Barrier (MPI_COMM_WORLD);    							//синхронизация
	}
}
