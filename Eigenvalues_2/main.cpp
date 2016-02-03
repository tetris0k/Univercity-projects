#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

int main(int argc, char *argv[]){
    int mode=0;
	int n = 1;
   	int MAX_OUTPUT_SIZE = 5;
	double *A;
	double *values;
	double eps = 1e-10;
	double t1,t2,p;
	double inv1 = 0.;
	double inv2 = 0.;
	FILE* input;
	int ret = 0;
	opterr = 0;
	while ((ret = getopt( argc, argv, "f::g:N:e:")) != -1) {
		switch (ret) {
		case 'e':
			eps = atof(optarg);
			break;
		case 'N':
			MAX_OUTPUT_SIZE = atoi(optarg);
			printf("Максимальный выходной размер: %d\n", MAX_OUTPUT_SIZE);
			break;
		case 'f':
			mode = 1;
			if (optarg != NULL){
				printf ("Ввод из файла %s\n", optarg);
				if ((input = fopen (optarg, "r")) == NULL) {
					perror ("Невозможно открыть файл\n");
					return -1;
				}
			} else {
				printf ("Ввод из файла input.txt\n");
				if ((input = fopen ("input.txt", "r")) == NULL) {
					perror ("Невозможно открыть файл\n");
					return -1;
				}
			}
			break;
		case 'g':
			if (mode == 1){
				break;
			}
			printf("Ввод из функции\n");
			mode = 2;
			n = atoi(optarg);
			printf("Матрица размера: %d\n", n);
			break;
		case '?':
			printf("Неизвестная опция -%c\n", optopt);
			printf("Поддерживаемые опции:\n");
			printf("f - ввод матрицы из файла, являющегося аргументом, или, если аргумента нет, ввод матрицы из файла input.txt;\n");
			printf("g - ввод матрицы из функции. Аргумент функции (обязателен) является размером матрицы;\n");
			printf("N - опция, аргумент которой (обязателен) задает максимальный выходной размер.\n");
			printf("e - точность. Аргумент обязателен.\n");
			return -1;
		}
	}
	if (mode == 1){
		if(fscanf (input, "%d", &n) == 0){
			printf ("не получилось сосканировать размер матрицы из файла\n");
			return -1;
		}
		A = new double [n*n];
		if (A == NULL){
			printf("не удалось выделить память под матрицу А\n");
			return -1;
		}
		if (file_input(input, n, A) != 0){
			return -1;
		}
		fclose(input);
	} else if (mode == 2){
		A = new double [n*n];
		if (A == NULL){
			printf("не удалось выделить память под матрицу А\n");
			return -1;
		}				
		if (func_input(A, n) != 0){
			return -1;
		}
	} else {
		printf("Требуется запустить программу с какими-либо параметрами\n");
		printf("Поддерживаемые параметры запуска:\n");
		printf("f - ввод матрицы из файла, являющегося аргументом, или, если аргумента нет, ввод матрицы из файла input.txt;\n");
		printf("g - ввод матрицы из функции. Аргумент функции (обязателен) является размером матрицы;\n");
		printf("N - опция, аргумент которой (обязателен) задает максимальный выходной размер.\n");
		printf("e - точность. Аргумент обязателен.\n");
		return -1;
	}
	for (int i = 0; i < n; i++){
		inv1 += A[i*n+i];
		for (int j = 0; j < n; j++){
			inv2 += A[i*n+j]*A[i*n+j];
		}
	}
	inv2 = sqrt(inv2);

   	PrintMatrix(n, A, MAX_OUTPUT_SIZE);
	printf("\n");
	FILE *out;
	out = fopen("output.txt", "wr");
	if(out == NULL){
		printf("не могу открыть выходной файл\n");
		return -1;
	}

	values = new double [n];
	if (values==NULL){
		printf("Не удалoсь выделить память под ответ\n");
		delete [] A;
		return -1;
	}

	
	double* sin = new double [n];
	if (sin==NULL){
		printf("Не удалoсь выделить память под синусы\n");
		delete [] A;
		delete [] values;
		return -1;
	}
	double* cos = new double [n];	
	if (cos==NULL){
		printf("Не удалoсь выделить память под косинусы\n");
		delete [] A;
		delete [] values;
		delete [] sin;
		return -1;
	}


	t1 = get_time();
	FindValues(n, A, values, eps, sin, cos);
	t2 = get_time();

	printf("Решение:\n");
	PrintVector(n, values, MAX_OUTPUT_SIZE);
	
	p = 0.;
	for (int i = 0; i < n; i++){
		inv1 -= values[i];
		p += values[i]*values[i];
	}
	inv2 -= sqrt(p);
	printf("\n");
		
	file_output(out,n, values);
	printf("время = %f\n", t2-t1);
	printf("невязка следа матрицы = %e\n", inv1);
	printf("норма вектора собственных значений = %g\n", inv2);
	delete []A;
	delete []values;
	delete []cos;
	delete []sin;
	fclose(out);
	return 0;
}


double get_time(){
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_usec/1000000.0 + (double)tv.tv_sec;
}
