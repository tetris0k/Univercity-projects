#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    int mode=0;
	int n = 1;
	int max = 10;
	double *A;
    double *b;
	double t1,t2;
	FILE* input;
	int ret = 0;
	opterr = 0;
	while ((ret = getopt( argc, argv, "f::g:N:")) != -1) {
		switch (ret) {
		case 'N':
			max = atoi(optarg);
			printf("Максимальный выходной размер: %d\n", max);
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
			if (mode == 1) break;
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
			return -1;
		}
	}
	if (mode == 1){
		if(fscanf (input, "%d", &n) == 0){
			printf ("не получилось сосканировать размер матрицы из файла\n");	
			return -1;
		}
		A = new double [n*(n+1)];
		if (A == NULL){
			printf("не удалось выделить память под матрицу А\n");
			return -1;
		}
		if (file_input(input, n, A) != 0){
			return -1;
		}			
		fclose(input);
	} else if (mode == 2){
		A = new double [n*(n+1)];
		if (func_input(A, n) != 0){
			return -1;
		}
	} else {
		printf("Требуется запустить программу с какими-либо параметрами\n");
		printf("Поддерживаемые параметры запуска:\n");
		printf("f - ввод матрицы из файла, являющегося аргументом, или, если аргумента нет, ввод матрицы из файла input.txt;\n");
		printf("g - ввод матрицы из функции. Аргумент функции (обязателен) является размером матрицы;\n");
		printf("N - опция, аргумент которой (обязателен) задает максимальный выходной размер.\n");
		return -1;
	}
	b = new double [n];
	if (b == NULL){
		printf("не удалось выделить память под вектор b\n");
		return -1;
	}
    for (int z=0;z<n;z++){
        b[z] = -A[z*(n+1)+n];
    }
   	int MAX_OUTPUT_SIZE = 10;

   	PrintMatrix(n, A, b, MAX_OUTPUT_SIZE);
	delete []b;
	
	FILE *out;
	out = fopen("output.txt", "wr");
	if(out == NULL){
		printf("не могу открыть выходной файл\n");
		return -1;
	}
	
	double *x = new double [n];
	if (x == NULL){
		printf("не удалось выделить память под вектор x\n");
		return -1;
	}
	t1 = get_full_time();
	Solve(n,A,x);
	t2 = get_full_time();
	file_output(out,n,max,x);
	printf("время = %f\n", t2-t1);
	printf("невязка = %e\n", nev(n,A,x));
	if (mode == 2){
		printf("норма погрешности = %e\n", error_rate(n,x));
	}	
	delete []A;
	delete []x;
	fclose(out); 
	return 0;
}


