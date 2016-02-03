#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "functions.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    int mode=0;
	int n = 1;
	int t = 3;
	int max = 10;
	double *A;
    double *b;
	double t1;
	FILE* input;
	int ret = 0;
	opterr = 0;
	while ((ret = getopt( argc, argv, "f::g:N:t:")) != -1) {
		switch (ret) {
		case 'N':
			if (optarg == NULL) return -1;
			max = atoi(optarg);
			printf("максимальный выходной размер: %d\n", max);
			break;
		case 't':
			if (optarg == NULL) return -1;
			t = atoi(optarg);
			break;		
		case 'f':
			mode = 1;
			if (optarg != NULL){
				printf ("ввод из файла %s\n", optarg);
				if ((input = fopen (optarg, "r")) == NULL) {
					perror ("невозможно открыть файл\n");
					return -1;
				}
			} else {
				printf ("ввод из файла input.txt\n");
				if ((input = fopen ("input.txt", "r")) == NULL) {
					perror ("невозможно открыть файл\n");
					return -1;
				}
			}
			break;
		case 'g':
			if (mode == 1) break;
			if (optarg == NULL) return -1;
			printf("ввод из функции\n");
			mode = 2;
			n = atoi(optarg);
			printf("матрица размера: %d\n", n);
			break;
		case '?':
			printf("неизвестная опция -%c\n", optopt);
			printf("поддерживаемые опции:\n");
			printf("f - ввод матрицы из файла, являющегося аргументом, или, если аргумента нет, ввод матрицы из файла input.txt;\n");
			printf("g - ввод матрицы из функции. Аргумент функции (обязателен) является размером матрицы;\n");
			printf("t - количество нитей;\n");
			printf("N - максимальный выходной размер.\n");
			return -1;
		}
	}
	printf("Используется %d нитей.\n", t); 	
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
		printf("t - количество нитей;\n");
		printf("N - опция, аргумент которой (обязателен) задает максимальный выходной размер.\n");
		return -1;
	}
	pthread_t *threads = new pthread_t [t];                     //инициализируем нити, выделяем под них память
	if (threads == NULL){
		printf("не удалось выделить память под массив нитей\n");
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
	printf("\n");	
	
	FILE *out;
	out = fopen("output.txt", "wr");
	if(out == NULL){
		printf("не могу открыть выходной файл\n");
		return -1;
	}
	
	double *x = new double [n];                    //сюда пишем решение
	if (x == NULL){
		printf("не удалось выделить память под вектор x\n");
		return -1;
	}
	double *E = new double [(n+1)*(n+1)];          //здесь будут храниться базисы
	if (E == NULL){
		printf("не удалось выделить память под матрицу Е\n");
		return -1;
	}
	solve *args = new solve [t];                   //массив структур для нитей
	if (args == NULL){
		printf("не удалось выделить память под аргументы нитей\n");
		return -1;
	}
	double *v = new double [n+1];                  //сюда запоминается вектор для master
	if (v == NULL){
		printf("не удалось выделить память под вектор v\n");
		return -1;
	}
	printf("\nрешение системы...\n");
	t1 = get_full_time();                         //get_full_time для счета астрономического времени
	Solve(n,A,x,t,threads,E,args,v);
	t1 = get_full_time() - t1;
	
	Nevyaska * arg = new Nevyaska [t];
	if (args == NULL){
		printf("не удалось выделить память под структуру значений для невязки\n");
		return -1;
	}
	printf("\n");
	file_output(out,n,max,x);
	printf("время решения системы = %f\n", t1);
	printf("\nподсчет невязки...\n");
	printf("невязка = %e\n", nev(n,A,x,t,threads,arg));
	if (mode == 2){
		printf("норма погрешности = %e\n", error_rate(n,x));
	}	
	delete []A;
	delete []x;
	delete []E;
	delete []v;
	delete []args;
	delete []threads;
	delete []arg;
	fclose(out); 
	return 0;
}


