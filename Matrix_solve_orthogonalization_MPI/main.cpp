#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include <unistd.h>

int main(int argc, char **argv){
    int mode = 0;									//изначально объявляем
	int n = 1;										//все переменные
	int max = 10;
	double *A;
	double *x;
	double *c;
	double t;
	int my_rank, count;
	FILE* input;
	FILE* output;
 
    MPI_Init (&argc, &argv);                      	//инициализируем запуск процессов
    MPI_Comm_size(MPI_COMM_WORLD, &count);        	//считываем их общее количество в count (в каждом процессе)
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);      	//считываем в каждом процессе его номер в my_rank
	MPI_Barrier (MPI_COMM_WORLD);                 	//синхронизация 
	int ret = 0;
	opterr = 0;	

	while ((ret = getopt( argc, argv, "f::g:N:")) != -1) {
		switch (ret) {
		case 'N':
			max = atoi(optarg);
			if (my_rank == 0)
				printf("Максимальный выходной размер: %d\n", max);
			break;
		case 'f':
			mode = 1;
			if (my_rank == 0){
				if (optarg != NULL){                	  			//открываем файлы только в 0 процессе.
					printf ("Ввод из файла %s\n", optarg);			//он все считывает, и раздает остальным
					if ((input = fopen (optarg, "r")) == NULL) {        
						perror ("Невозможно открыть файл\n");
						MPI_Abort (MPI_COMM_WORLD, 1);       		//обрывает работу всех процессов
						return -1;
					}			
				} else {
					printf ("Ввод из файла input.txt\n");
					if ((input = fopen ("input.txt", "r")) == NULL) {
						perror ("Невозможно открыть файл\n");
						MPI_Abort (MPI_COMM_WORLD, 1);
						return -1;
					}
				}
			}
			break;
		case 'g':
			if (mode == 1) break;
			if (my_rank == 0)                   					//печатаем все на экран только в 0 процессе
				printf("Ввод из функции\n");
			mode = 2;
			n = atoi(optarg);
			if (my_rank == 0)
				printf("Матрица размера: %d\n", n);
			break;
		case '?':
			if (my_rank == 0){                            
				printf("Неизвестная опция -%c\n", optopt);
				printf("Поддерживаемые опции:\n");
				printf("f - ввод матрицы из файла, являющегося аргументом, или, если аргумента нет, ввод матрицы из файла input.txt;\n");
				printf("g - ввод матрицы из функции. Аргумент функции (обязателен) является размером матрицы;\n");
				printf("N - опция, аргумент которой (обязателен) задает максимальный выходной размер.\n");
			}
			MPI_Finalize ();	                           			//завершение работы процессов (не аварийное). вызывается во всех процессах
			return -1;
		}
	}
	if (mode == 1){
		if (my_rank == 0){
			if(fscanf (input, "%d", &n) == 0){
				printf ("не получилось сосканировать размер матрицы из файла\n");		
				MPI_Abort (MPI_COMM_WORLD, 1);
				return -1;
			}
		}	
		MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);            	//рассылаем n всем процессам из 0 процесса
		MPI_Barrier (MPI_COMM_WORLD);
		A = new double [(int)((n/count) + 1) * (n+1)];              //матрицу А храним по строкам, длины n+1,
		c = new double [n+1];										//так как на последнем месте строки стоит b[i]
		if (A == NULL || c == NULL ){								//количество столбцов в процессе = цел_часть(n/кол-во_процессов)+1
			printf("не удалось выделить память под матрицу А в процессе %d\n", my_rank);
			MPI_Abort (MPI_COMM_WORLD, 1);
			return -1;
		}
		if (file_input(input, n, A, my_rank, count, c) != 0){       //считывание из файла. вызывается везде, 
			printf("Не получилось считать информацию из файла\n"); 	//0 процесс считывает информацию, и раздает остальным			
			MPI_Abort (MPI_COMM_WORLD, 1);			
			return -1;
		}			
		if (my_rank == 0){		
			fclose(input);
		}
	} else if (mode == 2){
		A = new double [(int)((n/count) + 1) * (n+1)];
		c = new double [n+1];										//массив для промежуточных пересылок
		if (A == NULL || c == NULL){
			printf ("не получилось выделить память под матрицу А в процессе %d\n", my_rank);
			MPI_Abort (MPI_COMM_WORLD, 1);
			return -1;
		}
		if (func_input(A, n, count, my_rank) != 0){                //считывание из функции. вызывается везде. каждый процесс
			printf("не получилось считать матрицу из функции\n");  //имеет доступ к func_input, и сам считывает то, что ему нужно
			MPI_Abort (MPI_COMM_WORLD, 1);
			return -1;
		}
	} else {
		if (my_rank == 0){
			printf("Требуется запустить программу с какими-либо параметрами\n");
			printf("Поддерживаемые параметры запуска:\n");
			printf("f - ввод матрицы из файла, являющегося аргументом, или, если аргумента нет, ввод матрицы из файла input.txt;\n");
			printf("g - ввод матрицы из функции. Аргумент функции (обязателен) является размером матрицы;\n");
			printf("N - опция, аргумент которой (обязателен) задает максимальный выходной размер.\n");
		}
		MPI_Finalize ();
		return -1;
	}

	if (my_rank == 0){                              	//открываем выходной файл в процессе 0
		output = fopen("output.txt", "wr");	
		if(output == NULL){
			printf("не могу открыть выходной файл\n");
			MPI_Abort (MPI_COMM_WORLD, 1);
			return -1;
		}
	}    

   	PrintMatrix(n, A, max, my_rank, count, c);			//вызывается во всех, но в итоге печатает все только 0 процес

	
	x = new double [n];
	if (x == NULL){
		printf("не удалось выделить память под вектор x в процессе %d\n", my_rank);
		MPI_Abort (MPI_COMM_WORLD, 1);
		return -1;
	}
	
	MPI_Barrier (MPI_COMM_WORLD);    							//считаем начальное время (как в Богачеве):
	t = MPI_Wtime();											//синхронизируем все процессы, начинаем отсчет

	Solve(n, A, x, my_rank, count, c);  						//после работы решение будет храниться во всех процессах

	MPI_Barrier (MPI_COMM_WORLD);   							//считаем конечное время (как в Богачеве):
	t = MPI_Wtime() - t;										//ждем, пока все процессы закончат, фиксируем время
	
	if (my_rank == 0){
		file_output(output,n,max,x);  		 					//вывод ответа вызывается в 0 процессе
		printf("время = %f\n", t);
	}
	double nevyas = nev(n,A,x,count,my_rank);  					//вычисление невязки вызывается во всех процессах,
	if (my_rank == 0){											//но ответ будет храниться в процессе 0
		printf("невязка = %e\n", nevyas);
		printf("норма погрешности = %e\n", error_rate(n,x));	//норма погрешности считается в 0 процессе
		fclose(output);
	}	
	delete []A;					//чистим память
	delete []x;
	delete []c;
	MPI_Finalize();				//завершаем работу всех процессов
	return 0;
}


