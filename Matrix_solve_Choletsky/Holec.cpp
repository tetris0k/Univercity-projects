#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "functions.h"
double get_time();


int main(int argc, char *argv[]){
	int q = 0;
	int n = 1;
	int max = 10;
	double t1,t2;
	int max_matrix_size = 10;	
	double *A;
	double *b;
	double *A1;
	FILE* f;
	int ret = 0;
	opterr = 0;
	while ((ret = getopt( argc, argv, "f::g:N:")) != -1) {
		switch (ret) {
		case 'N':
			max = atoi(optarg);
			printf("Max output size: %d\n", max);
			break;
		case 'f':
			q = 1;
			if (optarg != NULL){
				printf ("Input from file %s\n", optarg);
				if ((f = fopen (optarg, "r")) == NULL) {
					perror ("Can't open file.\n");
					return -1;
				}
			} else {
				printf ("Input from file input.txt\n");
				if ((f = fopen ("input.txt", "r")) == NULL) {
					perror ("Can't open file.\n");
					return -1;
				}
			}
			break;
		case 'g':
			printf("Input from function in file input_function.cpp\n");
			q = 2;
			n = atoi(optarg);
			printf("Size of matrix: %d\n", n);
			break;
		case '?':
			printf("Unknown option -%c\n", optopt);
			printf("Supported options:\n");
			printf("-f -- Input matrix from the argument file, or, if there is no argument, input matrix from the file input.txt;\n");
			printf("-g -- Input matrix from the function. Argument of this option means size of matrix and it is required;\n");
			printf("-N -- Maximum size of the solution, displayed on the screen. Argument is required.\n");
			return -1;
		}
	}
	if (q == 1){
		if(fscanf (f, "%d", &n) == 0){
			printf ("Problem with size of matrix.\n");
			return -1;
		}
		A = new double [n*(n+1)];
		b = new double [n];
		A1 = new double [n*n];
		if (A == NULL || b == NULL || A1 == NULL){
			printf ("Problem with memory allocation.\n");
			return -1;
		}
		if(file_input(f, n, A, A1, b) != 0){
			return -1;
		}
		fclose(f);
	} else if (q == 2){
		A = new double [n*(n+1)];
		b = new double [n];
		if (A == NULL || b == NULL){
			printf ("Problem with memory allocation.\n");
			return -1;
		}
		if (func_input(A, b, n) != 0){
			return -1;
		}
	} else {
		printf("It's recommended to run the program with boot options\n"); 
		printf("Supported options:\n");
		printf("-f -- Input matrix from the argument file, or, if there is no argument, input matrix from the file input.txt;\n");
		printf("-g -- Input matrix from the function. Argument of this option means size of the matrix and it's required;\n");
		printf("-N -- Maximum size of the solution, displayed on the screen. Argument is required.\n");
		printf ("\n\n");
		printf("Choose, where to take the matrix from:\n");	
		printf("1 -- from file input.txt\n");
		printf("2 -- from input function in file input_function.cpp\n");
		scanf("%d", &q);
		if (q == 1){
			f = fopen("input.txt","r");
			if(f == NULL){
				printf("Can't open input file.\n");
				return -1;
			}
			if(fscanf (f, "%d", &n) == 0){
				printf("Problem with size of matrix.\n");
				return -1;
			}
			A = new double [n*(n+1)/2];	
			b = new double [n];
			A1 = new double [n*n];
			if (A == NULL || b == NULL || A1 == NULL){
				printf ("Problem with memory allocation.\n");
				return -1;
			}
			if(file_input(f, n, A,A1, b) != 0){
				return -1;
			}
			fclose(f);
		} else if (q == 2){
			printf("Size of matrix: ");
			scanf("%d", &n);
			A = new double [n*(n+1)/2];
			b = new double [n];
			if (A == NULL || b == NULL){
				printf ("Problem with memory allocation.\n");
				return -1;
			}
			if (func_input(A, b, n) != 0){
				return -1;
			}	
		} else {	
			printf("Incorrect input way.\n");
			return -1;
		}
	}	
	FILE *g;
	g = fopen("output.txt", "wr");
	if(g == NULL){
		printf("Can't open output file.\n");
		return -1;
	}

	PrintMatrix(n,A,b,max_matrix_size);

	double *Y = new double [n];
	t1 = get_time();

	Solve(A,b,n,Y);

	t2 = get_time();

	file_output(g,n,max,Y);
	printf("Time = %f\n", t2-t1);
	if (q == 1){
		printf("Residual norm = %e\n", residual_norm1(n,A1,b,Y));
		delete []A1;
	} else if(q == 2){
		printf("Residual norm = %e\n", residual_norm(n,b,Y));
		printf("Error rate = %e\n", error_rate(n,Y));
	}
	delete []A;
	delete []b;
	delete []Y;
	fclose(g);
	return 0;
}			
	
double get_time(){
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_usec/1000000.0 + (double)tv.tv_sec;
}
