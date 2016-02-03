#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "functions.h"
#define EPS 10e-25

int main(int argc, char *argv[]){
	int q = 0;
	int n = 1;
	int max = 10;
	double t1;
	double *A;
	double *b;
	double *A1;
	FILE* f;
	int ret = 0;
	int t = 2;
	opterr = 0;
	while ((ret = getopt( argc, argv, "f::g:N:t:")) != -1) {
		switch (ret) {
		case 'N':
			if (optarg == NULL){
				perror("Please input max output size or do not use this option.\n");
				return -1;
			}
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
			if (q == 1) break;
			printf("Input from function in file input_function.cpp\n");
			q = 2;
			if (optarg == NULL){
				perror("Please input size of matrix to be generated.\n");
				return -1;
			}
			n = atoi(optarg);
			printf("Size of matrix: %d\n", n);
			break;
		case 't':
			if (optarg == NULL){
				perror("Please input count of threads, or do not use this option (in this case count of threads will instantly be equal to 2).\n");
				return -1;
			}
			t = atoi(optarg);
			break;
		case '?':
			printf("Unknown option -%c\n", optopt);
			printf("Supported options:\n");
			printf("-f -- Input matrix from the argument file, or, if there is no argument, input matrix from the file input.txt;\n");
			printf("-g -- Input matrix from the function. Argument of this option means size of matrix and it is required;\n");
			printf("-N -- Maximum size of the solution, displayed on the screen. Argument is required.\n");
			printf("-t -- Count of threads. Initial one is 2. \n");
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
		printf("Need to run the program with boot options\n"); 
		printf("Supported options:\n");
		printf("-f -- Input matrix from the argument file, or, if there is no argument, input matrix from the file input.txt;\n");
		printf("-g -- Input matrix from the function. Argument of this option means size of the matrix and it's required;\n");
		printf("-N -- Maximum size of the solution, displayed on the screen. Argument is required.\n");
		printf("-t -- Count of threads. Initial one is 2. \n");	
		printf ("\n\n");
		return -1;
	}
	FILE *g;
	g = fopen("output.txt", "wr");
	if(g == NULL){
		printf("Can't open output file.\n");
		return -1;
	}
	printf("Count of threads: %d\n", t);
	PrintMatrix(n,A,b,max);
	printf("\n");
	pthread_t *threads = new pthread_t [t];
	double *Y = new double [n];
	if(Y == NULL){
		printf("Can't allocate memory.\n");
		return -1;
	}
	double *X = new double [n];
	if(X == NULL){
		printf("Can't allocate memory.\n");
		return -1;
	}	
	double *d = new double [n];
	if(d == NULL){
		printf("Can't allocate memory.\n");
		return -1;
	}	
	ARGS *args = new ARGS [t];
	if(args == NULL){
		printf("Can't allocate memory.\n");
		return -1;
	}
	double * bb = new double [n];
	for (int i = 0; i < n; i++){
		bb[i] = b[i];
	}
	t1 = get_full_time();

	Solve(A,b,n,Y,t,threads,X,d,args,bb);

	t1 = get_full_time() - t1;

	file_output(g,n,max,Y);
	printf("\n");
	if (t1 < EPS)
		t1 =  0.000001;
	printf("Total time of solving = %f\n\n", t1);
	if (q == 1){
		printf("Residual norm = %e\n", residual_norm1(n,A1,b,Y,t, threads, args));
		delete []A1;
	} else if(q == 2){
		printf("Residual norm = %e\n", residual_norm(n,b,Y,t, threads, args));
		printf("Error rate = %e\n", error_rate(n,Y));
	}
	delete []X;
	delete []d;
	delete []args;
	delete []threads;
	delete []A;
	delete []b;
	delete []Y;
	delete []bb;
	fclose(g);
	return 0;
}			

