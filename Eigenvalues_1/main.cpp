#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "functions.h"
double get_time();


int main(int argc, char *argv[]){
	int z,q,n,max=5;
	double t1,t2;
	double inv1, inv2;
	double eps = 10e-10;
	double *A;
	FILE* f;
	int ret = 0;
	opterr = 0;
	q = 0;
	while ((ret = getopt( argc, argv, "f::g:N:e:")) != -1) {
		switch (ret) {
		case 'e':
			eps = atof(optarg);
			break;
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
			if (q == 1){
				break;
			}
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
			printf("-e -- Precision. Argument is required.\n"); 
			return -1;
		}
	}
	if (q == 1){
		if(fscanf (f, "%d", &n) == 0){
			printf("Problem with size of matrix.\n");
			return -1;
		}
		A = new double [n*n];
		if (A==NULL){
			printf ("Can't allocate memory\n");
			return -1;
		}	
		if(file_input(f, n, A) != 0){
			delete [] A;
			return -1;
		}
		fclose(f);
		if (issym(A,n) != 1){
			printf ("Matrix is not symmetric.\n");
			delete [] A;
			return -1;
		}
	} else if (q == 2){
		A = new double [n*n];
		if (A==NULL){
			printf ("Can't allocate memory\n");
			return -1;
		}
		if (func_input(A, n) != 0){
			delete [] A;
			return -1;
		}
	} else {	
		printf("Need to run the program with boot options\n"); 
		printf("Supported options:\n");
		printf("-f -- Input matrix from the argument file, or, if there is no argument, input matrix from the file input.txt;\n");
		printf("-g -- Input matrix from the function. Argument of this option means size of the matrix and it's required;\n");
		printf("-N -- Maximum size of the solution, displayed on the screen. Argument is required.\n");
		printf("-e -- Precision. Argument is required.\n");
		return -1;
	}
	printf("Precision: %e\n", eps);
	FILE *g;
	g = fopen("output.txt", "wr");
	if(g == NULL){	
		printf("Can't open output file.\n");
		delete [] A;
		return -1;
	}
	inv1 = inv_trace(A,n);
	inv2 = inv_length(A,n);
	
	PrintMatrix(n,A,max);
	printf("\n");

	double *Q1 = new double [n-1];
	if (Q1==NULL){
		printf ("Can't allocate memory\n");
		delete [] A;
		fclose(g);
		return -1;
	}

	double *Q2 = new double [n-1];
	if (Q2==NULL){
		printf ("Can't allocate memory\n");
		delete [] A;
		delete [] Q1;
		fclose(g);
		return -1;
	}

	for (int i = 0; i < n-1; i++){
		Q1[i] = 0.;
		Q2[i] = 0.;
	}

	double* v = new double [n];
	if (v == NULL){
		printf("Can't allocate memory\n");
		delete [] A;
		delete [] Q1;
		delete [] Q2;
		fclose(g);
		return -1;
	}
	t1 = get_time();
	z = values(A,n,Q1,Q2,v,eps);
	if (z == -1){
		delete [] A;
		delete [] Q1;
		delete [] Q2;
		delete [] v;
		fclose(g);
		return -1;
	}
	while (z==3){
		eps *= 10;
//		printf("Need to change precision. Changing to %e\n", eps);
		z = values(A,n,Q1,Q2,v,eps);
		if (z == -1){
			delete [] A;
		    delete [] Q1;
		    delete [] Q2;
			delete [] v;
			fclose(g);
			return -1;
		}
		if (eps > 1e-1){
		    printf("Too big numbers. Cannot find values.\n");
		    delete [] A;
		    delete [] Q1;
		    delete [] Q2;
		    delete [] v;
			fclose(g);
			return -1;
		}
	}
		
		
	t2 = get_time();
	
	file_output(g,n,max,v);
	for (int i = 0; i < n; i++){
		inv1 -= v[i];
	}
	inv2 -= vector_norm(v,n);

	printf("Time = %f\n", t2-t1);
	printf("First invariant (trace) = %e\n", fabs(inv1));
	printf("Second invariant (length) = %e\n", fabs(inv2));
	delete [] A;
	delete [] Q1;
	delete [] Q2;
	delete [] v;
	fclose(g);
	return 0;
}			
	
double get_time(){
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_usec/1000000.0 + (double)tv.tv_sec;
}
