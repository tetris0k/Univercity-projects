#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include "functions.h"

using namespace std;

int main(int argc, char **argv){
    int my_rank;
    int p,q;
    int n;
    double t;
    FILE* f;
	FILE* g;
    double *A;
    double *c;
    double *d;
    double *b;
    double *A1;
    double *x;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Barrier (MPI_COMM_WORLD);

    int max = 10;
    int ret = 0;
    opterr = 0;
    while ((ret = getopt( argc, argv, "f::g:N:")) != -1) {
		switch (ret) {
		case 'f':
			q = 1;
			if (my_rank == 0){
			    if (optarg != NULL){
				    printf ("Input from file %s\n", optarg);
				    if ((f = fopen (optarg, "r")) == NULL) {
					    perror ("Can't open file.\n");
					    MPI_Abort(MPI_COMM_WORLD, 1);
                        return -1;
				    }
			     
			    } else {
				    printf ("Input from file input.txt\n");
				    if ((f = fopen ("input.txt", "r")) == NULL) {
					    perror ("Can't open file.\n");
					    MPI_Abort(MPI_COMM_WORLD, 1);
                        return -1;
				    } 
			    }
			}
			break;
		case 'g':
			if (my_rank == 0) printf("Input from function in file input_function.cpp\n");
			q = 2;
			n = atoi(optarg);
			if (my_rank == 0) printf("Size of matrix: %d\n", n);
			break;
		case '?':
			if (my_rank == 0){
				printf("Unknown option -%c\n", optopt);
				printf("Supported options:\n");
				printf("-f -- Input matrix from the argument file, or, if there is no argument, input matrix from the file input.txt;\n");
				printf("-g -- Input matrix from the function. Argument of this option means size of matrix and it is required;\n");
			}
			MPI_Finalize ();
			return -1;
		}
	}
	if (q == 1){
		if (my_rank == 0) {
			if(fscanf (f, "%d", &n) == 0){
				printf ("Problem with size of matrix.\n");
				MPI_Abort(MPI_COMM_WORLD, 1);
				return -1;
			}
		}
		MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier (MPI_COMM_WORLD);
		A = new double [((int)(n/p)+1)*n];
		b = new double [n];
		A1 = new double [((int)(n/p)+1)*n];
		c = new double [n];
		if (A == NULL || b == NULL || A1 == NULL || c == NULL){
			printf ("Problem with memory allocation in %d process.\n", my_rank);
			MPI_Abort(MPI_COMM_WORLD, 1);
			return -1;
		}
		if (file_input(f, n, A, A1, b, p, my_rank, c) != 0){
			cout<<"file_input is not working"<<endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
			return -1;
		}
		if (my_rank == 0) fclose(f);
	} else if (q == 2){
		A = new double [((int)(n/p)+1)*n];
		b = new double [n];
		c = new double [n];
		if (A == NULL || b == NULL|| c == NULL){
			printf ("Problem with memory allocation.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return -1;
		}
		if (func_input(A, b, n, p, my_rank) != 0){
			cout<<"func_input is not working"<<endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
			return -1;
		}
	} else {
		if (my_rank == 0){
			printf("You need to run the program with boot options\n"); 
			printf("Supported options:\n");
			printf("-f -- Input matrix from the argument file, or, if there is no argument, input matrix from the file input.txt;\n");
			printf("-g -- Input matrix from the function. Argument of this option means size of the matrix and it's required;\n");
			printf ("\n\n");
		}
		MPI_Finalize ();
		return -1;
	}
	if (my_rank == 0){
		g = fopen("output.txt", "w");
		if(g == NULL){
			printf("Can't open output file.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return -1;
		}
	}
	MPI_Barrier (MPI_COMM_WORLD);

	PrintMatrix(n, A, b, max, p, my_rank, c);
	MPI_Barrier(MPI_COMM_WORLD);
	x = new double [n];
	d = new double [n];
	if (x == NULL || d == NULL){
		printf ("Problem with memory allocation in %d process.\n", my_rank);
		MPI_Abort(MPI_COMM_WORLD, 1);
		return -1;
	}
    
	MPI_Barrier (MPI_COMM_WORLD);
	t = MPI_Wtime();
    
	Solve(A, b, n, x, c, d, p, my_rank);                            //solvation
    
	MPI_Barrier (MPI_COMM_WORLD);
	t = MPI_Wtime() - t;
    
	if (my_rank == 0){
		file_output(g, n, max, x);
		cout << "Time: " << t << endl;
	}
	if (q == 1){                                         //residual norm
		double resid = residual_norm1(n, A1, b, x, p, my_rank);
		if (my_rank == 0)
			cout << "Residual norm: " << resid << endl;
		delete []A1; 	
	} else if (q == 2){
		if (my_rank == 0)
			cout << "Residual norm: " << residual_norm(n, b, x) <<endl;	
	}
	if (my_rank == 0)
		cout << "Error rate: " << error_rate(n, x) << endl;       //error rate   
	delete []A;
	delete []b;
	delete []c;
	delete []d;
	delete []x;
	if (my_rank == 0)
		fclose(g);
	MPI_Finalize();
	return 0;
}
