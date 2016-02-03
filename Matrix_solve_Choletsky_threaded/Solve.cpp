#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "functions.h"

static double threads_total_time;
static pthread_mutex_t threads_total_time_mutex = PTHREAD_MUTEX_INITIALIZER;

void Solve(double* A , double* b, int n, double* Y, int t, pthread_t* threads, double * X, double * d, ARGS *args, double * bb){
	threads_total_time = 0.;	
	rrij(A,b,d,n,X,Y,t,threads,args,bb);
//	R_solve(A,b,n,Y);
//	Dsolve(d,Y,n,X);
//	Rsolve(A,X,n,Y);
}

void * rrij_threaded(void * Ptr){
	ARGS *pargs = (ARGS*)Ptr;
	int n = pargs->n;
	double * A = pargs->A;
	double * d = pargs->d;
	double * b = pargs->b;
	double * Y = pargs->Y;
	double * X = pargs->X;
	double * bb = pargs->bb;
	int total_threads = pargs->total_threads;
	int thread_num = pargs->thread_num;
	int cs,ce;
	int perthread = 0;
	int Rest,tt;		
	int i,j,k;
	double t;
	printf("Thread %d started.\n", thread_num);
	t = get_time();

	if (thread_num == 0){
		d[0] = A[0] / fabs(A[0]);
		A[0] = sqrt(fabs(A[0]));
	}

	synchronize (total_threads);

	double v = 1./(A[0]*d[0]);	
	Rest = (n-1) % total_threads;
	if (!Rest) tt = total_threads; else tt = total_threads-1;
	if (tt)
		perthread = (n-1) / tt;
	if (thread_num < tt){
		cs = thread_num*perthread + 1;
		ce = (thread_num + 1)*perthread + 1;
	} else if (Rest) {
		cs = tt*perthread + 1;
		ce = n;
	}
	for (j = cs; j < ce; j++){
		A[j] *= v;
	}

	synchronize (total_threads);

	for (i = 1; i < n; i++){

		if (thread_num == 0){ 
			drii(A,d,i,n);
		}
	
		synchronize (total_threads);

		Rest = (n-i-1) % total_threads;
		if (!Rest) tt = total_threads; else tt = total_threads-1;
		if (tt)
			perthread = (n-i-1) / tt;
		if (thread_num < tt){
			cs = thread_num*perthread + 1;
			ce = (thread_num + 1)*perthread + 1;
		} else if (Rest) {
			cs = tt*perthread + 1;
			ce = n-i;
		}		
		for (j = cs; j < ce; j++){
			A[i*(n+1) + j - i*(i+1)/2] = rij(A, d, i, i+j, n);
		}
	
		synchronize (total_threads);	
	}
//R_Solve

/*	for (k = 0; k < n; k++){
		if (thread_num == 0) Y[k] = bb[k] / A[k*n + k - k*(k+1)/2];
		
		synchronize (total_threads);

		Rest = (n-k-1) % total_threads;
		if (!Rest) tt = total_threads; else tt = total_threads-1;
		if (tt)
			perthread = (n-k-1) / tt;
		if (thread_num < tt){
			cs = thread_num*perthread + k + 1;
			ce = (thread_num + 1)*perthread + k + 1;
		} else if (Rest) {
			cs = tt*perthread + k + 1;
			ce = n;
		}
		for (i = cs; i < ce; i++){
			bb[i] -= Y[k]*A[i + k*n - k*(k+1)/2];
		}

		synchronize (total_threads);
	}
*/
	if (thread_num == 0) R_solve(A,b,n,Y);
	synchronize(total_threads);

//D_Solve
	Rest = n % total_threads;
	if (!Rest) tt = total_threads; else tt = total_threads-1;
	if (tt)
		perthread = n / tt;
	if (thread_num < tt){
		cs = thread_num*perthread;
		ce = (thread_num + 1)*perthread;
	} else if (Rest) {
		cs = tt*perthread;
		ce = n;
	}
	for (i = cs; i < ce; i++){
		X[i] = Y[i]*d[i];
	}
	synchronize (total_threads);
//RSolve

/*	for (k = 1; k < n+1; k++){
		if (thread_num == 0) Y[n-k] = X[n-k] / A[n*(n+1)/2 - k  - k*(k-1)/2];
		
		synchronize (total_threads);

		Rest = (n-k) % total_threads;
		if (!Rest) tt = total_threads; else tt = total_threads-1;
		if (tt)
			perthread = (n-k) / tt;
		if (thread_num < tt){
			cs = thread_num*perthread;
			ce = (thread_num + 1)*perthread;
		} else if (Rest) {
			cs = tt*perthread;
			ce = n-k;
		}
		for (i = cs; i < ce; i++){
			X[n-k-i-1] -= Y[n-k]*A[n*(n+1)/2 - k - (k+i+1)*(k+i)/2];
		}

		synchronize (total_threads);
	}
	synchronize (total_threads);
*/	if (thread_num == 0) Rsolve(A,X,n,Y);
	synchronize (total_threads);

	t = get_time() - t;

	pthread_mutex_lock (&threads_total_time_mutex);
  	threads_total_time += t;
  	pthread_mutex_unlock (&threads_total_time_mutex);
	printf("Thread %d stopped. Time: %f\n", thread_num, t);
	return 0;
}	
	
void rrij(double * A, double * b, double * d, int n, double * X, double * Y,  int t, pthread_t* threads, ARGS * args, double * bb) {
	int j;
	for (j = 0; j < t; j++){
		args[j].A = A;
		args[j].b = b;
		args[j].d = d;
		args[j].X = X;
		args[j].Y = Y;
		args[j].n = n;
		args[j].thread_num = j;
		args[j].total_threads = t;
		args[j].bb = bb;
	}
	for (j = 0; j < t; j++){
    	if (pthread_create (threads + j, 0, rrij_threaded, args + j)){
			fprintf (stderr, "cannot create thread #%d!\n", j);
        }
    }
	for (j = 0; j < t; j++){
    	if (pthread_join (threads[j], 0))
        	fprintf (stderr, "cannot wait thread #%d!\n", j);
    }	
//	printf("Threads total time: %f\n", threads_total_time);
}
	
	
void drii(double *A, double* d, int i, int n){
	double l = 0.;
	for (int k=0; k<i; k++){
		l -= A[k*n + i - k*(k+1)/2]*A[k*n + i - k*(k+1)/2]*d[k];
	}
	
	l+=A[i*(n+1) - i*(i+1)/2];
	d[i] = l/fabs(l);
	A[i*(n+1) - i*(i+1)/2] = sqrt(fabs(l));
}

double rij(double *A, double *d, int i, int j, int n){
	double l = 0;
	int b;
	for (int k = 0; k < i; k++){
		b = k*(k+1)/2;
		l -= A[k*n + i - b]*A[k*n + j - b]*d[k];
	}
	b = i*(i+1)/2;
	l += A[i*n + j - b];
	l /= A[i*(n+1) - b]*d[i];
	return l;
}

void Rsolve(double *A, double *b, int n, double *x){
	x[n-1] = b[n-1] / A[n*(n+1)/2 - 1];
	for (int k = 2; k < n+1; k++){
		x[n-k] = b[n-k];
		for (int i = 1; i < k; i++){
			x[n-k] -= x[n-i]*A[n*(n+1)/2 - i - k*(k-1)/2];
		}
		x[n-k] /=  A[n*(n+1)/2 - k  - k*(k-1)/2];
	}
}
/*
void Dsolve(double *d, double *b, int n, double *x){
	for (int i = 0; i < n; i++){
		x[i] = b[i]*d[i];
	}
}
*/
void R_solve(double *A, double *b, int n, double *x){
	x[0] = b[0]/A[0];
	for (int k = 1; k < n; k++){
		x[k] = b[k];
		for (int i = 0; i <	k; i++){
			x[k] -= x[i]*A[k + i*n - i*(i+1)/2];
		}
		x[k] /= A[k*n + k - k*(k+1)/2];
	}
}

