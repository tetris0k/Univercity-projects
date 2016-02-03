#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "functions.h"

static double flg = 0.;
static double threads_total_time = 0.;

static pthread_mutex_t threads_total_time_mutex = PTHREAD_MUTEX_INITIALIZER;

double residual_norm1(int n, double *A, double *b, double *x, int t, pthread_t * threads, ARGS * args){
	flg = 0.;
	threads_total_time = 0.;
	int j;
	for (j = 0; j < t; j++){
		args[j].A = A;
		args[j].X = x;
		args[j].b = b;
		args[j].n = n;
		args[j].thread_num = j;
		args[j].total_threads = t;
	}
	for (j = 0; j < t; j++){
    	if (pthread_create (threads + j, 0, resid1_threaded, args + j)){
			fprintf (stderr, "cannot create thread #%d!\n", j);
        }
    }
	for (j = 0; j < t; j++){
    	if (pthread_join (threads[j], 0))
        	fprintf (stderr, "cannot wait thread #%d!\n", j);
    }	
	return sqrt(flg);
}

void * resid1_threaded(void * arg){
	ARGS * pargs = (ARGS*)arg;
	double * A = pargs->A;
	double * b = pargs->b;
	double * X = pargs->X;
	int n = pargs->n;
	int thread_num = pargs->thread_num;
	int total_threads = pargs->total_threads;
	int cs,ce;
	int perthread = 0;
	int Rest,tt;		
	int i,j;
	double t,tmp;
	printf("Thread %d started.\n", thread_num);
	t = get_time();	
	
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
		tmp = 0.;
		for (j = 0; j < n; j++){
			tmp += A[i*n+j]*X[j];
		}
		tmp -= b[i];
		pthread_mutex_lock (&threads_total_time_mutex);
		flg += tmp*tmp;
		pthread_mutex_unlock (&threads_total_time_mutex);
	}
	t = get_time() - t;
	pthread_mutex_lock (&threads_total_time_mutex);
	threads_total_time += t;
	pthread_mutex_unlock (&threads_total_time_mutex);
	printf("Thread %d stopped. Time: %f\n", thread_num, t);
	return 0;
}

double residual_norm(int n, double *b, double *x, int t, pthread_t * threads, ARGS * args){
	flg = 0.;
	threads_total_time = 0.;
	int j;
	for (j = 0; j < t; j++){
		args[j].X = x;
		args[j].b = b;
		args[j].n = n;
		args[j].thread_num = j;
		args[j].total_threads = t;
	}
	for (j = 0; j < t; j++){
    	if (pthread_create (threads + j, 0, resid_threaded, args + j)){
			fprintf (stderr, "cannot create thread #%d!\n", j);
        }
    }
	for (j = 0; j < t; j++){
    	if (pthread_join (threads[j], 0))
        	fprintf (stderr, "cannot wait thread #%d!\n", j);
    }	
	return sqrt(flg);
}

void * resid_threaded(void * arg){
	ARGS * pargs = (ARGS*)arg;
	double * b = pargs->b;
	double * X = pargs->X;
	int n = pargs->n;
	int thread_num = pargs->thread_num;
	int total_threads = pargs->total_threads;
	int cs,ce;
	int perthread = 0;
	int Rest,tt;		
	int i,j;
	double t,tmp;
	printf("Thread %d started.\n", thread_num);
	t = get_time();	
	
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
		tmp = 0.;
		for (j = 0; j < n; j++){
			tmp += Aij(i,j)*X[j];
		}
		tmp -= b[i];
		pthread_mutex_lock (&threads_total_time_mutex);
		flg += tmp*tmp;
		pthread_mutex_unlock (&threads_total_time_mutex);
	}
	t = get_time() - t;
	pthread_mutex_lock (&threads_total_time_mutex);
	threads_total_time += t;
	pthread_mutex_unlock (&threads_total_time_mutex);
	printf("Thread %d stopped. Time: %f\n", thread_num, t);
	return 0;
}

double error_rate(int n, double* x){
	double l = 0.;
	for (int i = 0; i < n; i += 2){
		l += (x[i]-1)*(x[i]-1);
	}
	for (int i = 1; i < n; i += 2){
		l += x[i]*x[i];
	}
	return sqrt(l);
}
	

void synchronize (int total_threads){
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;
	pthread_mutex_lock (&mutex);
	threads_in++;
	if (threads_in >= total_threads){
		threads_out = 0;
    	pthread_cond_broadcast (&condvar_in);
    } else {
   		while (threads_in < total_threads){
        	pthread_cond_wait (&condvar_in, &mutex);
        }
    }
	threads_out++;
	if (threads_out >= total_threads){
	    threads_in = 0;
	    pthread_cond_broadcast (&condvar_out);
    } else {
    	while (threads_out < total_threads){
	        pthread_cond_wait (&condvar_out, &mutex);
        }
    }
	pthread_mutex_unlock (&mutex);
}
