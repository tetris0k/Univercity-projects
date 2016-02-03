#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#include <pthread.h>


static double nevv = 0.;
static double threads_total_time = 0.;

static pthread_mutex_t threads_total_time_mutex = PTHREAD_MUTEX_INITIALIZER;


double nev(int n, double *A, double *x, int t, pthread_t* threads, Nevyaska * args){
	nevv = 0.;
	threads_total_time = 0.;
	int i;	
	for (i = 0; i < t; i++){                                	//заполняем структуры для нитей
		args[i].A = A;
		args[i].x = x;
		args[i].n = n;
		args[i].thread_num = i;
		args[i].total_threads = t;
	}
	for (i = 0; i < t; i++){                                	//запускаем каждую из нитей
    	if (pthread_create (threads + i, 0, Nevyaska_threaded, args + i)){      
			fprintf (stderr, "не получилось создать нить %d!\n", i);
        }
    }
	for (i = 0; i < t; i++){                                	//ждем, пока нити закончат работу
    	if (pthread_join (threads[i], 0))
        	fprintf (stderr, "не получилось ждать нить %d!\n", i);
    }	
	printf("общее время работы нитей невязки = %f\n", threads_total_time);
	return sqrt(nevv);
}

void * Nevyaska_threaded(void * Ptr){       					//получает ссылку на структуру
	Nevyaska *pargs = (Nevyaska*)Ptr;                 			//вытаскиваем структуру
	double * A = pargs->A;             	              			//считываем информацию из структуры
	double * x = pargs->x;
	int n = pargs->n;
	int thread_num = pargs->thread_num;
	int total_threads = pargs->total_threads;
	int i,j;
 	double t,m;
//	printf("Нить невязки %d начала работу.\n", thread_num);
	t = get_time();                                  			//get_time для счета времени внутри нити
	
	for (i = thread_num; i < n; i+=total_threads){
		m = 0.;
		for (j = 0; j < n; j++){
			m = m + A[i*(n+1) + j]*x[j];
		}
		m = m + A[i*(n+1)+n];
		pthread_mutex_lock (&threads_total_time_mutex);      	//блокируем нить для подсчета невязки (глобальная переменная)
  		nevv = nevv + m*m;
  		pthread_mutex_unlock (&threads_total_time_mutex);    	//разблокируем нить
	}
	
	t = get_time() - t;
//	printf("Нить невязки %d закончила работу. Время: %f\n", thread_num, t);
	pthread_mutex_lock (&threads_total_time_mutex);         	//блокируем нить для подсчета общего времени (глобальная переменная)
	threads_total_time = threads_total_time + t;	
	pthread_mutex_unlock (&threads_total_time_mutex);       	//разблокируем нить
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

