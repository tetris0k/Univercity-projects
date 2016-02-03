#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include <pthread.h>


static double threads_total_time = 0.;

static pthread_mutex_t threads_total_time_mutex = PTHREAD_MUTEX_INITIALIZER;

void Solve(int n, double* A, double* x, int t, pthread_t * threads, double * E, solve * args, double * v){
	threads_total_time = 0.;		
	reshenie(A,E,v,x,n,t,threads,args);
	printf("общее время работы нитей решения: %f\n", threads_total_time);	
}

void reshenie(double * A, double * E, double * v, double * X, int n, int t, pthread_t* threads, solve * args) {
	int i;	
	for (i = 0; i < t; i++){                       									//заполняем структуры для каждой из нитей
		args[i].A = A;
		args[i].E = E;
		args[i].v = v;
		args[i].X = X;
		args[i].n = n;
		args[i].thread_num = i;
		args[i].total_threads = t;
	}
	for (i = 0; i < t; i++){                       									//запускаем нити
    	if (pthread_create (threads + i, 0, reshenie_threaded, args + i)){
			fprintf (stderr, "не получилось создать нить %d\n", i);
        }
    } 
	for (i = 0; i < t; i++){                        								//ждем, пока они закончат работу
    	if (pthread_join (threads[i], 0))
        	fprintf (stderr, "не получилось ждать нить %d\n", i);
    }	
}

	
void * reshenie_threaded(void * arg){                   							//принимает ссылку на структуру
	solve *pargs = (solve*)arg;                  									//вытаскиваем входную структуру
	double * A = pargs->A;                     	 									//вытаскиваем информацию из структуры
	double * E = pargs->E;
	double * v = pargs->v;
	double * X = pargs->X;
	int n = pargs->n;
	int thread_num = pargs->thread_num;          									//номер нити
	int total_threads = pargs->total_threads;    									//всего нитей
	int i,j,k;
	double t,tmp;
	printf("нить решения номер %d начала работу\n", thread_num);
	t = get_time();                                     							//get_time для счета времени внутри нити
															
	for (i = thread_num; i < n+1; i+=total_threads){                    			//заполнение E. распараллелен цикл от 0 до n+1
		for (j = 0; j < n+1; j++){
			if (i == j){
				E[i*(n+1) + j] = 1; 
			} else {
				E[i*(n+1)+j] = 0;
			}
		}
	}

	synchronize (total_threads);

	for (k = 1; k < n+1; k++){
		for (i = thread_num; i < n+1; i+=total_threads){
			v[i] = E[i*(n+1) + k-1];												//nowv. распараллелен цикл от 0 до n+1
		}
																		
		synchronize (total_threads);

		for (i = thread_num + k; i < n+1; i+=total_threads){						//master.распараллелен цикл от k до n+1
			master(n,v,i,k,A,E);
		}
		
		synchronize (total_threads);
	}
	
	synchronize(total_threads);

	for (j = thread_num; j < n; j+=total_threads){                            		//запись решения в x. распараллелен цикл от 0 до n
		X[j] = E[(n+1)*j + n];
	}	

	synchronize (total_threads);

	t = get_time() - t;
	printf("нить решения номер %d закончила работу. время: %f\n", thread_num, t);
	pthread_mutex_lock (&threads_total_time_mutex);           						//блокируем нить для подсчета суммарного времени (глоб переменная)
	threads_total_time += t;
	pthread_mutex_unlock (&threads_total_time_mutex);         						//разблокируем нить
	return 0;
}



void master(int n, double *v, int i , int k, double *A, double *E){
	double x,y;
	x = smult(n,i,k,A,E);
	y = smult2(n,v,k,A);
	double a = -x/y;
	for (int j = 0; j < n+1; j++){
		E[j*(n+1)+i] = E[j*(n+1)+i] + a*v[j];
	}
}

double smult2(int n, double *v , int k, double *A){
	double res=0.;
	for (int j = 0; j < n+1; j++){
		res = res + A[(k-1)*(n+1)+j]*v[j];
	}	    
	return res;
}

double smult(int n, int i , int k, double *A, double *E){
	double res=0.;
	for (int j = 0; j < k; j++){
		res = res + A[(k-1)*(n+1)+j]*E[j*(n+1)+i];
	}	
	res = res + A[(k-1)*(n+1)+i]*E[i*(n+1)+i];    
	return res;
}
/*
void nowv(int n, double *E, double *v, int k){
	for(int i=0;i<n+1;i++){
		v[i]=E[i*(n+1)+k-1];
	}
}*/
void PrintMatrix(int n, double * A, double * b, int MAX_OUTPUT_SIZE)
{
    int i;
    int j;
    int nPrint;

    nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;

    for (i = 0; i < nPrint; ++i) 	{
        printf("| ");
        for (j = 0; j < nPrint; ++j)
            printf("%5.3g ", A[i * (n+1) + j]);
        if (nPrint < n)
            printf("\t...\t");
        printf("|");
        printf("%5.3g\n", b[i]);
    }
}




////////////////////////////////////////////////синхронизация от богачева


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
	
