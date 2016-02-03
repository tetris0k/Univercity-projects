#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "functions.h"

/* Вернуть процессорное время, затраченное на текущий
   поток. Берется время только
   самого потока, время работы системных вызовов не
   прибавляется. */
double get_time (void){
	struct rusage buf;
	getrusage (RUSAGE_THREAD, &buf);
	return   buf.ru_utime.tv_usec / 10000000.0 + (double)buf.ru_utime.tv_sec;
}
/*Астрономическое время*/
double get_full_time (void){
	struct timeval buf;
	gettimeofday (&buf, 0);
	return  buf.tv_usec / 10000000.0 + (double)buf.tv_sec; 
}
