#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "functions.h"

double get_full_time (void){
	struct timeval buf;
	gettimeofday (&buf, 0);
	return  buf.tv_usec / 10000000.0 + (double)buf.tv_sec; 
}
