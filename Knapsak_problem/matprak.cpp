#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <cmath>
using namespace std;
#define EPS 1e-20
static long int it=0;
int Compare(const void * , const void * );
double max(double , double );
double get_time();

struct elements {
	double weight;
	double cost;
	double frac;
};
double perebor (int , double , elements* );
double doit(int , double , elements* );
double maximumCost (elements* , int );
int main(){	
	int n;
	double res,res1,W;
	double t1,t2;
	elements *Elems;
	FILE* f;
	FILE* g;
	f = fopen("input.txt", "r");
	if (f==NULL) return -1;
	g = fopen("output.txt", "wr");
	if (g==NULL) return -1;
	if (fscanf(f,"%d ",&n)==0) return -1;
	Elems = new elements[n];
	for (int i=0;i<n;i++) {
		if (fscanf(f,"%lf ",&(Elems[i]).cost)==0) return -1;
	}	
	for (int i=0;i<n;i++) {
		if (fscanf(f,"%lf ",&Elems[i].weight)==0) return -1;	
	}
	if (fscanf(f,"%lf ", &W)==0) return -1;
//	res1 = perebor(n,W,Elems);
//	printf("ответ перебора = %lf\n", res1);
	t1 = get_time();	
	res = doit(n,W,Elems);
	t2 = get_time();
	printf("ответ = %lf \n", res);
	fprintf(g,"%f",res);
	printf("время = %e\n", t2-t1);
	printf("отношение = %e\n", (t2-t1)/it);
	delete[]Elems;
	fclose(f);
	fclose(g);
	return 0;
}

double doit(int n, double W, elements* Elems){
	it =  0;
	double za = maximumCost(Elems,n);		
	for (int i=0;i<n;i++){
		Elems[i].frac=Elems[i].cost / Elems[i].weight;
		it ++;
	}
	qsort (Elems, n, sizeof(elements), Compare);
	double intpart;
	double fracpart = modf(n*log10(n), &intpart);
	it += intpart;
	int i = n-1;
	double q,otv;
	q = 0.;
	otv = 0.;
	while ((q < W)&&(i >= 0)) {
		q+=Elems[i].weight;
		otv+=Elems[i].cost;
		i--;
		it += 3;
	}
	if (za > otv) otv = za;
	it++;
	printf("итерации = %ld\n", it);
	return(otv); 
}

double perebor (int k, double W, elements* Elems){
	if (k == 0){
		return 0.;
	} else if (W < EPS){
		return 0.;
	} else {
		return max( perebor(k-1,W,Elems), perebor(k-1,W - Elems[k-1].weight, Elems) + Elems[k-1].cost );
	}
}
	

double get_time(){
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_usec/1000000.0 + (double)tv.tv_sec;
}

int Compare(const void * a, const void * b) {
	elements x = *(elements*)a;
	elements y = *(elements*)b;
	if ((x.frac - y.frac) > EPS){
		return 1;
	} else if ((y.frac - x.frac)>EPS){
		return -1;
	} else {
		return 0;
	}
}
double maximumCost (elements* a, int n){
	double m = a[0].cost;
	for (int i = 1; i < n; i++){
		if (a[i].cost > m) m = a[i].cost;
		it++;
	}
	return m;
}

double max(double a, double b){
	if (a - b > EPS){
		return a;
	} else {
		return b;
	}
}
