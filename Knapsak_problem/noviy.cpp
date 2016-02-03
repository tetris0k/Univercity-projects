#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <cmath>
using namespace std;
#define EPS 1e-10
static long int it=0;
static double intpart = 0.;
int Compare(const void * , const void * );
double get_time();
double min(double , double );

struct elements {
	double weight;
	double cost;
	double frac;
};

double doit(int , double , elements* );
double maximumCost (elements* , int );
void elemsCost(elements* , int , double* );

int main(){	
	int n;
	double res,W;
	double t1;
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
	t1 = get_time();	
	res = doit(n,W,Elems);
	t1 = get_time()-t1;
	printf("ответ = %lf \n", res);
	fprintf(g,"%f",res);
	printf("время = %e\n", t1);
	printf("отношение = %e \n",t1/it);
	printf("предполагаемое отношение = %e \n", t1/intpart);
	delete[]Elems;
	fclose(f);
	fclose(g);
	return 0;
}
/*
double doit(int n, double W, elements* Elems){
	it =  0;
	double zb,q,otv,Wn,costN;
	int i,j,k;
	double za = maximumCost(Elems,n);	
	for (i = 0; i < n; i++){
		Elems[i].frac = Elems[i].cost / Elems[i].weight;
		it ++;
	}
	qsort (Elems, n, sizeof(elements), Compare);
	double intpart;
	double fracpart = modf(n*log10(n), &intpart);
	it += intpart;
	for (j = 0; j < n-1; j++){
		for (k = j + 1; k < n; k++){
			if (Elems[j].weight + Elems[k].weight - W < EPS){	
				it++;			
				Wn = W - Elems[j].weight - Elems[k].weight;				
				costN = min(Elems[j].cost, Elems[k].cost);	
				it += 2;			
				i = n-1;
				while ((Elems[i].cost - costN > EPS)&&(i > 0)){
					i--;
					it++;
				}
				q = 0.;
				otv = 0.;
				while ((Wn - q > EPS)&&(i >= 0)) {
					if ( (i != j) && (i != k) && (!(Elems[i].cost - costN > EPS)) ){
						q += Elems[i].weight;
						otv += Elems[i].cost;
						it += 2;
					}
					it++;
					i--;
				}
				zb = otv + Elems[j].cost + Elems[k].cost;
				if (zb > za)
					za = zb;
				it += 3;
			}
		}
	}
	printf("итерации = %ld\n", it);
	return(za); 
}
*/
double doit(int n, double W, elements* Elems){
	it =  0;
	intpart = 0.;
	double zb,q,otv,Wn,costN;
	int i,j,k;
	double za = maximumCost(Elems,n);
	for (i = 0; i < n; i++){
		Elems[i].frac = Elems[i].cost / Elems[i].weight;
		it ++;
	}
	qsort (Elems, n, sizeof(elements), Compare);
	double fracpart = modf(n*log10(n), &intpart);
	it += intpart;
	for (j = 0; j < n; j++){
		Wn = W - Elems[j].weight;
		if (j == n-1){
			i = n-2;
		} else {
			i = n-1;
		}
		q = 0.;
		otv = 0.;
		while ((Wn - q > EPS)&&(i >= 0)) {
			if (i != j){
				q += Elems[i].weight;
				otv += Elems[i].cost;
				it += 2;
			}
			it++;
			i--;
		}
		zb = otv + Elems[j].cost;
		if (zb > za) za = zb;
		it += 2;
	}
	cout<<"Итерации: "<<it<<endl;
	return za;
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
void elemsCost(elements* Elems, int n, double* dest){
	for (int i = 0; i < n; i++)
		dest[i] = Elems[i].cost;
}
	
double min(double a, double b){
	if (a - b > EPS){
		return b;
	} else {
		return a;
	}
}
