#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
using namespace std;
#define EPS 1e-10

int Compare(const void * , const void * );
double min(double , double );

struct elements {
	double weight;
	double cost;
	double frac;
};

double doit_old(int , double , elements* );
double maximumCost (elements* , int );
double doit_new(int , double , elements* );

int main(){	
	FILE* f;
	FILE* g1;
	FILE* g2;
	int n,n1;
	g1 = fopen("out1.txt", "wr");
	g2 = fopen("out2.txt", "wr");
	printf("количество предметов:\n");
	scanf("%d",&n);
	printf("длина цикла:\n");
	scanf("%d",&n1);
	double t,res1,res2,W;
	W = 3*n;
	elements *Elems = new elements[n];
	srand(time(NULL)); 
	for (int o = 0; o < n1; o++){
		f = fopen("input.txt", "wr");	
		for (int i=0;i<n;i++){
			t =(rand()%10000)/10000. + 10.;//стоимости от 10 до 11
			fprintf(f,"%f ",t);	
		}
		for (int i=0;i<n;i++){
			t =(rand()%(3000*n))/1000.+1.;//веса от 1 до размера рюкзака
			fprintf(f,"%f ",t);	
		}
		fclose(f);	

		f = fopen("input.txt", "r");
		for (int i=0;i<n;i++) {
			if (fscanf(f,"%lf ",&(Elems[i]).cost)==0) return -1;
		}	
		for (int i=0;i<n;i++) {
			if (fscanf(f,"%lf ",&Elems[i].weight)==0) return -1;	
		}
		res1 = doit_old(n,W,Elems);
		res2 = doit_new(n,W,Elems);
		fprintf(g1, "%lf\n", res1);
		fprintf(g2, "%lf\n", res2);
		fclose(f);
	}
	delete[]Elems;
	fclose(g1);
	fclose(g2);
	return 0;
}
/*
double doit_new(int n, double W, elements* Elems){
	double zb,q,otv,Wn,costN;
	int i,j,k;
	double za = maximumCost(Elems,n);	
	for (i = 0; i < n; i++){
		Elems[i].frac = Elems[i].cost / Elems[i].weight;
	}
	qsort (Elems, n, sizeof(elements), Compare);
	for (j = 0; j < n-1; j++){
		for (k = j + 1; k < n; k++){
			if (Elems[j].weight + Elems[k].weight - W < EPS){	
				Wn = W - Elems[j].weight - Elems[k].weight;				
				costN = min(Elems[j].cost, Elems[k].cost);	
				i = n-1;
				while ((Elems[i].cost - costN > EPS)&&(i > 0)){
					i--;
				}
				q = 0.;
				otv = 0.;
				while ((Wn - q > EPS)&&(i >= 0)) {
					if ( (i != j) && (i != k) && (!(Elems[i].cost - costN > EPS)) ){
						q += Elems[i].weight;
						otv += Elems[i].cost;
					}
					i--;
				}
				zb = otv + Elems[j].cost + Elems[k].cost;
				if (zb > za)
					za = zb;
			}
		}
	}
	return(za); 
}
*/
double doit_old(int n, double W, elements* Elems){
	double za = maximumCost(Elems,n);	
	for (int i=0;i<n;i++){
		Elems[i].frac=Elems[i].cost / Elems[i].weight;
	}
	qsort (Elems, n, sizeof(elements), Compare);
	int i = n-1;
	double q,otv;
	q = 0.;
	otv = 0.;
	while ((q < W)&&(i >= 0)) {
		q+=Elems[i].weight;
		otv+=Elems[i].cost;
		i--;
	}
	if (za > otv) otv = za;
	return(otv); 
}

double doit_new(int n, double W, elements* Elems){
	double zb,q,otv,Wn,costN;
	int i,j,k;
	double za = 0.;
	for (i = 0; i < n; i++){
		Elems[i].frac = Elems[i].cost / Elems[i].weight;
	}
	qsort (Elems, n, sizeof(elements), Compare);
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
			}
			i--;
		}
		zb = otv + Elems[j].cost;
		if (zb > za) za = zb;
	}
	return za;
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
	}
	return m;
}
	
double min(double a, double b){
	if (a - b > EPS){
		return b;
	} else {
		return a;
	}
}
