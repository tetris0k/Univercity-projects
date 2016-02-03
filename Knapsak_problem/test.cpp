#include <stdio.h>
#include <stdlib.h>
int main(){
	FILE* f;
	int a;
	double t,w;
	f = fopen("input.txt", "wr");
	int n;
	printf("Cколько предметов в наличии?\n");
	scanf("%d", &n);
	fprintf(f,"%d ",n);
	for (int i=0;i<n;i++){
		t =(rand()%100000)/100000. + 10.;//стоимости от 10 до 11
		fprintf(f,"%f ",t);	
	}
		for (int i=0;i<n;i++){
		t =(rand()%(30000*n))/10000.+1.;//веса от нуля до 100
		fprintf(f,"%f ",t);	
	}
	w=3*n;
	fprintf(f,"%f ",w);
	fclose(f); 
	return 0;
}
