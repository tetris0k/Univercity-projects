#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;

int main(){
	ofstream fout("input_words.txt");
	int l1, l2;
	cout<<"length of words:"<<endl;
	cin>>l1;
	cin>>l2;
	cout<<"your words:"<<endl;
	int t;
	while(l1>0){
		t = rand()%10;
		if (t > 1) t = 1; else t = 0;
		cout << t;
		fout << t;
		l1--;
	}
	fout << " ";
	cout << " ";
	while(l2>0){
		t = rand()%10;
		if (t > 1) t = 1; else t = 0;
		cout << t;
		fout << t;
		l2--;
	}
	fout<<" ";
	cout<<endl;
	fout.close();
	return 0;
}
