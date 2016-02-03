#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <cstring>
#define min(a,b) ((a<b) ? a : b)

static int tmp;

using namespace std;

double get_time();
void preKMP(string , int* );
bool KMP(string , string );

int main(void){
	ifstream fin("input_words.txt");
	ofstream fout("output_word.txt");
	int l,t = 0;
    string v1;
    string v2;
	string flag;
	string max;
    getline(fin,v1);

    char *s1 = new char [v1.size()+1];     //space cleaning
	strcpy(s1, v1.c_str());
	int i = 0;
	while(isblank(s1[i]) && (i < v1.size()+1)) i++;
	if (i > 0) v1.erase(0,i);
	strcpy(s1, v1.c_str());
    int j = 0;
	while ((!(isblank(s1[j])))&&(j < v1.size())) j++;
	int l1 = v1.size();
	v2 = v1.substr(j, l1-j);
	v1.erase(j,l1-j);
	strcpy(s1, v2.c_str());
	i = 0;
	while(isblank(s1[i])) i++;
	if (i > 0) v2.erase(0,i);
	while(!(isblank(s1[i]))) i++;
	if (i < v2.size()) v2.erase(i, v2.size()-i);
	delete []s1;


	l1 = v1.size();             //after that v1 will be smaller then v2
	if (l1 > v2.size()){
		flag = v1;
		v1 = v2;
		v2 = flag;
		l1 = v1.size();
	}

	int l2 = v2.size();
	bool ok = false;
	tmp = 0;
	i = l1;

	double t1 =  get_time();

	while ((i != 0)&&(t != 1)){
		for (j = 0; j <= l1-i; j++){
			flag = v1.substr(j,i);
			if (KMP(flag,v2)){
				max = flag;
				t = 1;
				break;
			}
		}
		i--;
	}

	t1 = get_time() - t1;

	if (t == 0){
		cout<<"there is no same subwords"<<endl;
		fout<<"there is no same subwords"<<endl;
	} else {
		cout<<"Max same subword: "<<max<<endl;
		fout<<max<<endl;
	}
	cout<<"time: "<<t1<<endl;
	cout<<"iterations: "<<tmp<<endl;
	cout<<"fraction: "<< t1/tmp <<endl;
	fin.close();
	fout.close();
    return 0;
}
double get_time(){
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_usec/1000000.0 + (double)tv.tv_sec;
}

void preKMP(string pattern, int* f){
    int m = pattern.length(), k;
    f[0] = -1;
    for (int i = 1; i < m; i++){
        k = f[i - 1];
        while (k >= 0){
			tmp++;
            if (pattern[k] == pattern[i - 1])
                break;
            else
                k = f[k];
        }
        f[i] = k + 1;
    }
}


bool KMP(string pattern, string target){   //Knuth-Morris-Pratt algorithm
    int m = pattern.length();
    int n = target.length();
    int f[m];
    preKMP(pattern, f);
    int i = 0;
    int k = 0;
    while (i < n){
        if (k == -1){
            i++;
            k = 0;
        } else if (target[i] == pattern[k]){
            i++;
            k++;
            tmp++;
            if (k == m)
                return 1;
        } else {
            k = f[k];
		}
    }
    return 0;
}

