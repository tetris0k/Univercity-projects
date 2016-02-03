typedef struct solve{					//структура для параллельного решения задачи
	double* A;							//матрица A
	double* E;							//матрица базисов
	double* v;							//вектор, который запоминаем (nowv)
	double* X;							//решение
	int n;								//размерность A
	int thread_num;						//номер нити
	int total_threads;					//всего нитей
} solve;	


typedef struct Nevyaska{				//структура для параллельного вычисления невязки
	double * A;							//матрица A
	double * x;							//решение
	int n; 								//размерность A
	int thread_num;						//номер нити
	int total_threads;					//всего нитей
} Nevyaska;

double smult(int ,int ,int ,double* , double*);
double smult2(int ,double * , int , double *);
void master(int , double *, int  , int , double *, double *);
//void nowv(int ,double *, double *, int );
void Solve(int , double* , double* ,int , pthread_t * , double *, solve * , double * );
int file_input(FILE * , int , double* );
void file_output(FILE * , int , int , double* );
double Aij(int , int );
int func_input(double* ,int );
double nev(int , double *, double *, int , pthread_t* , Nevyaska *);
void * Nevyaska_threaded(void * );
double error_rate(int , double* );
void PrintMatrix(int , double * , double * , int);
double get_time();
double get_full_time();
void reshenie(double * , double * , double * , double * , int , int , pthread_t* , solve * );
void * reshenie_threaded(void * );
void synchronize (int );
