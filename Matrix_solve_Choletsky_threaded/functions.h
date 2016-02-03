typedef struct ARGS{
	double* A;
	double* b;
	double* d;
	double* X;
	double* Y;
	int n;
	int thread_num;
	int total_threads;
	double * bb;
} ARGS;

double Aij(int , int );
int func_input(double* , double* , int );
int file_input(FILE* , int , double*, double* , double* );
void file_output(FILE* , int , int , double* );
void drii(double* , double* , int , int );
double rij(double* , double* , int , int , int );
void Rsolve(double* , double* , int , double* );
void Dsolve(double* , double* , int , double* );
void R_solve(double* , double* , int , double* );
void Solve(double* , double* , int , double* ,int , pthread_t * , double * , double * , ARGS * , double *);
double residual_norm(int , double *, double *, int , pthread_t * , ARGS * );
double error_rate(int ,double* );
double residual_norm1(int , double *, double *, double *, int , pthread_t * , ARGS * );
void * resid_threaded(void * );
void * resid1_threaded(void * );
void PrintMatrix(int , double * , double * , int );
double get_time();
double get_full_time();
void * rrij_threaded(void * );
void rrij(double * , double *, double * , int ,double *, double *, int , pthread_t * , ARGS *, double * );
void synchronize (int );


