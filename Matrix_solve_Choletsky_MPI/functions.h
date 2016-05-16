double Aij(int , int );
int func_input(double* , double* , int , int , int );
int file_input(FILE* , int , double*, double* , double* , int , int , double* );
void file_output(FILE* , int , int , double* );
/*double dii(double* , double* , int , int );
double rii(double* , double* , int , int );
double rij(double* , double* , int , int , int );
void Rsolve(double* , double* , int , double* );
void Dsolve(double* , double* , int , double* );
void R_solve(double* , double* , int , double* );
*/
void Solve(double* , double* , int , double* , double* , double* , int , int );
double residual_norm(int , double* ,double* );
double error_rate(int ,double* );
double residual_norm1(int , double* ,double* ,double* , int , int );
void PrintMatrix(int , double * , double * , int , int , int , double* );

