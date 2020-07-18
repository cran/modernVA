double* R_Vector(int);
double* R_VectorInit(int,double);
double** R_Matrix(int,int);
double** R_MatrixInit(int,int,double);
double quform(double*,double*,int);
void Rprintvec(char*,double*,int);
void RprintVecAsMat(char*,double*,int,int);
void RprintIvec(char*,int*,int);
void RprintIVecAsMat(char*,int*,int,int);
void mat_transpose(double*, double*, int, int);
int factorial(int);
void ran_mvnorm(double*,double*,int,double*,double*);
void ran_wish(int,double*,int,double*,double*,double*,double*);
void ran_dirich(double*, int, double*, double*);
double rinvgauss(double, double);
double rtnorm(double, double, double, double);
double dmvnorm(double*,double*,double*,int,double,double*,int);
double dinvwish(double*, double, double, int, int);
double dinvgamma(double, double, double, int);
double dtnorm(double, double, double, double, double, int);