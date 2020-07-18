/*************************************************************
 * Copyright (c) 2009 Steven McKay Curtis and Garritt L. Page
 *
 * I give permission for anyone to use and/or modify these
 * programs provided the following conditions are met:
 *
 * 1) This copyright message is retained.
 * 2) Any changes to the code are noted.
 *
 * These programs are provided WITHOUT ANY WARRNTY.
 *
 *************************************************************/

#include <R.h>
#include <Rmath.h>
#include "matrix.h"
#include <math.h>

#include <R_ext/Lapack.h>


//============================================================
// Allocating and Freeing Vectors and Matrices
//
// Note that the Matrix allocation functions allocate the
// memory for the entire matrix in a contiguous block.

//==================================================
// Vectors
double* R_Vector(int n){
    double* v;
    v = (double *) R_alloc(n,sizeof(double));
    return(v);
}

double* R_VectorInit(int n,double init){
    double* v;
    v = (double *) R_alloc(n,sizeof(double));
    for(int i=0; i<n; i++)
	v[i]=init;
    return(v);
}


//==================================================
// Matrices
double** R_Matrix(int nr, int nc){
    double** m;
    m = (double **) R_alloc(nr,sizeof(double*));
    *m = (double *) R_alloc(nr*nc,sizeof(double));
    for (int i=1; i<nr; i++)
	*(m+i) = *m + i*nc;
    return(m);
}

double** R_MatrixInit(int nr, int nc, double init){
    double** m;
    int i, j;
    m = (double **) R_alloc(nr,sizeof(double*));
    *m = (double *) R_alloc(nr*nc,sizeof(double));
    for (i=1; i<nr; i++)
	*(m+i) = *m + i*nc;
    for (i=0; i<nr; i++)
	for (j=0; j<nc; j++)
	    m[i][j] = init;
    return(m);
}



//============================================================
// Matrix Operations
//


double quform(double *x, double* A, int dim){
    int i,j;
    double sm=0.0;
    for (i=1; i<dim; i++)
	for (j=0; j<i; j++)
	    sm += x[i]*x[j]*A[i*dim+j];
    sm*=2;
    for(i=0; i<dim; i++)
	sm += x[i]*x[i]*A[i*dim+i];
    return(sm);
}




//************************************************************
// Printing Matrices and Vectors
//
// The following functions
//
// Rprintvec, Rprintmat, RprintIvec, RprintImat,
//
// are modified versions of functions
// provided by Howard Seltman at the following web page:
// http://www.stat.cmu.edu/~hseltman/files/vmr.c
//
// I have modified the functions to work with R and to
// provide slightly modified output to suit my tastes.
//

// for doubles
void Rprintvec(char* title, double *v, int l){
    if (title!=NULL)
	Rprintf("%s\n",title);
    for (int i=0; i<l; i++)
	Rprintf("%f\n", v[i]);
    Rprintf("\n");
    return;
}



// for integers
void RprintIvec(char* title, int* v, int n){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<n; i++)
	Rprintf("%i\n", v[i]);
    Rprintf("\n");
    return;
}



// Print a vector as a matrix
void RprintVecAsMat(char* title, double *v, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%f ", v[i*nc + j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Print a vector as a matrix
void RprintIVecAsMat(char* title, int *v, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%d ", v[i*nc + j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Matrix tranpose
void mat_transpose(double *mat, double *tmat, int nr, int nc)
	{

//		RprintVecAsMat("Uk", mat, 1, (nr)*(nc));
		int i, j;
		for(i = 0; i < nr; i++)
			{
				for(j = 0; j < nc; j++)
					{

						tmat[j*nr + i] = mat[i*nc + j];


					}

			}
//		RprintVecAsMat("Uk", mat, 1, (nr)*(nc));

	}




// factorial function using recursive arguments
int factorial(int n){
//	Rprintf("n = %d\n", n);
	int fact = 1;
	if(n < 0){ Rprintf("Cannot compute factorial of negative number "); return(0);}
	if(n <= 1 ){
		return(1);
	}else{
		fact = n*factorial(n-1);
	}
	return(fact);
}






//============================================================
// Random Number Generation
//

void ran_mvnorm(double* m, double* cholV, int dim, double* z, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a draw from a multivariate normal distribution.
 *
 * INPUTS:
 * m     = an array for the mean
 * cholV = an array for the cholesky decomposition
 *         of the variance (note this is a 1D array that
 *         that will be treated as 2D).
 * dim   = dimension of the multivariate normal distribution
 * z     = a scratch array of length dim
 *         to be filled with N(0,1) draws
 *
 * OUTPUTS:
 * out   = final output array to be filled with a random
 *         multivariate normal draw
 *
 *************************************************************/
    int i,j;
    for (i=0; i<dim; i++){
	z[i] = rnorm(0,1);
	out[i] = m[i];
	for (j=0; j<=i; j++)
	    out[i] += cholV[i*dim + j]*z[j];
    }
}


void ran_wish(int nu, double* cholS, int dim, double* z, double* x, double* zeros, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a random draw from Wishart distribution with
 * degrees of freedom nu and parameter matrix S
 *
 * INPUTS:
 * nu    = degrees of freedom
 * cholS = cholesky decomposition of matrix parameter S
 *         This is a 1D array but is accessed
 *         as a 2D array as in cholS[i*dim + j]
 * dim   = the dimension of the Wishart distribution
 * z     = scratch vector of length dim to be passed
 *         to the function ran_mvnorm
 * x     = scratch vector of length dim to be passed
 *         to the function ran_mvnorm
 * zeros = vector of zeros of length dim to be passed
 *         to the function ran_mvnorm as the mean
 * out   = matrix to contain the random draw from the
 *         Wishart distribution
 *
 *************************************************************/

    int i, j, k;

    /* Zero the "out" matrix */
    for (j=0; j<dim; j++)
	for (k=0; k<dim; k++)
	    out[j*dim + k] = 0.0;

    for (i=0; i<nu; i++){
	ran_mvnorm(zeros,cholS,dim,z,x);
	for (j=0; j<dim; j++)
	    for (k=0; k<=j; k++)
		out[j*dim + k] += x[j]*x[k];
    }

    /* fill the upper triangular part with lower triangular part */
    for (j=0; j<dim; j++)
	for (k=0; k<j; k++)
	    out[k*dim + j] = out[j*dim + k];
}


void ran_dirich(double *alpha, int dim, double* scratch, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a random draw from a Dirichlet distribution with
 * parameters alpha
 *
 * INPUTS:
 * alpha   = parameter vector associated with dirichlet
 * dim     = dimension of alpha vector
 * scratch = vector of length dim to be passed
 * out     = array that holds the random values
 *
 *************************************************************/

    int h;
	double sg = 0;
    /* Zero the "out" matrix */
//	RprintVecAsMat("alpha", alpha, 1, dim);
	for(h = 0; h < dim; h++)
		{

			scratch[h] = rgamma(alpha[h], 1);
//			Rprintf("rgamma = %f\n", scratch[h]);
			sg = sg + scratch[h];
//			Rprintf("sg = %f\n", sg);
		}

//	Rprintf("sg = %f\n", sg);
	for(h = 0; h < dim; h++) out[h] = scratch[h]/sg;

//	RprintVecAsMat("out", out, 1, dim);
}



double rinvgauss(double mu, double lambda){
/***************************************************************
 * PURPOSE:
 * Generate a random draw from an inverse gaussian distribution
 *
 * INPUTS:
 * mu   = mean (or location parameter) mu >= 0
 * lambda = variance (or scale parameter)
 *
 ***************************************************************/

	double y2, u, r2, r1, out;

	y2 = rchisq(1.0);
	u = runif(0,1);

	r2 = mu/(2*lambda)*(2*lambda + mu*y2 +
	                    sqrt(4*lambda*mu*y2 + mu*mu * y2*y2));
//	Rprintf("r2 = %f\n", r2);
	r1 = mu*mu/r2;
//	Rprintf("r1 = %f\n", r1);
	out = r2;
	if(u < mu/(mu+r1)) out = r1;
	return(out);
}



/* The following provides a function that allows me to sample from a
   truncated normal distribution with the following arguments.  This
   function relies heavly on the Rmath library.  As an example here
   are the arguments for the pnorm function
   double pnorm(double x, double mu, double sigma, int lower_tail,
				int give_log );


  	m - mean of the normal prior to truncation;
  	s - standard deviation of normal distribution prior to truncation;
  	a - Lower bound of the truncation
  	b - Upper bound of the truncation */

double rtnorm(double m, double s, double a, double b){

  int jj;
  double a_term, b_term, Un, p, rtdraw, tmp, z;

//	Rprintf("a = %f\n", a);
//	Rprintf("b = %f\n", b);

//	Rprintf("(a-m)/s = %f\n", (a-m)/s);
	a_term = pnorm((a-m)/s, 0.0, 1.0, 1, 0);

//	Rprintf("(b-m)/s = %f\n", (b-m)/s);
	b_term = pnorm((b-m)/s, 0.0, 1.0, 1, 0);

//	Rprintf("a_term = %f\n", a_term);
//	Rprintf("b_term = %f\n", b_term);

	Un = runif(0,1);

	p = a_term + Un*(b_term - a_term);

	rtdraw = m + s*qnorm(p, 0.0, 1.0, 1, 0);

//	Rprintf("rtdraw = %f\n", rtdraw);
//	Rprintf("p = %f\n", p);

	if(p == 1.0){
		jj = 1;

		while(jj != 0){

			a = (a - m)/s;
    		b = (b - m)/s;

//			Rprintf("a = %f\n", a);
			tmp = (a + sqrt(a*a + 4))/2;
//			Rprintf("tmp = %f\n", tmp);
        	z = rexp(1/tmp) + a;
//			Rprintf("z = %f\n", z);
        	Un = runif(0,1);
//			Rprintf("Un = %f\n", Un);
//			Rprintf("exp(-(z - tmp)*(z - tmp)/2) = %f\n", exp(-(z - tmp)*(z - tmp)/2));

        	if((Un <= exp(-(z - tmp)*(z - tmp)/2)) & (z <= b)) jj = 0;
//			Rprintf("jj = %d\n", jj);

		}

//		Rprintf("z = %f\n", z);

		rtdraw = z*s + m;


	}


//	Rprintf("rtdraw = %f\n", rtdraw);

	return(rtdraw);

}


//============================================================
// Multivariate density functions
//
/* The following provides a density function for the multivariate
   normal distribution.  This function relies heavily on matrix.c
   To use it I must create an matrix.o(object file)

	*y -  is the observation vector for which the density will be computed
  	*mu - mean vector
  	*iSig - the inverse covariance matrix as a contiguos vector in row-major form
  	dim - is the dimension of the multivariate distribution
	ld - the log of the determinant of Sigma
	scr - scratch vector to hold
	logout - logical determines if log density is returned
*/

//Multivariate normal
double dmvnorm(double* y, double* mu, double* iSig, int dim, double ld, double* scr, int logout){
    int i;
    double qf, out;
    for(i=0; i<dim; i++) scr[i] = y[i] - mu[i];
    qf = quform(scr,iSig,dim);
    out = -(double) dim*M_LN_SQRT_2PI - 0.5*(ld + qf);
    if (logout)	return out;
    return exp(out);
}






/* The following provides a density function for the Inverse Wishart
   distribution.  This function relies heavily on matrix.c
   To use it I must create an matrix.o(object file)

	*ASigInv - is ASigma^{-1} that is found in the trace of the density.  It is a
	           contiguos vector of memory that is in row-major form
	detSig - is the determinant of the argument matrix
	*A - is the dimxdim scale matrix that is in contiguos vector row-major form
	detA - is the determinant of A
  	df - degrees of freedom of the inverse-wishart function
  	dim - is the dimension of the scale matrix distribution */

double dinvwish(double *ASigInv, double detSig, double detA, int df, int dim)
	{

		int i;
		double C, density;
		double gamprod;
		double trace;
		double p1, p2, p3, p4;


		gamprod = 1;
		for(i = 0; i < dim; i++)
			{

				gamprod = gamprod*gammafn(0.5*(df + 1 - (i+1)));

			}
//		Rprintf("gamprod = %f\n", gamprod);
		trace = 0.0;
		for(i = 0; i < dim*dim; i++)
			{

				if(i % (dim+1) == 0) trace = trace + ASigInv[i];

			}
//		Rprintf("trace = %f\n", trace);
		p1 = 0.5*df*dim;
		p2 = 0.25*dim*(dim - 1);
		p3 = 0.5*df;
		p4 = -0.5*(df + dim + 1);
//		Rprintf("p1 = %f\n", p1);Rprintf("p2 = %f\n", p2);Rprintf("p3 = %f\n", p3);Rprintf("p4 = %f\n", p4);

		C = 1/(pow(2, p1)*pow(M_PI, p2)*gamprod);
//		Rprintf("detA = %f\n", detA);	Rprintf("detSig = %f\n", detSig);
//		Rprintf("C = %f\n", C);
		density = C*pow(detA, p3)*pow(detSig, p4)*exp(-0.5*trace);
//		Rprintf("density = %f\n", density);
		return(density);


	}




// Inverse Gamma density
// the parameterization here provides a mean of beta/(alpha - 1)
double dinvgamma(double y, double alpha, double beta, int logout){

//	Rprintf("alpha = %f\n", alpha);
//	Rprintf("beta = %f\n", beta);

	double ldens;

	ldens = alpha*log(beta) - lgamma(alpha) - (alpha + 1)*log(y) - beta*(1/y);

	if(logout) return ldens;
	return exp(ldens);

}




// density for univariate truncated normal
double dtnorm(double x, double mu, double sigma, double l, double u, int logout){

	double den, num, ldensity;

	den = pnorm5(u,mu,sigma, 1, 0) - pnorm5(l,mu,sigma, 1, 0);
	num = dnorm(x, mu, sigma,0);

	ldensity = log(num) - log(den);

	if(logout) return ldensity;
	return exp(ldensity);
}









