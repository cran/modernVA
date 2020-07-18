/*************************************************************
 * Copyright (c) 2012 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit the quantile value added model in Page Orellana y San Martin
 * It employs quantile mixed model regression as described in  ()
 * Through this model the condicional quantile is readily available, but work was
 * necessary to estimate the marginal quantile.
 *
 *
 *************************************************************/
#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double fy(double y0, double lambda, double tau, double Xbeta, double sig2a){
	double out;
//	Rprintf("y0 = %f\n", y0);
//	Rprintf("lambda = %f\n", lambda);
//	Rprintf("tau = %f\n", tau);
//	Rprintf("Xbeta = %f\n", Xbeta);
//	Rprintf("sig2a = %f\n", sig2a);

	out = (1/lambda)*tau*(1-tau)*dnorm(y0, Xbeta, sqrt(sig2a), 0)*
	(
		1/(dnorm(y0, Xbeta + (sig2a/lambda)*(tau-1), sqrt(sig2a), 0) + 1e-308)*
		pnorm(0, y0-Xbeta - (sig2a/lambda)*(tau-1), sqrt(sig2a), 1, 0) +

		1/(dnorm(y0, Xbeta + (sig2a/lambda)*tau, sqrt(sig2a), 0) + 1e-308)*
		pnorm(0, y0-Xbeta - (sig2a/lambda)*tau, sqrt(sig2a), 0, 0)
	);
	return(out);
}



/***************************************************************************************************
* The following are the inputs of the function that are read from R
*
* nobsi = number of students per school
* nschool = number of schools
* ncov = number of covariates
* school = sum(nobsi) vector indicating to which school student belongs
* tau = the quantile for which quantile value added is being calculated
* y = vector containing response values (sum(nobsi))
* xmat = sum(nobsi) x ncov array containing covariate values
* y0 = vector that contains values that discretize the range of y
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned
* verbose = logical where information is printed to screen
***************************************************************************************************/


void mcmcloop_qva(int *draws, int *burn, int *thin, int *nobsi, int *nschool, int *ncov,
			  int *school, int *ny0,
			  double *tau, double *y, double *xmat, double *y0, double *priorVal,
			  int *verbose,
			  double *beta, double *alpha, double *v, double *a,
			  double *sig2a, double *lambda, double *cVA,
			  double *mVA, double *qVA, double *Q){


	// i - MCMC iterate;
	// ii - MCMC save iterate;
	// j - nschool iterate (through schools) ;
	// k - covariate iterate
	// kk - second covariate iterate
	// t - iterate that cycles through all obs
	// tt - second iterate that cycles through all obs

	int i, j, k, kk, t, tt;
  	int ii = 0;
	int N=0;
	for(j = 0; j < *nschool; j++){N = N + nobsi[j];}


	double *X = R_VectorInit((*ncov)*N, 0.0);
	double *tX = R_VectorInit((*ncov)*N, 0.0);

	for(t=0; t<N*(*ncov); t++) X[t] = xmat[t];

	mat_transpose(X, tX, N, *ncov);

//	RprintVecAsMat("X", X, N, *ncov);





	// =============================================================================================
	//
	// Memory vectors to hold a single MCMC iterate
	//
	// =============================================================================================

	double *beta_iter = R_VectorInit((*ncov), 0.5);
	double *alpha_iter = R_VectorInit(*nschool, 0.0);
	alpha_iter[0] = -10.0, alpha_iter[1] = -3.3333, alpha_iter[2] = 3.3333;
	alpha_iter[3] = 10.0;

	double *v_iter = R_VectorInit(N, 1.0);

	double lambda_iter = 0.1;
	double sig2a_iter = 62.5;
	double a_iter=0.0;

	double *cVA_iter = R_Vector(*nschool);
	double *mVA_iter = R_Vector(*nschool);
	double *mVAmc_iter = R_Vector(*nschool);
	double *qVA_iter = R_Vector(*nschool);
	double *qVAmc_iter = R_Vector(*nschool);

	double *Q_iter = R_Vector((N));
	double *Qmc_iter = R_Vector((N));

//	RprintVecAsMat("beta_iter", beta_iter, *nobs, *ncov);
//	RprintVecAsMat("alpha_iter", alpha_iter, 1, *nschool);




	// =============================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// =============================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector((*ncov)*(N));
	double *scr2 = R_Vector((*ncov)*(N));
	double *scr3 = R_Vector((*ncov)*N);


	//stuff I need to update alpha and a
	double sumv, sumvy_xb, xb, mstar, s2star, suma;

	//stuff I need to update lambda and sig2a
	double astar, bstar, ssq, mu;

	//stuff I need to update beta
	double ldo;
	double *u = R_Vector(N);
	double *Sstar = R_Vector((*ncov)*(*ncov));
	double *Mstar = R_Vector(N);
//	double *iV = R_Vector(N*N);
	double *tXiVX = R_Vector((*ncov)*(*ncov));
	double *tXiV = R_Vector((*ncov)*(N));


	//stuff I need to update cVA, mVA, qVA;
	double smargy, spmargy , qmargy;
	double *margy = R_Vector(*ny0);
	double *pmargy = R_Vector(*ny0);

	int ns=100, kp;
	double Z, xi, ald, qmargymc;
	double *margymc = R_Vector(ns);

	// =============================================================================================
	//
	// Prior distribution values;
	//
	// =============================================================================================

	// priorVal s2b=100, al=1, bl=1, as=1, bs=1, ma=0, s2a=100^2
	// priors for beta (mean = 0);
	double s2b = priorVal[0];
	double *iB0 = R_Vector((*ncov)*(*ncov));

	for(k=0; k<*ncov; k++){
		for(kk=0; kk<*ncov; kk++){
			iB0[k*(*ncov)+kk] = 0.0;
			if(k==kk)iB0[k*(*ncov)+kk] = (1/s2b);
		}
	}

	// priors for lambda
	double al=priorVal[1], bl = priorVal[2];

	// priors for sig2
	double as=priorVal[3], bs = priorVal[4];

	// prior for a
	double ma=priorVal[5], s2a = priorVal[6];

	if(*verbose){
	  Rprintf("number of schools = %d\n", *nschool);
    Rprintf("N = %d\n", N);

	  Rprintf("tau = %.2f\n", *tau);
	  Rprintf("ncov = %d\n\n", *ncov);

    Rprintf("Prior values being used are: \n s2b = %.1f \n al = %.1f, bl = %.1f,\n as = %.1f, bs = %.1f\n ma = %.1f, s2a = %.1f\n\n", s2b, al, bl, as, bs, ma, s2a);

	}



	GetRNGstate();


	// =============================================================================================
	//
	// start of the mcmc algorithm;
	//
	// =============================================================================================

	double calc_time = 0.0;
	clock_t  begin = clock();

	for(i = 0; i < *draws; i++){

//		if(*verbose & ((i+1) % 1000 == 0)){
		if(*verbose){
//			time_t now;
//			time(&now);

//		  Rprintf("mcmc iter = %d ===================================================== \n", i+1);
//      Rprintf("%s", ctime(&now));
      clock_t ith_iterate = clock();
		  calc_time = (ith_iterate - begin)/CLOCKS_PER_SEC;

      Rprintf("  Progress:%.1f%%, Time:%.1f seconds\r", ((double) (i+1) / (double) (*draws))*100.0, calc_time);
//    fflush(stdout);

		}



		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update alpha.  There are nschool of them (using gibbs sampler)
		//

		//////////////////////////////////////////////////////////////////////////////////
		for(j = 0; j < *nschool; j++){

			sumv = 0.0, sumvy_xb = 0.0;
			for(t = 0; t < N; t++){

				if(school[t] == j+1){
					xb = 0.0;
					for(k = 0; k < *ncov; k++){
						xb = xb + beta_iter[k]*xmat[t*(*ncov)+k];
					}

					sumv = sumv + v_iter[t];
					sumvy_xb = sumvy_xb + v_iter[t]*(y[t] - xb);
				}

			}


			s2star = 1.0/(0.5*(1.0/lambda_iter)*(*tau)*(1-(*tau))*sumv + (1/sig2a_iter));
			mstar = s2star*(0.5*(1.0/lambda_iter)*(*tau)*(1-(*tau))*sumvy_xb -
			         		0.5*(1.0/lambda_iter)*nobsi[j]*(1-2*(*tau)) +
			         		(1/sig2a_iter)*a_iter);

			alpha_iter[j] = rnorm(mstar, sqrt(s2star));


		}



		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update lambda (using gibbs-sampler)
		//
		//////////////////////////////////////////////////////////////////////////////////

		ssq = 0.0;
		sumv = 0.0;
		for(t = 0; t < N; t++){

			sumv = sumv + (1/v_iter[t]);

 			xb = 0.0;
			for(k = 0; k < *ncov; k++){
				xb = xb + beta_iter[k]*xmat[t*(*ncov)+k];
			}

			mu = (1-2*(*tau))/((*tau)*(1-(*tau))*v_iter[t]) + xb + alpha_iter[school[t]-1];
			ssq = ssq + v_iter[t]*(y[t]-mu)*(y[t]-mu);
		}


		astar = 1.5*(N) + al;
		bstar = 0.25*(*tau)*(1-(*tau))*ssq + (1/bl) + sumv;

		lambda_iter = 1/rgamma(astar, 1/bstar);


    //////////////////////////////////////////////////////////////////////////////////
    //
		// update sig2a. The variability of random effects.
		//
		//////////////////////////////////////////////////////////////////////////////////
		ssq = 0.0;
		for(j = 0; j < *nschool; j++){
			ssq = ssq + (alpha_iter[j]-a_iter)*(alpha_iter[j]-a_iter);
		}


		astar = as + 0.5*(*nschool);
		bstar = 0.5*ssq + 1/bs;

		sig2a_iter = 1/rgamma(astar, 1/bstar);

    //////////////////////////////////////////////////////////////////////////////////
    //
		// update a. The variability of random effects.
		//
		//////////////////////////////////////////////////////////////////////////////////
		suma = 0.0;
		for(j = 0; j < *nschool; j++){
			suma = suma + (alpha_iter[j]);
		}


		s2star = 1/(*nschool/sig2a_iter + 1/s2a);
		mstar = s2star*((1/sig2a_iter)*suma + (1/s2a)*ma);


		a_iter = rnorm(mstar, sqrt(s2star));

    ////////////////////////////////////////////////////////////////////////
    //
		// update v. This requires sampling from inverse gaussian.
		//
		////////////////////////////////////////////////////////////////////////

		for(t = 0; t < N; t++){
 			xb = 0.0;
			for(k = 0; k < *ncov; k++){
				xb = xb + beta_iter[k]*xmat[t*(*ncov)+k];
			}


			mstar = 1.0/((*tau)*(1-(*tau))*fabs(y[t]-xb-alpha_iter[school[t]-1]));
			s2star = 1.0/(lambda_iter*2*(*tau)*(1-(*tau)));

			v_iter[t] = rinvgauss(mstar, s2star);


		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// update beta
		//
		//////////////////////////////////////////////////////////////////////////////////

		for(t = 0; t < N; t++){
			u[t] = y[t] - (1-2*(*tau))/((*tau)*(1-*tau))*(1/v_iter[t]) - alpha_iter[school[t]-1];

			for(k = 0; k < *ncov; k++){

				tXiV[k*N+t] = X[t*(*ncov) + k]*((*tau)*(1-(*tau)))/(2*lambda_iter) * v_iter[t];

			}
		}

		// So that I don't have to create a huge diagonal matrix I am performing the matrix
		// product of X'V^{-1}X here.
		for(k = 0; k < (*ncov)*(*ncov); k++) tXiVX[k]=0.0;

		for(k = 0; k < *ncov; k++){
			for(kk = 0; kk < *ncov; kk++){
				for(t = 0; t < N; t++){

					tXiVX[k*(*ncov)+kk] = tXiVX[k*(*ncov)+kk] +
											X[t*(*ncov) + k]*X[t*(*ncov) + kk]*
											((*tau)*(1-(*tau)))/(2*lambda_iter) * v_iter[t];
				}

			}
		}

		matrix_product(tXiV, u, scr3, *ncov, 1, N);


		for(k=0; k < *ncov; k++){
			for(kk = 0; kk <*ncov; kk++){
				Sstar[k*(*ncov)+kk] = tXiVX[k*(*ncov)+kk] + iB0[k*(*ncov)+kk];
			}
		}


		cholesky(Sstar, (*ncov), &ldo);
		inverse_from_cholesky(Sstar, scr1, scr2, (*ncov));

		matrix_product(Sstar, scr3, Mstar, *ncov, 1.0, *ncov);

		cholesky(Sstar, *ncov, &ldo);

 		ran_mvnorm(Mstar, Sstar, *ncov, scr2, beta_iter);






		//////////////////////////////////////////////////////////////////////////////////
		//
		// calculate cVA, mVA and qVA and write to file
		//
		//////////////////////////////////////////////////////////////////////////////////
		for(j = 0; j < *nschool; j++){

			cVA_iter[j] = 0.0;
			mVA_iter[j] = 0.0;
			mVAmc_iter[j] = 0.0;

		}

		// conditional value added is simply Xbeta + alpha_j
		for(t = 0; t < N; t++){

			// Calculate Xbeta
 			xb = 0.0;
			for(k = 0; k < *ncov; k++){
				xb = xb + beta_iter[k]*xmat[t*(*ncov)+k];
			}

			// Calculate f(y0) where f is maringal of y;
			smargy = 0.0;
			for(tt = 0; tt < *ny0; tt++){

				margy[tt] = fy(y0[tt], lambda_iter, *tau, xb+a_iter, sig2a_iter);
				smargy = smargy + margy[tt];
			}

			// Calculate the densities of (is this correct)
			for(tt = 0; tt < *ny0; tt++){
				pmargy[tt] = margy[tt]/smargy;
			}


			// find the tau percentile (is this correct)
			spmargy = 0.0;
			for(tt = 0; tt < *ny0; tt++){

				spmargy = spmargy + pmargy[tt];

				if(spmargy >= *tau){

					qmargy = y0[tt+1];
					break;
				}
			}


			Q_iter[t] = qmargy;

			cVA_iter[school[t]-1] = cVA_iter[school[t]-1] +
									(xb + alpha_iter[school[t]-1])/(double) nobsi[school[t]-1];
			mVA_iter[school[t]-1] = mVA_iter[school[t]-1] +
									(qmargy)/(double) nobsi[school[t]-1];



			// Calculate marginal quantile using MC type arguments
			for(tt=0; tt < ns; tt++){

				Z = rnorm(0,1);
				xi = rexp(lambda_iter);
				ald = sqrt((2*xi*lambda_iter)/((*tau)*(1-(*tau))))*Z +
				       ((1-2*(*tau))/((*tau)*(1-(*tau))))*xi;
				margymc[tt] = rnorm(xb + ald, sqrt(sig2a_iter));

			}


			R_rsort(margymc,  ns);
			kp = (int) ns*(*tau)-1;
			qmargymc = margymc[kp];


			Qmc_iter[t] = qmargymc;

			mVAmc_iter[school[t]-1] = mVAmc_iter[school[t]-1] +
									 (qmargymc)/(double) nobsi[school[t]-1];



		}

		for(j = 0; j < *nschool; j++){

			qVA_iter[j] = cVA_iter[j] - mVA_iter[j];
			qVAmc_iter[j] = cVA_iter[j] - mVAmc_iter[j];

		}


		if((i > (*burn-1)) & (i % (*thin) == 0)){


			for(t = 0; t < N; t++){
				Q[ii*(N) + t] = Q_iter[t];
				v[ii*(N) + t] = v_iter[t];
		  	}

		  	for(j = 0; j < *nschool; j++){
		    	cVA[ii*(*nschool) + j] = cVA_iter[j];
		    	mVA[ii*(*nschool) + j] = mVA_iter[j];
		    	qVA[ii*(*nschool) + j] = qVA_iter[j];
		    	alpha[ii*(*nschool) + j] = alpha_iter[j];
		  	}

			for(k=0; k<*ncov; k++){
				beta[ii*(*ncov) + k] = beta_iter[k];
			}

			a[ii] = a_iter;
			sig2a[ii] = sig2a_iter;
			lambda[ii] = lambda_iter;

      		ii = ii+1;
		}


	}




	PutRNGstate();


}


