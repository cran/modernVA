/*************************************************************
 * Copyright (c) 2012 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit the time dependent value-added model.
 *
 *
 * We are not bounding the presistence parameter between (-1,1)
 * therefore the "time series" model is no stationary.
 * For more details look to the write up.
 *
 *************************************************************/
#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>




/***************************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned
* nobsi1 = number of students per school in time period 1
* nobsi2 = number of students per school in time period 2
* nschool = number of schools
* ncov = number of covariates
* school1 = sum(nobsi1) vector indicating to which school student belongs in 2nd time period
* school2 = sum(nobsi2) vector indicating to which school student belongs in 2nd time period
* y1 = vector containing response values (sum(nobsi1)) from time period 1
* xmat1 = sum(nobsi1) x ncov array containing covariate values
* y2 = vector containing response values (sum(nobsi2)) from time period 2
* xmat2 = sum(nobsi2) x ncov array containing covariate values
*
* ngroup = number of groups considered so as to borrow strength
* nci = number of schools in each group
* groupID = vector indicating to which group each school belongs
*
* priors = vector containing prior parameter values
* 			mub, s2b, at, bt, as, bs
* model = scalar indicating which model to use
* 		0 - Model0, phi1 = phi0 = gamma = 0
* 		1 - Model1, gamma = 0
*       2 - Model2, non parameters are equal to 0
*
*
* Output
* beta1 = nout x ncov matrix of MCMC iterates
* beta2 = nout x ncov matrix of MCMC iterates
* alpha1 = nout x nschool matrix of MCMC iterates
* alpha2 = nout x nschool matrix of MCMC iterates
* sig21 = nout vector of MCMC iterates
* sig22 = nout vector of MCMC iterates
* tau21 = nout vector of MCMC iterates
* tau22 = nout vector of MCMC iterates
* phi02 = nout x ngroup matrix  of MCMC iterates
* phi12 = nout x ngroup matrix  of MCMC iterates
* phi01 = nout x ngroup matrix  of MCMC iterates
* gamma = nout x ngroup matrix of MCMC iterates
* lpml = scalar with lpml value
* waic = scalar containing waic value
*
***************************************************************************************************/

void mcmcloopM2b1(int *draws, int *burn, int *thin, int *nobsi1, int *nobsi2,
			  int *nschool, int *ncov,
			  int *school1, int *school2,
			  double *y1, double *xmat1, double *y2, double *xmat2,
			  int *ngroup, int *groupID,
			  double *priors,	int *model, double *MHsd, int *verbose,
			  double *beta1, double *beta2, double *alpha1, double *alpha2,
			  double *sig21, double *sig22, double *tau21, double *tau22,
			  double *phi02, double *phi12, double *gamma, double *phi01,
			  double *lpml, double *waic){


	// i - MCMC iterate;
	// j - nschool iterate (through schools) ;
	// k - covariate iterate
	// kk - second covariate iterate
	// t - iterate that cycles through all obs
	// c - iterate through groupIDs to group phi



	int i, j, k, kk, t, c;
	int ii = 0;

	int nout = (*draws - *burn)/(*thin);

	int N1=0, N2=0;
  	for(j = 0; j < *nschool; j++){
		N1 = N1 + nobsi1[j];
		N2 = N2 + nobsi2[j];
	}



	double *X1 = R_VectorInit((*ncov)*N1, 0.0);
	double *X2 = R_VectorInit((*ncov)*N2, 0.0);

	double *tX1 = R_VectorInit((*ncov)*N1, 0.0);
	double *tX2 = R_VectorInit((*ncov)*N2, 0.0);

	for(t=0; t<N1*(*ncov); t++) X1[t] = xmat1[t];
	for(t=0; t<N2*(*ncov); t++) X2[t] = xmat2[t];

	mat_transpose(X1, tX1, N1, *ncov);
	mat_transpose(X2, tX2, N2, *ncov);



	// Create a vector of school specific means for cohort one (ybari.1)
	double y2_1=0.0;
	double *sY1X1 = R_VectorInit((*ncov), 0.0);
	double *tX1X1 = R_VectorInit((*ncov)*(*ncov), 0.0);
	double *ybari_1 = R_VectorInit(*nschool, 0.0);
	double *xbari_1 = R_VectorInit((*nschool)*(*ncov), 0.0);
	for(t = 0; t < N1; t++){
	  ybari_1[school1[t]-1] = ybari_1[school1[t]-1] + y1[t]/nobsi1[school1[t]-1];
		for(k = 0; k < *ncov; k++){
		  sY1X1[k] = sY1X1[k] + y1[t]*X1[t*(*ncov) + k];
			xbari_1[(school1[t]-1)*(*ncov) + k] = xbari_1[(school1[t]-1)*(*ncov) + k] + xmat1[t*(*ncov)+k]/nobsi1[(school1[t]-1)];
      for(kk=0; kk < *ncov; kk++){
        tX1X1[k*(*ncov) + kk] = tX1X1[k*(*ncov) + kk] + X1[t*(*ncov) + k]*X1[t*(*ncov) + kk];
      }
		}
		y2_1 = y2_1 + y1[t]*y1[t];
	}
	double sybar2 = 0.0;
	double sybar1=0.0;
	for(j = 0; j < *nschool; j++){
		sybar2 = sybar2 + ybari_1[j]*ybari_1[j];
		sybar1 = sybar1 + ybari_1[j];
	}


	// Create a vector of school specific means for cohort two (ybari.2)
	double y2_2=0.0;
	double *sY2X2 = R_VectorInit((*ncov), 0.0);
	double *tX2X2 = R_VectorInit((*ncov)*(*ncov), 0.0);
	double *ybari_2 = R_VectorInit(*nschool, 0.0);
	double *xbari_2 = R_VectorInit((*nschool)*(*ncov), 0.0);
	for(t = 0; t < N2; t++){
	  ybari_2[school2[t]-1] = ybari_2[school2[t]-1] + y2[t]/nobsi2[school2[t]-1];
		for(k = 0; k < *ncov; k++){
		  sY2X2[k] = sY2X2[k] + y2[t]*X2[t*(*ncov) + k];
			xbari_2[(school2[t]-1)*(*ncov) + k] = xbari_2[(school2[t]-1)*(*ncov) + k] + xmat2[t*(*ncov)+k]/nobsi2[(school2[t]-1)];
      for(kk=0; kk < *ncov; kk++){
        tX2X2[k*(*ncov) + kk] = tX2X2[k*(*ncov) + kk] + X2[t*(*ncov) + k]*X2[t*(*ncov) + kk];
      }
		}
		y2_2 = y2_2 + y2[t]*y2[t];
	}





	// =============================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// =============================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector((*ncov)*(N1));
	double *scr2 = R_Vector((*ncov)*(N1));
	double *scr3 = R_Vector((*ncov)*(N1));

	// stuff that I need to update gamma
	double sumy2, phi1sq, mstar, s2star, a1_mn, ma, va, xb;

	// stuff that I need to update phi
	double phio, phin, lln, llo, llr, uu, mno, mnn, varo, varn,suma;

	//stuff I need to update tau2 and sig2 time period 1 and 2
	double ssq1, ssq2, a2_1, a2_2;
	double asy_1, asy_2, asxb_1, asxb_2;
	double syxb_1, syxb_2, bxxb_1, bxxb_2;

	//stuff I need to update tau2 and sig2 time period 1 and 2
	double astar, bstar, ssq;


	//stuff I need to update beta1 and beta2
	double ldo;
	double *Sstar = R_Vector((*ncov)*(*ncov));
	double *Mstar = R_Vector(N1);


	// =============================================================================================
	//
	// Prior distribution values;
	//
	// =============================================================================================


	// priors for beta (mean = 0) and theta;
	double s2b = priors[1];
	double *iB0 = R_Vector((*ncov)*(*ncov));
	double *mub = R_VectorInit((*ncov), priors[0]);
	double *iB0mub = R_Vector((*ncov));

	for(k=0; k<*ncov; k++){
		for(kk=0; kk<*ncov; kk++){
			iB0[k*(*ncov)+kk] = 0.0;
			if(k==kk)iB0[k*(*ncov)+kk] = (1/s2b);
		}
	}

	matrix_product(iB0, mub, iB0mub, *ncov, 1, *ncov);

	//	RprintVecAsMat("iB0", iB0, *ncov, *ncov);

	// priors for both tau2's
	double at=priors[2], bt = priors[3];

	// priors for both sig2's
	double as=priors[4], bs = priors[5];

	// priors for phi0
	double mp = priors[10], s2p = priors[11];

	// priors for phi00
	double mp00 = priors[12], s2p00 = priors[13];

	// priors for phi1
	double lp1 = priors[8], up1 = priors[9];

	// priors for gamma
	double mg=priors[6], s2g=priors[7];



	// =============================================================================================
	//
	// Memory vectors to hold a single MCMC iterate
	//
	// =============================================================================================


	GetRNGstate();


	double *beta1_iter = R_VectorInit((*ncov), rnorm(priors[0], sqrt(priors[1])));
	double *beta2_iter = R_VectorInit((*ncov), rnorm(priors[0], sqrt(priors[1])));

	double *alpha1_iter = R_VectorInit(*nschool, rnorm(mp, sqrt(s2p)));
	double *alpha2_iter = R_VectorInit(*nschool, rnorm(mp, sqrt(s2p)));

	double sig21_iter = 1/rgamma(as, 1/bs);
	double sig22_iter = 1/rgamma(as, 1/bs);

	double tau21_iter = 1/rgamma(at, 1/bt);
	double tau22_iter = 1/rgamma(at, 1/bt);

	double *phi0_iter = R_VectorInit(*ngroup, rnorm(mp, sqrt(s2p)));
	double phi00_iter = rnorm(mp00, sqrt(s2p00));

	double *phi1_iter = R_VectorInit(*ngroup, runif(lp1, up1));
	if(*model != 1){
		for(c=0; c<*ngroup; c++) phi1_iter[c] = 0.0;
	}

	double *gamma_iter = R_VectorInit(*ngroup, rnorm(mg, sqrt(s2g)));
	if(*model != 2){
		for(c=0; c<*ngroup; c++) gamma_iter[c] = 0.0;
	}




	double *like_iter = R_VectorInit((N1+N2), 0.0);
	double *mnllike = R_VectorInit((N1+N2), 0.0);
	double *mnlike = R_VectorInit((N1+N2), 0.0);
	double *CPO = R_VectorInit((N1+N2), 0.0);

	double lpml_iter;
	double elppdWAIC;



	// Stuff that I need for MH step
	double csigPHI1 = MHsd[0];


	if(*verbose){
		Rprintf("N1 = %d\n", N1);
		Rprintf("N2 = %d\n", N2);

		Rprintf("nschool = %d\n", *nschool);
		Rprintf("ncov = %d\n", *ncov);
		Rprintf("ngroup = %d\n", *ngroup);

		Rprintf("sybar1 = %f\n", sybar1);
		Rprintf("sybar1^2 = %f\n", sybar2);

		RprintVecAsMat("iB0", iB0, *ncov, *ncov);
		Rprintf("csigPHI1 = %f\n\n", csigPHI1);

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
		// Update gamma.  quality of school parameter
		//
		//////////////////////////////////////////////////////////////////////////////////
		if((*model == 2) | (*model==3)){
			for(c = 0; c < *ngroup; c++){

				phi1sq = phi1_iter[c]*phi1_iter[c];
				sumy2 = 0.0;
				a1_mn = 0.0;
				for(j = 0; j < (*nschool); j++){

					sumy2 = sumy2 + ybari_1[j]*ybari_1[j];
					a1_mn = a1_mn +  ybari_1[j]*(alpha2_iter[j] - phi0_iter[c] -
					                                    phi1_iter[c]*alpha1_iter[j]);
				}


				s2star = 1.0/(1.0/(tau22_iter*(1-phi1sq))*sumy2  +  1/s2g);
				mstar = s2star*(1.0/(tau22_iter*(1-phi1sq))*a1_mn + mg*(1/s2g));

//				Rprintf("mstar = %f\n", mstar);
//				Rprintf("s2star = %f\n", s2star);

				gamma_iter[c] = rnorm(mstar, sqrt(s2star));

			}
		}
//		RprintVecAsMat("gamma_iter = %f\n", gamma_iter);


		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update phi1.  time dependence parameter
		//
		//////////////////////////////////////////////////////////////////////////////////
		if((*model == 1) | (*model==3)){
			for(c = 0; c < (*ngroup); c++){
				phio = phi1_iter[c];
				phin = rnorm(phio, csigPHI1);
				if((phin < up1) & (phin > lp1)){
					lln = 0.0, llo = 0.0;
					for(j = 0; j < (*nschool); j++){
						if(groupID[j] == c+1){

							varo = tau22_iter*(1-phio*phio);
							varn = tau22_iter*(1-phin*phin);

							mno = phi0_iter[c] + phio*alpha1_iter[j] + gamma_iter[c]*(ybari_1[j]);
							mnn = phi0_iter[c] + phin*alpha1_iter[j] + gamma_iter[c]*(ybari_1[j]);

							llo = llo + dnorm(alpha2_iter[j], mno, sqrt(varo), 1);
							lln = lln + dnorm(alpha2_iter[j], mnn, sqrt(varn), 1);
						}
					}

					llo = llo + dunif(phio, lp1,up1,1);
					lln = lln + dunif(phin, lp1,up1,1);

//					llo = llo + dtnorm(phio, theta_iter,sqrt(nu2_iter),-1.0,1.0,1);
//					lln = lln + dtnorm(phin, theta_iter,sqrt(nu2_iter),-1.0,1.0,1);

//					Rprintf("llo = %f\n", llo);
//					Rprintf("lln = %f\n", lln);

					llr = lln - llo;

					uu = runif(0.0,1.0);

					if(llr > log(uu)) phi1_iter[c] = phin;

				}
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update phi0.  time intercept
		//
		//////////////////////////////////////////////////////////////////////////////////

		for(c = 0; c < (*ngroup); c++){

			a1_mn = 0.0;
			for(j = 0; j < (*nschool); j++){
				if(groupID[j] == c+1){

					a1_mn = a1_mn +  (alpha2_iter[j] - phi1_iter[c]*alpha1_iter[j]-
						                                       gamma_iter[c]*ybari_1[j]);
				}
			}

			phi1sq = phi1_iter[c]*phi1_iter[c];


			s2star = 1.0/(*nschool/(tau22_iter*(1-phi1sq))  +  1/s2p);
			mstar = s2star*(1.0/(tau22_iter*(1-phi1sq))*a1_mn + mp*(1/s2p));

			phi0_iter[c] = rnorm(mstar, sqrt(s2star));
		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update phi00.  mean of alpha1i
		//
		//////////////////////////////////////////////////////////////////////////////////

		suma = 0.0;
		for(j = 0; j < (*nschool); j++){

			suma = suma +  alpha1_iter[j];
		}

		s2star = 1.0/(*nschool/tau21_iter  +  1/s2p00);
		mstar = s2star*((1.0/tau21_iter)*suma + mp00*(1/s2p00));

//			Rprintf("mstar = %f\n", mstar);
//			Rprintf("s2star = %f\n", s2star);

		phi00_iter = rnorm(mstar, sqrt(s2star));

/*
    //////////////////////////////////////////////////////////////////////////////////
    //
		// update sig21. The variability associated with y1.
		//
		//////////////////////////////////////////////////////////////////////////////////
		ssq = 0.0;
		for(t = 0; t < N1; t++){

	 		xb = 0.0;
			for(k = 0; k < *ncov; k++){
				xb = xb + beta1_iter[k]*xmat1[t*(*ncov)+k];
			}
			ssq = ssq + (y1[t]-xb-alpha1_iter[school1[t]-1])*
			            (y1[t]-xb-alpha1_iter[school1[t]-1]);
		}


		astar = as + 0.5*(N1);
		bstar = 0.5*ssq + 1/bs;

		sig21_iter = 1/rgamma(astar, 1/bstar);


    //////////////////////////////////////////////////////////////////////////////////
    //
		// update sig22. The variability associated with y2.
		//
		//////////////////////////////////////////////////////////////////////////////////
		ssq = 0.0;
		for(t = 0; t < N2; t++){

 			xb = 0.0;
			for(k = 0; k < *ncov; k++){
				xb = xb + beta2_iter[k]*xmat2[t*(*ncov)+k];
			}

			ssq = ssq + (y2[t]-xb-alpha2_iter[school2[t]-1])*
			            (y2[t]-xb-alpha2_iter[school2[t]-1]);
		}

		astar = as + 0.5*(N2);
		bstar = 0.5*ssq + 1/bs;

		sig22_iter = 1/rgamma(astar, 1/bstar);
*/
    //////////////////////////////////////////////////////////////////////////////////
    //
		// update sig21 and sig22 at the same time.
		//
		//////////////////////////////////////////////////////////////////////////////////

		a2_1 = 0.0, a2_2=0.0;
		asy_1 = 0.0, asy_2=0.0;
		asxb_1 =0.0, asxb_2=0.0;
		for(j=0; j<(*nschool); j++){

		  a2_1 = a2_1 + nobsi1[j]*alpha1_iter[j]*alpha1_iter[j];
		  asy_1 = asy_1 + nobsi1[j]*ybari_1[j]*alpha1_iter[j];

		  a2_2 = a2_2 + nobsi2[j]*alpha2_iter[j]*alpha2_iter[j];
		  asy_2 = asy_2 + nobsi2[j]*ybari_2[j]*alpha2_iter[j];

		  for(k=0; k< *ncov; k++){
		    asxb_1 = asxb_1 + nobsi1[j]*xbari_1[j*(*ncov)+k]*alpha1_iter[j]*beta1_iter[k];
		    asxb_2 = asxb_2 + nobsi2[j]*xbari_2[j*(*ncov)+k]*alpha2_iter[j]*beta2_iter[k];
		  }
		}
		syxb_1=0.0, syxb_2=0.0;
		bxxb_1=0.0, bxxb_2=0.0;
		for(k=0; k<*ncov; k++){
		  syxb_1 = syxb_1 + sY1X1[k]*beta1_iter[k];
		  syxb_2 = syxb_2 + sY2X2[k]*beta2_iter[k];
		  for(kk=0; kk<*ncov; kk++){
		    bxxb_1 = bxxb_1 + beta1_iter[k]*tX1X1[k*(*ncov) + kk]*beta1_iter[kk];
		    bxxb_2 = bxxb_2 + beta2_iter[k]*tX2X2[k*(*ncov) + kk]*beta2_iter[kk];
		  }
		}

		ssq1 = y2_1 + bxxb_1 + a2_1 - 2*asy_1 - 2*syxb_1 + 2*asxb_1;
		ssq2 = y2_2 + bxxb_2 + a2_2 - 2*asy_2 - 2*syxb_2 + 2*asxb_2;

//    Rprintf("ssq = %f\n", ssq);

		astar = as + 0.5*(N1);
		bstar = 0.5*ssq1 + 1/bs;

		sig21_iter = 1/rgamma(astar, 1/bstar);


		astar = as + 0.5*(N2);
		bstar = 0.5*ssq2 + 1/bs;

		sig22_iter = 1/rgamma(astar, 1/bstar);



    //////////////////////////////////////////////////////////////////////////////////
    //
		// update tau2. The variability associated with random effects.
		//
		//////////////////////////////////////////////////////////////////////////////////
/*
		ssq = 0.0;

		for(j = 0; j < *nschool; j++){
			ssq = ssq + (alpha1_iter[j] - phi00_iter)*(alpha1_iter[j] - phi00_iter) +
				           (1/(1-phi1_iter[groupID[j]-1]*phi1_iter[groupID[j]-1]))*
			    	       (alpha2_iter[j] - phi0_iter[groupID[j]-1] - phi1_iter[groupID[j]-1]*alpha1_iter[j] - gamma_iter[groupID[j]-1]*ybari_1[j])*
			       		   (alpha2_iter[j] - phi0_iter[groupID[j]-1] - phi1_iter[groupID[j]-1]*alpha1_iter[j] - gamma_iter[groupID[j]-1]*ybari_1[j]);
		}

//			Rprintf("ssq = %f\n", ssq);

		astar = at + (*nschool);
		bstar = 0.5*ssq + 1/bt;

//			Rprintf("astar = %f\n", astar);
//			Rprintf("bstar = %f\n", bstar);

		tau2_iter = 1/rgamma(astar, 1/bstar);
*/
    //////////////////////////////////////////////////////////////////////////////////
    //
		// update tau21. The variability associated with alpha1.
		//
		//////////////////////////////////////////////////////////////////////////////////

		ssq = 0.0;

		for(j = 0; j < *nschool; j++){

		  ssq = ssq + (alpha1_iter[j] - phi00_iter)*(alpha1_iter[j] - phi00_iter) ;
		}


		astar = at + 0.5*(*nschool);
		bstar = 0.5*ssq + 1/bt;

		tau21_iter = 1/rgamma(astar, 1/bstar);


    //////////////////////////////////////////////////////////////////////////////////
    //
		// update tau22. The variability associated with alpha1.
		//
		//////////////////////////////////////////////////////////////////////////////////

		ssq = 0.0;

		for(j = 0; j < *nschool; j++){


		  ssq = ssq +   (1.0/(1.0-phi1_iter[groupID[j]-1]*phi1_iter[groupID[j]-1]))*
				    	       (alpha2_iter[j] - phi0_iter[groupID[j]-1] - phi1_iter[groupID[j]-1]*alpha1_iter[j] - gamma_iter[groupID[j]-1]*ybari_1[j])*
				        	   (alpha2_iter[j] - phi0_iter[groupID[j]-1] - phi1_iter[groupID[j]-1]*alpha1_iter[j] - gamma_iter[groupID[j]-1]*ybari_1[j]);
		}


		astar = at + 0.5*(*nschool);
		bstar = 0.5*ssq + 1/bt;

		tau22_iter = 1/rgamma(astar, 1/bstar);



		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update alpha1.  School random effect in time period 1
		//
		//////////////////////////////////////////////////////////////////////////////////
		for(j = 0; j < *nschool; j++){

//			Rprintf("j = %d\n", j);

			phi1sq = phi1_iter[groupID[j]-1]*phi1_iter[groupID[j]-1];
			va = 1.0/tau21_iter + phi1sq/(tau22_iter*(1-phi1sq));

			s2star = 1.0/((nobsi1[j]/sig21_iter) +  va);

			ma = (1.0/tau21_iter)*phi00_iter +
			     (1.0/(tau22_iter*(1-phi1sq)))*(phi1_iter[groupID[j]-1]*alpha2_iter[j] -
			                                   phi1_iter[groupID[j]-1]*phi0_iter[groupID[j]-1] -
			                                   phi1_iter[groupID[j]-1]*gamma_iter[groupID[j]-1]*ybari_1[j]);
			xb = 0.0;
			for(k=0; k<*ncov; k++){
			  xb = xb + xbari_1[j*(*ncov) + k]*beta1_iter[k];
			}

			mstar = s2star*((nobsi1[j]/sig21_iter)*(ybari_1[j] - xb) + ma);

			alpha1_iter[j] = rnorm(mstar, sqrt(s2star));

		}

//		RprintVecAsMat("alpha1_iter", alpha1_iter, 1, *nschool) ;


		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update alpha2.  School random effect in time period 2 that depends on alpha1
		//
		//////////////////////////////////////////////////////////////////////////////////
		for(j = 0; j < *nschool; j++){

//			Rprintf("j = %d ============ \n", j);

			va = tau22_iter*(1-phi1_iter[groupID[j]-1]*phi1_iter[groupID[j]-1]);

			s2star = 1.0/((nobsi2[j]/sig22_iter) + (1.0/va));

			xb = 0.0;
			for(k=0; k<*ncov; k++){
			  xb = xb + xbari_2[j*(*ncov) + k]*beta2_iter[k];
			}

			mstar = s2star*((nobsi2[j]/sig22_iter)*(ybari_2[j] - xb) +
			                (1.0/va)*
			                (phi0_iter[groupID[j]-1] + phi1_iter[groupID[j]-1]*alpha1_iter[j] + gamma_iter[groupID[j]-1]*ybari_1[j]));


//			Rprintf("mstar = %f\n", mstar);
//			Rprintf("s2star = %f\n", s2star);

			alpha2_iter[j] = rnorm(mstar, sqrt(s2star));


//			Rprintf("alpha2_iterj = %f\n", alpha2_iter[j]);


		}

//		RprintVecAsMat("alpha2_iter", alpha2_iter, 1, *nschool) ;




		//////////////////////////////////////////////////////////////////////////////////
		//
		// update beta1.  Slope associated with first time period
		//
		//////////////////////////////////////////////////////////////////////////////////
		for(k=0; k < *ncov; k++){
		  scr3[k] = 0.0;
		  for(j=0; j<*nschool; j++){
		    scr3[k] = scr3[k] + (1.0/sig21_iter)*(xbari_1[j*(*ncov) + k]*nobsi1[j]*alpha1_iter[j]);
		  }
		  scr3[k] = (1.0/sig21_iter)*sY1X1[k] - scr3[k] + iB0mub[k];
			for(kk = 0; kk <*ncov; kk++){
				Sstar[k*(*ncov)+kk] = (1/sig21_iter)*tX1X1[k*(*ncov)+kk] +  iB0[k*(*ncov)+kk];
			}
		}

		cholesky(Sstar, (*ncov), &ldo);
		inverse_from_cholesky(Sstar, scr1, scr2, (*ncov));

		matrix_product(Sstar, scr3, Mstar, *ncov, 1.0, *ncov);

//		RprintVecAsMat("Mstar", Mstar, 1, *ncov);
//		RprintVecAsMat("Sstar", Sstar, *ncov, *ncov);

		cholesky(Sstar, *ncov, &ldo);

 		ran_mvnorm(Mstar, Sstar, *ncov, scr2, beta1_iter);

//		RprintVecAsMat("beta1_iter", beta1_iter, 1, *ncov);




		//////////////////////////////////////////////////////////////////////////////////
		//
		// update beta2.  Slope associated with first time period
		//
		//////////////////////////////////////////////////////////////////////////////////
		for(k=0; k < *ncov; k++){
		  scr3[k] = 0.0;
		  for(j=0; j<*nschool; j++){
		    scr3[k] = scr3[k] + (1.0/sig22_iter)*(xbari_2[j*(*ncov) + k]*nobsi2[j]*alpha2_iter[j]);
		  }
		  scr3[k] = (1.0/sig22_iter)*sY2X2[k] - scr3[k] + iB0mub[k];
			for(kk = 0; kk <*ncov; kk++){
				Sstar[k*(*ncov)+kk] = (1/sig22_iter)*tX2X2[k*(*ncov)+kk] +  iB0[k*(*ncov)+kk];
			}
		}

//		RprintVecAsMat("scr3", scr3, 1, *ncov);
//		RprintVecAsMat("Sstar", Sstar, *ncov, *ncov);

		cholesky(Sstar, (*ncov), &ldo);
		inverse_from_cholesky(Sstar, scr1, scr2, (*ncov));

		matrix_product(Sstar, scr3, Mstar, *ncov, 1.0, *ncov);

//		RprintVecAsMat("Mstar", Mstar, 1, *ncov);
//		RprintVecAsMat("Sstar", Sstar, *ncov, *ncov);

		cholesky(Sstar, *ncov, &ldo);

 		ran_mvnorm(Mstar, Sstar, *ncov, scr2, beta2_iter);


//		RprintVecAsMat("beta2_iter", beta2_iter, 1, *ncov);

		//////////////////////////////////////////////////////////////////////////////////
		//
		// update likelihood to compute LPML and WAIC and
		//
		//////////////////////////////////////////////////////////////////////////////////

		lpml_iter=0.0;
		if((i > (*burn-1)) & (i % (*thin) == 0)){
			for(t = 0; t < N1; t++){

 				xb = 0.0;
				for(k = 0; k < *ncov; k++){
					xb = xb + beta1_iter[k]*xmat1[t*(*ncov)+k];
				}

//				Rprintf("xb = %f\n", xb);

				like_iter[t] = dnorm(y1[t], xb + alpha1_iter[school1[t]-1], sqrt(sig21_iter), 0);
//				Rprintf("like_iter = %f\n", like_iter[j]);
//				Rprintf("like_iter = %40.9f\n", like_iter[j]);

				// These are needed for WAIC
				mnlike[t] = mnlike[t] + (like_iter[t])/(double) nout;
				mnllike[t] = mnllike[t] + log(like_iter[t])/(double) nout;

				CPO[t] = CPO[t] + (1/(double) nout)*
				                  (1/like_iter[t]);

//				Rprintf("CPO = %f\n", CPO[j]);


			}

			for(t = 0; t < N2; t++){

 				xb = 0.0;
				for(k = 0; k < *ncov; k++){
					xb = xb + beta2_iter[k]*xmat2[t*(*ncov)+k];
				}

//				Rprintf("xb = %f\n", xb);

				like_iter[N1+t] = dnorm(y2[t], xb + alpha2_iter[school2[t]-1], sqrt(sig22_iter), 0);
//				Rprintf("like_iter = %f\n", like_iter[j]);
//				Rprintf("like_iter = %40.9f\n", like_iter[j]);

				// These are needed for WAIC
				mnlike[N1+t] = mnlike[N1+t] + (like_iter[N1+t])/(double) nout;
				mnllike[N1+t] = mnllike[N1+t] + log(like_iter[N1+t])/(double) nout;

				CPO[N1+t] = CPO[N1+t] + (1/(double) nout)*
				                  (1/like_iter[N1+t]);

//				Rprintf("CPO = %f\n", CPO[j]);


			}

//			if(i == (*draws-1)) Rprintf("xb = %f\n", xb);

		}



		//////////////////////////////////////////////////////////////////////////////////
		//
		// Save MCMC iterates
		//
		//////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin ==0)){

//			Rprintf("ii = %d\n", ii);


			tau21[ii] = tau21_iter;
			tau22[ii] = tau22_iter;


			sig21[ii] = sig21_iter;
			sig22[ii] = sig22_iter;

			phi01[ii] = phi00_iter;

			for(c = 0; c < *ngroup; c++){
				phi12[ii*(*ngroup) + c] = phi1_iter[c];
				phi02[ii*(*ngroup) + c] = phi0_iter[c];
				gamma[ii*(*ngroup) + c] = gamma_iter[c];
			}

			for(j = 0; j < *nschool; j ++){

				alpha1[ii*(*nschool) + j] = alpha1_iter[j];
				alpha2[ii*(*nschool) + j] = alpha2_iter[j];


			}
			for(k = 0; k < *ncov; k++){
				beta1[ii*(*ncov) + k] = beta1_iter[k];
				beta2[ii*(*ncov) + k] = beta2_iter[k];

			}


			ii = ii+1;

		}



	}


	//////////////////////////////////////////////////////////////////////////////////
	//
	// Computing LPML and WAIC  (see Gelman article in lit review folder)
	// An issue with WAIC is that considering spatially structure data
	//
	//////////////////////////////////////////////////////////////////////////////////

	lpml_iter = 0.0;
	elppdWAIC = 0.0;


	for(t = 0; t < (N1+N2); t++){

		lpml_iter = lpml_iter + log(1/CPO[t]);
		elppdWAIC = elppdWAIC + (2*mnllike[t] - log(mnlike[t]));

	}

	lpml[0] = lpml_iter;
	waic[0] = -2*elppdWAIC;



	PutRNGstate();

}


