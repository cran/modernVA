% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDCurves.R
\name{tdVA}
\alias{tdVA}
\title{R wrapper that accesses C code to fit temporal dependent value-added models for two cohorts}
\description{
\code{tdVA} is the main function used to fit the temporally dependent value-added model for two cohorts
}

\usage{
tdVA(y1, xmat1,
     y2, xmat2,
     school1, school2,
     groupID=NULL,
     model=0,
     priors=c(0, 100^2, 1, 1, 1, 1, 0, 100^2, -1, 1, 0, 100^2, 0, 100^2),
     var.global=TRUE,
     MHsd=c(0.2),
     nchains=1,
     draws=50000, burn=40000, thin=10,
     verbose=FALSE)
}

\arguments{
\item{y1}{numeric vector (response variable) of length N1 for cohort 1.  Must be in long format.}

\item{y2}{numeric vector (response variable) of Length N2 for cohort 2.  Must be in long format.}

\item{xmat1}{N1 x p matrix of covariates for cohort 1 (column of 1's must NOT be included).}

\item{xmat2}{N2 x p matrix of covariates for cohort 2 (column of 1's must NOT be included).}

\item{school1}{numeric vector indicating to which school each student belongs for cohort 1.  These labels must be contiguous labels and start with 1}

\item{school2}{numeric vector indicating to which school each student belongs for cohort 2.  These labels must be contiguous labels and start with 1}

\item{groupID}{Optional vector that identifies to which group a school belongs.  If NULL there is no grouping}

\item{model}{Integer indicating which value-added model is to be fit

              0 - Independent school effects between the two cohorts.

              1 - Temporally dependent school effects between two cohorts based on a non-statinary AR(1) process,

              2 - Temporally dependent school effects based on previous cohorts post-test performance.

              3 - Full model that includes both an AR(1) type correlation and one based on previous cohorts post-test performance.}

\item{priors}{Vector of prior distribution parameter values.

  mb - prior mean for beta1 and beta2, default is 0.

  s2b - prior variance for beta1 and beta2, default is 100^2.

  at - prior shape for tau22 and tau21, default is 1.

  bt - prior rate for tau22 and tau21, default is 1.

  as - prior shape for sigma2, default is 1.

  bs - prior rate for sigma2, default is 1.

  mg - prior mean for gamma2, default is 0.  (only used if model = 2)

  s2g - prior variance for gamma2, default is 100^2.   (only used if model = 2)

  lp12 - prior lower bound for for phi12, default is -1.  (only used if model = 1)

  up12 - prior upper bound for for phi12, default is 1. (only used if model = 1)

  mp02 - prior mean for phi02, default is 0.

  s202 - prior variance for phi02, default is 100^2.

  mp01 - prior mean for phi01, default is 0.

  s201 - prior variance for phi01, default is 100^2.

}


\item{var.global}{Logical argument. If true, then a model with common sigma21 and sigma22 among schools is fit.  If false, then a model with school-specific sigma21i and sigma22i is fit.}

\item{MHsd}{Tuning parameter associated with M-H step of phi12.  Default is 0.2}

\item{nchains}{number of MCMC chains to run. Default is 1}

\item{draws}{number of MCMC iterates to be collected. default is 50,000}

\item{burn}{number of MCMC iterates discared as burn-in. default is 40,000}

\item{thin}{number by which the MCMC chain is thinne. default is 10}

\item{verbose}{Logical indicating if progress of MCMC algorithm should be printed to screen along with other data summaries}
}

\value{
This function returns a list that contains MCMC iterates for all the model parameters in addition to the value-added estimates and intervals for each of the two cohorts
}

\examples{
# Generate data from model 1 of San Martin el al.
m <- 25 # number of schools
ni <- 20 # number of students per school
N <- m*ni

# specify parameter values to generate data
beta1 <- 0.6
beta2 <- 0.75;
sig21 <- 100;
sig22 <- 100;
tau2 <- 100
phi02 <- 0;
phi12 <- 0.75;
phi01 <- 0

X1 <- rnorm(N, 0, sqrt(200))
X2 <- rnorm(N, 0, sqrt(200))

alpha1 <- rnorm(m, phi01, sqrt(tau2))
alpha2 <- rnorm(m, phi02 + phi12*alpha1, sqrt(tau2*(1-phi12^2)))

Y1 <- rep(alpha1, each=ni) + X1*beta1 + rnorm(N, 0, sqrt(sig21))
Y2 <- rep(alpha2, each=ni) + X2*beta2 + rnorm(N, 0, sqrt(sig22))

# Create school vector indicating to which school each observation belongs
school1 <- rep(1:m, each=ni)
school2 <- rep(1:m, each=ni)


# design matrix only one covariate and no intercept
X1i <- cbind(X1)
X2i <- cbind(X2)


fit <- tdVA(y1=Y1,xmat1=X1i,y2=Y2,xmat2=X2i,
              school1=school1,school2=school2, groupID=NULL, model=2,
              var.global=TRUE, nchains=1)

# Value-added estimates of cohort 1 and 2 with 95% credible intervals.  See paper for details
cbind(fit$VA1.estimate,t(fit$VA1.intervals))
cbind(fit$VA2.estimate,t(fit$VA2.intervals))


}
