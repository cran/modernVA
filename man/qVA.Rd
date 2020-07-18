\name{qVA}
\alias{qVA}
\title{R wrapper that accesses C code to fit quantile value-added model}
\description{
\code{qVA} is the main function used to fit hierarchical model that produces quantile value-added estimates.
}
\usage{
qVA(y,
    xmat,
    school,
    tau=0.5,
    draws=1100, burn=100, thin=1,
    priorVal=c(100^2, 1, 1, 1, 1, 0, 100^2),
    verbose=FALSE)
}

\arguments{
\item{y}{numeric vector response variable.  Must be in long format.}

\item{xmat}{N x p matrix of covariates (column of 1's must NOT be included) where N is total number of observations.}

\item{school}{vector indicating to which school each student belongs.  School labels must be contiguous and start with 1.}

\item{tau}{quantile specification.  The median is used as a default (tau=0.5)}

\item{priorVal}{vector of prior distribution parameter values.

  s2b - prior variance for beta, default is 100^2.

  al - prior shape for lambda, default is 1.

  bl - prior rate for lambda, default is 1.

  as - prior shape for sigma2, default is 1.

  bs - prior rate for sigma2, default is 1.

  ma - mean for a, default is 0.

  s2a - variance for a, default is 100^2
}

\item{draws}{total number of MCMC iterates to be collected. default is 1100}

\item{burn}{number of total MCMC iterates discared as burn-in. default is 100}

\item{thin}{number by which the MCMC chain is thinned. default is 1.  Note that the number of MCMC iterates provided is (draws - burn)/thin.}

\item{verbose}{Logical indicating if MCMC progress and other data summaries should be printed to screen}
}
\value{
This function returns a list containing MCMC iterates that correspond to model parameters and school-specific quantile value-added estimates.  In order to provide more detail, in what follows let

"T" - be the number of MCMC iterates collected (draws - burn)/thin,

"M" - be the number of schools,

"N" - be the total number of observations.

"p" - be the number of covariates

The output list contains the following

beta - an matrix of dimension (T, p) containing MCMC iterates associated with quantile regression covariate estimates.

alpha - an matrix of dimension (T, M) containing MCMC iterates assocated with school-specific random effects.

v - matrix of dimension (T, N) containing MCMC iterates of auxiliary variable.

a - matrix of dimension (T, 1) contaning MCMC iterates of mean of the school-specific random effects (alpha).

sig2a - a matrix of dimension (T, 1) containing MCMC iterates of variance of the school-specific random effects (alpha).

lambda - a matrix of dimension (T, M) containing MCMC interates associated with the lambda parameter of asymmetric laplace distribution.

cVA - a matrix of dimension (T, M) containing MCMC interates associated with conditional value-added for each school.

mVA - a matrix of dimension (T, M) containing MCMC iterates associated with the marginal value-added for each school.

qVA - a matrix of dimension (T, M) containing MCMC iterates associated with each schools quantile value-added.

Q - a matrix of dimension (T, N) containing MCMC iterates associated with the marginal quantile valued-added regression value for each student (i.e., averaging over school).

}
\references{
Page, Garritt L.; San Martín, Ernesto; Orellana, Javiera; Gonzalez, Jorge. (2017) ``Exploring Complete School Effectiveness via Quantile Value-Added'' Journal of the Royal Statistical Society: Series A 180(1) 315-340
}

\examples{
# Example with synthetic data
tau <- 0.75

m <- 4 # number of schools
n <- 25 # number of students
N <- m*n
p <- 1 # number of covariates
betaT <- 0.5
alphaT <- seq(-10,10, length=m)

# Generate from the asymmetric Laplace
# using a mixture of a Normal and an Exponential
lambdaT <- 0.1;
xi <- rexp(N, 1/lambdaT)
epsilon <- (sqrt((lambdaT*2*xi)/(tau*(1-tau)))*rnorm(N,0,1) +
           (1-2*tau)/(tau*(1-tau))*xi)
epsilon <- rnorm(N,0,1)
alphavec <-  rep(alphaT, each=n)

x <- rnorm(N,250,1)
y <- x*betaT  + alphavec + epsilon

X <- cbind(x)


school <- rep(1:m, each=n)


\donttest{
fitQ3 <- qVA(y=y, xmat=X, school=school, tau=0.75, verbose=FALSE)

# quantile value-added estimates with 95% credible intervals for each school
qVA.est <- apply(fitQ3$qVA,2,mean)
qVA.int <- apply(fitQ3$qVA,2,function(x) quantile(x, c(0.025, 0.975)))

beta <- fitQ3$beta
alpha <- fitQ3$alpha

mVA <- fitQ3$mVA
cVA <- fitQ3$cVA
Q <- fitQ3$Q

# Plot results.
plot(x,y, col=rep(c("red","blue","green","orange"), each=n), pch=19)

# Plot Q3 quantile regression line for each school
lines(X[school==1,],
       (X[school==1,])*mean(beta) + apply(alpha,2,mean)[1], col='red', lwd=3)
lines(X[school==2,],
       (X[school==2,])*mean(beta) + apply(alpha,2,mean)[2], col='blue', lwd=3)
lines(X[school==3,],
       (X[school==3,])*mean(beta) + apply(alpha,2,mean)[3], col='green', lwd=3)
lines(X[school==4,],
       (X[school==4,])*mean(beta) + apply(alpha,2,mean)[4], col='orange', lwd=3)

# Plot the marginal VA for each school
points(tapply(X, school,mean), apply(mVA,2,mean),
       col=c("red","blue","green","orange"), pch=4, cex=2, lwd=2)

# Plot the conditional VA for each school
points(tapply(X, school,mean), apply(cVA,2,mean),
      col=c("red","blue","green","orange"),pch=10,cex=2, lwd=2)

# Plot the "global" Q3 quantile regression line.
points(X, apply(Q,2,mean), type='l', lwd=2)
}
}
