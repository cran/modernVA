

# Wrapper function that access C-code that fits the quantile value-added model

qVA <- function(y,xmat,school,tau=0.5,draws=1100,burn=100,thin=1,
                priorVal=c(100^2, 1, 1, 1, 1, 0, 100^2), verbose=FALSE){

  # priorVal is vector of prior values as follows
  # s2b - prior variance for beta
  # al - prior shape for lambda
  # bl - prior rate for lambda
  # as - prior shape for sigma2
  # bs - prior rate for sigma2
  # ma - mean for a
  # s2a - variance for a

  out <- NULL

  nobsi = table(school)
  nschool = length(nobsi)
  ncov = dim(xmat)[2]
  if(is.null(dim(xmat))) ncov = 1

  y0 <- seq(min(y),max(y),length=100)
  ny0 = length(y0)

  nout <- (draws - burn)/thin

  # Empty vectors that permit storing MCMC iterates
  beta <- matrix(0, nrow=nout, ncol=ncov)
  alpha <- cVA <- mVA <- qVA <- matrix(0, nrow=nout, ncol=nschool)
  v <- Q<-matrix(0, nrow=nout, ncol=length(school))
  a<-sig2a<-lambda <- rep(0, length=nout)

  C.out <- .C("mcmcloop_qva",
              as.integer(draws), as.integer(burn),
              as.integer(thin),as.integer(nobsi),
              as.integer(nschool), as.integer(ncov),
              as.integer(school),as.integer(ny0),
              as.double(tau),
              as.double(y), as.double(t(xmat)),
              as.double(y0), as.double(priorVal),
              as.integer(verbose),
              beta.out=as.double(beta),alpha.out=as.double(alpha),
              v.out=as.double(v), a.out=as.double(a),
              sig2a.out=as.double(sig2a), lambda.out=as.double(lambda),
              cVA.out=as.double(cVA),mVA.out=as.double(mVA),
              qVA.out=as.double(qVA),Q.out=as.double(Q)
  )

  out$beta <- matrix(C.out$beta.out, nrow=nout, byrow=TRUE)
  out$alpha <- matrix(C.out$alpha.out, nrow=nout, byrow=TRUE)
  out$v <- matrix(C.out$v.out, nrow=nout, byrow=TRUE)
  out$a <- matrix(C.out$a.out, nrow=nout, byrow=TRUE)
  out$sig2a <- matrix(C.out$sig2a.out, nrow=nout, byrow=TRUE)
  out$lambda <- matrix(C.out$lambda.out, nrow=nout, byrow=TRUE)
  out$cVA <- matrix(C.out$cVA.out, nrow=nout, byrow=TRUE)
  out$mVA <- matrix(C.out$mVA.out, nrow=nout, byrow=TRUE)
  out$qVA <- matrix(C.out$qVA.out, nrow=nout, byrow=TRUE)
  out$Q <- matrix(C.out$Q.out, nrow=nout, byrow=TRUE)
  out

}


# Wrapper function that access C-code that fits the time-dependent value-added model

tdVA <- function(y1,xmat1,y2,xmat2,school1,school2,groupID=NULL,model=0,
				 priors=c(0,100^2,1, 1, 1, 1, 0,100^2,-1, 1, 0,100^2,0,100^2),
				 var.global=TRUE,MHsd=c(0.2),
				 nchains=1,draws=50000, burn=40000, thin=10,verbose=FALSE){
	# priors - vector of size 12 that contains hyper-prior values for following parameters
	#	beta - slot 1 - mean,  slot 2 - variance
	#	tau2 - slot 3 - shape, slot 4 - rate
	#	sig2 - slot 5 - shape, slot 6 - rate
	#	gamma - slot 7 - mean, slot 8 - variance
	#	phi12 - slot 9 - lower bound, slot 10 - upper bound
	# phi02 - slot 11 - mean, slot 12 - variance
	# phi01 - slot 13 - mean, slot 14 - variance

	# model - 0 independent (phi12 = gamma = 0)
	#         1 AR(1) (gamma = 0)
	#         2 shock (phi12 = 0)

	# MHvar are the MH tunning parameters.  The
	# one entry is for phi12 csigPHI1,


  out <- NULL
	uc1 <- unique(school1)
	uc2 <- unique(school2)
	diff1 <- diff(school1)
	diff2 <- diff(school2)

	if(sum(as.numeric(diff1) > 1) > 0) stop("School labels for cohort 1 are not contiguous")
	if(sum(as.numeric(diff2) > 1) > 0) stop("School labels for cohort 2 are not contiguous")

	nschool = length(uc1)

	nobsi1 = tapply(school1,school1,length)
	nobsi2 = tapply(school2,school2,length)

	ncov = dim(cbind(xmat1))[2]
	ngroup = length(unique(groupID))

	nout <- (draws - burn)/thin

	mcmc_out <- list()

	for(cc in 1:nchains){

		if(verbose) message("chain = ", cc)

		if(is.null(groupID)){
			if(var.global){
				beta1 <- beta2 <- matrix(0, nrow=nout, ncol=ncov)
				alpha1 <- alpha2 <- matrix(0, nrow=nout, ncol=nschool)
				sig21 <- sig22 <- tau2 <- tau21 <- tau22 <- phi02 <- phi12 <- gamma <- phi01 <- rep(0,nout)
				lpml <- waic <- rep(0,1)
				if(verbose) message("Fitting model with constant variance for each school and no grouping")
				C.out <- .C("mcmcloopM2a1",
    	          			as.integer(draws), as.integer(burn), as.integer(thin),
        	      			as.integer(nobsi1), as.integer(nobsi2),
            	  			as.integer(nschool), as.integer(ncov),
              				as.integer(school1),as.integer(school2),
              				as.double(y1), as.double(t(xmat1)),
              				as.double(y2), as.double(t(xmat2)),
              				as.double(priors), as.integer(model),
              				as.double(MHsd), as.integer(verbose),
              				beta1.out=as.double(beta1),beta2.out=as.double(beta2),
              				alpha1.out=as.double(alpha1), alpha2.out=as.double(alpha2),
              				sig21.out=as.double(sig21), sig22.out=as.double(sig22),
              				tau21.out=as.double(tau21), tau22.out=as.double(tau22),
              				phi02.out=as.double(phi02),phi12.out=as.double(phi12),
              				gamma.out=as.double(gamma),phi01.out=as.double(phi01),
              				lpml.out=as.double(lpml), waic.out=as.double(waic))
			} else {

				beta1 <- beta2 <- matrix(0, nrow=nout, ncol=ncov)
				alpha1 <- alpha2 <- sig21 <- sig22 <- matrix(0, nrow=nout, ncol=nschool)
				tau21 <- tau22 <- phi02 <- phi12 <- gamma <- phi01 <- rep(0,nout)
				lpml <- waic <- rep(0,1)
				if(verbose) message("Fitting model with school-specific variance and no grouping")
				C.out <- .C("mcmcloopM2a2",
    	          			as.integer(draws), as.integer(burn), as.integer(thin),
        	      			as.integer(nobsi1), as.integer(nobsi2),
            	  			as.integer(nschool), as.integer(ncov),
              				as.integer(school1),as.integer(school2),
              				as.double(y1), as.double(t(xmat1)),
              				as.double(y2), as.double(t(xmat2)),
              				as.double(priors), as.integer(model),
               				as.double(MHsd),as.integer(verbose),
              				beta1.out=as.double(beta1),beta2.out=as.double(beta2),
              				alpha1.out=as.double(alpha1), alpha2.out=as.double(alpha2),
              				sig21.out=as.double(sig21), sig22.out=as.double(sig22),
              				tau21.out=as.double(tau21), tau22.out=as.double(tau22),
              				phi02.out=as.double(phi02),phi12.out=as.double(phi12),
              				gamma.out=as.double(gamma),phi01.out=as.double(phi01),
             				  lpml.out=as.double(lpml), waic.out=as.double(waic))


			}

		}

		if(!is.null(groupID)){
			if(var.global){
				beta1 <- beta2 <- matrix(0, nrow=nout, ncol=ncov)
				alpha1 <- alpha2 <- matrix(0, nrow=nout, ncol=nschool)
				phi02 <- phi12 <- gamma <-  matrix(0, nrow=nout, ncol=ngroup)
				sig21 <- sig22 <- tau21 <- tau22  <-  phi01 <- rep(0,nout)
				lpml <- waic <- rep(0,1)
				if(verbose) message("Fitting model with constant variance for each school and with grouping")
				C.out <- .C("mcmcloopM2b1",
    	          		as.integer(draws), as.integer(burn), as.integer(thin),
        	      		as.integer(nobsi1), as.integer(nobsi2),
            	  		as.integer(nschool), as.integer(ncov),
              			as.integer(school1),as.integer(school2),
              			as.double(y1), as.double(t(xmat1)),
              			as.double(y2), as.double(t(xmat2)),
              			as.integer(ngroup), as.integer(groupID),
              			as.double(priors), as.integer(model),
               			as.double(MHsd),as.integer(verbose),
             			  beta1.out=as.double(beta1),beta2.out=as.double(beta2),
              			alpha1.out=as.double(alpha1), alpha2.out=as.double(alpha2),
              			sig21.out=as.double(sig21), sig22.out=as.double(sig22),
              			tau21.out=as.double(tau21), tau22.out=as.double(tau22),
              			phi02.out=as.double(phi02),phi12.out=as.double(phi12),
              			gamma.out=as.double(gamma),phi01.out=as.double(phi01),
              			lpml.out=as.double(lpml), waic.out=as.double(waic))
			} else {

				beta1 <- beta2 <- matrix(0, nrow=nout, ncol=ncov)
				alpha1 <- alpha2 <- sig21 <- sig22 <- preda1 <- preda2 <-  matrix(0, nrow=nout, ncol=nschool)
				phi02 <- phi12 <- gamma <- matrix(0, nrow=nout, ncol=ngroup)
				tau21  <- tau22 <-  phi01 <- rep(0,nout)
				lpml <- waic <- rep(0,1)
				if(verbose) message("Fitting model with school-specific variance and with grouping")
				C.out <- .C("mcmcloopM2b2",
    	          		as.integer(draws), as.integer(burn), as.integer(thin),
        	      		as.integer(nobsi1), as.integer(nobsi2),
            	  		as.integer(nschool), as.integer(ncov),
              			as.integer(school1),as.integer(school2),
              			as.double(y1), as.double(t(xmat1)),
              			as.double(y2), as.double(t(xmat2)),
              			as.integer(ngroup), as.integer(groupID),
              			as.double(priors), as.integer(model),
               			as.double(MHsd),as.integer(verbose),
             			  beta1.out=as.double(beta1),beta2.out=as.double(beta2),
              			alpha1.out=as.double(alpha1), alpha2.out=as.double(alpha2),
              			sig21.out=as.double(sig21), sig22.out=as.double(sig22),
              			tau21.out=as.double(tau21), tau22.out=as.double(tau22),
              			phi02.out=as.double(phi02),phi12.out=as.double(phi12),
              			gamma.out=as.double(gamma),phi01.out=as.double(phi01),
             			  lpml.out=as.double(lpml), waic.out=as.double(waic))

			}

		}

		out$beta1 <- matrix(C.out$beta1.out, nrow=nout, byrow=TRUE)
		out$beta2 <- matrix(C.out$beta2.out, nrow=nout, byrow=TRUE)
		out$alpha1 <- matrix(C.out$alpha1.out, nrow=nout, byrow=TRUE)
		out$alpha2 <- matrix(C.out$alpha2.out, nrow=nout, byrow=TRUE)
		out$sig21 <- matrix(C.out$sig21.out, nrow=nout, byrow=TRUE)
		out$sig22 <- matrix(C.out$sig22.out, nrow=nout, byrow=TRUE)
		out$phi02 <- matrix(C.out$phi02.out, nrow=nout, byrow=TRUE)
		if(model==1|model==3) out$phi12 <- matrix(C.out$phi12.out, nrow=nout, byrow=TRUE)
		out$phi01 <- matrix(C.out$phi01.out, nrow=nout, byrow=TRUE)
		out$tau21 <- matrix(C.out$tau21.out, nrow=nout, byrow=TRUE)
		out$tau22 <- matrix(C.out$tau22.out, nrow=nout, byrow=TRUE)
		if(model==2|model==3) out$gamma <- matrix(C.out$gamma.out, nrow=nout, byrow=TRUE)
		out$lpml <- C.out$lpml.out
		out$waic <- C.out$waic.out


   # HPD calculator (see TeachingDemos)
   emp.hpd <- function(x, conf=0.95){
     conf <- min(conf, 1-conf)
     n <- length(x)
     nn <- round( n*conf )
	   x <- sort(x)
	   xx <- x[ (n-nn+1):n ] - x[1:nn]
	   m <- min(xx)
	   nnn <- which(xx==m)[1]
	   return( c( x[ nnn ], x[ n-nn+nnn ] ) )
   }


    groupIDt <- groupID
		if(is.null(groupID)) {groupIDt <- rep(1, nschool)}
		# Value-Added estimates based on posterior mean of VA estimators
    #
		# Value-added for first cohort is the same for all models

		out$VA1.draws <- out$alpha1 - c(out$phi01)
		out$VA1.estimate <- apply(out$VA1.draws,2,mean)
		out$VA1.intervals <- apply(out$VA1.draws,2,emp.hpd)

		# value-added for second cohorts depend on model
		VA2.draws <- matrix(NA, nout, nschool)

		if(model == 0){

			for(j in 1: nschool){
				VA2.draws[,j] <- out$alpha2[,j] - c(out$phi02[,groupIDt[j]])
			}

			out$VA2.draws <- VA2.draws
			out$VA2.estimate <- apply(VA2.draws,2,mean)
			out$VA2.intervals <- apply(VA2.draws,2,emp.hpd)
		}
		if(model==1){

			for(j in 1: nschool){
				VA2.draws[,j] <- out$alpha2[,j] - c(out$phi02[,groupIDt[j]]) -
				                                  c(out$phi12[,groupIDt[j]])*c(out$phi01)
			}

			out$VA2.draws <- VA2.draws
			out$VA2.estimate <- apply(VA2.draws,2,mean)
			out$VA2.intervals <- apply(VA2.draws,2,emp.hpd)
		}
		if(model==2){

			mcmcSUM <- matrix(0,nout, nschool)
			for(k in 1:ncol(xmat1)){
				X1i.mn <- cbind(tapply(xmat1[,k],school1,mean))
				X2i.mn <- cbind(tapply(xmat2[,k],school2,mean))
				EX1igivenX2i <- lm(X1i.mn ~ X2i.mn)$fitted
				for(j in 1: nschool){
					mcmcSUM[,j] <- mcmcSUM[,j] + (out$gamma[,groupIDt[j]])*out$beta1[,k]*EX1igivenX2i[j]
				}


			}

			for(j in 1: nschool){
				VA2.draws[,j] <- out$alpha2[,j] -
				                 c(out$phi02[,groupIDt[j]]) -
				                 c(out$gamma[,groupIDt[j]])*c(out$phi01) -
				                 mcmcSUM[,j]
			}

			out$VA2.draws <- VA2.draws
			out$VA2.estimate <- apply(out$VA2.draws,2,mean)
			out$VA2.intervals <- apply(out$VA2.draws,2,emp.hpd)

		}


		if(model==3){

		  mcmcSUM <- matrix(0,nout, nschool)
		  for(k in 1:ncol(xmat1)){
		    X1i.mn <- cbind(tapply(xmat1[,k],school1,mean))
		    X2i.mn <- cbind(tapply(xmat2[,k],school2,mean))
		    EX1igivenX2i <- lm(X1i.mn ~ X2i.mn)$fitted
		    for(j in 1: nschool){
		      mcmcSUM[,j] <- mcmcSUM[,j] + (out$gamma[,groupIDt[j]])*out$beta1[,k]*EX1igivenX2i[j]
		    }


		  }

		  for(j in 1: nschool){
		    VA2.draws[,j] <- out$alpha2[,j] -
		                     c(out$phi02[,groupIDt[j]]) -
		                     c(out$phi01)*(c(out$phi12[,groupIDt[j]]) +
		                                   c(out$gamma[,groupIDt[j]])) -
		                     mcmcSUM[,j]
		  }

		  out$VA2.draws <- VA2.draws
		  out$VA2.estimate <- apply(out$VA2.draws,2,mean)
		  out$VA2.intervals <- apply(out$VA2.draws,2,emp.hpd)

		}


		mcmc_out[[cc]] <- out
	}

	if(nchains==1) allout <- mcmc_out[[1]]
	if(nchains>1) allout <- mcmc_out
	allout

}

