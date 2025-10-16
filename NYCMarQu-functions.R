################################################################################
# This file contains the functions for fitting the quantile curves reported in:#
# 'Quantile regression for longitudinal within-race running data: the 2022 New #
# York City Marathon' - (submitted 2025) - D'Haen, Flórez, Molenberghs,        #
# Van Keilegom, Delecluse, Verhasselt.                                         #
#                                                                              #
# Code authors: Alvaro Flórez & Myrthe D'Haen                                  #
# Last revised on: 16/10/2025                                                  #
#                                                                              #
################################################################################


#### Part I: data preparation --------------------------------------------------
IntervalSpeed = function(data){ 
  c(data$dis[1]/(data$time[1]/60), diff(data$dis)/(diff(data$time)/60))
}

CumulativeSpeed = function(data){ 
  data$dis/(data$time/60)
}

Acceleration <- function(data){
  c(NA, diff(data$intSpeed)/(diff(data$time)/60))
}

SpeedDecline = function(data){
  c(NA, tail(data$intSpeed, -1)/head(data$cumSpeed, -1))
}

SplittingVector = function(data){
  data$cumSpeed/tail(data$cumSpeed,1)
}
# ------------------------------------------------------------------------------


#### Part II: Model fitting ----------------------------------------------------

# Note: the code for the PWE estimator was mostly written by Alvaro Flórez in
# the context of the paper
# Verhasselt, Flórez, Molenberghs, and Van Keilegom (2025) - Copula-based
# pairwise estimator for quantile regression with hierarchical missing data.
# Statistical Modelling, 25(2):129-149.

### Function for the pairwise estimator (Gaussian copula) ...
pwe.ALDcop.parallel = function(ID,X,y,tau=0.5,EPS=0.05,parm.ini=NULL,
                               weighted=F,weights=NULL,Rstr='AR',mtimes=NULL,
                               miss.data.method='AC',var.estimate=T,ncores=8){
  
  # ID: id of each individual
  # X: covariates matrix
  # y: response vector
  # tau: quantile level
  # EPS: smoothing parameter
  # parm.ini: initial parameters
  # weighted: IPW estimator
  # weights: matrix of weights for each pair of observation
  # Rstr: correlation structure (UN: unstructured, AR: autorregresive, CS: compound-symmetry) 
  # mtimes: measurement times
  # miss.data.method: CC: only considering complete cases, CS: only considering complete pair, AC: considering also incomplete pair (only apply for PWE)
  # var.estimate: compute the variance matrix for beta estimated if TRUE, does not compute the variance matrix otherwise
  # ncores: number of cores used to compute the score vector and Hessian matrix
  
  pALaD = function (q, mu, phi, alpha){
      ifelse(q > mu, (1 - (1 - alpha) * exp(-alpha * (q - mu)/phi)), 
             (alpha * exp((1 - alpha) * (q - mu)/phi)))
    }
  is.PD = function(M){
      if(anyNA(M)){return(FALSE)}
      if (!is.matrix(M)) 
        stop("x is not a matrix.")
      eigs <- eigen(M,symmetric = TRUE,only.values=T)$values
      if (any(is.complex(eigs))) 
        return(FALSE)
      if (min(eigs) > 1e-8) 
        pd <- TRUE
      else pd <- FALSE
      return(pd)
    }
  Theta2Corr = function(d,Thetaval){
    # p: dimension the correlation matrix
    # Thetaval: vector of the Theta matix (Vec Theta)
    ThetaMat = matrix(0,d,d)
    ThetaMat[lower.tri(ThetaMat)] =Thetaval
    cos.Theta = cos(ThetaMat)
    sin.Theta = sin(ThetaMat)
    U = matrix(0,d,d)
    for(j in 2:d){
      U[j,1:j]=cumprod(c(1,sin.Theta[j,1:(j-1)]))*cos.Theta[j,1:j]
    }
    U[1,1] = 1
    R = tcrossprod(U)
    return(R)
  }
  Corr2Theta = function(R){
    # p: dimension the correlation matrix
    # Thetaval: vector of the Theta matix (Vec Theta)
    B = chol(R)
    B = t(B)
    d = ncol(R)
    U = matrix(0,d,d)
    
    U[,1] = acos(B[,1])
    
    if(d==2){U[2,1] = acos(B[2,1])}else{
      for(i in 3:(d)){
        for(j in 2:(i-1)){
          #U[i,j]=acos(B[i,j]*(1/prod(sin(acos(B[i,1:(j-1)])))))
          U[i,j]=acos(B[i,j]*(1/prod(sin(U[i,1:(j-1)]))))
        }
      }
      
    }
    Theta = U[lower.tri(U)]
    return(Theta)
  }
  ARstr = function(mtimes,rho){
    delta <- mtimes
    rho^abs(outer(delta,delta,"-"))
  }
  
  ALC.pair.loglike = function(parm,Y,Xmat,Alp,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL,ncores=NULL){
    if(is.matrix(Y)){
      N = ncol(Y)    
      d = nrow(Y)
    }else{
      N=1
      d=length(Y)
    }
    
    p = ncol(Xmat)
    
    Beta=parm[1:p]
    phi=exp(parm[(p+1):(p+d)])
    crho = parm[-c(1:(p+d))]
    if(any(phi <0)){
      return(NA)
    }
    
    if(length(crho)==1){
      if(Rstr=='AR'){
        rho = exp(crho)/(exp(crho)+1)
        Rho=ARstr(mtimes,rho)
      }else{
        if(Rstr=='CS'){
          rho = (exp(crho) - 1/(d-1))/(1 + exp(crho))
          Rho = matrix(rho,d,d)
          diag(Rho)  = 1          
        }
      }
    }else{
      Thetaval = pi*exp(rho)/(1+exp(rho))
      Rho = Theta2Corr(d,Thetaval)
    }
    if(!is.PD(Rho)){return(NA)}
    
    pairs = combn(d,2)
    if(N > 1){    
      Mu = matrix(Xmat%*%Beta,d,N)
      
      Yc = Y - Mu
      Yc.s = Yc/phi
      log.f = mapply(function(x){
        u = Yc.s[x,]
        ll = log(Alp*(1-Alp)/phi[x]) + ifelse(u>0,-(Alp *u) ,(1 - Alp) * u) 
        if(EPS > 0){
          ll = ll  + EPS/2*log(EPS + abs(u))
        }
        return(ll)
      },x=1:d)
      
      U = mapply(function(x){
        u = pALaD(q=Y[x,],mu=Mu[x,],phi=phi[x],alpha=Alp)
        u[u>0.9999]=0.9999
        u
      },x=1:d)
      
      Uq = qnorm(U)
      ll.pairs = mapply(function(x){
        rho.x = c(Rho[pairs[1,x],pairs[2,x]])
        U.x = Uq[,pairs[,x]]
        log.f.x = log.f[,pairs[,x]]
        log.fC.x = (2*rho.x*rowProds(U.x,method='expSumLog') - rowSums2(U.x^2)*rho.x^2)/(2*(1-rho.x^2)) - 0.5*log(1-rho.x^2)
        if(anyNA(log.fC.x)){log.fC.x[is.na(log.fC.x)] = 0}
        ll.x = rowSums2(log.f.x,na.rm=T) + log.fC.x
        if(miss.data.method=='CS'){
          CS.x = apply(log.f.x,1,anyNA)
          ll.x[CS.x] = 0
        }
        if(!is.null(weights)){
          weights.x = weights[,x]
          ll.x = ll.x*weights.x
        }
        return(sum(ll.x))
      },x=1:ncol(pairs))
      
      return(sum(ll.pairs,na.rm = T))
    }else{
      Mu = c(Xmat%*%Beta)
      Yc = Y - Mu
      Yc.s = Yc/phi
      log.f = log(Alp*(1-Alp)/phi) + ifelse(Yc.s>0,-(Alp *Yc.s) ,(1 - Alp) * Yc.s)
      if(EPS > 0){
        log.f = log.f  + EPS/2*log(EPS + abs(Yc.s))
      }
      
      U = pALaD(q=Y,mu=Mu,phi=phi,alpha=Alp)
      U[U>0.9999]=0.9999
      
      
      Uq = qnorm(U)
      ll.pairs = mapply(function(x){
        rho.x = c(Rho[pairs[1,x],pairs[2,x]])
        U.x = Uq[pairs[,x]]
        log.f.x = sum(log.f[pairs[,x]],na.rm = T)
        log.fC.x = (2*rho.x*prod(U.x) - sum(U.x^2)*rho.x^2)/(2*(1-rho.x^2)) - 0.5*log(1-rho.x^2)
        if(anyNA(log.fC.x)){log.fC.x[is.na(log.fC.x)] = 0}
        
        ll.x = log.f.x + log.fC.x
        
        if(miss.data.method=='CS'){
          CS.x = anyNA(log.f.x)
          ll.x = 0
        }
        
        if(!is.null(weights)){
          weights.x = weights[x]
          ll.x = ll.x*weights.x
        }
        return(ll.x)
      },x=1:ncol(pairs))
      
      return(sum(ll.pairs,na.rm = T))
      
    }
  }
  grad.ALC.pair.loglike = function(parm,Y,Xmat,Alp,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL){
    
    k=1e-6
    if(is.matrix(Y)){
      d = nrow(Y)    
    }else{
      d = length(Y)
    }
    Grad = mapply(function(x){
      parm.eval.1 = parm.eval.2 = parm
      parm.eval.1[x] = parm[x] + k/2
      parm.eval.2[x] = parm[x] - k/2
      ll.1 = ALC.pair.loglike(parm.eval.1,Y,Xmat,Alp,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      ll.2 = ALC.pair.loglike(parm.eval.2,Y,Xmat,Alp,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      deriv = (ll.1-ll.2)/(k)
      return(deriv)
    },x=1:length(parm))
    return(Grad)
  }
  
  
  grad.ALC.pair.loglike.par = function(parm,Y,Xmat,Alp,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL,ncores=ncores){
    
    k=1e-6
    if(is.matrix(Y)){
      d = nrow(Y)    
    }else{
      d = length(Y)
    }
    varlist =c('ALC.pair.loglike','Y','Xmat','Alp','EPS','weights','miss.data.method','Rstr','mtimes','k','parm')
    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl=cl, varlist=varlist, envir=environment())
    clusterEvalQ(cl=cl,library('matrixStats'))    
    
    
    Grad = pblapply(X=1:length(parm),function(X){
      parm.eval.1 = parm.eval.2 = parm
      parm.eval.1[X] = parm[X] + k/2
      parm.eval.2[X] = parm[X] - k/2
      ll.1 = ALC.pair.loglike(parm.eval.1,Y,Xmat,Alp,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      ll.2 = ALC.pair.loglike(parm.eval.2,Y,Xmat,Alp,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      deriv = (ll.1-ll.2)/(k)
      return(deriv)
    },cl=cl)
    stopCluster(cl)
    return(unlist(Grad))
  }
  hessian.ALC.pair.loglike = function(parm,Y,Xmat,Alp,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL,ncores=ncores){
    k=1e-6
    evals =   combinations(length(parm),2,repeats.allowed = T)
    d = nrow(Y)
    ll = ALC.pair.loglike(parm,Y,X,Alp,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)    
    
    varlist =c('ALC.pair.loglike','Y','Xmat','Alp','EPS','weights','miss.data.method','Rstr','mtimes','k','evals','evals',
               'parm','ll')
    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl=cl, varlist=varlist, envir=environment())
    clusterEvalQ(cl=cl,library('matrixStats'))    
    
    Hessian.vals = pblapply(X=1:nrow(evals),function(X){
      eval = evals[X,]
      if(length(unique(eval))==1){
        parm.pos = parm.neg = parm
        parm.pos[unique(eval)] = parm[unique(eval)] + k
        parm.neg[unique(eval)] = parm[unique(eval)] - k
        #     if(d > 2){
        ll.pos = ALC.pair.loglike(parm.pos,Y,Xmat,Alp,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.neg = ALC.pair.loglike(parm.neg,Y,Xmat,Alp,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        #     }else{
        #     }
        fit = (ll.pos  + ll.neg - 2*ll)/k^2
      }else{
        parm.pos.pos = parm.pos.neg = parm.neg.neg = parm.neg.pos = parm
        parm.pos.pos[eval] = parm[eval] + k
        parm.pos.neg[eval] = parm[eval] + c(k,-k)
        parm.neg.neg[eval] = parm[eval] -k
        parm.neg.pos[eval] = parm[eval] + c(-k,k)
        #    if(d > 2){
        ll.pos.pos = ALC.pair.loglike(parm.pos.pos,Y,Xmat,Alp,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.pos.neg = ALC.pair.loglike(parm.pos.neg,Y,Xmat,Alp,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.neg.neg = ALC.pair.loglike(parm.neg.neg,Y,Xmat,Alp,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.neg.pos = ALC.pair.loglike(parm.neg.pos,Y,Xmat,Alp,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        fit = (ll.pos.pos - ll.pos.neg - ll.neg.pos + ll.neg.neg)/(4*k^2)
      }
      return(fit)
    },cl=cl)
    stopCluster(cl)
    p = length(parm)
    Hessian = invvech(unlist(Hessian.vals))
    return(Hessian)
  }
  
  d= unique(table(ID))
  N = length(y)/d
  p = ncol(X)
  Y = matrix(y,d,N,byrow = F)
  alpha = tau
  
  
  # Initial values
  if(is.null(parm.ini)){
    mod.ini = lm(y~X-1,na.action="na.exclude")
    res = y  - predict(mod.ini,as.data.frame(X[,-1]))
    phi.ini = apply(matrix(res,d,N),1,sd,na.rm=T) * sqrt((alpha^2*(1-alpha)^2)/(1-2*alpha+2*alpha^2))
    est.ini = rq(y~X-1,tau=tau)
    Beta.ini = est.ini$coef
    res.ini = matrix(y - predict.rq(est.ini,as.data.frame(X[,-1])),d,N)
    pseudo.cor =cor(t(res.ini),use='complete.obs')
    if(Rstr=='AR'){
      rho.ini = (mean(diag(pseudo.cor[-1,-d])))^(1/(mtimes[2]-mtimes[1]))
      crho.ini = log(rho.ini/(1 - rho.ini))
      parm.ini = c(Beta.ini,log(phi.ini),crho.ini)      
      
    }else{
      rho.ini = Corr2Theta(pseudo.cor)
      crho.ini = log(rho.ini/(pi - rho.ini))
      parm.ini = c(Beta.ini,log(phi.ini),crho.ini)      
    }
    
  }
  n.par=length(parm.ini)
  
  
  if(!weighted){
    weights=NULL
  }
  
  

  fit = maxNR(ALC.pair.loglike, grad = grad.ALC.pair.loglike.par, hess = hessian.ALC.pair.loglike, start=parm.ini,
              constraints = NULL, finalHessian = T, bhhhHessian=FALSE,
              fixed = NULL, activePar = NULL, control = list(printLevel=2),
              Y=Y,Xmat=X,Alp=alpha,weights=weights, EPS=EPS,miss.data.method=miss.data.method,Rstr=Rstr,
              mtimes=mtimes,ncores=ncores)

  parmest = fit$estimate
  ll = fit$maximum
  Grad = fit$gradient
  if(var.estimate){
    A.inv = fit$hessian
    
    A=tryCatch(solve(-A.inv),error=function(e){NULL})
    if(!is.null(A)){
      
      B = pbmapply(function(x){
        X.i = X[ID==x,]
        y.i = Y[,x]
        b.i =  grad.ALC.pair.loglike(parmest, Y=y.i,Xmat=X.i,Alp=alpha,weights=weights[x,], 
                                     EPS=EPS,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        b.i = matrix(b.i,nrow = 1)
        B.i = crossprod(b.i,b.i)
        return(vech(B.i))
      },x=1:N)
      B = invvech(apply(B,1,sum))
      V = A%*%B%*%A
    }else{
      V = matrix(NA,nrow(A.inv),ncol(A.inv))
      A.inv = matrix(NA,length(parmest),length(parmest))
    }
  }else{
    V = matrix(NA,length(parmest),length(parmest))
    A.inv = matrix(NA,length(parmest),length(parmest))
  }
  
  Est.tr = parmest
  
  Beta = Est.tr[1:p]
  Phi = exp(Est.tr[(p+1):(p+d)])
  if(Rstr=='AR'){
    rho.ar = exp(Est.tr[p+d+1])/(exp(Est.tr[p+d+1])+1)
    Rho = list(rho = rho.ar,Rho=ARstr(mtimes,rho.ar))
  }else{
    rho.tr = Est.tr[-c(1:(p+d))]
    rho = pi*exp(rho.tr)/(1+exp(rho.tr))
    Rho = Theta2Corr(d,rho)
  }
  
  output = list(Beta=Beta,Phi=Phi,Rho=Rho,Grad=Grad,ll=ll,V=V,Hessian=A.inv,parm = fit$estimate,code=fit$code,
                EPS=EPS,tau=tau)
  # Beta: estimated quantile regression coefficients
  # Phi : estimated scale parameters for the ALD
  # Rho: estimated correlation matrix for the copula
  # Grad: score function evaluated at solution
  # ll: calculated log-likelihood (or log-pseudolikelihood) at solution
  # V: estimated variance matrix for estimated parameters
  
  return(output)
}

# ... and its variant with a t-copula
pwe.ALDcop.t.parallel = function(ID,X,y,tau=0.5,v=4,EPS=0.05,parm.ini=NULL,
                                 weighted=F,weights=NULL,Rstr='AR',mtimes=NULL,
                                 miss.data.method='AC',var.estimate=T,ncores=8){
  
  # ID: id of each individual
  # X: covariates matrix
  # y: response vector
  # tau: quantile level
  # v: degrees of freedom for the Student-t
  # EPS: smoothing parameter
  # parm.ini: initial parameters
  # weighted: IPW estimator
  # weights: matrix of weights for each pair of observation
  # Rstr: correlation structure (UN: unstructured, AR: autorregresive, CS: compound-symmetry) 
  # mtimes: measurement times
  # miss.data.method: CC: only considering complete cases, CS: only considering complete pair, AC: considering also incomplete pair (only apply for PWE)
  # var.estimate: compute the variance matrix for beta estimated if TRUE, does not compute the variance matrix otherwise
  # ncores: number of cores used to compute the score vector and Hessian matrix
  
  pALaD = function (q, mu, phi, alpha){
    ifelse(q > mu, (1 - (1 - alpha) * exp(-alpha * (q - mu)/phi)), 
           (alpha * exp((1 - alpha) * (q - mu)/phi)))
  }
  is.PD = function(M){
    if(anyNA(M)){return(FALSE)}
    if (!is.matrix(M)) 
      stop("x is not a matrix.")
    eigs <- eigen(M,symmetric = TRUE,only.values=T)$values
    if (any(is.complex(eigs))) 
      return(FALSE)
    if (min(eigs) > 1e-8) 
      pd <- TRUE
    else pd <- FALSE
    return(pd)
  }
  Theta2Corr = function(d,Thetaval){
    # p: dimension the correlation matrix
    # Thetaval: vector of the Theta matix (Vec Theta)
    ThetaMat = matrix(0,d,d)
    ThetaMat[lower.tri(ThetaMat)] =Thetaval
    cos.Theta = cos(ThetaMat)
    sin.Theta = sin(ThetaMat)
    U = matrix(0,d,d)
    for(j in 2:d){
      U[j,1:j]=cumprod(c(1,sin.Theta[j,1:(j-1)]))*cos.Theta[j,1:j]
    }
    U[1,1] = 1
    R = tcrossprod(U)
    return(R)
  }
  Corr2Theta = function(R){
    # p: dimension the correlation matrix
    # Thetaval: vector of the Theta matix (Vec Theta)
    B = chol(R)
    B = t(B)
    d = ncol(R)
    U = matrix(0,d,d)
    
    U[,1] = acos(B[,1])
    
    if(d==2){U[2,1] = acos(B[2,1])}else{
      for(i in 3:(d)){
        for(j in 2:(i-1)){
          #U[i,j]=acos(B[i,j]*(1/prod(sin(acos(B[i,1:(j-1)])))))
          U[i,j]=acos(B[i,j]*(1/prod(sin(U[i,1:(j-1)]))))
        }
      }
      
    }
    Theta = U[lower.tri(U)]
    return(Theta)
  }
  ARstr = function(mtimes,rho){
    delta <- mtimes
    rho^abs(outer(delta,delta,"-"))
  }
  
  ALC.pair.loglike = function(parm,Y,Xmat,Alp,v,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL,ncores=NULL){
    if(is.matrix(Y)){
      N = ncol(Y)    
      d = nrow(Y)
    }else{
      N=1
      d=length(Y)
    }
    
    p = ncol(Xmat)
    
    Beta=parm[1:p]
    phi=exp(parm[(p+1):(p+d)])
    crho = parm[-c(1:(p+d))]
    if(any(phi <0)){
      return(NA)
    }
    
    if(length(crho)==1){
      if(Rstr=='AR'){
        rho = exp(crho)/(exp(crho)+1)
        Rho=ARstr(mtimes,rho)
      }else{
        if(Rstr=='CS'){
          rho = (exp(crho) - 1/(d-1))/(1 + exp(crho))
          Rho = matrix(rho,d,d)
          diag(Rho)  = 1          
        }
      }
    }else{
      Thetaval = pi*exp(rho)/(1+exp(rho))
      Rho = Theta2Corr(d,Thetaval)
    }
    if(!is.PD(Rho)){return(NA)}
    
    pairs = combn(d,2)
    if(N > 1){    
      Mu = matrix(Xmat%*%Beta,d,N)
      
      Yc = Y - Mu
      Yc.s = Yc/phi
      log.f = mapply(function(x){
        u = Yc.s[x,]
        ll = log(Alp*(1-Alp)/phi[x]) + ifelse(u>0,-(Alp *u) ,(1 - Alp) * u) 
        if(EPS > 0){
          ll = ll  + EPS/2*log(EPS + abs(u))
        }
        return(ll)
      },x=1:d)
      
      U = mapply(function(x){
        u = pALaD(q=Y[x,],mu=Mu[x,],phi=phi[x],alpha=Alp)
        u[u>0.9999]=0.9999
        u
      },x=1:d)
      
      
      Uq = qt(U,v)
      logKv = lgamma(v/2) - 2*lgamma((v+1)/2) + lgamma(v/2+1) 
      
      ll.pairs = mapply(function(x){
        rho.x = c(Rho[pairs[1,x],pairs[2,x]])
        U.x = Uq[,pairs[,x]]
        log.f.x = log.f[,pairs[,x]]
        log.fC.x = logKv - 0.5*log(1-rho.x^2) -0.5*(v+2)*log(1+(rowSums2(U.x^2) - 2*rho.x*rowProds(U.x))/(v*(1-rho.x^2))) + 
          0.5*(v+1)* rowSums2(log(1+U.x^2/v))
        if(anyNA(log.fC.x)){log.fC.x[is.na(log.fC.x)] = 0}
        ll.x = rowSums2(log.f.x,na.rm=T) + log.fC.x
        if(miss.data.method=='CS'){
          CS.x = apply(log.f.x,1,anyNA)
          ll.x[CS.x] = 0
        }
        if(!is.null(weights)){
          weights.x = weights[,x]
          ll.x = ll.x*weights.x
        }
        return(sum(ll.x))
      },x=1:ncol(pairs))
      
      return(sum(ll.pairs,na.rm = T))
    }else{
      Mu = c(Xmat%*%Beta)
      Yc = Y - Mu
      Yc.s = Yc/phi
      log.f = log(Alp*(1-Alp)/phi) + ifelse(Yc.s>0,-(Alp *Yc.s) ,(1 - Alp) * Yc.s)
      if(EPS > 0){
        log.f = log.f  + EPS/2*log(EPS + abs(Yc.s))
      }
      
      U = pALaD(q=Y,mu=Mu,phi=phi,alpha=Alp)
      U[U>0.9999]=0.9999
      
      
      Uq = qt(U,v)
      logKv = lgamma(v/2) - 2*lgamma((v+1)/2) + lgamma(v/2+1) # log K(v)
      
      ll.pairs = mapply(function(x){
        rho.x = c(Rho[pairs[1,x],pairs[2,x]])
        U.x = Uq[pairs[,x]]
        log.f.x = sum(log.f[pairs[,x]],na.rm = T)
        log.fC.x = logKv - 0.5*log(1-rho.x^2) -0.5*(v+2)*log(1+(sum(U.x^2) - 2*rho.x*prod(U.x))/(v*(1-rho.x^2))) + 
          0.5*(v+1)* sum(log(1+U.x^2/v))
        if(anyNA(log.fC.x)){log.fC.x[is.na(log.fC.x)] = 0}
        
        ll.x = log.f.x + log.fC.x
        
        if(miss.data.method=='CS'){
          CS.x = anyNA(log.f.x)
          ll.x = 0
        }
        
        if(!is.null(weights)){
          weights.x = weights[x]
          ll.x = ll.x*weights.x
        }
        return(ll.x)
      },x=1:ncol(pairs))
      
      return(sum(ll.pairs,na.rm = T))
      
    }
  }
  grad.ALC.pair.loglike = function(parm,Y,Xmat,Alp,v,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL){
    
    k=1e-6
    if(is.matrix(Y)){
      d = nrow(Y)    
    }else{
      d = length(Y)
    }
    Grad = mapply(function(x){
      parm.eval.1 = parm.eval.2 = parm
      parm.eval.1[x] = parm[x] + k/2
      parm.eval.2[x] = parm[x] - k/2
      ll.1 = ALC.pair.loglike(parm.eval.1,Y,Xmat,Alp,v,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      ll.2 = ALC.pair.loglike(parm.eval.2,Y,Xmat,Alp,v,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      deriv = (ll.1-ll.2)/(k)
      return(deriv)
    },x=1:length(parm))
    return(Grad)
  }
  
  
  grad.ALC.pair.loglike.par = function(parm,Y,Xmat,Alp,v,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL,ncores=ncores){
    
    k=1e-6
    if(is.matrix(Y)){
      d = nrow(Y)    
    }else{
      d = length(Y)
    }
    varlist =c('ALC.pair.loglike','Y','Xmat','Alp','v','EPS','weights','miss.data.method','Rstr','mtimes','k','parm')
    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl=cl, varlist=varlist, envir=environment())
    clusterEvalQ(cl=cl,library('matrixStats'))    
    
    
    Grad = pblapply(X=1:length(parm),function(X){
      parm.eval.1 = parm.eval.2 = parm
      parm.eval.1[X] = parm[X] + k/2
      parm.eval.2[X] = parm[X] - k/2
      ll.1 = ALC.pair.loglike(parm.eval.1,Y,Xmat,Alp,v,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      ll.2 = ALC.pair.loglike(parm.eval.2,Y,Xmat,Alp,v,EPS=EPS,weights=weights,miss.data.method='AC',Rstr=Rstr,mtimes=mtimes)
      deriv = (ll.1-ll.2)/(k)
      return(deriv)
    },cl=cl)
    stopCluster(cl)
    return(unlist(Grad))
  }
  hessian.ALC.pair.loglike = function(parm,Y,Xmat,Alp,v,EPS,weights=NULL,miss.data.method='AC',Rstr = 'UN',mtimes = NULL,ncores=ncores){
    k=1e-6
    evals =   combinations(length(parm),2,repeats.allowed = T)
    d = nrow(Y)
    ll = ALC.pair.loglike(parm,Y,X,Alp,v,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)    
    
    varlist =c('ALC.pair.loglike','Y','Xmat','Alp', 'v','EPS','weights','miss.data.method','Rstr','mtimes','k','evals','evals',
               'parm','ll')
    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl=cl, varlist=varlist, envir=environment())
    clusterEvalQ(cl=cl,library('matrixStats'))    
    
    Hessian.vals = pblapply(X=1:nrow(evals),function(X){
      eval = evals[X,]
      if(length(unique(eval))==1){
        parm.pos = parm.neg = parm
        parm.pos[unique(eval)] = parm[unique(eval)] + k
        parm.neg[unique(eval)] = parm[unique(eval)] - k
        #     if(d > 2){
        ll.pos = ALC.pair.loglike(parm.pos,Y,Xmat,Alp,v,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.neg = ALC.pair.loglike(parm.neg,Y,Xmat,Alp,v,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        #     }else{
        #     }
        fit = (ll.pos  + ll.neg - 2*ll)/k^2
      }else{
        parm.pos.pos = parm.pos.neg = parm.neg.neg = parm.neg.pos = parm
        parm.pos.pos[eval] = parm[eval] + k
        parm.pos.neg[eval] = parm[eval] + c(k,-k)
        parm.neg.neg[eval] = parm[eval] -k
        parm.neg.pos[eval] = parm[eval] + c(-k,k)
        #    if(d > 2){
        ll.pos.pos = ALC.pair.loglike(parm.pos.pos,Y,Xmat,Alp,v,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.pos.neg = ALC.pair.loglike(parm.pos.neg,Y,Xmat,Alp,v,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.neg.neg = ALC.pair.loglike(parm.neg.neg,Y,Xmat,Alp,v,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        ll.neg.pos = ALC.pair.loglike(parm.neg.pos,Y,Xmat,Alp,v,EPS=EPS,weights = weights,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        fit = (ll.pos.pos - ll.pos.neg - ll.neg.pos + ll.neg.neg)/(4*k^2)
      }
      return(fit)
    },cl=cl)
    stopCluster(cl)
    p = length(parm)
    Hessian = invvech(unlist(Hessian.vals))
    return(Hessian)
  }
  
  d= unique(table(ID))
  N = length(y)/d
  p = ncol(X)
  Y = matrix(y,d,N,byrow = F)
  alpha = tau
  
  
  # Initial values
  if(is.null(parm.ini)){
    mod.ini = lm(y~X-1,na.action="na.exclude")
    res = y  - predict(mod.ini,as.data.frame(X[,-1]))
    phi.ini = apply(matrix(res,d,N),1,sd,na.rm=T) * sqrt((alpha^2*(1-alpha)^2)/(1-2*alpha+2*alpha^2))
    est.ini = rq(y~X-1,tau=tau)
    Beta.ini = est.ini$coef
    res.ini = matrix(y - predict.rq(est.ini,as.data.frame(X[,-1])),d,N)
    pseudo.cor =cor(t(res.ini),use='complete.obs')
    if(Rstr=='AR'){
      rho.ini = (mean(diag(pseudo.cor[-1,-d])))^(1/(mtimes[2]-mtimes[1]))
      crho.ini = log(rho.ini/(1 - rho.ini))
      parm.ini = c(Beta.ini,log(phi.ini),crho.ini)      
      
    }else{
      rho.ini = Corr2Theta(pseudo.cor)
      crho.ini = log(rho.ini/(pi - rho.ini))
      parm.ini = c(Beta.ini,log(phi.ini),crho.ini)      
    }
    
  }
  n.par=length(parm.ini)
  
  
  if(!weighted){
    weights=NULL
  }
  
  
  
  fit = maxNR(ALC.pair.loglike, grad = grad.ALC.pair.loglike.par, hess = hessian.ALC.pair.loglike, start=parm.ini,
              constraints = NULL, finalHessian = T, bhhhHessian=FALSE,
              fixed = NULL, activePar = NULL, control = list(printLevel=2),
              Y=Y,Xmat=X,Alp=alpha,v=v,weights=weights, EPS=EPS,miss.data.method=miss.data.method,Rstr=Rstr,
              mtimes=mtimes,ncores=ncores)
  
  parmest = fit$estimate
  ll = fit$maximum
  Grad = fit$gradient
  if(var.estimate){
    A.inv = fit$hessian
    
    A=tryCatch(solve(-A.inv),error=function(e){NULL})
    if(!is.null(A)){
      
      B = pbmapply(function(x){
        X.i = X[ID==x,]
        y.i = Y[,x]
        b.i =  grad.ALC.pair.loglike(parmest, Y=y.i,Xmat=X.i,Alp=alpha,v=v,weights=weights[x,], 
                                     EPS=EPS,miss.data.method=miss.data.method,Rstr=Rstr,mtimes=mtimes)
        b.i = matrix(b.i,nrow = 1)
        B.i = crossprod(b.i,b.i)
        return(vech(B.i))
      },x=1:N)
      B = invvech(apply(B,1,sum))
      Vmat = A%*%B%*%A
    }else{
      warning('hessian is not pd')
      Vmat = matrix(NA,nrow(A.inv),ncol(A.inv))
      A.inv = matrix(NA,length(parmest),length(parmest))
    }
  }else{
    Vmat = matrix(NA,length(parmest),length(parmest))
    A.inv = matrix(NA,length(parmest),length(parmest))
  }
  
  Est.tr = parmest
  
  Beta = Est.tr[1:p]
  Phi = exp(Est.tr[(p+1):(p+d)])
  if(Rstr=='AR'){
    rho.ar = exp(Est.tr[p+d+1])/(exp(Est.tr[p+d+1])+1)
    Rho = list(rho = rho.ar,Rho=ARstr(mtimes,rho.ar))
  }else{
    rho.tr = Est.tr[-c(1:(p+d))]
    rho = pi*exp(rho.tr)/(1+exp(rho.tr))
    Rho = Theta2Corr(d,rho)
  }
  
  output = list(Beta=Beta,Phi=Phi,Rho=Rho,Grad=Grad,ll=ll,V=Vmat,Hessian=A.inv,parm = fit$estimate,code=fit$code,
                EPS=EPS,tau=tau)
  # Beta: estimated quantile regression coefficients
  # Phi : estimated scale parameters for the ALD
  # Rho: estimated correlation matrix for the copula
  # Grad: score function evaluated at solution
  # ll: calculated log-likelihood (or log-pseudolikelihood) at solution
  # V: estimated variance matrix for estimated parameters
  
  return(output)
}
# ------------------------------------------------------------------------------


#### Part III: Model evaluation and visualisation ------------------------------

### Function for obtaining tables with the estimates

summary_LQC = function(varName, qu_lvls, UM_eps = 0.1, CM_eps = 0.1,
                       UM_path, CM_path, UM_formula, CM_formula, 
                       UM_infix, CM_infix, UM_suffix, CM_suffix,
                       noUM = FALSE, noCM = FALSE){
  # use noUM and noCM in case output is wanted for one of the models only. In
  # that case, corresponding path, formula, infix and suffix can be anything.
  
  if (! noUM){
    est.tau.uns = lapply(X = seq_along(qu_lvls),function(x){
      tau.val = qu_lvls[x]*100
      name.x = paste(paste0(UM_path,varName),'pwe',tau.val, UM_infix,
                     'eps', UM_eps, UM_suffix,
                     sep='.')
      est.x = tryCatch(list.load(name.x),error=function(e){
        message(paste('no unstr results for tau=',tau.val))})
      
      if(is.null(est.x)){return(NA)}
      
      ## Column names for the covariates
      if (UM_formula == "SexAge"){
        if (varName %in% c("intSpeed", "cumSpeed")){
          Label = c('intercept','Sex(M)','K10','K15','K20','Half','K25','K30',
                    'M20','K35','K40','MAR','AgeSc',
                    'K10*Sex(M)','K15*Sex(M)','K20*Sex(M)','Half*Sex(M)',
                    'K25*Sex(M)','K30*Sex(M)','M20*Sex(M)','K35*Sex(M)',
                    'K40*Sex(M)','MAR*Sex(M)')
          
        }
        if (varName %in% c("accel", "speedDec")){
          Label = c('intercept','Sex(M)','K15','K20','Half','K25','K30',
                    'M20','K35','K40','MAR','AgeSc',
                    'K15*Sex(M)','K20*Sex(M)','Half*Sex(M)',
                    'K25*Sex(M)','K30*Sex(M)','M20*Sex(M)','K35*Sex(M)',
                    'K40*Sex(M)','MAR*Sex(M)')
        }
        if (varName == "splitVec"){
          Label = c('intercept','Sex(M)','K10','K15','K20','Half','K25','K30',
                    'M20','K35','K40','AgeSc',
                    'K10*Sex(M)','K15*Sex(M)','K20*Sex(M)','Half*Sex(M)',
                    'K25*Sex(M)','K30*Sex(M)','M20*Sex(M)','K35*Sex(M)','K40*Sex(M)')
        }
      }
      
      if (UM_formula == "SexAgeInter"){
        if (varName %in% c("intSpeed", "cumSpeed")){
          Label = c('intercept','Sex(M)','AgeSc', 'K10','K15','K20','Half','K25','K30',
                    'M20','K35','K40','MAR','AgeSc*Sex(M)',
                    'K10*Sex(M)','K15*Sex(M)','K20*Sex(M)','Half*Sex(M)',
                    'K25*Sex(M)','K30*Sex(M)','M20*Sex(M)','K35*Sex(M)',
                    'K40*Sex(M)','MAR*Sex(M)',
                    'K10*AgeSc','K15*AgeSc','K20*AgeSc','Half*AgeSc',
                    'K25*AgeSc','K30*AgeSc','M20*AgeSc','K35*AgeSc',
                    'K40*AgeSc','MAR*AgeSc')
        }
        if (varName %in% c("accel", "speedDec")){
          Label = c('intercept','Sex(M)','AgeSc', 'K15','K20','Half','K25','K30',
                    'M20','K35','K40','MAR','AgeSc*Sex(M)',
                    'K15*Sex(M)','K20*Sex(M)','Half*Sex(M)',
                    'K25*Sex(M)','K30*Sex(M)','M20*Sex(M)','K35*Sex(M)',
                    'K40*Sex(M)','MAR*Sex(M)',
                    'K15*AgeSc','K20*AgeSc','Half*AgeSc',
                    'K25*AgeSc','K30*AgeSc','M20*AgeSc','K35*AgeSc',
                    'K40*AgeSc','MAR*AgeSc')
        }
        if (varName == "splitVec"){
          Label = c('intercept','Sex(M)','AgeSc', 'K10','K15','K20','Half','K25','K30',
                    'M20','K35','K40','AgeSc*Sex(M)',
                    'K10*Sex(M)','K15*Sex(M)','K20*Sex(M)','Half*Sex(M)',
                    'K25*Sex(M)','K30*Sex(M)','M20*Sex(M)','K35*Sex(M)',
                    'K40*Sex(M)',
                    'K10*AgeSc','K15*AgeSc','K20*AgeSc','Half*AgeSc',
                    'K25*AgeSc','K30*AgeSc','M20*AgeSc','K35*AgeSc',
                    'K40*AgeSc')
        }
      }
      
      est = est.x$Beta
      std.err = sqrt(diag(est.x$V))[1:length(est)]
      z.value = est/std.err
      grad = est.x$Grad[1:length(est)]
      p.value = pnorm(abs(z.value),lower.tail = F)*2
      output = data.frame(tau=tau.val,effect=Label,est = round(est,4), SE = round(std.err,6),
                          z.value=round(z.value,2),p.value=round(p.value,6),grad=round(grad,6))
      
    })
  } else {
    est.tau.uns <- NA
    }
  
  if (! noCM){
    est.tau.cubic = lapply(X = seq_along(qu_lvls),function(x){
      tau.val = qu_lvls[x]*100
      name.x = paste(paste0(CM_path,varName),'pwe',tau.val, CM_infix,
                     'eps', CM_eps, CM_suffix,
                     sep='.')
      est.x = tryCatch(list.load(name.x),error=function(e){
        message(paste('no cubic results for tau=',tau.val))})
      
      if(is.null(est.x)){return(NA)}
      
      ## Column names for the covariates
      if (CM_formula == "SexAge"){
        
        if (varName == 'splitVec'){
          Label = c('rSlope','rd','rd2','rd3','rd*Sex(M)','rd2*Sex(M)','rd3*Sex(M)','rd*AgeSc')
          
        } else {
          Label = c('intercept','Sex(M)','d','d2','d3','d*Sex(M)','d2*Sex(M)','d3*Sex(M)','AgeSc','Slope')
        }
      }
      
      if (CM_formula == "SexAgeInter"){
        
        if (varName == 'splitVec'){
          Label = c('rSlope','rd','rd2','rd3','rd*Sex(M)','rd2*Sex(M)',
                    'rd3*Sex(M)','rd*AgeSc', 'rd2*AgeSc', 'rd3*AgeSc',
                    'rd*Sex(M)*AgeSc')
        } else {
          Label = c('intercept','Sex(M)','d','d2','d3','d*Sex(M)','d2*Sex(M)',
                    'd3*Sex(M)','AgeSc','Slope', 'Sex(M)*AgeSc', 'd*AgeSc',
                    'd2*AgeSc','d3*AgeSc')
        }
      }
      
      est = est.x$Beta
      std.err = sqrt(diag(est.x$V))[1:length(est)]
      z.value = est/std.err
      grad = est.x$Grad[1:length(est)]
      p.value = pnorm(abs(z.value),lower.tail = F)*2
      output = data.frame(tau=tau.val,effect=Label,est = round(est,4) , SE = round(std.err,6),
                          z.value=round(z.value,2),p.value=round(p.value,6),grad=round(grad,6))
    })
  } else {
    est.tau.cubic <- NA
    }
  
  return(list(outUns = est.tau.uns, outCub = est.tau.cubic))
}


### Functions for visualisation of the fitted models

# for variables other than SplitVec
plot_LQC<- function(Data, y_var_name, y_label, y_range, y_ticks,
                    female_box_location, male_box_location,
                    qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                    plot_title, storage_path, storage_suffix,
                    UM_path, CM_path, UQR_path,
                    UM_formula, CM_formula, UQR_formula,
                    UM_infix, CM_infix, UM_suffix, CM_suffix,
                    UM_eps = 0.1, CM_eps = 0.1,
                    Age = 35, SampleSize = 50, Seed = 123,
                    omit_indices = c(),
                    save_separate_panels = FALSE){
  
  # UM_eps = smoothing parameter used in UM model to be plotted
  # UM_path = location (folder) with result of relevant fitted UM model
  # UM_formula = categorical variable representing the model formula used in the
  #     fitting, with currently implemented options "SexAge" and "SexAgeInter"
  # (with or without interactions, respectively, between sex and scaled age)
  
  # idem for CM_{eps,path,formula} and UQR_{path,formula} (for UQR, no smoothing)
  
  # {UM,CM}_{infix,suffix}: a trick for compatibility with different storage
  # names, see below.
  
  
  # omit any indices if requested, e.g. where y_var_name is undefined
  if (length(omit_indices) == 0){
    remaining_indices <- 1:11
    distances_remaining <- Dist.vec
  } else {
    # omit the variables on the indicated positions
    remaining_indices <- (1:11)[-omit_indices]
    distances_remaining <- Dist.vec[remaining_indices]
  }
  dist_num <- length(distances_remaining)
  
  # select people from the correct age group + only remaining distances
  Age.G = levels(Data$AgeGroup)[floor(Age/10)]
  Data <- Data[(Data$dis %in% distances_remaining &
                  Data$AgeGroup==Age.G),]
  
  
  # Results first for female, then for male runners
  sexes <- c('F', 'M')
  sexes_full <- c('Female', 'Male')
  
  plot_panels <-
    lapply(X = seq_along(sexes),
           function(sex_index){
             ### get appropriate (sub)set of data
             Sex = sexes[sex_index]
             Data.x = Data[Data$Sex == Sex,]
             
             if (sex_index == 1){
               box_location <- female_box_location
             } else {
               box_location <- male_box_location
             }
             
             if(nrow(Data.x)/dist_num > SampleSize){
               set.seed(Seed)
               bibs = sample(unique(Data.x$bib),SampleSize)
               Data.x = Data.x[Data.x$bib %in% bibs,]
             }
             
             ### get one plot panel with estimates according to all three methods
             sex = as.double(Sex == 'M')
             AgeSc = (Age-18)/10
             
             
             ## (1) PWE - cubic
             if (CM_formula == "SexAge"){
               cp.cubic = cbind(Data.x$dist,Data.x$dist^2,Data.x$dist^3)
               slope = Data.x$Slope
               x.cubic = cbind(rep(1,dist_num),rep(sex,dist_num),cp.cubic,
                               cp.cubic*sex,rep(AgeSc,dist_num),slope)
             }
             
             if (CM_formula == "SexAgeInter"){
               cp.cubic = cbind(Data.x$dist,Data.x$dist^2,Data.x$dist^3)
               slope = Data.x$Slope
               x.cubic = cbind(rep(1,dist_num),rep(sex,dist_num),cp.cubic,
                               cp.cubic*sex,rep(AgeSc,dist_num),slope,
                               rep(sex*AgeSc, dist_num), cp.cubic*AgeSc)
             }
             # get estimates per quantile level
             est.tau.cubic = lapply(X = seq_along(qu_lvls),function(x){
               tau.val = qu_lvls[x]*100
               name.x = paste(paste0(CM_path,y_var_name),'pwe',tau.val, CM_infix,
                              'eps', CM_eps, CM_suffix,
                              sep='.')
               est.x = tryCatch(list.load(name.x),error=function(e){
                 message(paste('no cubic results for tau=',tau.val))})
               
               betaest.x = est.x$Beta
               output = data.frame(predQu =  x.cubic%*%betaest.x,
                                   dis = Data.x$dis, bib = x, artificial = 'dash')
               output$tau = qu_lvls[x]
               output$tau = as.factor(output$tau)
               return(output)
             })
             est.tau.cubic = do.call('rbind',est.tau.cubic)
             est.tau.cubic$tau = as.factor(est.tau.cubic$tau)
             
             
             ## (2) PWE - unstructured
             if (UM_formula == "SexAge"){
               cp = rbind(rep(0,dist_num-1),diag(dist_num-1))
               x.uns = cbind(rep(1,dist_num),rep(sex,dist_num),cp,
                             rep(AgeSc,dist_num),cp*sex)
             }
             if (UM_formula == "SexAgeInter"){
               cp = rbind(rep(0,dist_num-1),diag(dist_num-1))
               x.uns = cbind(rep(1,dist_num),rep(sex,dist_num),
                             rep(AgeSc,dist_num),cp,rep(sex*AgeSc,dist_num),
                             cp*sex,cp*AgeSc)
             }
             
             est.tau.uns = lapply(X = seq_along(qu_lvls),function(x){
               tau.val = qu_lvls[x]*100
               name.x = paste(paste0(UM_path,y_var_name),'pwe',tau.val, UM_infix,
                              'eps',UM_eps,UM_suffix,
                              sep='.')
               est.x = tryCatch(list.load(name.x),error=function(e){
                 message(paste('no unstr. results for tau=',tau.val))})
               
               betaest.x = est.x$Beta
               output = data.frame(predQu =  x.uns%*%betaest.x,
                                   dis = Data.x$dis, bib = x, artificial = 'solid')
               output$tau = qu_lvls[x]
               #output$tau =as.factor(output$tau)
               return(output)
             })
             est.tau.uns = do.call('rbind',est.tau.uns)
             est.tau.uns$tau = as.factor(est.tau.uns$tau)
             
             
             ## (3) UQR
             # getting estimates over all tau, all distances (no SE)
             uqr.results <- tryCatch(list.load(
               paste0(UQR_path,'UQR.', UQR_formula,".", y_var_name, '.Rdata')),error=function(e){
                 message(paste('no UQR results for ',y_var_name))})
             all_UQR_df <- get_all_est_SE(summary_list = uqr.results,
                                          taus = qu_lvls)
             
             
             if (UQR_formula == "SexAge"){
               covariate <- c(1, sex, AgeSc) 
             }
             
             if (UQR_formula == "SexAgeInter"){
               covariate <- c(1, sex, AgeSc, sex*AgeSc) 
             } 
             
             est.indices <- 1 + 3*(seq_along(qu_lvls)-1)
             beta.estimates <- all_UQR_df[,est.indices]
             
             # will recycle covariate over all distances:
             comp.wise.prods <- covariate*beta.estimates
             
             # create data frame with (# distances) rows and (# taus)*4 columns
             # (each time predicted quantile, dis, bib + additional value to
             # force the legends)
             qu.estimates.list <- lapply(
               X=seq_along(qu_lvls),
               function(j){
                 output <- sum.by.k.and.append(
                   comp.wise.prods[,j],to.be.appended =
                     data.frame(dis = Data.x$dis, bib = j, artificial = 'dot'),
                   k = length(covariate))
                 output$tau = qu_lvls[j]
                 output$tau = as.factor(output$tau)
                 return(output)
               }
             )
             
             all.uqr.out <- do.call('rbind',qu.estimates.list)
             
             est.tau.cubic <- cbind(distinct(est.tau.cubic), method = 'CM')
             est.tau.uns <- cbind(distinct(est.tau.uns), method = 'UM')
             est.tau.uqr <- cbind(distinct(all.uqr.out), method = 'UQR')
             
             
             ## (4) plotting part: very ad hoc, workaround solution for having
             # 2 legends!
             # artificially repeat data three times, to introduce 3 linetypes
             # in order to make them occur for the legend (for the curves, the
             # dotted and dashed line will be overlaid by solid anyway)
             repData.x = cbind(rbind(Data.x, Data.x, Data.x),
                               artificial = as.factor(rep(c('dot', 'dash', 'solid'),
                                                          each = nrow(Data.x))))
             # print(head(repData.x))
             # print(tail(repData.x))
             
             # basis plot layer
             # https://stackoverflow.com/questions/53581787/ggplot2-combining-group-color-and-linetype
             # https://stackoverflow.com/questions/68723238/changing-the-color-of-the-legend-linetype-in-ggplot
             # https://stackoverflow.com/questions/21734557/r-ggplot-ylim-doesnt-work
             p = ggplot(repData.x,
                        aes(x = dis, y = .data[[y_var_name]], group = interaction(bib, artificial))) +
               geom_line(aes(x = dis, y = .data[[y_var_name]], linetype = artificial),
                         color='grey') + xlab('Distance (km)') +
               theme(text = element_text(size=15)) +
               scale_linetype_manual(name = "Methods",
                                     values = c(1, 2, 3),
                                     labels = c("UM", "CM", "UQR")) +
               guides(linetype = guide_legend(override.aes = list(color = "black")))
             
             p.pwe.uqr = p + geom_line(data=est.tau.cubic,
                                       aes(x = dis, y = predQu, color = tau),
                                       linetype = "dashed", linewidth = 1) +
               geom_line(data=est.tau.uns,
                         aes(x = dis, y = predQu, color = tau),
                         linetype = "solid", linewidth = 1) +
               geom_line(data=est.tau.uqr,
                         aes(x = dis, y = predQu, color = tau),
                         linetype = "dotted", linewidth = 1) +
               scale_colour_manual(name = "Levels", values =
                                     brewer.pal(length(qu_lvls), "YlGnBu")) +
               ylab(y_label) + coord_cartesian(ylim = y_range) +
               scale_y_continuous(breaks = y_ticks) +
               geom_label(
                 label= sexes_full[sex_index],
                 x = box_location[1],
                 y = box_location[2],
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = 0.35,
                 color = "black",
                 fill="white"
               )
             
             ### Separate saving if requested
             if (save_separate_panels){
               ggsave(path = storage_path,
                      filename = paste0("LQC_", y_var_name,"_", storage_suffix,'_',Sex,".png"),
                      width = 6, height = 4, device='png')
             }
             
             return(p.pwe.uqr)
           })
  
  # Combine results in one plot
  pCombi <- wrap_plots(plot_panels, ncol = 2) +
    plot_layout(guides = 'collect', axes = 'collect_x',
                axis_titles = 'collect') +
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.47),
                                  text = element_text(size=15)))
  ggsave(path = storage_path,
         filename = paste0("LQC_", y_var_name,"_", storage_suffix,"_Combi.png"),
         width = 10, height = 4, device='png')
  
}

# for SplitVec (different cubic model there)
plot_LQC_splitVec<- function(Data, y_var_name, y_label,y_range, y_ticks,
                             female_box_location, male_box_location,
                             qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                             plot_title, storage_path, storage_suffix,
                             UM_path, CCM_path, UQR_path,
                             UM_formula, CCM_formula, UQR_formula,
                             UM_infix, CCM_infix, UM_suffix, CCM_suffix,
                             UM_eps = 0.1, CCM_eps = 0.1,
                             Age = 35, SampleSize = 50, Seed = 123,
                             omit_indices = c(11),
                             save_separate_panels = FALSE){
  
  # omit any indices if requested, e.g. where y_var_name is undefined
  if (length(omit_indices) == 0){
    remaining_indices <- 1:11
    distances_remaining <- Dist.vec
  } else {
    # omit the variables on the indicated positions
    remaining_indices <- (1:11)[-omit_indices]
    distances_remaining <- Dist.vec[remaining_indices]
  }
  dist_num <- length(distances_remaining)
  
  # select people from the correct age group + only remaining distances
  Age.G = levels(Data$AgeGroup)[floor(Age/10)]
  Data <- Data[(Data$dis %in% distances_remaining &
                  Data$AgeGroup==Age.G),]
  
  
  # Results first for female, then for male runners
  sexes <- c('F', 'M')
  sexes_full <- c('Female', 'Male')
  
  plot_panels <-
    lapply(X = seq_along(sexes),
           function(sex_index){
             ### get appropriate (sub)set of data
             Sex = sexes[sex_index]
             Data.x = Data[Data$Sex == Sex,]
             
             if (sex_index == 1){
               box_location <- female_box_location
             } else {
               box_location <- male_box_location
             }
             
             if(nrow(Data.x)/dist_num > SampleSize){
               set.seed(Seed)
               bibs = sample(unique(Data.x$bib),SampleSize)
               Data.x = Data.x[Data.x$bib %in% bibs,]
             }
             
             ### get one plot panel with estimates according to all three methods
             sex = as.double(Sex == 'M')
             AgeSc = (Age-18)/10
             
             ## (1) PWE - CONSTRAINED cubic
             if (CCM_formula == "SexAge"){
               rd = Data.x$revdist
               cp.cubic = cbind(rd,rd^2,rd^3)
               rSlope = Data.x$revSlope
               x.cubconstr = cbind(rSlope, cp.cubic,cp.cubic*sex,rd*AgeSc)
             }
             
             if (CCM_formula == "SexAgeInter"){
               rd = Data.x$revdist
               cp.cubic = cbind(rd, rd^2, rd^3)
               rSlope = Data.x$revSlope
               x.cubconstr = cbind(rSlope,cp.cubic,cp.cubic*sex,cp.cubic*AgeSc,rd*sex*AgeSc)
             }
             
             
             
             # get estimates per quantile level
             est.tau.cubconstr = lapply(X = seq_along(qu_lvls),function(x){
               tau.val = qu_lvls[x]*100
               name.x = paste(paste0(CCM_path,y_var_name),'pwe',tau.val, CCM_infix,
                              'eps', CCM_eps, CCM_suffix,
                              sep='.')
               
               
               est.x = tryCatch(list.load(name.x),error=function(e){
                 message(paste('no constrained cubic results for tau=',tau.val))})
               
               betaest.x = est.x$Beta
               output = data.frame(predQu =  1 + x.cubconstr%*%betaest.x,
                                   dis = Data.x$dis, bib = x, artificial = 'dash')
               output$tau = qu_lvls[x]
               output$tau = as.factor(output$tau)
               return(output)
             })
             
             est.tau.cubic = do.call('rbind',est.tau.cubconstr)
             est.tau.cubic$tau = as.factor(est.tau.cubic$tau)
             
             ## (2) PWE - unstructured
             if (UM_formula == "SexAge"){
               cp = rbind(rep(0,dist_num-1),diag(dist_num-1))
               x.uns = cbind(rep(1,dist_num),rep(sex,dist_num),cp,
                             rep(AgeSc,dist_num),cp*sex)
             }
             if (UM_formula == "SexAgeInter"){
               cp = rbind(rep(0,dist_num-1),diag(dist_num-1))
               x.uns = cbind(rep(1,dist_num),rep(sex,dist_num),
                             rep(AgeSc,dist_num),cp,rep(sex*AgeSc,dist_num),
                             cp*sex,cp*AgeSc)
             }
             
             est.tau.uns = lapply(X = seq_along(qu_lvls),function(x){
               tau.val = qu_lvls[x]*100
               name.x = paste(paste0(UM_path,y_var_name),'pwe',tau.val,UM_infix,
                              'eps',UM_eps,UM_suffix,
                              sep='.')
               est.x = tryCatch(list.load(name.x),error=function(e){
                 message(paste('no unstr. results for tau=',tau.val))})
               
               betaest.x = est.x$Beta 
               output = data.frame(predQu =  x.uns%*%betaest.x,
                                   dis = Data.x$dis, bib = x, artificial = 'solid')
               output$tau = qu_lvls[x]
               #output$tau =as.factor(output$tau)
               return(output)
             })
             est.tau.uns = do.call('rbind',est.tau.uns)
             est.tau.uns$tau = as.factor(est.tau.uns$tau)
             
             ## (3) UQR
             # getting estimates over all tau, all distances (no SE)
             uqr.results <- tryCatch(list.load(
               paste0(UQR_path,'UQR.', UQR_formula, ".", y_var_name, '.Rdata')),error=function(e){
                 message(paste('no UQR results for ',y_var_name))})
             all_UQR_df <- get_all_est_SE(summary_list = uqr.results,
                                          taus = qu_lvls)
             
             if (UQR_formula == "SexAge"){
               covariate <- c(1, sex, AgeSc) 
             }
             
             if (UQR_formula == "SexAgeInter"){
               covariate <- c(1, sex, AgeSc, sex*AgeSc) 
             } 
             
             est.indices <- 1 + 3*(seq_along(qu_lvls)-1)
             beta.estimates <- all_UQR_df[,est.indices]
             
             # will recycle covariate over all distances:
             comp.wise.prods <- covariate*beta.estimates
             
             # create data frame with (# distances) rows and (# taus)*4 columns
             # (each time predicted quantile, dis, bib + additional value to
             # force the legends)
             qu.estimates.list <- lapply(
               X=seq_along(qu_lvls),
               function(j){
                 output <- sum.by.k.and.append(
                   comp.wise.prods[,j],to.be.appended =
                     data.frame(dis = Data.x$dis, bib = j, artificial = 'dot'),
                   k = length(covariate))
                 output$tau = qu_lvls[j]
                 output$tau = as.factor(output$tau)
                 return(output)
               }
             )                             
             
             all.uqr.out <- do.call('rbind',qu.estimates.list)
             
             est.tau.cubic <- cbind(distinct(est.tau.cubic), method = 'CCM')
             est.tau.uns <- cbind(distinct(est.tau.uns), method = 'UM')
             est.tau.uqr <- cbind(distinct(all.uqr.out), method = 'UQR')
             
             
             ## (4) plotting part: very ad hoc, workaround solution for having
             # 2 legends!
             # artificially repeat data three times, to introduce 3 linetypes
             # in order to make them occur for the legend (for the curves, the
             # dotted and dashed line will be overlaid by solid anyway)
             repData.x = cbind(rbind(Data.x, Data.x, Data.x),
                               artificial = as.factor(rep(c('dot', 'dash', 'solid'),
                                                          each = nrow(Data.x))))
             # print(head(repData.x))
             # print(tail(repData.x))
             
             # basis plot layer
             # https://stackoverflow.com/questions/53581787/ggplot2-combining-group-color-and-linetype
             # https://stackoverflow.com/questions/68723238/changing-the-color-of-the-legend-linetype-in-ggplot
             # https://stackoverflow.com/questions/21734557/r-ggplot-ylim-doesnt-work
             p = ggplot(repData.x,
                        aes(x = dis, y = .data[[y_var_name]], group = interaction(bib, artificial))) +
               geom_line(aes(x = dis, y = .data[[y_var_name]], linetype = artificial),
                         color='grey') + xlab('Distance (km)') +
               theme(text = element_text(size=15)) +
               scale_linetype_manual(name = "Methods",
                                     values = c(1, 2, 3),
                                     labels = c("UM", "CCM", "UQR")) +
               guides(linetype = guide_legend(override.aes = list(color = "black")))
             
             p.pwe.uqr = p + geom_line(data=est.tau.cubic,
                                       aes(x = dis, y = predQu, color = tau),
                                       linetype = "dashed", linewidth = 1) +
               geom_line(data=est.tau.uns,
                         aes(x = dis, y = predQu, color = tau),
                         linetype = "solid", linewidth = 1) +
               geom_line(data=est.tau.uqr,
                         aes(x = dis, y = predQu, color = tau),
                         linetype = "dotted", linewidth = 1) +
               scale_colour_manual(name = "Levels", values =
                                     brewer.pal(length(qu_lvls), "YlGnBu")) +
               ylab(y_label) + coord_cartesian(ylim = y_range) +
               scale_y_continuous(breaks = y_ticks) +
               geom_label(
                 label= sexes_full[sex_index],
                 x = box_location[1],
                 y = box_location[2],
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = 0.35,
                 color = "black",
                 fill="white"
               )
             
             ### Separate saving if requested
             if (save_separate_panels){
               ggsave(path = storage_path,
                      filename = paste0("LQC_", y_var_name,'_',storage_suffix,'_',Sex,".png"),
                      width = 6, height = 4, device='png')
             }
             
             return(p.pwe.uqr)
           })
  
  # Combine results in one plot
  pCombi <- wrap_plots(plot_panels, ncol = 2) +
    plot_layout(guides = 'collect', axes = 'collect_x',
                axis_titles = 'collect') +
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.47),
                                  text = element_text(size=15)))
  ggsave(path = storage_path,
         filename = paste0("LQC_", y_var_name,'_',storage_suffix,"_Combi.png"),
         width = 10, height = 4, device='png')
  
}


### Functions for comparison SexAge variant (with less interaction terms) to
# SexAgeInter variant (including sex-age interaction, as well as age-distance
# interaction): visual comparison + AIC value. Comparison is done for all models
# (UM, CM, UQR)
compare_LQC<- function(Data, y_var_name, y_label, y_range, y_ticks,
                       female_box_location, male_box_location,
                       qu_lvls = c(0.1,0.25,0.5,0.75,0.9),
                       plot_title, storage_path,
                       modelType = "UM",
                       res_formula_vec, res_path_vec, res_eps_vec,
                       res_infix_vec, res_suffix_vec,
                       Age = 35, SampleSize = 50, Seed = 123,
                       omit_indices = c()){
  # modelType = UM / CM / UQR (type of models to be compared)
  # res_formula_vec: different models within the specified type, that are to be
  # compared, e.g. c(AgeSex, AgeSexInter) compares the model with and without
  # extra interaction terms
  
  # res_path_vec, res_eps_vec: vectors with length corresponding to the number
  # of models to be compared, where each entry corresponds to an UM_path/UM_eps
  # if modelType = UM, to CM_path/CM_eps if modelType = CM, etc. (cf. plot_LQC)
  
  dataSize <- nrow(Data)/11 # for BIC computation, all runners matter, not just
  # those of the age group or sex of interest
  
  # omit any indices if requested, e.g. where y_var_name is undefined
  if (length(omit_indices) == 0){
    remaining_indices <- 1:11
    distances_remaining <- Dist.vec
  } else {
    # omit the variables on the indicated positions
    remaining_indices <- (1:11)[-omit_indices]
    distances_remaining <- Dist.vec[remaining_indices]
  }
  dist_num <- length(distances_remaining)
  
  # select people from the correct age group + only remaining distances
  Age.G = levels(Data$AgeGroup)[floor(Age/10)]
  Data <- Data[(Data$dis %in% distances_remaining &
                  Data$AgeGroup==Age.G),]
  
  
  # Results first for female, then for male runners
  sexes <- c('F', 'M')
  sexes_full <- c('Female', 'Male')
  
  
  plot_panels <- vector(mode = "list", length = 2)
  AIC.vals <- matrix(rep(0, length(res_formula_vec)*length(qu_lvls)),
                     nrow = length(qu_lvls))
  BIC.vals <- matrix(rep(0, length(res_formula_vec)*length(qu_lvls)),
                     nrow = length(qu_lvls))
  est.list <- vector(mode = "list", length = length(res_formula_vec))
  
  for (sex_index in seq_along(sexes)){
    ### get appropriate (sub)set of data
    Sex = sexes[sex_index]
    Data.x = Data[Data$Sex == Sex,]
    
    if (sex_index == 1){
      box_location <- female_box_location
    } else {
      box_location <- male_box_location
    }
    
    if(nrow(Data.x)/dist_num > SampleSize){
      set.seed(Seed)
      bibs = sample(unique(Data.x$bib),SampleSize)
      Data.x = Data.x[Data.x$bib %in% bibs,]
    }
    
    sex = as.double(Sex == 'M')
    AgeSc = (Age-18)/10
    
    ### one panel with estimates according to all models to be compared
    # + compute (pseudo-)AIC values in case of PWE
    
    ## PWE - cubic
    if (modelType %in% c("CM", "CCM")){ # CCM ~ splitVec, CM ~ other variables
      
      for (model.index in seq_along(res_formula_vec)){
        CM_formula <- res_formula_vec[model.index]
        CM_path <- res_path_vec[model.index]
        CM_eps <- res_eps_vec[model.index]
        CM_infix <- res_infix_vec[model.index]
        CM_suffix <- res_suffix_vec[model.index]
        
        if (y_var_name == "splitVec"){ # reversed distance and slope
          if (CM_formula == "SexAge"){
            rd = Data.x$revdist
            cp.cubic = cbind(rd, rd^2, rd^3)
            rslope = Data.x$revSlope
            x.cubic = cbind(rslope,cp.cubic,cp.cubic*sex,rd*AgeSc)
          }
          
          if (CM_formula == "SexAgeInter"){
            rd = Data.x$revdist
            cp.cubic = cbind(rd, rd^2, rd^3)
            rslope = Data.x$revSlope
            x.cubic = cbind(rslope,cp.cubic,cp.cubic*sex,cp.cubic*AgeSc,rd*sex*AgeSc)
          }
          
        } else { # normal distance and slope
          if (CM_formula == "SexAge"){
            cp.cubic = cbind(Data.x$dist,Data.x$dist^2,Data.x$dist^3)
            slope = Data.x$Slope
            x.cubic = cbind(rep(1,dist_num),rep(sex,dist_num),cp.cubic,
                            cp.cubic*sex,rep(AgeSc,dist_num),slope)
          }
          
          if (CM_formula == "SexAgeInter"){
            cp.cubic = cbind(Data.x$dist,Data.x$dist^2,Data.x$dist^3)
            slope = Data.x$Slope
            x.cubic = cbind(rep(1,dist_num),rep(sex,dist_num),cp.cubic,
                            cp.cubic*sex,rep(AgeSc,dist_num),slope,
                            rep(sex*AgeSc, dist_num), cp.cubic*AgeSc)
          }
        }
        
        # get estimates per quantile level + compute AIC/BIC
        est.tau.cubic = lapply(X = seq_along(qu_lvls),function(x){
          tau.val = qu_lvls[x]*100
          name.x = paste(paste0(CM_path,y_var_name),'pwe',tau.val, CM_infix,
                         'eps', CM_eps, CM_suffix,
                         sep='.')
          est.x = tryCatch(list.load(name.x),error=function(e){
            message(paste('no cubic results for tau=',tau.val))})
          
          betaest.x = est.x$Beta
          
          output = data.frame(predQu =  1 + x.cubic%*%betaest.x,
                              dis = Data.x$dis, bib = x,
                              artificial = model.index)
          output$tau = qu_lvls[x]
          output$aic = 2*length(betaest.x) - 2*est.x$ll
          output$bic = log(dataSize)*length(betaest.x) - 2*est.x$ll
          return(output)
        })
        est.tau.cubic = do.call('rbind',est.tau.cubic)
        est.tau.cubic$tau = as.factor(est.tau.cubic$tau)
        est.tau.cubic$artificial = as.factor(est.tau.cubic$artificial)
        
        est.tau.cubic <- cbind(distinct(est.tau.cubic), method = 'CM')
        
        est.list[[model.index]] <- est.tau.cubic
        if (sex_index == 1){ # identical for both sexes, so compute only once
          AIC.vals[, model.index] <- est.tau.cubic$aic[1 + dist_num*(seq_along(qu_lvls) - 1)]
          BIC.vals[, model.index] <- est.tau.cubic$bic[1 + dist_num*(seq_along(qu_lvls) - 1)]
        }
      }
      
    }
    
    
    ## PWE - unstructured
    if (modelType == "UM"){
      
      for (model.index in seq_along(res_formula_vec)){
        UM_formula <- res_formula_vec[model.index]
        UM_path <- res_path_vec[model.index]
        UM_eps <- res_eps_vec[model.index]
        UM_infix <- res_infix_vec[model.index]
        UM_suffix <- res_suffix_vec[model.index]
        
        if (UM_formula == "SexAge"){
          cp = rbind(rep(0,dist_num-1),diag(dist_num-1))
          x.uns = cbind(rep(1,dist_num),rep(sex,dist_num),cp,
                        rep(AgeSc,dist_num),cp*sex)
        }
        if (UM_formula == "SexAgeInter"){
          cp = rbind(rep(0,dist_num-1),diag(dist_num-1))
          x.uns = cbind(rep(1,dist_num),rep(sex,dist_num),
                        rep(AgeSc,dist_num),cp,rep(sex*AgeSc,dist_num),
                        cp*sex,cp*AgeSc)
        }
        
        # get estimates per quantile level + compute AIC/BIC
        est.tau.uns = lapply(X = seq_along(qu_lvls),function(x){
          tau.val = qu_lvls[x]*100
          name.x = paste(paste0(UM_path,y_var_name),'pwe',tau.val, UM_infix,
                         'eps', UM_eps, UM_suffix,
                         sep='.')
          est.x = tryCatch(list.load(name.x),error=function(e){
            message(paste('no unstr. results for tau=',tau.val))})
          
          betaest.x = est.x$Beta 
          output = data.frame(predQu =  x.uns%*%betaest.x,
                              dis = Data.x$dis, bib = x,
                              artificial = model.index)
          output$tau = qu_lvls[x]
          output$aic = 2*length(betaest.x) - 2*est.x$ll
          output$bic = log(dataSize)*length(betaest.x) - 2*est.x$ll
          return(output)
        })
        est.tau.uns = do.call('rbind',est.tau.uns)
        est.tau.uns$tau = as.factor(est.tau.uns$tau)
        est.tau.uns$artificial = as.factor(est.tau.uns$artificial)
        
        est.tau.uns <- cbind(distinct(est.tau.uns), method = 'UM')
        
        est.list[[model.index]] <- est.tau.uns
        if (sex_index == 1){ # identical for both sexes, so compute only once
          AIC.vals[, model.index] <- est.tau.uns$aic[1 + dist_num*(seq_along(qu_lvls) - 1)]
          BIC.vals[, model.index] <- est.tau.uns$bic[1 + dist_num*(seq_along(qu_lvls) - 1)]
        }
      }
      
    }
    
    
    if (modelType == "UQR"){
      
      for (model.index in seq_along(res_formula_vec)){
        UQR_formula <- res_formula_vec[model.index]
        UQR_path <- res_path_vec[model.index]
        UQR_eps <- res_eps_vec[model.index]
        UQR_infix <- res_infix_vec[model.index]
        UQR_suffix <- res_suffix_vec[model.index]
        
        # getting estimates over all tau, all distances (no SE)
        uqr.results <- tryCatch(list.load(
          paste0(UQR_path,'UQR.', UQR_formula,".", y_var_name, '.Rdata')),error=function(e){
            message(paste('no UQR results for ',y_var_name))})
        all_UQR_df <- get_all_est_SE(summary_list = uqr.results,
                                     taus = qu_lvls)
        
        
        if (UQR_formula == "SexAge"){
          covariate <- c(1, sex, AgeSc) 
        }
        
        if (UQR_formula == "SexAgeInter"){
          covariate <- c(1, sex, AgeSc, sex*AgeSc) 
        } 
        
        est.indices <- 1 + 3*(seq_along(qu_lvls)-1)
        beta.estimates <- all_UQR_df[,est.indices]
        
        # will recycle covariate over all distances:
        comp.wise.prods <- covariate*beta.estimates
        
        # create data frame with (# distances) rows and (# taus)*4 columns
        # (each time predicted quantile, dis, bib + additional value to
        # force the legends)
        qu.estimates.list <- lapply(
          X=seq_along(qu_lvls),
          function(j){
            output <- sum.by.k.and.append(
              comp.wise.prods[,j],to.be.appended =
                data.frame(dis = Data.x$dis, bib = j,
                           artificial = model.index),
              k = length(covariate))
            output$tau = qu_lvls[j]
            output$tau = as.factor(output$tau)
            return(output)
          }
        )
        
        all.uqr.out <- do.call('rbind',qu.estimates.list)
        
        est.tau.uqr <- cbind(distinct(all.uqr.out), method = 'UQR')
        
        
        est.list[[model.index]] <- all.uqr.out
      }
    }
    
    
    # print(est.list)
    # print(AIC.vals)
    
    
    ## (4) plotting part: very ad hoc, workaround solution for having
    # 2 legends!
    # artificially repeat data, to introduce different linetypes
    # in order to make them occur for the legend (for the curves, the
    # dotted and dashed line will be overlaid by solid anyway)
    
    
    if (length(res_formula_vec) == 3){
      repData.x = cbind(rbind(Data.x, Data.x, Data.x),
                        artificial = as.factor(rep(c('dot', 'dash', 'solid'),
                                                   each = nrow(Data.x))))
      # print(head(repData.x))
      # print(tail(repData.x))
      
      # basis plot layer
      # https://stackoverflow.com/questions/53581787/ggplot2-combining-group-color-and-linetype
      # https://stackoverflow.com/questions/68723238/changing-the-color-of-the-legend-linetype-in-ggplot
      # https://stackoverflow.com/questions/21734557/r-ggplot-ylim-doesnt-work
      p = ggplot(repData.x,
                 aes(x = dis, y = .data[[y_var_name]], group = interaction(bib, artificial))) +
        geom_line(aes(x = dis, y = .data[[y_var_name]], linetype = artificial),
                  color='grey') + xlab('Distance (km)') +
        theme(text = element_text(size=15)) +
        scale_linetype_manual(name = "Methods",
                              values = c(1, 2, 3),
                              labels = res_formula_vec) +
        guides(linetype = guide_legend(override.aes = list(color = "black")))
      
      p.pwe.uqr = p + geom_line(data=est.list[[1]],
                                aes(x = dis, y = predQu, color = tau),
                                linetype = "dashed", size = 1) +
        geom_line(data=est.list[[2]],
                  aes(x = dis, y = predQu, color = tau),
                  linetype = "solid", size = 1) +
        geom_line(data=est.list[[3]],
                  aes(x = dis, y = predQu, color = tau),
                  linetype = "dotted", size = 1) +
        scale_colour_manual(name = "Levels", values =
                              brewer.pal(length(qu_lvls), "YlGnBu")) +
        ylab(y_label) + coord_cartesian(ylim = y_range) +
        scale_y_continuous(breaks = y_ticks) +
        geom_label(
          label= sexes_full[sex_index],
          x = box_location[1],
          y = box_location[2],
          label.padding = unit(0.40, "lines"), # Rectangle size around label
          label.size = 0.35,
          color = "black",
          fill="white"
        )
    }
    
    if (length(res_formula_vec) == 2){
      repData.x = cbind(rbind(Data.x, Data.x),
                        artificial = as.factor(rep(c('dash', 'solid'),
                                                   each = nrow(Data.x))))
      # print(head(repData.x))
      # print(tail(repData.x))
      
      # basis plot layer
      # https://stackoverflow.com/questions/53581787/ggplot2-combining-group-color-and-linetype
      # https://stackoverflow.com/questions/68723238/changing-the-color-of-the-legend-linetype-in-ggplot
      # https://stackoverflow.com/questions/21734557/r-ggplot-ylim-doesnt-work
      p = ggplot(repData.x,
                 aes(x = dis, y = .data[[y_var_name]], group = interaction(bib, artificial))) +
        geom_line(aes(x = dis, y = .data[[y_var_name]], linetype = artificial),
                  color='grey') + xlab('Distance (km)') +
        theme(text = element_text(size=15)) +
        scale_linetype_manual(name = "Methods",
                              values = c(1, 2),
                              labels = res_formula_vec) +
        guides(linetype = guide_legend(override.aes = list(color = "black")))
      
      p.pwe.uqr = p + geom_line(data=est.list[[1]],
                                aes(x = dis, y = predQu, color = tau),
                                linetype = "dashed", size = 1) +
        geom_line(data=est.list[[2]],
                  aes(x = dis, y = predQu, color = tau),
                  linetype = "solid", size = 1) +
        scale_colour_manual(name = "Levels", values =
                              brewer.pal(length(qu_lvls), "YlGnBu")) +
        ylab(y_label) + coord_cartesian(ylim = y_range) +
        scale_y_continuous(breaks = y_ticks) +
        geom_label(
          label= sexes_full[sex_index],
          x = box_location[1],
          y = box_location[2],
          label.padding = unit(0.40, "lines"), # Rectangle size around label
          label.size = 0.35,
          color = "black",
          fill="white"
        )
    }
    
    
    plot_panels[[sex_index]] <- p.pwe.uqr
  }
  
  
  # Combine results in one plot
  pCombi <- wrap_plots(plot_panels, ncol = 2) +
    plot_layout(guides = 'collect', axes = 'collect_x',
                axis_titles = 'collect') +
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.47),
                                  text = element_text(size=15)))
  ggsave(path = storage_path,
         filename = paste0("Comp_LQC_", y_var_name,"_Combi_", modelType, ".png"),
         width = 10, height = 4, device='png')
  
  # print(est.list)
  # print(AIC.vals)
  df.AIC <- as.data.frame(AIC.vals)
  colnames(df.AIC) <- res_formula_vec
  rownames(df.AIC) <- qu_lvls
  
  df.BIC <- as.data.frame(BIC.vals)
  colnames(df.BIC) <- res_formula_vec
  rownames(df.BIC) <- qu_lvls
  return(list(AIC = df.AIC, BIC = df.BIC))
}
# ------------------------------------------------------------------------------


#### Part IV: Comparison to UQR ------------------------------------------------

# get estimates and SE for one distance, put in T pairs of columns (for T taus)
get_est_SE <- function(summary_list, dist_index, taus){
  this_dist_summary <- summary_list[[dist_index]]
  
  if (length(taus) > 1){
    coeff_df <- coefficients(this_dist_summary[[1]])[,c(1,2)]
    for (j in 2:length(taus)){
      coeff_df <- cbind(coeff_df, NA, coefficients(this_dist_summary[[j]])[,c(1,2)])
    }
  } else {
    coeff_df <- coefficients(this_dist_summary)[,c(1,2)]
  }
  return(coeff_df)
}

# get estimates and SE for all distances, put in T pairs of columns (for T taus)
# and in #(model terms)*D rows (for D distances)
get_all_est_SE <- function(summary_list, taus){
  # distance 1: table of size (# model terms) x 2T
  all_coeff_df <- get_est_SE(summary_list = summary_list,
                             dist_index = 1,
                             taus = taus)
  # all other distances 2, ..., D                           
  if (length(summary_list) > 1){
    for (d in 2:length(summary_list)){
      all_coeff_df <- rbind(all_coeff_df, get_est_SE(summary_list = summary_list,
                                                     dist_index = d,
                                                     taus = taus))
    }
  }
  return(all_coeff_df)
}

# assuming model with k coefficients (including intercepts)
sum.by.k.and.append <- function(inp.vec, to.be.appended, k){
  res.number <- length(inp.vec)/k
  predQu <- sapply(seq_len(res.number),
                 function(j){sum(inp.vec[((j-1)*k + 1):((j-1)*k + k)])})
  return(cbind(to.be.appended, predQu))
}


## CI construction
# distIndex = index corresponding to the distance of interest, one of D1, ... D11
# note: theoretical formulas in reporting <-> order of variables for fitting!
getCIstring.fun <- function(effectName, distIndex, est, varMat, varName, UM_formula){
  # UM_formula = "SexAge" / "SexAgeInter"
  # effectName = "interc" / "sex" / "Age" / (only for SexAgeInter option) "SexAge"
  # varName = "intSpeed" / "splitVec"
  
  SE <- sqrt(diag(varMat))
  
  ### Depending on the UM_formula, effectName, and varName, determine basisIndex
  # and distShift s.t. est and SE for distIndex 1 can be found at [indexBasis],
  # and extra term for distIndex > 1 at [indexShifted]
  
  # model before review (less interaction terms)
  if (UM_formula == "SexAge"){
    
    ## for intercepts
    if (effectName == "interc"){
      indexBasis <- 1
      indexShifted <- distIndex + 1
    }
    
    ## for sex effect
    if (effectName == "Sex"){
      indexBasis <- 2
      if (varName == "intSpeed"){
        indexShifted <- distIndex + 12
      }
      if (varName == "splitVec"){
        indexShifted <- distIndex + 11
      }
    }
    
    ## for age effect (common over all distances -> restricting to distIndex = 1 ok)
    if (effectName == "Age"){
      distIndex <- 1 # hack for consistency with other effects
      if (varName == "intSpeed"){
        indexBasis <- 13
      }
      if (varName == "splitVec"){
        indexBasis <- 12
      }
    }
    
    ## age x sex effect: never present in the SexAge model
  }
  
  # model after review (more interaction terms)
  if (UM_formula == "SexAgeInter"){
    ## for intercepts
    if (effectName == "interc"){
      indexBasis <- 1
      indexShifted <- distIndex + 2
    }
    
    ## for sex effect
    if (effectName == "Sex"){
      indexBasis <- 2
      if (varName == "intSpeed"){
        indexShifted <- distIndex + 13
      }
      if (varName == "splitVec"){
        indexShifted <- distIndex + 12
      }
    }
    
    ## for age effect
    if (effectName == "Age"){
      indexBasis <- 3
      if (varName == "intSpeed"){
        indexShifted <- distIndex + 23
      }
      if (varName == "splitVec"){
        indexShifted <- distIndex + 21
      }
    }
    
    ## for age x sex effect (common over all distances -> restricting to distIndex = 1 ok)
    if (effectName == "SexAge"){
      distIndex <- 1 # hack for consistency with other effects
      if (varName == "intSpeed"){
        indexBasis <- 14
      }
      if (varName == "splitVec"){
        indexBasis <- 13
      }
    }
  }
  
  ### getting the corresponding estimates and SE
  
  if (distIndex == 1){
    this.est <- est[indexBasis]
    this.SE <- SE[indexBasis]
  } else { # distance index 2, ..., 10/11
    this.est <- est[indexBasis] + est[indexShifted]
    this.var <- varMat[indexBasis,indexBasis] +
      varMat[indexShifted, indexShifted] +
      2*varMat[indexBasis, indexShifted] #var of sum equals sum of variances + 2*Cov
    this.SE <- sqrt(this.var)
  }
  
  
  ### output: CI
  return(paste0("(", round(this.est - 1.96*this.SE, digits = 4), ", ",
                round(this.est + 1.96*this.SE, digits = 4),")"))
}
get.uqr.CI <- function(df.uqr, this.index){
  this.CI <- paste0("(",round(df.uqr$est[this.index] -
                                1.96*df.uqr$SE[this.index], digits = 4),
                    ", ",
                    round(df.uqr$est[this.index] +
                            1.96*df.uqr$SE[this.index], digits = 4),
                    ")")
  return(this.CI)
}

get.UM.quartiles.CI.df <- function(varName, UM_eps, UM_path, UM_infix, UM_suffix, UM_formula){
  res.25 <- list.load(paste(paste0(UM_path,varName),'pwe',25, UM_infix,
                            'eps', UM_eps, UM_suffix,
                            sep='.'))
  res.50 <- list.load(paste(paste0(UM_path,varName),'pwe',50, UM_infix,
                            'eps', UM_eps, UM_suffix,
                            sep='.'))
  res.75 <- list.load(paste(paste0(UM_path,varName),'pwe',75, UM_infix,
                            'eps', UM_eps, UM_suffix,
                            sep='.'))
  est.25 <- res.25$Beta
  varMat.25 <- res.25$V[1:length(est.25),1:length(est.25)]
  est.50 <- res.50$Beta
  varMat.50 <- res.50$V[1:length(est.50),1:length(est.50)]
  est.75 <- res.75$Beta
  varMat.75 <- res.75$V[1:length(est.75),1:length(est.75)]
  
  if (varName == 'intSpeed'){
    max.distance <- 11
  }
  if (varName == "splitVec"){
    max.distance <- 10
  }
  
  if (UM_formula == "SexAge"){
    tab.res <- do.call('rbind', lapply(1:max.distance, function(j){
      c(getCIstring.fun(effectName = "interc", distIndex = j, est = est.25,
                        varMat = varMat.25, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Sex", distIndex = j, est = est.25,
                        varMat = varMat.25, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Age", distIndex = 1, est = est.25,
                        varMat = varMat.25, varName = varName, UM_formula = UM_formula),
        NA,
        getCIstring.fun(effectName = "interc", distIndex = j, est = est.50,
                        varMat = varMat.50, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Sex", distIndex = j, est = est.50,
                        varMat = varMat.50, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Age", distIndex = 1, est = est.50,
                        varMat = varMat.50, varName = varName, UM_formula = UM_formula),
        NA,
        getCIstring.fun(effectName = "interc", distIndex = j, est = est.75,
                        varMat = varMat.75, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Sex", distIndex = j, est = est.75,
                        varMat = varMat.75, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Age", distIndex = 1, est = est.75,
                        varMat = varMat.75, varName = varName, UM_formula = UM_formula))
    }))
  }
  
  if (UM_formula == "SexAgeInter"){
    tab.res <- do.call('rbind', lapply(1:max.distance, function(j){
      c(getCIstring.fun(effectName = "interc", distIndex = j, est = est.25,
                        varMat = varMat.25, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Sex", distIndex = j, est = est.25,
                        varMat = varMat.25, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Age", distIndex = j, est = est.25,
                        varMat = varMat.25, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "SexAge", distIndex = 1, est = est.25,
                        varMat = varMat.25, varName = varName, UM_formula = UM_formula),
        NA,
        getCIstring.fun(effectName = "interc", distIndex = j, est = est.50,
                        varMat = varMat.50, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Sex", distIndex = j, est = est.50,
                        varMat = varMat.50, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Age", distIndex = j, est = est.50,
                        varMat = varMat.50, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "SexAge", distIndex = 1, est = est.50,
                        varMat = varMat.50, varName = varName, UM_formula = UM_formula),
        NA,
        getCIstring.fun(effectName = "interc", distIndex = j, est = est.75,
                        varMat = varMat.75, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Sex", distIndex = j, est = est.75,
                        varMat = varMat.75, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "Age", distIndex = j, est = est.75,
                        varMat = varMat.75, varName = varName, UM_formula = UM_formula),
        getCIstring.fun(effectName = "SexAge", distIndex = 1, est = est.75,
                        varMat = varMat.75, varName = varName, UM_formula = UM_formula))
    }))
  }
  
  
  return(tab.res)
}

get.UQR.quartiles.CI.df <- function(varName,UM_formula){
  uqr.res <- list.load(paste0(uqr.dir,'/uqr.',UM_formula, '.',varName,'.Rdata'))
  
  uqr.res.df <- get_all_est_SE(summary_list = uqr.res,
                               taus = c(0.1,0.25,0.5,0.75,0.9))
  
  uqr.res.25 <- data.frame(est = uqr.res.df[,4],
                           SE = uqr.res.df[,5])
  
  uqr.res.25$lbound <- round(uqr.res.25$est - 1.96*uqr.res.25$SE, digits = 4)
  uqr.res.25$rbound <- round(uqr.res.25$est + 1.96*uqr.res.25$SE, digits = 4)
  
  uqr.res.50 <- data.frame(est = uqr.res.df[,7],
                           SE = uqr.res.df[,8])
  
  uqr.res.50$lbound <- round(uqr.res.50$est - 1.96*uqr.res.50$SE, digits = 4)
  uqr.res.50$rbound <- round(uqr.res.50$est + 1.96*uqr.res.50$SE, digits = 4)
  
  uqr.res.75 <- data.frame(est = uqr.res.df[,10],
                           SE = uqr.res.df[,11])
  
  uqr.res.75$lbound <- round(uqr.res.75$est - 1.96*uqr.res.75$SE, digits = 4)
  uqr.res.75$rbound <- round(uqr.res.75$est + 1.96*uqr.res.75$SE, digits = 4)
  
  if (varName == 'intSpeed'){
    max.distance <- 11
  }
  if (varName == 'splitVec'){
    max.distance <- 10
  }
  
  
  if (UM_formula == 'SexAge'){
    NumOfCoeff <- 3 # number of UQR model coefficients, including intercept
    
    interc.index <- 1 + NumOfCoeff*(0:(max.distance - 1))
    sex.index <- 2 + NumOfCoeff*(0:(max.distance - 1))
    age.index <- 3 + NumOfCoeff*(0:(max.distance - 1))
    
    tab.res.uqr <- do.call('rbind', lapply(1:max.distance, function(j){
      c(get.uqr.CI(df.uqr = uqr.res.25, this.index = interc.index[j]),
        get.uqr.CI(df.uqr = uqr.res.25, this.index = sex.index[j]),
        get.uqr.CI(df.uqr = uqr.res.25, this.index = age.index[j]),
        NA,
        get.uqr.CI(df.uqr = uqr.res.50, this.index = interc.index[j]),
        get.uqr.CI(df.uqr = uqr.res.50, this.index = sex.index[j]),
        get.uqr.CI(df.uqr = uqr.res.50, this.index = age.index[j]),
        NA,
        get.uqr.CI(df.uqr = uqr.res.75, this.index = interc.index[j]),
        get.uqr.CI(df.uqr = uqr.res.75, this.index = sex.index[j]),
        get.uqr.CI(df.uqr = uqr.res.75, this.index = age.index[j])
      )
    }))
  }
  
  if (UM_formula == 'SexAgeInter'){
    NumOfCoeff <- 4 # number of UQR model coefficients, including intercept
    
    interc.index <- 1 + NumOfCoeff*(0:(max.distance - 1))
    sex.index <- 2 + NumOfCoeff*(0:(max.distance - 1))
    age.index <- 3 + NumOfCoeff*(0:(max.distance - 1))
    sexAge.index <- 4 + NumOfCoeff*(0:(max.distance - 1))
    
    tab.res.uqr <- do.call('rbind', lapply(1:max.distance, function(j){
      c(get.uqr.CI(df.uqr = uqr.res.25, this.index = interc.index[j]),
        get.uqr.CI(df.uqr = uqr.res.25, this.index = sex.index[j]),
        get.uqr.CI(df.uqr = uqr.res.25, this.index = age.index[j]),
        get.uqr.CI(df.uqr = uqr.res.25, this.index = sexAge.index[j]),
        NA,
        get.uqr.CI(df.uqr = uqr.res.50, this.index = interc.index[j]),
        get.uqr.CI(df.uqr = uqr.res.50, this.index = sex.index[j]),
        get.uqr.CI(df.uqr = uqr.res.50, this.index = age.index[j]),
        get.uqr.CI(df.uqr = uqr.res.50, this.index = sexAge.index[j]),
        NA,
        get.uqr.CI(df.uqr = uqr.res.75, this.index = interc.index[j]),
        get.uqr.CI(df.uqr = uqr.res.75, this.index = sex.index[j]),
        get.uqr.CI(df.uqr = uqr.res.75, this.index = age.index[j]),
        get.uqr.CI(df.uqr = uqr.res.75, this.index = sexAge.index[j])
      )
    }))
  }
  
  return(tab.res.uqr)
}
# ------------------------------------------------------------------------------


#### Part V: Comparison to classical methods -----------------------------------

# auxiliary function for ANOVA output, always assumes 11 distance indices
get.signif.mat <- function(df.diff){
  mat.signif <- matrix(NA,11,11)
  for (i in 1:10){
    for (j in (i+1):11){
      iString <- ifelse(i < 10,  paste0("0", i), toString(i))
      jString <- ifelse(j < 10,  paste0("0", j), toString(j))
      diffString <- paste0("D", iString, " - D", jString)
      mat.signif[i,j] <- df.diff$p.value[df.diff$`1` == diffString] < 0.05
    }
  }
  return(mat.signif)
}

# auxiliary function for lm output; each time row of est, of SE and of p-values
get.summ <- function(myModel){
  allCoeff <- summary(myModel)$coefficients
  rbind(as.numeric(allCoeff[, "Estimate"]),
        as.numeric(allCoeff[, "Std. Error"]),
        as.numeric(allCoeff[, "Pr(>|t|)"]))
}

# ------------------------------------------------------------------------------


#### Part VI: Exploratory (quantile-based) methods -----------------------------

# creating the "quantiles" (at each distance point) for men and women.
# put append_col_dist = FALSE if no additional column is to be appended
# put append_col_dist = j (between 1 and 11) to append an additional column with
# repetition of the quantile value at the specified distance (to ease logical
# indexing)

allQuantiles <- function(data, var_name, append_col_dist = 0){
  data$qu.lvl <- rep(0, dim(data)[1])
  data$qu.short <- rep(0, dim(data)[1])
  for (dist_index in 1:11){
    current_dist_data <- data[data$dis == Dist.vec[dist_index], ]
    current_dist_ecdf <- ecdf(x = current_dist_data[[var_name]])
    current_qu_long <- current_dist_ecdf(current_dist_data[[var_name]])
    data[data$dis == Dist.vec[dist_index], ]$qu.lvl <- current_qu_long
    
    # also append rounded version, that can be used for easier stratification
    # ensure that 0.020...01 up to 0.03 get assigned to 0.03 -> avoid 'round'
    current_qu_short <- 0.01*ceiling(100*current_qu_long)
    data[data$dis == Dist.vec[dist_index], ]$qu.short <- current_qu_short  
  }
  
  if (append_col_dist %in% (1:11)){
    data$qu.det <- rep(data$qu.short[data$dis == Dist.vec[append_col_dist]],
                       each = 11)
  }
  
  return(data)
}
# ----


### For summary curves ###

# qu_determ_pos = VECTOR of determining positions (each will correspond to a
# row in the plotting; columns are determined by F (left), M (right))
plot_EQCC <- function(data_female, data_male, y_var_name, y_label,
                      qu_var_name, qu_lvls = c(0.25,0.5,0.75),
                      qu_determ_pos = c(2,5,11),
                      plot_title, file_path,
                      female_box_location,
                      male_box_location,
                      omit_indices = c(),
                      save_separate_panels = FALSE){
  # omit any indices if requested, e.g. where y_var_name is undefined
  if (length(omit_indices) == 0){
    remaining_indices <- 1:11
    distances_remaining <- Dist.vec
  } else {
    # omit the variables on the indicated positions
    remaining_indices <- (1:11)[-omit_indices]
    distances_remaining <- Dist.vec[remaining_indices]
  }
  
  # Compute curves to be plotted for female runners
  data_female_subs <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             # extend data with all quantile levels w.r.t. the given variable
             data_ext <- allQuantiles(data = data_female, var_name = qu_var_name,
                                      append_col_dist = qu_determ_pos[determ_index])
             
             # compute pointwise averages of response variable y_var_name
             # summarise in matrix, col ~ distance point, row ~ quantile level
             y_summ_mat <- matrix(rep(0,(length(remaining_indices))*length(qu_lvls)),
                                  ncol = length(remaining_indices), byrow = TRUE)
             
             for (qu_index in seq_along(qu_lvls)){
               this.qu = data_ext[((data_ext$qu.det > (qu_lvls[qu_index] - 0.001)) &
                                     (data_ext$qu.det < (qu_lvls[qu_index] + 0.001))),]
               
               # case distinction: pointwise average within this "quantile"
               for (dist_index in (seq_along(remaining_indices))){
                 y_summ_mat[qu_index, dist_index] =
                   mean(this.qu[[y_var_name]][this.qu$dis ==
                                                (distances_remaining[dist_index])])
               }
             }
             
             return(y_summ_mat)
           })
  
  # Compute curves to be plotted for male runners
  data_male_subs <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             # extend data with all quantile levels w.r.t. the given variable
             data_ext <- allQuantiles(data = data_male, var_name = qu_var_name,
                                      append_col_dist = qu_determ_pos[determ_index])
             
             # compute pointwise averages of response variable y_var_name
             # summarise in matrix, col ~ distance point, row ~ quantile level
             y_summ_mat <- matrix(rep(0,(length(remaining_indices))*length(qu_lvls)),
                                  ncol = length(remaining_indices), byrow = TRUE)
             
             for (qu_index in seq_along(qu_lvls)){
               this.qu = data_ext[((data_ext$qu.det > (qu_lvls[qu_index] - 0.001)) &
                                     (data_ext$qu.det < (qu_lvls[qu_index] + 0.001))),]
               
               # case distinction: pointwise average within this "quantile"
               for (dist_index in (seq_along(remaining_indices))){
                 y_summ_mat[qu_index, dist_index] =
                   mean(this.qu[[y_var_name]][this.qu$dis ==
                                                (distances_remaining[dist_index])])
               }
             }
             
             return(y_summ_mat)
           })
  
  
  # use common range over all data subsets to ensure plotting on the same scale
  y_range <- range(data_female_subs, data_male_subs)
  
  
  # Compute all panels of the plot for female runners (+ save separately)
  female_plots <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             y_summ_mat <- data_female_subs[[determ_index]]
             df_quantiles <- data.frame(cbind(rep(qu_lvls, each =
                                                    length(remaining_indices)),
                                              rep(distances_remaining,length(qu_lvls)),
                                              as.vector(t(y_summ_mat))))
             colnames(df_quantiles) <- c("qu_lvl", "distance","summary")
             
             
             # https://stackoverflow.com/questions/20301922/r-ggplot2-legend-should-be-discrete-and-not-continuous
             p = ggplot(df_quantiles,
                        aes(x = distance, y = .data[["summary"]],
                            group = factor(qu_lvl), color = factor(qu_lvl))
             ) + geom_line(size = 1.2)
             p1 = p + geom_point(size=1.7) +
               xlab('Distance (km)') +
               ylab(y_label) +
               ylim(y_range) +
               scale_colour_manual(name = "Levels", values =
                                     brewer.pal(length(qu_lvls), "YlGnBu")) +
               theme(legend.position ="right") +
               theme(text = element_text(size=15)) +
               theme(panel.background = element_rect(fill = "grey90", colour = "grey90",
                                                     size = 2, linetype = "solid")) +
               geom_label(
                 label=paste0("Female, ", format(round(Dist.vec[qu_determ_pos[determ_index]], 1), nsmall = 1), " km"),
                 x = female_box_location[1],
                 y = female_box_location[2],
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = 0.35,
                 color = "black",
                 fill="white"
               )
             # https://r-graph-gallery.com/275-add-text-labels-with-ggplot2.html#:~:text=Adding%20text%20with%20geom_text()&text=It%20works%20pretty%20much%20the,along%20X%20and%20Y%20axis 
             
             
             if (save_separate_panels){
               # put plot title only on separately saved plot
               p2 = p1 + ggtitle(plot_title)
               
               ggsave(path = file_path,
                      filename = paste0("QCC_d", qu_determ_pos[determ_index],
                                        "_", y_var_name, "_F.png"),
                      width = 6, height = 4, device='png')
             }
             
             return(p1)
           })
  
  # Compute all panels of the plot for female runners (+ save separately)
  male_plots <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             y_summ_mat <- data_male_subs[[determ_index]]
             df_quantiles <- data.frame(cbind(rep(qu_lvls, each =
                                                    length(remaining_indices)),
                                              rep(distances_remaining,length(qu_lvls)),
                                              as.vector(t(y_summ_mat))))
             colnames(df_quantiles) <- c("qu_lvl", "distance","summary")
             
             
             # https://stackoverflow.com/questions/20301922/r-ggplot2-legend-should-be-discrete-and-not-continuous
             p = ggplot(df_quantiles,
                        aes(x = distance, y = .data[["summary"]],
                            group = factor(qu_lvl), color = factor(qu_lvl))
             ) + geom_line(size = 1.2)
             p1 = p + geom_point(size=1.7) +
               xlab('Distance (km)') +
               ylab(y_label) +
               ylim(y_range) +
               scale_colour_manual(name = "Levels", values =
                                     brewer.pal(length(qu_lvls), "YlGnBu")) +
               theme(legend.position ="right") +
               theme(text = element_text(size=15)) +
               theme(panel.background = element_rect(fill = "grey90", colour = "grey90",
                                                     size = 2, linetype = "solid")) +
               geom_label(
                 label=paste0("Male, ", format(round(Dist.vec[qu_determ_pos[determ_index]], 1), nsmall = 1), " km"),
                 x = male_box_location[1],
                 y = male_box_location[2],
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = 0.35,
                 color = "black",
                 fill="white"
               )
             # https://r-graph-gallery.com/275-add-text-labels-with-ggplot2.html#:~:text=Adding%20text%20with%20geom_text()&text=It%20works%20pretty%20much%20the,along%20X%20and%20Y%20axis 
             
             if (save_separate_panels){
               # put plot title only on separately saved plot
               p2 = p1 + ggtitle(plot_title)
               
               ggsave(path = file_path,
                      filename = paste0("QCC_d", qu_determ_pos[determ_index],
                                        "_", y_var_name, "_M.png"),
                      width = 6, height = 4, device='png')
             }
             
             return(p1)
           })
  
  
  all_plots <- c(female_plots, male_plots)
  
  # Combine all female (left) and male (right) plots on one ggplot
  # https://patchwork.data-imaginist.com/reference/plot_layout.html
  # https://patchwork.data-imaginist.com/articles/guides/annotation.html
  # https://patchwork.data-imaginist.com/articles/guides/layout.html
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange/51352933#51352933
  
  
  pCombi <- wrap_plots(all_plots, ncol = 2, byrow = FALSE) +
    plot_layout(guides = 'collect', axes = 'collect_x',
                axis_titles = 'collect') +
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.47),
                                  text = element_text(size=15)))
  ggsave(path = file_path,
         filename = paste0("EQCC_", y_var_name,"_Combi.png"),
         width = 10, height = 4*length(qu_determ_pos), device='png')
}


# Nearly the same as plot_EQCC, but just with some modifications for the
# plotting: specified ticks + extra top space to enable placement of box in the
# corner with Male/Female, without overlaying curves.

# Note: this could be done more efficiently, e.g. y_var_name will always be
# qu.short for EMQC, just performed minimal modifications for convenience
plot_EMQC <- function(data_female, data_male, y_var_name, y_label,
                      qu_var_name, qu_lvls = c(0.25,0.5,0.75),
                      qu_determ_pos = c(2,5,11),
                      plot_title, file_path,
                      female_box_location,
                      male_box_location,
                      omit_indices = c(),
                      save_separate_panels = FALSE){
  # omit any indices if requested, e.g. where y_var_name is undefined
  if (length(omit_indices) == 0){
    remaining_indices <- 1:11
    distances_remaining <- Dist.vec
  } else {
    # omit the variables on the indicated positions
    remaining_indices <- (1:11)[-omit_indices]
    distances_remaining <- Dist.vec[remaining_indices]
  }
  
  # Compute curves to be plotted for female runners
  data_female_subs <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             # extend data with all quantile levels w.r.t. the given variable
             data_ext <- allQuantiles(data = data_female, var_name = qu_var_name,
                                      append_col_dist = qu_determ_pos[determ_index])
             
             # compute pointwise averages of response variable y_var_name
             # summarise in matrix, col ~ distance point, row ~ quantile level
             y_summ_mat <- matrix(rep(0,(length(remaining_indices))*length(qu_lvls)),
                                  ncol = length(remaining_indices), byrow = TRUE)
             
             for (qu_index in seq_along(qu_lvls)){
               this.qu = data_ext[((data_ext$qu.det > (qu_lvls[qu_index] - 0.001)) &
                                     (data_ext$qu.det < (qu_lvls[qu_index] + 0.001))),]
               
               # case distinction: pointwise average within this "quantile"
               for (dist_index in (seq_along(remaining_indices))){
                 y_summ_mat[qu_index, dist_index] =
                   mean(this.qu[[y_var_name]][this.qu$dis ==
                                                (distances_remaining[dist_index])])
               }
             }
             
             return(y_summ_mat)
           })
  
  # Compute curves to be plotted for male runners
  data_male_subs <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             # extend data with all quantile levels w.r.t. the given variable
             data_ext <- allQuantiles(data = data_male, var_name = qu_var_name,
                                      append_col_dist = qu_determ_pos[determ_index])
             
             # compute pointwise averages of response variable y_var_name
             # summarise in matrix, col ~ distance point, row ~ quantile level
             y_summ_mat <- matrix(rep(0,(length(remaining_indices))*length(qu_lvls)),
                                  ncol = length(remaining_indices), byrow = TRUE)
             
             for (qu_index in seq_along(qu_lvls)){
               this.qu = data_ext[((data_ext$qu.det > (qu_lvls[qu_index] - 0.001)) &
                                     (data_ext$qu.det < (qu_lvls[qu_index] + 0.001))),]
               
               # case distinction: pointwise average within this "quantile"
               for (dist_index in (seq_along(remaining_indices))){
                 y_summ_mat[qu_index, dist_index] =
                   mean(this.qu[[y_var_name]][this.qu$dis ==
                                                (distances_remaining[dist_index])])
               }
             }
             
             return(y_summ_mat)
           })
  
  
  # use common range over all data subsets to ensure plotting on the same scale
  y_range <- c(0, 1.1) # enable F/M corner box for meta-quantiles, too
  
  # Compute all panels of the plot for female runners (+ save separately)
  female_plots <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             y_summ_mat <- data_female_subs[[determ_index]]
             df_quantiles <- data.frame(cbind(rep(qu_lvls, each =
                                                    length(remaining_indices)),
                                              rep(distances_remaining,length(qu_lvls)),
                                              as.vector(t(y_summ_mat))))
             colnames(df_quantiles) <- c("qu_lvl", "distance","summary")
             
             
             # https://stackoverflow.com/questions/20301922/r-ggplot2-legend-should-be-discrete-and-not-continuous
             p = ggplot(df_quantiles,
                        aes(x = distance, y = .data[["summary"]],
                            group = factor(qu_lvl), color = factor(qu_lvl))
             ) + geom_line(size = 1.2)
             p1 = p + geom_point(size=1.7) +
               xlab('Distance (km)') +
               ylab(y_label) +
               coord_cartesian(ylim = y_range) +
               scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
               scale_colour_manual(name = "Levels", values =
                                     brewer.pal(length(qu_lvls), "YlGnBu")) +
               theme(legend.position ="right") +
               theme(text = element_text(size=15)) +
               theme(panel.background = element_rect(fill = "grey90", colour = "grey90",
                                                     size = 2, linetype = "solid")) +
               geom_label(
                 label=paste0("Female, ", format(round(Dist.vec[qu_determ_pos[determ_index]], 1), nsmall = 1), " km"),
                 x = female_box_location[1],
                 y = female_box_location[2],
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = 0.35,
                 color = "black",
                 fill="white"
               )
             # https://r-graph-gallery.com/275-add-text-labels-with-ggplot2.html#:~:text=Adding%20text%20with%20geom_text()&text=It%20works%20pretty%20much%20the,along%20X%20and%20Y%20axis 
             
             
             if (save_separate_panels){
               # put plot title only on separately saved plot
               p2 = p1 + ggtitle(plot_title)
               
               ggsave(path = file_path,
                      filename = paste0("QCC_d", qu_determ_pos[determ_index],
                                        "_", y_var_name, "_F.png"),
                      width = 6, height = 4, device='png')
             }
             
             return(p1)
           })
  
  # Compute all panels of the plot for female runners (+ save separately)
  male_plots <-
    lapply(X = seq_along(qu_determ_pos),
           function(determ_index){
             y_summ_mat <- data_male_subs[[determ_index]]
             df_quantiles <- data.frame(cbind(rep(qu_lvls, each =
                                                    length(remaining_indices)),
                                              rep(distances_remaining,length(qu_lvls)),
                                              as.vector(t(y_summ_mat))))
             colnames(df_quantiles) <- c("qu_lvl", "distance","summary")
             
             
             # https://stackoverflow.com/questions/20301922/r-ggplot2-legend-should-be-discrete-and-not-continuous
             p = ggplot(df_quantiles,
                        aes(x = distance, y = .data[["summary"]],
                            group = factor(qu_lvl), color = factor(qu_lvl))
             ) + geom_line(size = 1.2)
             p1 = p + geom_point(size=1.7) +
               xlab('Distance (km)') +
               ylab(y_label) +
               coord_cartesian(ylim = y_range) +
               scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)) +
               scale_colour_manual(name = "Levels", values =
                                     brewer.pal(length(qu_lvls), "YlGnBu")) +
               theme(legend.position ="right") +
               theme(text = element_text(size=15)) +
               theme(panel.background = element_rect(fill = "grey90", colour = "grey90",
                                                     size = 2, linetype = "solid")) +
               geom_label(
                 label=paste0("Male, ", format(round(Dist.vec[qu_determ_pos[determ_index]], 1), nsmall = 1), " km"),
                 x = male_box_location[1],
                 y = male_box_location[2],
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = 0.35,
                 color = "black",
                 fill="white"
               )
             # https://r-graph-gallery.com/275-add-text-labels-with-ggplot2.html#:~:text=Adding%20text%20with%20geom_text()&text=It%20works%20pretty%20much%20the,along%20X%20and%20Y%20axis 
             
             if (save_separate_panels){
               # put plot title only on separately saved plot
               p2 = p1 + ggtitle(plot_title)
               
               ggsave(path = file_path,
                      filename = paste0("QCC_d", qu_determ_pos[determ_index],
                                        "_", y_var_name, "_M.png"),
                      width = 6, height = 4, device='png')
             }
             
             return(p1)
           })
  
  
  all_plots <- c(female_plots, male_plots)
  
  # Combine all female (left) and male (right) plots on one ggplot
  # https://patchwork.data-imaginist.com/reference/plot_layout.html
  # https://patchwork.data-imaginist.com/articles/guides/annotation.html
  # https://patchwork.data-imaginist.com/articles/guides/layout.html
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange/51352933#51352933
  
  
  pCombi <- wrap_plots(all_plots, ncol = 2, byrow = FALSE) +
    plot_layout(guides = 'collect', axes = 'collect_x',
                axis_titles = 'collect') +
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.47),
                                  text = element_text(size=15)))
  ggsave(path = file_path,
         filename = 'EMQC_Combi.png',
         width = 10, height = 4*length(qu_determ_pos), device='png')
}


### For class curves ###

get_panels_EQCC_class <- function(data_sub, y_var_name, y_label,
                                  qu_var_name, qu_lvls = c(0.05,0.1,0.25,0.5,0.75,0.9,0.95),
                                  qu_determ_pos = 0,
                                  box_location,
                                  omit_indices = c(),
                                  overrule_yrange = c()){
  
  # extend data with all quantile levels w.r.t. the given variable
  data_ext <- allQuantiles(data = data_sub, var_name = qu_var_name,
                           append_col_dist = qu_determ_pos)
  
  y_range = range(data_ext[[y_var_name]])
  if (length(overrule_yrange) > 0){
    y_range <- overrule_yrange
  }
  
  plot_panels <-
    lapply(X = seq_along(qu_lvls),
           function(qu_index){
             this.qu = data_ext[((data_ext$qu.det > (qu_lvls[qu_index] - 0.001)) &
                                   (data_ext$qu.det < (qu_lvls[qu_index] + 0.001))),]
             
             size.qu = nrow(this.qu)/11
             
             # avoid too many lines: take random subset (without repetition) of size n
             if(size.qu > 100){ # will always be the case, since only 100 quantile levels
               # and more than 100*100 runners
               set.seed(123)
               bib.x = sample(unique(this.qu$bib),100)
               this.qu = this.qu[this.qu$bib %in% bib.x,]
             }
             
             # determine position at the finish (within this quantile batch)
             this.qu$pos.x = rep(rank(this.qu$Mtime[this.qu$distance=="MAR"]),each=11)
             mid = median(this.qu$pos.x)
             
             ## plotting: points, connecting lines, correct colours, ...
             
             if (length(omit_indices) > 0){
               # omit the variables on the indicated positions
               this.qu <- this.qu[! this.qu$distance %in%
                                    this.qu$distance[omit_indices],]
             }
             
             # actual plotting
             p = ggplot(this.qu,
                        aes(x = dis, y = .data[[y_var_name]], group = bib,
                            color=.data[["pos.x"]])
             ) + geom_line() 
             p1fig = p + geom_point() +
               xlab('Distance (km)') +
               ylab(y_label) +
               coord_cartesian(ylim = y_range) +
               scale_color_gradient2(midpoint = mid,
                                     low = "blue",
                                     high = "red",
                                     mid ='yellow') +
               theme(legend.position ="none") +
               theme(text = element_text(size=15)) +
               geom_label(
                 label=paste0(format(round(Dist.vec[qu_determ_pos], 1),
                                     nsmall = 1), " km, ",
                              format(qu_lvls[qu_index], nsmall = 2)),
                 x = box_location[1],
                 y = box_location[2],
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = 0.35,
                 color = "black",
                 fill="white"
               )
             
             return(p1fig)
           })
  
  return(plot_panels)
}


# ------------------------------------------------------------------------------
