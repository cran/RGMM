###### Reconstruction des valeurs propres a partir de la methode de gradient iterative
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
#### niter = nombre iteration max pour la descente de gradient
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations  pour differentes iteration
#### pas =pas de la descente de gradient


GradMC=function(mc_sample_size=1000,niter=10,vp,epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=vp)
{
  p=length(vp)
  vp2=init
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  vplist=c()
  for (l in 1:niter)
  {
    ####mi se a zero des esperance
    E1=rep(0,p)
    E2=0
    for (i in 1:mc_sample_size)
    {
      Z=Y[i,]
      E1= E1 +  Z^2*(sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
      E2= E2 +  (sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
    }
    vp0=vp2
    vp2=vp2 - pas[l]*mc_sample_size^(-1)*vp2*E1 + pas[l]*mc_sample_size^(-1) * (vp)*E2
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==l)) >0)
    {
      vplist=cbind(vplist,vp2)
    }
    #   if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=l,vplist=vplist))
}


###### Reconstruction des valeurs propres a partir de la methode de point fixe
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
#### niter = nombre iteration max pour l'algo
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations de point fixe pour differentes iteration


FixMC=function(mc_sample_size=1000,niter=10,vp,epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=vp)
{
  p=length(vp)
  vplist=c()
  vp2=init
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  for (l in 1:niter)
  {
    ####mi se a zero des esperance
    E1=rep(0,p)
    E2=0
    for (i in 1:mc_sample_size)
    {
      Z=Y[i,]
      E1= E1 +  Z^2*(sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
      E2= E2 +  (sum(( (vp)-vp2*(Z^2))^2) + sum((vp2 * Z^2)%*%t((vp2 * Z^2))) - sum(((vp2 * Z^2)^2))  )^(-0.5)
    }
    vp0=vp2
    vp2 = E2* (vp)*(E1)^(-1)
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==l)) >0)
    {
      vplist=cbind(vplist,vp2)
    }
    #    if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=l))
}


###### Reconstruction des valeurs propres a partir de la methode de gradient recursive
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations pour differentes tailles d'?chantillon



RobbinsMC=function(mc_sample_size=1000,vp,epsilon=10^(-8),alpha=0.75,c=2,w=2,samp=mc_sample_size,init=vp)
{
  p=length(vp)
  vp2=init
  lambda=init
  lambdalist=c()
  vplist=c()
  Y=matrix(rnorm(mc_sample_size*p),ncol=p)
  # X2=matrix(rnorm(mc_sample_size*p),ncol=p)
  slog=1
  for (i in 1:mc_sample_size)
  {
    Z=Y[i,]
    # Z2=X2[i,]
    E1=    Z^2*(sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    vp0=vp2
    E2=    (sum(( (vp)-lambda*(Z^2))^2) + sum((lambda * Z^2)%*%t((lambda * Z^2))) - sum(((lambda * Z^2)^2))  )^(-0.5)
    lambda =lambda  - c*i^(-alpha)*lambda*E1 + c*i^(-alpha)* (vp)*E2
    slog=slog+log(i+1)^w
    vp2=vp2+log(i+1)^w *((slog)^(-1)) *(lambda - vp2)
    eps=sqrt(sum((vp2-vp0)^2))
    if (length(which(samp==i) > 0))
    {
      lambdalist=cbind(lambdalist,lambda)
      vplist=cbind(vplist,vp2)
    }
    #   if ( eps < epsilon ) break;
  }
  return(list(vp=vp2,niter=i, lambdalist=lambdalist, vplist=vplist))
}




Weiszfeld_init <- function (X,init=rep(0,ncol(X)), weights = NULL, epsilon=1e-08, nitermax = 100)
{
  X <- as.matrix(X)
  if (is.null(weights)) weights <- rep(1,nrow(X))

  return(Weiszfeld_init_rcpp(X,init=init, weights= weights, epsilon=epsilon, nitermax=nitermax))
}

###fonction R pour estimer la matrice de covariance mediane via Weiszfeld en reglant l'initialisation


WeiszfeldCov_init <- function(X, init=rep(0,ncol(X)),init_cov=(diag(ncol(X))), weights=NULL, scores=2, epsilon=1e-08, nitermax = 100){
  ### Computation of the Geometric covariation matrix
  ### output : (geometric) median (1 x p numeric vector) and (geometric) median covariation matrix (p x p)
  ### require library(rARPACK)
  X <- as.matrix(X)
  n <- nrow(X)
  if (is.null(weights)){ weights <- rep(1,n)}
  Wmed.est <- Weiszfeld_init_rcpp(X,init=init,weights = weights,epsilon=epsilon,nitermax=nitermax)
  WMCM.est <- WeiszfeldCovMat_init_rcpp(X,init_cov=init_cov, median_est = Wmed.est$median,weights=weights,epsilon=epsilon,nitermax=nitermax)
  if (scores==FALSE){
    return(list(median = Wmed.est$median, covmedian=WMCM.est$median, iterm = Wmed.est$iter, itercov = WMCM.est$iter))
  }
  else {
    ### Computation of the eigenvectors and scores
    vectors <- RSpectra::eigs_sym(WMCM.est$median, scores)$vectors
    vscores = sweep(X,2,Wmed.est$median)%*%vectors
    return(list(median=Wmed.est$median, covmedian=WMCM.est$median, scores=vscores, vectors=vectors, iterm = Wmed.est$iter, itercov = WMCM.est$iter))
  }
}

###fonction RCPP pour estimer la median via la descente de gradient en reglant l'initialisation





Gmedian_init <- function (X, init = rep(0,ncol(X)), gamma = 2, alpha = 0.75, nstart=1, epsilon=1e-08,weights=NULL)
{
  X <- as.matrix(X)
  if (is.null(weights)) weights <- rep(1,nrow(X))
  med.X = Gmedianrowvec_init_rcpp(X,init=init,weights=weights, gamma=gamma, alpha=alpha, nstart=nstart, epsilon=epsilon)
  return(med.X)
}


GmedianCov_init <- function(X, init=rep(0,ncol(X)),init_cov=diag(ncol(X)),weights=NULL, scores=2, gamma=2, gc=2, alpha=0.75, nstart=1){
  ### Computation of the Geometric covariation matrix
  ### with averaged stochastic gradient algorithms
  ### input : X   (n x p matrix, n observations, in dimension p)
  ### output : (geometric) median (1 x p numeric vector) and (geometric) median covariation matrix (p x p)
  ### require library(rARPACK)
  if (is.null(weights)) weights <- rep(1,nrow(X))
  Gmed.est = Gmedian_init(X,init=init,weights=weights,gamma=gamma,alpha=alpha,nstart=nstart)
  GMCM.est = MedianCovMatRow_init_rcpp(X,init_cov=init_cov,Gmedian=Gmed.est,weights=weights,gamma=gc,alpha=alpha,nstart=nstart)
  if (scores==FALSE){
    return(list(median = Gmed.est,covmedian=GMCM.est))
  }
  else {
    ### Computation of the eigenvectors and scores
    vectors <- RSpectra::eigs_sym(GMCM.est, scores)$vectors
    scores = sweep(X,2,Gmed.est)%*%vectors
    return(list(median=Gmed.est,covmedian=GMCM.est,scores=scores,vectors=vectors))
  }
}

gen_ech=function(n=500,d=5,pcont=0,df=1,cont="Student",min=-5,max=5)
{
  Sigma= (diag(d)*2)
  X=c()
  Tclassif=c()
  mean=rep(0,d)

  X=mvtnorm::rmvnorm(n,mean=mean,sigma=Sigma)

  Tclassif=c(Tclassif,rep(1,n))
  Sigma=diag( (1:d))
  mean=rep(3,d)
  #np=rbinom(1,size = n,prob = 1-pcont)
  X=rbind(X,mvtnorm::rmvnorm(n ,mean=mean,sigma=Sigma))
  #if (n-np >0){
  #X=rbind(X,matrix(rep(mean,n-np),nrow=n-np,byrow = T)+ mvtnorm::rmvt(n = n-np,df = 2,sigma = Sigma))
  #X=rbind(X,mvtnorm::rmvnorm(n,mean=mean,sigma=Sigma))
  #}
  Tclassif=c(Tclassif,rep(2,n))

  Sigma=  (diag((1:d)^(-1)))
  mean=- rep(3,d)
  X=rbind(X,mvtnorm::rmvnorm(n ,mean=mean,sigma=Sigma))
  Tclassif=c(Tclassif,rep(3,n))
  #X=X+ rep(5,d)


  if (pcont>0)
  {
    if(cont=='Unif')
    {
      Z=matrix(runif(pcont*n*d,min,max),ncol=d)
    }
    if (cont=='Student')
    {
      Z=matrix(rt(pcont*n*d,df=df),ncol=d)
    }
    I=sample(1:(3*n),size=pcont*n)
    X[I,]=Z
    Tclassif[I]="outliers"
  }


  mel=sample.int(3*n)
  Tclassif=Tclassif[mel]
  X=X[mel,]
  return(list(X=X,classif=Tclassif))
}


gen_K=function(n=500,d=5,K=3,pcont=0,df=1,cont="Student",min=-5,max=5,radius=5)
{
  Sigma= (diag(d)*2)
  X=c()
  Tclassif=c()
  mean=rep(0,d)
  for (k in 1:K)
  {
    Z=rnorm(d)
    Z=radius*Z/sqrt(sum(Z^2))
    X=rbind(X,mvtnorm::rmvnorm(n ,mean=Z))
    Tclassif=c(Tclassif,rep(k,n))
  }

  if (pcont>0)
  {
    if(cont=='Unif')
    {
      Z=matrix(runif(pcont*n*d,min,max),ncol=d)
    }
    if (cont=='Student')
    {
      Z=matrix(rt(pcont*n*d,df=df),ncol=d)
    }
    I=sample(1:(3*n),size=pcont*n)
    X[I,]=Z
    Tclassif[I]="outliers"
  }


  mel=sample.int(3*n)
  Tclassif=Tclassif[mel]
  X=X[mel,]
  return(list(X=X,classif=Tclassif))
}

# Simul function for robust GMM

SimulGMMcontStudent <- function(nk, muG, sigmaG, delta, dfT, muT, sigmaT){
  # muG <- parm$muG; sigmaG <- parm$sigmaG; delta <- delta; dfT <- parm$dfT; muT <- parm$muT; sigmaT <- parm$sigmaT
  # nk = vecteur des effectifs des composantes
  # muG = matrice des moyennes par composantes
  # sigmaG = 'array' des variances par composantes
  # delta = taux de contamination
  # dfT = ddl de la student contaminantes
  # muT = matrice des moyennes des students contaminantes
  # sigmaT = 'matrice'array' des variances pardes students contaminantes
  K <- length(nk); p <- ncol(muG); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkG <- round(nk*(1-delta)); nkA <- nk - nkG
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvnorm(n=nkG[k], mean=muG[k, ], sigma=sigmaG[k, , ]),
              mvtnorm::rmvt(n=nkA[k], delta=muT[k, ], sigma=sigmaT[k, , ], df=dfT))
      C[nkCum[k]+nkG[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvnorm(n=nk[k], mean=muG[k, ], sigma=sigmaG[k, , ])
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}

SimulGMMcontUniform <- function(nk, muG, sigmaG, delta, lUnif, uUnif){
  # muG <- parm$muG; sigmaG <- parm$sigmaG; lUnif <- parm$lUnif; uUnif <- parm$uUnif
  # nk = vecteur des effectifs des composantes
  # muG = matrice des moyennes par composantes
  # sigmaG = 'array' des variances par composantes
  # delta = taux de contamination
  # lUnif = vecteur des bornes inf des uniformes contaminantes
  # uUnif = vecteur des bornes sup des uniformes contaminantes
  K <- length(nk); p <- ncol(muG); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkG <- round(nk*(1-delta)); nkA <- nk - nkG
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvnorm(n=nkG[k], mean=muG[k, ], sigma=sigmaG[k, , ]),
              sapply(1:p, function(j){runif(n=nkA[k], min=lUnif[j], max=uUnif[j])}))
      C[nkCum[k]+nkG[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvnorm(n=nk[k], mean=muG[k, ], sigma=sigmaG[k, , ])
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}

