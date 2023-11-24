###### Reconstruction des valeurs propres a partir de la methode de gradient iterative
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
#### niter = nombre iteration max pour la descente de gradient
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations  pour differentes iteration
#### pas =pas de la descente de gradient


TGradMC=function(mc_sample_size=1000,niter=10,vp,df=3,
                 epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=vp)
{
  p=length(vp)
  vp2=init
  vplist=c()
  for (l in 1:niter)
  {
    ####mi se a zero des esperance
    E1=rep(0,p)
    E2=0
    for (i in 1:mc_sample_size)
    {
      Z=rnorm(p)/sqrt(rchisq(1,df=df)) *sqrt((df-2))
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


TFixMC=function(mc_sample_size=1000,niter=10,vp,epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=vp,df=3)
{
  p=length(vp)
  vplist=c()
  vp2=init
  for (l in 1:niter)
  {
    ####mi se a zero des esperance
    E1=rep(0,p)
    E2=0
    for (i in 1:mc_sample_size)
    {
      Z=rnorm(p)/sqrt(rchisq(1,df=df)) *sqrt((df-2) )
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



TRobbinsMC=function(mc_sample_size=1000,vp,epsilon=10^(-8),alpha=0.75,c=2,w=2,samp=mc_sample_size,init=vp,df=3)
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
    Z= rnorm(p)/sqrt(rchisq(1,df=df)) *sqrt((df-2) )
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




Tgen_ech=function(n=500,d=5,pcont=0,df=3,dfcont=1,cont="Student",min=-5,max=5)
{
  Sigma= (diag(d)*2)
  X=c()
  Tclassif=c()
  mean=rep(0,d)

  X=mvtnorm::rmvt(n,delta=mean,sigma=(df-2)/df*Sigma,df=df)

  Tclassif=c(Tclassif,rep(1,n))
  Sigma=diag( (1:d))
  mean=rep(3,d)
  #np=rbinom(1,size = n,prob = 1-pcont)
  X=rbind(X,mvtnorm::rmvt(n,delta=mean,sigma=(df-2)/df*Sigma,df=df))
  #if (n-np >0){
  #X=rbind(X,matrix(rep(mean,n-np),nrow=n-np,byrow = T)+ mvtnorm::rmvt(n = n-np,df = 2,sigma = Sigma))
  #X=rbind(X,mvtnorm::rmvnorm(n,mean=mean,sigma=Sigma))
  #}
  Tclassif=c(Tclassif,rep(2,n))

  Sigma=  (diag((1:d)^(-1)))
  mean=- rep(3,d)
  X=rbind(X,mvtnorm::rmvt(n,delta=mean,sigma=(df-2)/df*Sigma,df=df))
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
      Z=matrix(rt(pcont*n*d,df=dfcont),ncol=d)
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


Tgen_K=function(n=500,d=5,K=3,pcont=0,df=3,dfcont=1,cont="Student",min=-5,max=5,radius=5)
{
  Sigma= (diag(d)*2)
  X=c()
  Tclassif=c()
  mean=rep(0,d)
  for (k in 1:K)
  {
    Z=rnorm(d)
    Z=radius*Z/sqrt(sum(Z^2))
    X=rbind(X,mvtnorm::rmvt(n,delta=Z,sigma=(df-2)/df*Sigma,df=df))
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
      Z=matrix(rt(pcont*n*d,df=dfcont),ncol=d)
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


SimulTMMcontStudent <- function(nk, dfT0, muT0, sigmaT0, delta, dfT1, muT1, sigmaT1){
  # muT0 <- parm$muT0; sigmaT0 <- parm$sigmaT0; delta <- delta; dfT1 <- parm$dfT1; muT1 <- parm$muT1; sigmaT1 <- parm$sigmaT1
  # nk = vecteur des effectifs des composantes
  # dfT0 = ddl de la student du mélange
  # muT0 = matrice des moyennes par composantes
  # sigmaT0 = 'array' des variances par composantes
  # delta = taux de contamination
  # dfT1 = ddl de la student contaminantes
  # muT1 = matrice des moyennes des students contaminantes
  # sigmaT1 = 'matrice'array' des variances pardes students contaminantes
  K <- length(nk); p <- ncol(muT0); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkT0 <- round(nk*(1-delta)); nkA <- nk - nkT0
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvt(n=nkT0[k], delta=muT0[k, ], sigma=sigmaT0[k, , ], df=dfT0),
              mvtnorm::rmvt(n=nkA[k], delta=muT1[k, ], sigma=sigmaT1[k, , ], df=dfT1))
      C[nkCum[k]+nkT0[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvt(n=nk[k], delta=muT0[k, ], sigma=sigmaT0[k, , ], df=dfT0)
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}

SimulTMMcontUniform <- function(nk, dfT, muT, sigmaT, delta, lUnif, uUnif){
  # muT <- parm$muT; sigmaT <- parm$sigmaT; lUnif <- parm$lUnif; uUnif <- parm$uUnif
  # nk = vecteur des effectifs des composantes
  # dfT = ddl de la student du mélange
  # muT = matrice des moyennes par composantes
  # sigmaT = 'array' des variances par composantes
  # delta = taux de contamination
  # lUnif = vecteur des bornes inf des uniformes contaminantes
  # uUnif = vecteur des bornes sup des uniformes contaminantes
  K <- length(nk); p <- ncol(muT); n <- sum(nk)
  nkCum <- c(0, cumsum(nk))
  nkT <- round(nk*(1-delta)); nkA <- nk - nkT
  Z <- C <- rep(0, n); X <- matrix(NA, n, p)
  for(k in 1:K){
    if(nkA[k] > 0){
      X[nkCum[k]+(1:nk[k]), ] <-
        rbind(mvtnorm::rmvt(n=nkT[k], delta=muT[k, ], sigma=sigmaT[k, , ], df=dfT),
              sapply(1:p, function(j){runif(n=nkA[k], min=lUnif[j], max=uUnif[j])}))
      C[nkCum[k]+nkT[k]+(1:nkA[k])] <- 1
    }else{
      X[nkCum[k]+(1:nk[k]), ] <-
        mvtnorm::rmvt(n=nk[k], delta=muT[k, ], sigma=sigmaT[k, , ], df=dfT)
    }
    Z[nkCum[k]+(1:nk[k])] <- k
  }
  randOrder <- rank(runif(n)); Z <- Z[randOrder]; C <- C[randOrder]; X <- X[randOrder, ];
  return(list(Z=Z, C=C, X=X))
}


