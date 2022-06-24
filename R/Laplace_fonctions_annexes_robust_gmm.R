###### Reconstruction des valeurs propres a partir de la methode de gradient iterative
##### mc_sample_size= Taille de l'echantillon pour la methode de monte  carlo
#### niter = nombre iteration max pour la descente de gradient
### epsilon= arret si diff entre deux iteration < epsilon
#### samp= si on veut ressortir les valeurs des estimations  pour differentes iteration
#### pas =pas de la descente de gradient


LGradMC=function(mc_sample_size=1000,niter=10,vp,
                 epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=rep(0,length(vp)))
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
      Z=as.vector(LaplacesDemon::rmvl(1,mu=rep(0,p),Sigma=diag(p)))
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


LFixMC=function(mc_sample_size=1000,niter=10,vp,epsilon=10^(-8),pas=rep(1,niter),samp=niter,init=rep(0,length(vp)),df=3)
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
      Z=as.vector(LaplacesDemon::rmvl(1,mu=rep(0,p),Sigma=diag(p)))
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



LRobbinsMC=function(mc_sample_size=1000,vp,epsilon=10^(-8),alpha=0.75,c=2,w=2,samp=mc_sample_size,init=rep(0,length(vp)))
{
  p=as.numeric(length(vp))
  vp2=init
  lambda=init
  lambdalist=c()
  vplist=c()
  # X2=matrix(rnorm(mc_sample_size*p),ncol=p)
  slog=1
  for (i in 1:mc_sample_size)
  {
    Z=as.vector(LaplacesDemon::rmvl(1,mu=rep(0,p),Sigma=diag(p)))
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






Lgen_ech=function(n=500,d=5,pcont=0,df=3,dfcont=1,cont="Student",min=-5,max=5)
{
  Sigma= (diag(d)*2)
  X=c()
  Tclassif=c()
  mean=rep(0,d)

  X=LaplacesDemon::rmvl(n,mu=mean,Sigma=Sigma)

  Tclassif=c(Tclassif,rep(1,n))
  Sigma=diag( (1:d))
  mean=rep(3,d)
  #np=rbinom(1,size = n,prob = 1-pcont)
  X=rbind(X,LaplacesDemon::rmvl(n,mu=mean,Sigma=Sigma))
  #if (n-np >0){
  #X=rbind(X,matrix(rep(mean,n-np),nrow=n-np,byrow = T)+ mvtnorm::rmvt(n = n-np,df = 2,sigma = Sigma))
  #X=rbind(X,mvtnorm::rmvnorm(n,mean=mean,sigma=Sigma))
  #}
  Tclassif=c(Tclassif,rep(2,n))

  Sigma=  (diag((1:d)^(-1)))
  mean=- rep(3,d)
  X=rbind(X,LaplacesDemon::rmvl(n,mu=mean,Sigma=Sigma))
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


Lgen_K=function(n=500,d=5,K=3,pcont=0,df=3,dfcont=1,cont="Student",min=-5,max=5,radius=5)
{
  Sigma= (diag(d)*2)
  X=c()
  Tclassif=c()
  mean=rep(0,d)
  for (k in 1:K)
  {
    Z=rnorm(d)
    Z=radius*Z/sqrt(sum(Z^2))
    X=rbind(X,LaplacesDemon::rmvl(n,mu=Z,Sigma=Sigma,df=df))
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


