RobMM=function(X,nclust=2:5,model="Gaussian",ninit=0,nitermax=50,niterEM=50,niterMC=50,df=3,epsvp=0,
               mc_sample_size=1000, LogLike=-10^10,init='genie',epsPi=10^-4,epsout=-100,
               alpha=0.75,c=ncol(X),w=2,epsilon=10^(-8),
               methodMC="RobbinsMC", par=TRUE,methodMCM="Weiszfeld")
{
  if (model=='Gaussian')
  {
    resultat=RGMM(X,nclust=nclust,ninit=ninit,nitermax=nitermax,niterEM=niterEM,
                  niterMC=niterMC,epsvp=epsvp,mc_sample_size=mc_sample_size,
                  LogLike=LogLike,init=init,epsPi=epsPi,epsout=epsout,
                  alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC, par=par)
  }
  if (model=='Student')
  {
    resultat=RTMM(X,nclust=nclust,ninit=ninit,nitermax=nitermax,niterEM=niterEM,
                  niterMC=niterMC,epsvp=epsvp,mc_sample_size=mc_sample_size,df=df,
                  LogLike=LogLike,init=init,epsPi=epsPi,epsout=epsout,
                  alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC, par=par)
  }
  if (model=='Laplace')
  {
    resultat=RLMM(X,nclust=nclust,ninit=ninit,nitermax=nitermax,niterEM=niterEM,
                  niterMC=niterMC,epsvp=epsvp,mc_sample_size=mc_sample_size,
                  LogLike=LogLike,init=init,epsPi=epsPi,epsout=epsout,
                  alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC, par=par)
  }
  return(resultat)

}

RobVar=function(X,c=2,alpha=0.75,model='Gaussian',methodMCM='Weiszfeld',
                methodMC='Robbins' ,mc_sample_size=1000,init=rep(0,ncol(X)),init_cov=diag(ncol(X)),
                epsilon=10^(-8),w=2,initvp=rep(0,length(vp)),df=3,niterMC=50,
                cgrad=2,niterWeisz=50,epsWeisz=10^-8,alphaMedian=0.75,cmedian=2)
{
  d=ncol(X)
  if(methodMCM=='Weiszfeld')
  {
    mcm=WeiszfeldCov_init(X=X,epsilon=epsWeisz,nitermax=niterWeisz,scores=F)
  }
  if(methodMCM=='Gmedian')
  {
    mcm=GmedianCov_init(X=X,init = init,init_cov =init_cov,scores=0,gamma=cgrad,gc=cmedian)
  }
  eig=eigen(mcm$covmedian)
  vec=eig$vectors
  vp=eig$values
  if (model=='Gaussian')
  {
    if (methodMC=='Robbins')
    {
      vpvar=RobbinsMC(mc_sample_size=mc_sample_size,vp=vp,
                      epsilon=epsilon,alpha=alpha,c=c,w=w,init=initvp)
    }
    if (methodMC=='Fix')
    {
      vpvar=FixMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                  epsilon=epsilon,init=initvp)
    }
    if (methodMC=='Grad')
    {
      if(length(c)> 1)
      {
        vpvar=GradMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                     epsilon=epsilon,init=initvp,pas=c)
      }
      if (length(c) == 1)
      {
        vpvar=GradMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                     epsilon=epsilon,init=initvp,pas=c*((1:(niterMC))^(-0.5)))
      }
    }
  }
  if (model=='Student')
  {
    if (methodMC=='Robbins')
    {
      vpvar=TRobbinsMC(mc_sample_size=mc_sample_size,vp=vp,df=df,
                       epsilon=epsilon,alpha=alpha,c=c,w=w,init=initvp)
    }
    if (methodMC=='Fix')
    {
      vpvar=TFixMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                   epsilon=epsilon,init=initvp,df=df)
    }
    if (methodMC=='Grad')
    {
      if(length(c)> 1)
      {
        vpvar=TGradMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                      epsilon=epsilon,init=initvp,pas=c,df=df)
      }
      if (length(c) == 1)
      {
        vpvar=TGradMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                      epsilon=epsilon,init=initvp,pas=c*((1:(niterMC))^(-0.5)))
      }
    }
  }
  if (model=='Laplace')
  {
    if (methodMC=='Robbins')
    {
      vpvar=LRobbinsMC(mc_sample_size=mc_sample_size,vp=vp,
                       epsilon=epsilon,alpha=alpha,c=c,w=w,init=initvp)
    }
    if (methodMC=='Fix')
    {
      vpvar=LFixMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                   epsilon=epsilon,init=initvp )
    }
    if (methodMC=='Grad')
    {
      if(length(c)> 1)
      {
        vpvar=LGradMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                      epsilon=epsilon,init=initvp,pas=c)
      }
      if (length(c) == 1)
      {
        vpvar=LGradMC(mc_sample_size=mc_sample_size,vp=vp,niter=niterMC,
                      epsilon=epsilon,init=initvp,pas=c*((1:(niterMC))^(-0.5)))
      }
    }
  }



  lambda=vpvar$vp
  variance=t(matrix(vec,ncol=d,byrow=TRUE))%*%diag(lambda)%*%(matrix(vec,ncol=d,byrow=TRUE))
  resultat=list(median=mcm$median,variance=variance,covmedian=mcm$covmedian)
  return(resultat)
}



Gen_MM=function(nk=NA, df=3, mu=NA, Sigma=FALSE, delta=0,cont="Student",
                model="Gaussian", dfcont=1, mucont=FALSE, Sigmacont=FALSE, minU=-20, maxU=20)
{
  if (length(is.na(nk)) ==1)
  {
  if (is.na(nk)==TRUE)
  {
    if (is.matrix(mu)==FALSE)
    {
      nk=rep(500,3)
    }
    if (is.matrix(mu) != FALSE)
    {
      nk=rep(500,nrow(mu))
    }
  }
  if (is.matrix(mu)==FALSE)
  {
    mu=c()
    for (i in 1:length(nk))
    {
      Z=rnorm(3)
      mu=rbind(mu,Z)
    }
  }
  }
  if (is.matrix(mucont)==FALSE)
  {
    mucont=mu
  }
  if (is.array(Sigma)==FALSE)
  {
    p=ncol(mu)
    Sigma <- array(dim=c(length(nk), p, p))
    for (i in 1:length(nk))
    {
      Sigma[i, ,]=diag(p)
    }
  }
  if (is.array(Sigmacont)==FALSE)
  {
    Sigmacont=Sigma
  }
  if (cont=="Student")
  {
    if (model=="Student")
    {
      resultat=SimulTMMcontStudent(nk=nk, dfT0=df, muT0=mu, sigmaT0=Sigma,
                                   delta=delta, dfT1=dfcont,
                                   muT1=mucont, sigmaT1=Sigmacont)
    }
    if (model=='Gaussian')
    {
      resultat=SimulGMMcontStudent(nk=nk, muG=mu, sigmaG=Sigma,
                                   delta=delta, dfT=dfcont,
                                   muT=mucont, sigmaT=Sigmacont)
    }
  }
  if (cont=="Unif")
  {
    if (model=="Student")
    {
      resultat=SimulTMMcontUniform(nk=nk, dfT=df, muT=mu, sigmaT=Sigma,
                                   delta=delta, lUnif=minU, uUnif=maxU)
    }
    if (model=="Gaussian")
    {
      resultat=SimulGMMcontUniform(nk=nk, muG=mu, sigmaG=Sigma,
                                   delta=delta, lUnif=minU, uUnif=maxU)
    }
  }
  return(resultat)
}
