Robust_GMM=function(X,K=2,ninit=10,nitermax=50,niterEM=50,niterMC=50,
                    mc_sample_size=1000, LogLike=-10^10,arret=10^(-4),epsvp=0,
                    alpha=0.75,c=ncol(X),w=2,epsilon=10^(-3),epsPi=10^-4,initprop=F,epsout=-100,
                    methodMC="RobbinsMC",methodMCM="Weiszfeld")
{
  d=ncol(X)
  n=nrow(X)
  finalcenters=matrix(0,nrow=K,ncol=ncol(X))
  finalcluster=matrix(0,nrow=nrow(X),ncol=K)
  finalvar=matrix(0,nrow=ncol(X),ncol=K*ncol(X))
  finalprop=rep(0,K)
  finalniter=0
  if (length(initprop) >0){
    classif=initprop
    centers=matrix(0,ncol=ncol(X),nrow=K)
    lambda=rep(0,d*K)
    Pi=matrix(0,nrow=n,ncol=K)
    prop=rep(1/K,K)
    for (k in 1:K)
    {
      I=which(classif==k)
      Pi[I,k]=1
    }
    Sigma=c()
    l=0
    for( k in 1:K)
    {
      Sigma=cbind(Sigma,diag(d))
    }
    dist=10^10
    while (l < niterEM && dist > K*ncol(X)*arret)
    {
      l=l+1
      centers2=centers
      #### Mise ? jour des centres et des variances et des poids
      for (k in 1:K)
      {
        if(length(which(Pi[,k] >0))!=0)
        {
          prop[k]=mean(Pi[,k])
          if (methodMCM=="Weiszfeld_init")
          {
            weisz=WeiszfeldCov_init(X,init = centers[k,],init_cov =Sigma[,((k-1)*d+1):(k*d)],weights=(Pi[,k])/sum(Pi[,k]),scores=F,nitermax = nitermax,epsilon = epsilon)
          }
          if (methodMCM=="Weiszfeld")
          {
            weisz=WeiszfeldCov_init(X,weights=(Pi[,k])/sum(Pi[,k]),scores=F,nitermax = nitermax,epsilon = epsilon)
          }
          #         weisz=Gmedian::WeiszfeldCov(X,weights=(Pi[,k]),scores=F,nitermax = nitermax,epsilon = epsilon)
          if (methodMCM=="Gmedian_init")
          {
            weisz=GmedianCov_init(X,init = centers[k,],init_cov =Sigma[,((k-1)*d+1):(k*d)],weights=(Pi[,k]),scores=F)
          }
          if (methodMCM=="Gmedian")
          {
            weisz=GmedianCov_init(X,weights=(Pi[,k]),scores=F)
          }
          if (length(which(is.na(weisz$covmedian)==T))==0){
            if (length(which(is.infinite(weisz$covmedian)==T))==0){
              centers[k,]=weisz$median
              eig=eigen(weisz$covmedian)
              vec=eig$vectors
              vp=eig$values
              lambdak=lambda[((k-1)*d+1):(k*d)]

              ####calculs des vrais valeurs propres
              if (methodMC=="RobbinsMC")
              {
                #        mcm=RobbinsMC(mc_sample_size = mc_sample_size,vp = vp,init=lambdak,
                mcm=RobbinsMC(mc_sample_size = mc_sample_size,vp = vp,
                              alpha = alpha,w = w,epsilon = epsilon,c = c)
              }
              if (methodMC=="GradMC")
              {
                #         mcm=GradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                mcm=GradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,
                           epsilon = epsilon,pas  = c*(1:niterMC)^(-0.5))
              }
              if (methodMC=="FixMC")
              {
                #         mcm=FixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                mcm=FixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,
                          epsilon = epsilon)
              }
              lambdak=apply(cbind(mcm$vp,rep(epsvp,length(mcm$vp))),1,FUN=max)
              lambda[((k-1)*d+1):(k*d)]=lambdak
              Sigma[,((k-1)*d+1):(k*d)]= t(matrix(vec,ncol=d,byrow=T))%*%diag(lambdak)%*%(matrix(vec,ncol=d,byrow=T))
            }
          }
        }
      }
      dist=sum((centers2-centers)^2)
      #### mise a jour des probas
      for (k in 1 : K)
      {
        var=Sigma[,((k-1)*d+1):(k*d)]
        cen=centers[k,]
        Pi[,k] = mvtnorm::dmvnorm(X,mean=cen,sigma = var,log=T)
      }
      Pi=Pi - apply(Pi,1,max)
      for (k in 1:K)
      {
        I=which(Pi[,k]< -100)
        Pi[I,k] = -100
        Pi[,k] = prop[k]*exp(Pi[,k])
      }
      Pi=Pi/rowSums(Pi)
      Pi=Pi+epsPi
      Pi=Pi/rowSums(Pi)
#      if(sum(is.na(Pi))>0)
#      {
#        cat('Pi fout la merde : i=',K,o,l,'\n')
#      }
    }
    LogLikeEM=0
    outliers=c()
    for (k in 1 : K)
    {
      var=Sigma[,((k-1)*d+1):(k*d)]
      cen=centers[k,]
      Pilog =  mvtnorm::dmvnorm(X,mean=cen,sigma = var, log = T)
      outliers=cbind(outliers,Pilog)
      I=which(Pilog < epsout)
      Pilog[I]=epsout
      LogLikeEM=LogLikeEM+prop[k] *  exp(Pilog)
    }
    I = apply(outliers,1,max)
    outliers=which(I< epsout)
    LogLikeEM=sum(log(LogLikeEM))
    if (length(which(is.na(LogLikeEM)==T))==0){
      if(LogLikeEM> LogLike)
      {
        finalcenters=centers
        finalcluster=Pi
        finalvar=Sigma
        LogLike=LogLikeEM
        finalniter=l
        finalprop=prop
        finaloutliers=outliers
      }
    }

  }
  if (ninit > 0){
    for (o in 1:ninit)
    {
      LogLikeEM=0
      prop=rep(1/K,K)
      lambda=rep(0,d*K)
      Pi=matrix(1,nrow=n,ncol=K)
      Sigma=c()
      centers=X[sample(1:n,K),]
      l=0
      dist=10^10
      for( k in 1:K)
      {
        Sigma=cbind(Sigma,diag(d))
      }

      while (l < niterEM && dist > K*ncol(X)*arret)
      {
        l=l+1
        #### mise a jour des probas
        for (k in 1 : K)
        {
          var=Sigma[,((k-1)*d+1):(k*d)]
          cen=centers[k,]
          Pi[,k] = mvtnorm::dmvnorm(X,mean=cen,sigma = var,log=T)
        }
        Pi=Pi - apply(Pi,1,max)
        for (k in 1:K)
        {
          I=which(Pi[,k]< -100)
          Pi[I,k] = -100
          Pi[,k] = prop[k]*exp(Pi[,k])
        }
        Pi=Pi/rowSums(Pi)
        Pi=Pi+epsPi
        Pi=Pi/rowSums(Pi)
#        if(sum(is.na(Pi))>0)
#        {
#          cat('Pi fout la merde : i=',K,o,l,'\n')
#        }
        centers2=centers
        #### Mise ? jour des centres et des variances et des poids
        for (k in 1:K)
        {
          if(length(which(Pi[,k] >0))!=0)
          {
            prop[k]=mean(Pi[,k])
            if (methodMCM=="Weiszfeld_init")
            {
              weisz=WeiszfeldCov_init(X,init = centers[k,],init_cov =Sigma[,((k-1)*d+1):(k*d)],weights=(Pi[,k])/sum(Pi[,k]),scores=F,nitermax = nitermax,epsilon = epsilon)
            }
            if (methodMCM=="Weiszfeld")
            {
              weisz=WeiszfeldCov_init(X,weights=(Pi[,k])/sum(Pi[,k]),scores=F,nitermax = nitermax,epsilon = epsilon)
            }
            #         weisz=Gmedian::WeiszfeldCov(X,weights=(Pi[,k]),scores=F,nitermax = nitermax,epsilon = epsilon)
            if (methodMCM=="Gmedian_init")
            {
              weisz=GmedianCov_init(X,init = centers[k,],init_cov =Sigma[,((k-1)*d+1):(k*d)],weights=(Pi[,k]),scores=F)
            }
            if (methodMCM=="Gmedian")
            {
              weisz=GmedianCov_init(X,weights=(Pi[,k]),scores=F)
            }
            if (length(which(is.na(weisz$covmedian)==T))==0){
              if (length(which(is.infinite(weisz$covmedian)==T))==0){
                centers[k,]=weisz$median
                eig=eigen(weisz$covmedian)
                vec=eig$vectors
                vp=eig$values
                lambdak=lambda[((k-1)*d+1):(k*d)]

                ####calculs des vrais valeurs propres
                if (methodMC=="RobbinsMC")
                {
                  #        mcm=RobbinsMC(mc_sample_size = mc_sample_size,vp = vp,init=lambdak,
                  mcm=RobbinsMC(mc_sample_size = mc_sample_size,vp = vp,
                                alpha = alpha,w = w,epsilon = epsilon,c = c)
                }
                if (methodMC=="GradMC")
                {
                  #         mcm=GradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                  mcm=GradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,
                             epsilon = epsilon,pas  = c*(1:niterMC)^(-0.5))
                }
                if (methodMC=="FixMC")
                {
                  #         mcm=FixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                  mcm=FixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,
                            epsilon = epsilon)
                }
                lambdak=apply(cbind(mcm$vp,rep(epsvp,length(mcm$vp))),1,FUN=max)
                lambda[((k-1)*d+1):(k*d)]=lambdak
                Sigma[,((k-1)*d+1):(k*d)]= t(matrix(vec,ncol=d,byrow=T))%*%diag(lambdak)%*%(matrix(vec,ncol=d,byrow=T))
              }
            }
          }
        }
        dist=sum((centers2-centers)^2)
      }
      LogLikeEM=0
      outliers=c()
      for (k in 1 : K)
      {
        var=Sigma[,((k-1)*d+1):(k*d)]
        cen=centers[k,]
        Pilog =  mvtnorm::dmvnorm(X,mean=cen,sigma = var, log = T)
        outliers=cbind(outliers,Pilog)
        I=which(Pilog < epsout)
        Pilog[I]=epsout
        LogLikeEM=LogLikeEM+prop[k] *  exp(Pilog)
      }
      I = apply(outliers,1,max)
      outliers=which(I< epsout)
      LogLikeEM=sum(log(LogLikeEM))
      if (length(which(is.na(LogLikeEM)==T))==0){
        if(LogLikeEM> LogLike)
        {
          finalcenters=centers
          finalcluster=Pi
          finalvar=Sigma
          LogLike=LogLikeEM
          finalniter=l
          finalprop=prop
          finaloutliers=outliers
        }
      }

    }
  }
  return(list(centers=finalcenters,Sigma=finalvar,Loglike=LogLike, Pi=finalcluster,niter=finalniter,initEM=initprop,prop=finalprop,outliers=finaloutliers))
}



RGMM=function(X,nclust=2:5,ninit=10,nitermax=50,niterEM=50,niterMC=50,epsvp=epsvp,
              mc_sample_size=1000, LogLike=-10^10,init=T,epsPi=10^-4,epsout=-100,
              alpha=0.75,c=ncol(X),w=2,epsilon=10^(-8),
              methodMC="RobbinsMC", par=T,methodMCM="Weiszfeld")
{
  initprop=F
  if (init==T)
  {
    clas=mclust::hcVVV(data=X)
  }
  if (length(nclust)==1)
  {
    K=nclust
    if (init==T)
    {
      initprop=  mclust::hclass(clas,nclust)
    }
    if (init=='Mclust')
    {
      initprop=  mclust::hclass(clas,nclust)
    }
    if (init=='genie')
    {
      initprop=  genieclust::genie(X,k=nclust)
    }
    resultat=Robust_GMM(X,K=nclust,ninit=ninit,nitermax=nitermax,niterEM=niterEM,epsPi=epsPi,epsout=epsout,epsvp=epsvp,
                        niterMC=niterMC,mc_sample_size=mc_sample_size, LogLike=LogLike,initprop=initprop,
                        alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,methodMCM=methodMCM)
    a=resultat$Pi*log(resultat$Pi)
    bestresult=resultat
    I=which(is.na(a))
    a[I]=0
    ICL= resultat$Loglike - log(nrow(X)/2)*(nclust*ncol(X) + nclust*ncol(X)*(ncol(X)+1)/2) + sum(a)
  }
  if (length(nclust)>1)
  {
    if (par ==T)
    {
      numCores = min(parallel::detectCores()-2,length(nclust))
      cl = parallel::makeCluster(numCores  )
      #  clusterExport(cl, varlist = c('Robust_GMM','WeiszfeldCov_init','Weiszfeld_init_rcpp',
      #                               'GradMC','FixMC','RobbinsMC','Gmedian_init','GmedianCov_init','Gmedianrowvec_init_rcpp',
      #                               'MedianCovMatW_init_rcpp','Weiszfeld_init','MedianCovMatRow_init_rcpp',
      #                               'X','nclust' ,'ninit' ,'nitermax' ,'niterEM','niterMC',
      #                               'mc_sample_size', 'LogLike' ,
      #                               'alpha' ,'c' ,'w' ,'epsilon' ,
      #                               'methodMC'), envir = environment())
      doParallel::registerDoParallel(cl)
      #  iterations <- nclust[length(nclust)]
      #  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
      #  progress <- function(n) setTxtProgressBar(pb, n)
      #  opts <- list(progress = progress)
      resultat=foreach::foreach(K=nclust,.multicombine = TRUE, .inorder = T,.packages = c("dplyr","mvtnorm","Gmedian","Rcpp","RcppArmadillo",'mixtools',"mclust"),.combine='list')  %dopar%
        {
          if (init==T)
          {
            initprop=  mclust::hclass(clas,K)
          }
          if (init=='genie')
          {
            initprop=  genieclust::genie(X,k=K)
          }
          #         cat('Running for : K=',K,'\n')
#          cat("Running K =", min(nclust), "...\n")
          resultatk=Robust_GMM(X,K=K,ninit=ninit,nitermax=nitermax,niterEM=niterEM,epsout=epsout,epsvp=epsvp,
                               niterMC=niterMC,mc_sample_size=mc_sample_size, LogLike=LogLike,initprop=initprop,epsPi=epsPi,
                               alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,methodMCM=methodMCM)
          return(resultatk)
        }

      parallel::stopCluster(cl)
    }

    if (par ==F)
    {
      resultat=foreach::foreach(K=nclust,.multicombine = TRUE, .inorder = TRUE,.packages = c("dplyr","mvtnorm","Gmedian","Rcpp","RcppArmadillo",'mixtools'),.combine='list')  %do%
        {
          if (init==T)
          {
            initprop=  mclust::hclass(clas,nclust)
          }
          if (init=='genie')
          {
            initprop=  genieclust::genie(X,k=K)
          }
#          cat('Running for : K=',K,'\n')
          resultatk=Robust_GMM(X,K=K,ninit=ninit,nitermax=nitermax,niterEM=niterEM,epsPi=epsPi,epsout=epsout,epsvp=epsvp,
                               niterMC=niterMC,mc_sample_size=mc_sample_size, LogLike=LogLike,initprop=initprop,
                               alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,methodMCM=methodMCM)
          return(resultatk)
        }

    }
    ICL=c()
    for (i in 1:length(nclust))
    {
      a=resultat[[i]]$Pi*log(resultat[[i]]$Pi)
      I=which(is.na(a))
      a[I]=0
      ICL=c(ICL,resultat[[i]]$Loglike - log(nrow(X)/2)*(nclust[i]*ncol(X) + nclust[i]*ncol(X)*(ncol(X)+1)/2) + sum(a) )
    }
    k=which.max(ICL)
    bestresult=resultat[[k]]
  }
  return(list(allresults=resultat,bestresult=bestresult,ICL=ICL))
}

