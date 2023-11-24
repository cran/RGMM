Robust_TMM=function(X,K=2,ninit=10,nitermax=50,niterEM=50,niterMC=50,df=3,scale='none',
                    mc_sample_size=1000, LogLike=-10^10,arret=10^(-4),epsvp=10^-4,
                    alpha=0.75,c=ncol(X),w=2,epsilon=10^(-3),epsPi=10^-4,initprop=F,epsout=-100,
                    methodMC="RobbinsMC",methodMCM="Weiszfeld")
{
  if (scale=='robust')
  {
    X=DescTools::RobScale(X)
    scales=attr(X,"scaled:scale")
    init_loc=attr(X,"scaled:center")
    c=1
  }
  d=ncol(X)
  n=nrow(X)
  finalcenters=matrix(0,nrow=K,ncol=ncol(X))
  finalcluster=matrix(0,nrow=nrow(X),ncol=K)
  finaldensities=matrix(0,nrow=nrow(X),ncol=K)
  densities=matrix(0,nrow=nrow(X),ncol=K)
  finalvar=matrix(0,nrow=ncol(X),ncol=K*ncol(X))
  finalprop=rep(0,K)
  finalniter=0
  finaloutliers=c()
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
      if (K > 1){
        centers2=centers
      }
      if (K == 1){
        centers2=as.vector(centers)
      }
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
              if (K>1){
                centers[k,]=weisz$median
              }
              if (K==1){
                centers=weisz$median
              }
              eig=eigen(weisz$covmedian)
              vec=eig$vectors
              vp=eig$values
              lambdak=lambda[((k-1)*d+1):(k*d)]

              ####calculs des vrais valeurs propres
              if (methodMC=="RobbinsMC")
              {
                #        mcm=RobbinsMC(mc_sample_size = mc_sample_size,vp = vp,init=lambdak,
                mcm=TRobbinsMC(mc_sample_size = mc_sample_size,vp = vp,df=df,
                               alpha = alpha,w = w,epsilon = epsilon,c = c)
              }
              if (methodMC=="GradMC")
              {
                #         mcm=GradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                mcm=TGradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,df=df,
                            epsilon = epsilon,pas  = c*(1:niterMC)^(-0.5))
              }
              if (methodMC=="FixMC")
              {
                #         mcm=FixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                mcm=TFixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,df=df,
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
        if (K>1){
          cen=(centers[k,])
        }
        if (K==1){
          cen=c(centers)
        }
        Pi[,k] = mvtnorm::dmvt(X,delta=cen,sigma = (df-2)/df*var,df=df,log=T)
      }
      Phiik=Pi
      Pi=Pi - apply(Pi,1,max)
      for (k in 1:K)
      {
        I=which(Pi[,k]< epsout)
        Pi[I,k] = (epsout)
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
      I=which(Phiik[,k] < epsout)
      Phiik[I,k]=epsout
    }
    LogLikeEM=sum(Pi*Phiik)+ sum(Pi*log(prop))- sum(Pi*log(Pi))
    I = apply(Phiik,1,max)
    outliers=which(I<= epsout)
    if (length(which(is.na(LogLikeEM)==T))==0){
      if(LogLikeEM> LogLike)
      {
        finalcenters=centers
        finalcluster=Pi
        finalvar=Sigma
        LogLike=LogLikeEM
        finaldensities=Phiik
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
      centers=as.matrix(X[sample(1:n,K),],ncol=d)
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
          if (K>1){
            cen=(centers[k,])
          }
          if (K==1){
            cen=c(centers)
          }
          Pi[,k] = mvtnorm::dmvt(X,delta=cen,sigma = (df-2)/df*var,df=df,log=T)
        }
        Phiik=Pi
        Pi=Pi - apply(Pi,1,max)
        for (k in 1:K)
        {
          I=which(Pi[,k]< (epsout))
          Pi[I,k] = epsout
          Pi[,k] = prop[k]*exp(Pi[,k])
        }
        Pi=Pi/rowSums(Pi)
        Pi=Pi+epsPi
        Pi=Pi/rowSums(Pi)
#        if(sum(is.na(Pi))>0)
#        {
#          cat('Pi fout la merde : i=',K,o,l,'\n')
#        }
        if (K > 1){
          centers2=centers
        }
        if (K == 1){
          centers2=as.vector(centers)
        }
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
                if (K>1){
                  centers[k,]=weisz$median
                }
                if (K==1){
                  centers=weisz$median
                }
                eig=eigen(weisz$covmedian)
                vec=eig$vectors
                vp=eig$values
                lambdak=lambda[((k-1)*d+1):(k*d)]

                ####calculs des vrais valeurs propres
                if (methodMC=="RobbinsMC")
                {
                  #        mcm=RobbinsMC(mc_sample_size = mc_sample_size,vp = vp,init=lambdak,
                  mcm=TRobbinsMC(mc_sample_size = mc_sample_size,vp = vp,df=df,
                                 alpha = alpha,w = w,epsilon = epsilon,c = c)
                }
                if (methodMC=="GradMC")
                {
                  #         mcm=GradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                  mcm=TGradMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,df=df,
                              epsilon = epsilon,pas  = c*(1:niterMC)^(-0.5))
                }
                if (methodMC=="FixMC")
                {
                  #         mcm=FixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,init=lambdak,
                  mcm=TFixMC(mc_sample_size = mc_sample_size,vp = vp,niter = niterMC,df=df,
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
        I=which(Phiik[,k] < epsout)
        Phiik[I,k]=epsout
      }
      LogLikeEM=sum(Pi*Phiik)+ sum(Pi*log(prop))- sum(Pi*log(Pi))
      I = apply(Phiik,1,max)
      outliers=which(I<= epsout)
      if (length(which(is.na(LogLikeEM)==T))==0){
        if(LogLikeEM> LogLike)
        {
          finalcenters=centers
          finalcluster=Pi
          finalvar=Sigma
          LogLike=LogLikeEM
          finalniter=l
          finaldensities=Phiik
          finalprop=prop
          finaloutliers=outliers
        }
      }

    }
  }
  if (scale=='robust')
  {
    for (k in 1:K){
      finalvar[,((k-1)*d+1):(k*d)]=diag(scales)%*%finalvar[,((k-1)*d+1):(k*d)]%*%diag(scales)
      finalcenters[k,]=finalcenters[k,]%*%diag(scales)+ init_loc
    }
  }
  return(list(centers=finalcenters,Sigma=finalvar,Loglike=LogLike, Pi=finalcluster,niter=finalniter,initEM=initprop,prop=finalprop,outliers=finaloutliers))
}


RTMM=function(X,nclust=2:5,ninit=10,nitermax=50,niterEM=50,niterMC=50,df=3,epsvp=10^-4,
              mc_sample_size=1000, LogLike=-10^10,init=T,epsPi=10^-4,epsout=-100,
              alpha=0.75,c=ncol(X),w=2,epsilon=10^(-8),criterion='ICL',scale='none',
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
    Kopt=nclust
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
    resultat=Robust_TMM(X,K=nclust,df=df,ninit=ninit,nitermax=nitermax,niterEM=niterEM,epsPi=epsPi,epsout=epsout,epsvp=epsvp,
                        niterMC=niterMC,mc_sample_size=mc_sample_size, LogLike=LogLike,initprop=initprop,
                        alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,methodMCM=methodMCM,scale=scale)
    a=resultat$Pi*log(resultat$Pi)
    bestresult=resultat
    I=which(is.na(a))
    a[I]=0
    ICL=bestresult$Loglike - 0.5*log(nrow(X))*(nclust-1+  nclust*ncol(X) + nclust*ncol(X)*(ncol(X)+1)/2) - sum(resultat$Pi*log(resultat$Pi))
    BIC=bestresult$Loglike - 0.5*log(nrow(X))*(nclust-1+ nclust*ncol(X) + nclust*ncol(X)*(ncol(X)+1)/2)
  }
  if (length(nclust)>1)
  {
    if (par ==T)
    {
      numCores = min(detectCores()-2,length(nclust))
      cl = makeCluster(numCores  )
      #  clusterExport(cl, varlist = c('Robust_GMM','WeiszfeldCov_init','Weiszfeld_init_rcpp',
      #                               'GradMC','FixMC','RobbinsMC','Gmedian_init','GmedianCov_init','Gmedianrowvec_init_rcpp',
      #                               'MedianCovMatW_init_rcpp','Weiszfeld_init','MedianCovMatRow_init_rcpp',
      #                               'X','nclust' ,'ninit' ,'nitermax' ,'niterEM','niterMC',
      #                               'mc_sample_size', 'LogLike' ,
      #                               'alpha' ,'c' ,'w' ,'epsilon' ,
      #                               'methodMC'), envir = environment())
      registerDoParallel(cl)
      #  iterations <- nclust[length(nclust)]
      #  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
      #  progress <- function(n) setTxtProgressBar(pb, n)
      #  opts <- list(progress = progress)
      resultat=foreach(K=nclust,.multicombine = TRUE, .inorder = T,.packages = c("dplyr","mvtnorm","Gmedian","Rcpp","RcppArmadillo",'mixtools',"mclust"),.combine='list')  %dopar%
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
          resultatk=Robust_TMM(X,K=K,df=df,ninit=ninit,nitermax=nitermax,niterEM=niterEM,epsout=epsout,epsvp=epsvp,
                               niterMC=niterMC,mc_sample_size=mc_sample_size, LogLike=LogLike,initprop=initprop,epsPi=epsPi,
                               alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,methodMCM=methodMCM,scale=scale)
          return(resultatk)
        }

      stopCluster(cl)
    }

    if (par ==F)
    {
      resultat=foreach(K=nclust,.multicombine = TRUE, .inorder = T,.packages = c("dplyr","mvtnorm","Gmedian","Rcpp","RcppArmadillo",'mixtools'),.combine='list')  %do%
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
          resultatk=Robust_TMM(X,K=K,df=df,ninit=ninit,nitermax=nitermax,niterEM=niterEM,epsPi=epsPi,epsout=epsout,epsvp=epsvp,
                               niterMC=niterMC,mc_sample_size=mc_sample_size, LogLike=LogLike,initprop=initprop,
                               alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,methodMCM=methodMCM,scale=scale)
          return(resultatk)
        }

    }
    ICL=c()
    BIC=c()
    for (i in 1:length(nclust))
    {
      ICL=c(ICL,resultat[[i]]$Loglike - 0.5*log(nrow(X))*(nclust[i]-1+  nclust[i]*ncol(X) + nclust[i]*ncol(X)*(ncol(X)+1)/2) - sum(resultat[[i]]$Pi*log(resultat[[i]]$Pi)))
      BIC=c(BIC,resultat[[i]]$Loglike - 0.5*log(nrow(X))*(nclust[i]-1+ nclust[i]*ncol(X) + nclust[i]*ncol(X)*(ncol(X)+1)/2) )
    }
    if (criterion=='ICL'){
    k=which.max(ICL)
    Kopt=nclust[k]
    }
    if (criterion=='BIC'){
      k=which.max(BIC)
      Kopt=nclust[k]
    }
    bestresult=resultat[[k]]
  }
  return(list(allresults=resultat,bestresult=bestresult,ICL=ICL,BIC=BIC,data=X,nclust=nclust,Kopt=Kopt))
}

