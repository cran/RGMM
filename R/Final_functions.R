utils::globalVariables(c( "x", "y", "coul", "couleur",'BIC','ICL'))


RobMM=function(X,nclust=2:5,model="Gaussian",ninit=10,nitermax=50,
               niterEM=50,niterMC=50,df=3,epsvp=10^(-4),
               mc_sample_size=1000, LogLike=-Inf,init='genie',
               epsPi=10^-4,epsout=-20,
               alpha=0.75,c=ncol(X),w=2,epsilon=10^(-8),criterion='ICL',
               methodMC="RobbinsMC", par=TRUE,methodMCM="Weiszfeld")
{
  if (model=='Gaussian')
  {
    resultat=RGMM(X,nclust=nclust,ninit=ninit,nitermax=nitermax,niterEM=niterEM,
                  niterMC=niterMC,epsvp=epsvp,mc_sample_size=mc_sample_size,
                  LogLike=LogLike,init=init,epsPi=epsPi,epsout=epsout,criterion=criterion,
                  alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,
                  methodMCM=methodMCM,par=par)
  }
  if (model=='Student')
  {
    resultat=RTMM(X,nclust=nclust,ninit=ninit,nitermax=nitermax,niterEM=niterEM,
                  niterMC=niterMC,epsvp=epsvp,mc_sample_size=mc_sample_size,df=df,
                  LogLike=LogLike,init=init,epsPi=epsPi,epsout=epsout,criterion=criterion,
                  alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,
                  methodMCM=methodMCM,par=par)
  }
  if (model=='Laplace')
  {
    resultat=RLMM(X,nclust=nclust,ninit=ninit,nitermax=nitermax,niterEM=niterEM,
                  niterMC=niterMC,epsvp=epsvp,mc_sample_size=mc_sample_size,
                  LogLike=LogLike,init=init,epsPi=epsPi,epsout=epsout,criterion=criterion,
                  alpha=alpha,c=c,w=w,epsilon=epsilon,methodMC=methodMC,
                  methodMCM=methodMCM,par=par)
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

RMMplot=function(a,outliers=TRUE,graph=c('Two_Dim','Two_Dim_Uncertainty','ICL','BIC','Profiles','Uncertainty'),bestresult=TRUE,K=FALSE)
{
  if (bestresult==TRUE){
  X=a$data
  d=ncol(X)
  nclust=a$nclust
  Ksel=a$Kopt
  if (d==1)
  {
    print('d must be larger than 1')
  }
  if (d >1)
  {
    cluster=apply(a$bestresult$Pi,1,which.max)
    if (d > 2){
      Gmed=WeiszfeldCov_init(X,scores = 2)
      vec=Gmed$vectors
    }
    if (d==2)
    {
      Gmed=WeiszfeldCov_init(X,scores = F)
      eig=eigen(Gmed$covmedian)
      vec=eig$vectors
    }
    med0=Gmed$median
    if (length(which(graph=='Two_Dim'))>0){
      med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
      Xtrans=matrix(0,nrow=nrow(X),ncol=2)
      for(i in 1:nrow(X))
      {
        Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
        Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
      }
      if (outliers==T){
        if (length(a$bestresult$outliers)>0){
          I=a$bestresult$outliers
          dataplot=data.frame(x=Xtrans[-I,1],y=Xtrans[-I,2],K=as.character(cluster[-I]))
        }
        if (length(a$bestresult$outliers)==0){
          I=a$bestresult$outliers
          dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster))
        }
      }
      if (outliers==F)
      {
        dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster))
      }
      plo=  ggplot(data=dataplot) +
        scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
        scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
        geom_point(aes(x = x, y = y,color = K))
      print(plo)
    }
    if (length(which(graph=='Two_Dim_Uncertainty'))>0){
      Pi=a$bestresult$Pi
      Pi=Pi*log(Pi)
      Uncertainty=-rowSums(Pi)/max(abs(-rowSums(Pi)))
      med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
      Xtrans=matrix(0,nrow=nrow(X),ncol=2)
      for(i in 1:nrow(X))
      {
        Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
        Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
      }
      if (outliers==T){
        if (length(a$bestresult$outliers)>0){
          I=a$bestresult$outliers
          dataplot=data.frame(x=Xtrans[-I,1],y=Xtrans[-I,2],K=as.character(cluster[-I]),Uncertainty=Uncertainty[-I])
        }
        if (length(a$bestresult$outliers)==0){
          I=a$bestresult$outliers
          dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster),Uncertainty=Uncertainty)
        }
      }
      if (outliers==F)
      {
        dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster),Uncertainty=Uncertainty)
      }
      plo2=  ggplot(data=dataplot) +
        scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
        scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
        scale_size_continuous(range=c(1,4))+
        geom_point(aes(x = x, y = y,color = K,size=  Uncertainty))
      print(plo2)
    }
    if (length(which(graph=='ICL'))>0){
      data_icl=data.frame(K=a$nclust,ICL=a$ICL)
      plo_icl=ggplot(data=data_icl)+
        geom_line(aes(x=K,y=ICL))
      print(plo_icl)
    }
    if (length(which(graph=='BIC'))>0){
      data_bic=data.frame(K=a$nclust,BIC=a$BIC)
      plo_bic=ggplot(data=data_bic)+
        geom_line(aes(x=K,y=BIC))
      print(plo_bic)
    }
    if (length(which(graph=='Profiles'))>0){
      if (outliers==T)
      {
        if (length(a$bestresult$outliers)>0){
          Xt=X[-(a$bestresult$outliers),]
          cluster2=cluster[-(a$bestresult$outliers)]
        }
        if (length(a$bestresult$outliers)==0)
        {
          Xt=X
          cluster2=cluster
        }
      }
      if (outliers==F)
      {
        Xt=X
        cluster2=cluster
      }
      if (d>2){
        dataclust=data.frame(x=1:d,
                             y=t(Xt))
        datacen=data.frame(xt=1:d,t(a$bestresult$centers))
        ind=c()
        Kclu=c()
        for (o in 1:nrow(Xt))
        {
          ind=c(ind,rep(o,d))
          Kclu=c(Kclu,rep(as.character(cluster2[o]),d))
        }
        Kcen=c()
        indK=c()
        col=c()
        for (o in 1:Ksel)
        {
          I=length(which(Kclu==as.character(o)))
          Kcen=c(Kcen,rep(as.character(o),d))
          indK=c(indK,rep((o+max(ind)),d))
        }
        data_new <- melt(dataclust, id = c("x"))
        data_cen <- melt(datacen, id = c("xt"))
        datanew=data.frame( x=c(data_new$x,data_cen$xt),y=c(data_new$value,data_cen$value),
                            couleur=c(ind,indK),Kclu=c(Kclu,Kcen),
                            coul=c(rep('grey',length(data_new$x)),rep('red',length(data_cen$xt))))

        prof_plot=  ggplot(datanew,aes(x=x , y=y,color=coul))+
          scale_color_manual(values=c('grey','red'))+
          scale_x_continuous('') +
          scale_y_continuous('Profiles',limits=c(min(Xt),max(Xt))) +
          geom_line(mapping=aes(x=x , y=y , group=couleur),show.legend = F,alpha=0.4)+
          #  geom_line(aes(color=(coul)))+
          facet_wrap(~as.factor(Kclu))
        print(prof_plot)
      }
      if (d==2){
        dataclust=data.frame(x=1:d,
                             y=t(Xt))
        Kclu=c()
        for (o in 1:nrow(Xt))
        {
          Kclu=c(Kclu,as.character(cluster2[o]))
        }
        Kcen=c()
        for (o in 1:Ksel)
        {
          Kcen=c(Kcen,as.character(o))
        }
        data_clust=data.frame(x=c(Xt[,1],(a$bestresult$centers)[,1]),y=c(Xt[,2],(a$bestresult$centers)[,2]),
                              K=c(Kclu,Kcen),coul=c(rep('grey',nrow(Xt)),rep('red',Ksel)))
        prof_plot=  ggplot(data_clust,aes(x=x , y=y,color=coul))+
          scale_color_manual(values=c('grey','red'))+
          scale_x_continuous('') +
          scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
          geom_point(mapping=aes(x=x , y=y),show.legend = F,alpha=0.4)+
          #  geom_line(aes(color=(coul)))+
          facet_wrap(~as.factor(K))
        print(prof_plot)
      }
    }
    if (length(which(graph=='Uncertainty'))>0)
    {
      Pi=a$bestresult$Pi
      Pi=Pi*log(Pi)
      Uncertainty=-rowSums(Pi)/max(abs(-rowSums(Pi)))
      ord=order(Uncertainty)
      Uncertainty=Uncertainty[ord]
      data_unc=data.frame(X=1:nrow(X),Uncertainty=Uncertainty)
      plot_unc=ggplot(data=data_unc,aes(x=X,y=Uncertainty))+
        geom_bar(stat='identity')
      print(plot_unc)
    }
  }
  }
  if (K!=FALSE){
    X=a$data
    d=ncol(X)
    nclust=a$nclust
    if (length(which(nclust==K))==0)
    {
      print('K does not belong to nclust')
    }
    if (length(which(nclust==K))>0){
    Ksel=K
    Knclust=which(nclust==Ksel)
    if (d==1)
    {
      print('d must be larger than 1')
    }
    if (d >1)
    {
      cluster=apply(a$allresults[[Knclust]]$Pi,1,which.max)
      if (d > 2){
        Gmed=WeiszfeldCov_init(X,scores = 2)
        vec=Gmed$vectors
      }
      if (d==2)
      {
        Gmed=WeiszfeldCov_init(X,scores = F)
        eig=eigen(Gmed$covmedian)
        vec=eig$vectors
      }
      med0=Gmed$median
      if (length(which(graph=='Two_Dim'))>0){
        med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
        Xtrans=matrix(0,nrow=nrow(X),ncol=2)
        for(i in 1:nrow(X))
        {
          Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
          Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
        }
        if (outliers==T){
          if (length(a$allresults[[Knclust]]$outliers)>0){
            I=a$allresults[[Knclust]]$outliers
            dataplot=data.frame(x=Xtrans[-I,1],y=Xtrans[-I,2],K=as.character(cluster[-I]))
          }
          if (length(a$allresults[[Knclust]]$outliers)==0){
            I=a$allresults[[Knclust]]$outliers
            dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster))
          }
        }
        if (outliers==F)
        {
          dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster))
        }
        plo=  ggplot(data=dataplot) +
          scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
          scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
          geom_point(aes(x = x, y = y,color = K))
        print(plo)
      }
      if (length(which(graph=='Two_Dim_Uncertainty'))>0){
        Pi=a$allresults[[Knclust]]$Pi
        Pi=Pi*log(Pi)
        Uncertainty=-rowSums(Pi)/max(abs(-rowSums(Pi)))
        med=c(sum(med0*vec[,1]),sum(med0*vec[,2]))
        Xtrans=matrix(0,nrow=nrow(X),ncol=2)
        for(i in 1:nrow(X))
        {
          Xtrans[i,1]=sum(X[i,]*vec[,1])-med[1]
          Xtrans[i,2]=sum(X[i,]*vec[,2])-med[2]
        }
        if (outliers==T){
          if (length(a$allresults[[Knclust]]$outliers)>0){
            I=a$allresults[[Knclust]]$outliers
            dataplot=data.frame(x=Xtrans[-I,1],y=Xtrans[-I,2],K=as.character(cluster[-I]),Uncertainty=Uncertainty[-I])
          }
          if (length(a$allresults[[Knclust]]$outliers)==0){
            I=a$allresults[[Knclust]]$outliers
            dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster),Uncertainty=Uncertainty)
          }
        }
        if (outliers==F)
        {
          dataplot=data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster),Uncertainty=Uncertainty)
        }
        plo2=  ggplot(data=dataplot) +
          scale_x_continuous('First Principal Componant',  limits = c(min(dataplot$x),max(dataplot$x))) +
          scale_y_continuous('Second Principal Componant',limits = c(min(dataplot$y),max(dataplot$y))) +
          scale_size_continuous(range=c(1,4))+
          geom_point(aes(x = x, y = y,color = K,size=  Uncertainty))
        print(plo2)
      }
      if (length(which(graph=='ICL'))>0){
        data_icl=data.frame(K=a$nclust,ICL=a$ICL)
        plo_icl=ggplot(data=data_icl)+
          geom_line(aes(x=K,y=ICL))
        print(plo_icl)
      }
      if (length(which(graph=='BIC'))>0){
        data_bic=data.frame(K=a$nclust,BIC=a$BIC)
        plo_bic=ggplot(data=data_bic)+
          geom_line(aes(x=K,y=BIC))
        print(plo_bic)
      }
      if (length(which(graph=='Profiles'))>0){
        if (outliers==T)
        {
          if (length(a$allresults[[Knclust]]$outliers)>0){
            Xt=X[-(a$allresults[[Knclust]]$outliers),]
            cluster2=cluster[-(a$allresults[[Knclust]]$outliers)]
          }
          if (length(a$allresults[[Knclust]]$outliers)==0)
          {
            Xt=X
            cluster2=cluster
          }
        }
        if (outliers==F)
        {
          Xt=X
          cluster2=cluster
        }
        if (d>2){
          dataclust=data.frame(x=1:d,
                               y=t(Xt))
          datacen=data.frame(xt=1:d,t(a$allresults[[Knclust]]$centers))
          ind=c()
          Kclu=c()
          for (o in 1:nrow(Xt))
          {
            ind=c(ind,rep(o,d))
            Kclu=c(Kclu,rep(as.character(cluster2[o]),d))
          }
          Kcen=c()
          indK=c()
          col=c()
          for (o in 1:Ksel)
          {
            I=length(which(Kclu==as.character(o)))
            Kcen=c(Kcen,rep(as.character(o),d))
            indK=c(indK,rep((o+max(ind)),d))
          }
          data_new <- melt(dataclust, id = c("x"))
          data_cen <- melt(datacen, id = c("xt"))
          datanew=data.frame( x=c(data_new$x,data_cen$xt),y=c(data_new$value,data_cen$value),
                              couleur=c(ind,indK),Kclu=c(Kclu,Kcen),
                              coul=c(rep('grey',length(data_new$x)),rep('red',length(data_cen$xt))))

          prof_plot=  ggplot(datanew,aes(x=x , y=y,color=coul))+
            scale_color_manual(values=c('grey','red'))+
            scale_x_continuous('') +
            scale_y_continuous('Profiles',limits=c(min(Xt),max(Xt))) +
            geom_line(mapping=aes(x=x , y=y , group=couleur),show.legend = F,alpha=0.4)+
            #  geom_line(aes(color=(coul)))+
            facet_wrap(~as.factor(Kclu))
          print(prof_plot)
        }
        if (d==2){
          dataclust=data.frame(x=1:d,
                               y=t(Xt))
          Kclu=c()
          for (o in 1:nrow(Xt))
          {
            Kclu=c(Kclu,as.character(cluster2[o]))
          }
          Kcen=c()
          for (o in 1:Ksel)
          {
            Kcen=c(Kcen,as.character(o))
          }
          data_clust=data.frame(x=c(Xt[,1],(a$allresults[[Knclust]]$centers)[,1]),y=c(Xt[,2],(a$allresults[[Knclust]]$centers)[,2]),
                                K=c(Kclu,Kcen),coul=c(rep('grey',nrow(Xt)),rep('red',Ksel)))
          prof_plot=  ggplot(data_clust,aes(x=x , y=y,color=coul))+
            scale_color_manual(values=c('grey','red'))+
            scale_x_continuous('') +
            scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
            geom_point(mapping=aes(x=x , y=y),show.legend = F,alpha=0.4)+
            #  geom_line(aes(color=(coul)))+
            facet_wrap(~as.factor(K))
          print(prof_plot)
        }
      }
      if (length(which(graph=='Uncertainty'))>0)
      {
        Pi=a$allresults[[Knclust]]$Pi
        Pi=Pi*log(Pi)
        Uncertainty=-rowSums(Pi)/max(abs(-rowSums(Pi)))
        ord=order(Uncertainty)
        Uncertainty=Uncertainty[ord]
        data_unc=data.frame(X=1:nrow(X),Uncertainty=Uncertainty)
        plot_unc=ggplot(data=data_unc,aes(x=X,y=Uncertainty))+
          geom_bar(stat='identity')
        print(plot_unc)
      }
    }
    }
  }
}
