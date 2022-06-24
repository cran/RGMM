#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

 using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List Weiszfeld_init_rcpp(const arma::mat X, const arma::rowvec init, const arma::rowvec weights, const double epsilon = 1e-08, int nitermax = 100)
{
  // X : n * p  matrix
  // Inputs
  const int n = X.n_rows ;
  const int p = X.n_cols ;
  // Containers and intialisation of the algorithm
  arma::rowvec meanvec(p);
  arma::rowvec poids(n);
  arma::rowvec medvec(p);
  double diffxn, normxm = 10;
  int iter = 0;
  medvec=init;
  meanvec=init;
  /* Boucle Weiszfeld */
    while (iter < nitermax and normxm > epsilon )
    {
      for (int it=0 ; it < n ; it++)
      {
        diffxn = norm(X.row(it)-meanvec);
        if (diffxn > 0) {poids(it) = weights(it)/diffxn;}
        else {poids(it)=0;}
      }
      poids = poids/sum(poids); /* normalisation */
        medvec = poids*X; /* mise a jour de la mediane */
        normxm = norm(medvec-meanvec)/sqrt(double(p));
      meanvec = medvec;
      iter++;
    }
  // Returns ;
  Rcpp::List ret ;
  ret["median"] = medvec ;
  ret["iter"] = iter ;
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
Rcpp::List WeiszfeldCovMat_init_rcpp(const arma::mat  X,  const arma::mat  init_cov, const arma::rowvec median_est, const arma::rowvec weights, const double epsilon = 1e-08, int nitermax = 100)
{
  // X : n * p  matrix
  // Inputs
  const int n = X.n_rows ;
  const int p = X.n_cols ;
  // Containers
  arma::mat Xcent(n,p);
  for (int it=0 ; it < n ; it++)
  {
    Xcent.row(it) = X.row(it)-median_est;
  }
  // arma::mat medinit = arma::trans(Xcent) * Xcent/n;
  arma::mat medinit = init_cov;
  arma::mat  medest(p,p);
  // medest = init_cov;
  arma::rowvec poids(n);
  double diffxn, normxm = 1;
  int iter = 0;


  while (iter < nitermax and normxm > epsilon )
  {
    medest.fill(0.0);
    for (int it=0 ; it < n ; it++)
    {
      diffxn = arma::norm(arma::trans(Xcent.row(it))*Xcent.row(it)-medinit,"fro");
      if (diffxn > 0) {poids(it) = weights(it)/diffxn;}
      else {poids(it)=0;}
    }
    poids = poids/sum(poids);
    for (int it=0 ; it < n ; it++)
    {
      medest += poids(it)*arma::trans(Xcent.row(it))*Xcent.row(it);
    }
  normxm = norm(medest-medinit,"fro")/p;
  medinit = medest;
  iter++;
  }
  // Returns ;
  Rcpp::List ret ;
  ret["median"] = medest ;
  ret["iter"] = iter ;
  ret["poids"] = poids;
  return Rcpp::wrap(ret);
}


// [[Rcpp::export]]
Rcpp::NumericVector Gmedianrowvec_init_rcpp(const arma::mat X, const arma::rowvec init, const arma::rowvec weights, const double gamma = 2, const double alpha = 0.75, const int nstart = 1, const double epsilon = 1e-8)
{
  // X : n * p  matrix
  // Inputs
  const int n = X.n_rows ;
  const int p = X.n_cols ;
  // Containers and intialisation of the algorithm
  arma::rowvec medvec = init ;
  arma::rowvec medrm = init;
  double poids, normxm ;
  // Number of replications of the algorithm
  for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
    // Stochastic gradient algorithms
    for (int it = 1 ; it < n ; it++)
    {
      normxm = arma::norm(X.row(it)-medrm);
      if (normxm > epsilon) {
        poids = weights(it)*sqrt(double(p))*gamma*pow(double(it+1),-alpha)/normxm;
        medrm += poids * (X.row(it)-medrm) ;
      }
      medvec += (medrm-medvec)/(it+1);
    }
  }
  // Returns ;
  return Rcpp::wrap(medvec);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix MedianCovMatRow_init_rcpp(const arma::mat  X,  const arma::mat  init_cov, const arma::rowvec  Gmedian, const arma::rowvec weights, const double gamma = 2, const double alpha = 0.75, const int nstart = 1){
  // X : n * p  matrix
  // Inputs
  const int n = X.n_rows ;
  const int p = X.n_cols ;
  // Containers
  arma::mat medav(p,p);
  medav=init_cov;
  arma::mat  diffmat(p,p), medrm(p,p),diffmed(p,p)  ;
  double nrmrm , poids ;


  // Initialization of the algorithm
  //    diffmed = X.row(0)-Gmedian ;
  //    medav =   arma::mat init_cov;
  medrm = medav ;


  // Number of replications of the algorithm
  for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
    // Stochastic gradient algorithms
    for (int it = 1 ; it < n ; it++)
    {
      diffmed = X.row(it)-Gmedian ;
      diffmat = arma::trans(diffmed)*diffmed;
      diffmat -= medrm ;
      nrmrm = arma::norm(diffmat,"fro") ; // Frobenius norm divided by the dimension
      //           nrmrm = arma::sqrt(arma::sum(arma::square(diffmat)))/p ;
      poids = p *weights(it) * gamma * pow(double(it+1),-alpha) * pow(nrmrm,-1);
      medrm +=   poids*diffmat ;
      medav += (medrm-medav)/(it+1);
    }
  }
  return Rcpp::wrap(medav);
}





