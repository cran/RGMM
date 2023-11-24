\name{RobVar}
\alias{RobVar}
\title{RobVar}

\description{Robust estimate of the variance
}
\usage{
RobVar(X, c=2, alpha=0.75, model='Gaussian', methodMCM='Weiszfeld',
                methodMC='Robbins' , mc_sample_size=1000, init=rep(0, ncol(X)),
                init_cov=diag(ncol(X)),
                epsilon=10^(-8), w=2, df=3, niterMC=50,
                cgrad=2, niterWeisz=50, epsWeisz=10^-8, alphaMedian=0.75, cmedian=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{A numeric matrix of whose rows correspond to observations.}
\item{c}{A positive scalar giving the constant in the stepsequence of the Robbins-Monro or Gradient method if \code{methodMC='RobbinsMC'} or \code{methodMC='GradMC'}. Default is \code{2}.}
\item{alpha}{A scalar between 1/2 and 1 giving the power in the stepsequence for the Robbins-Monro algorithm is \code{methodMC='RobbinsMC'}. Default is \code{0.75}.  }
\item{model}{A string character specifying the model: can be \code{'Gaussian'} (default), \code{'Student'} or \code{'Laplace'}. }
\item{methodMCM}{A string character specifying the method to estimate the Median Covariation Matrix. Can be \code{'Gmedian'} or \code{'Weiszfeld'} (defualt).}
\item{methodMC}{A string character specifying the method to estimate robustly the variance. Can be \code{'Robbins'} (default), \code{'Fix'} or \code{'Grad'}.}
\item{mc_sample_size}{A positive integer giving the number of data simulated for the Monte-Carlo method. Default is \code{1000}.}
\item{init}{A numeric vector giving the initialization for estimating the median.}
\item{init_cov}{A numeric matrix giving an initialization for estimating the Median Covariation Matrix.}
\item{epsilon}{A positive scalar giving a stoping condition for algorithm.}
\item{w}{A positive integer specifying the power for the weighted averaged Robbins-Monro algorithm if \code{methodMC='RobbinsMC'}.}
\item{df}{An integer larger (or equal) than \code{3} specifying the degrees of freedom for the Student law if \code{model='Student'}. See also \code{\link{Gen_MM}}. Default is \code{3}.}
\item{niterMC}{An integer giving the number of iterations for iterative algorithms if the selected method is \code{'Grad'} or \code{'Fix'}. Default is \code{50}.}
\item{cgrad}{A numeric vector with positive values giving the stepsequence of the gradient algorithm for estimating the variance if \code{methodMC='Grad'}. Its length has to be equal to \code{niter}.}
\item{niterWeisz}{A positive integer giving the maximum number of iterations for the Weiszfeld algorithms if \code{methodMCM='Weiszfeld'}. Default is \code{50}.}
\item{epsWeisz}{A stopping factor for the Weiszfeld algorithm.}
\item{alphaMedian}{A scalar betwwen 1/2 and 1  giving the power of the stepsequence of the gradient algorithm for estimating the Median Covariation Matrix if \code{methodMCM='Gmedian'}. Default is \code{0.75}.}
\item{cmedian}{A positive scalar giving the constant in the stepsequence of the gradient algorithm for estimating the Median Covariation Matrix if \code{methodMCM='Gmedian'}. Default is \code{2}.}
}


\value{
An object of class \code{list} with the following outputs:
\item{median}{The median of \code{X}.}
\item{variance}{The robust variance of \code{X}.}
\item{median}{The Median Covariation Matrix of \code{X}.}
}
\references{
Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.

Cardot, H. and Godichon-Baggioni, A. (2017). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Components Analysis.  \emph{Test}, 26(3), 461-480

Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and associated data depth. \emph{Proc. Natl. Acad. Sci. USA}, 97(4):1423-1426.

}
\seealso{
See also \code{\link{RobMM}} and \code{\link{Gen_MM}}.
}
\examples{

n <- 2000
d <- 5
Sigma <-diag(1:d)
mean <- rep(0,d)
X <- mvtnorm::rmvnorm(n,mean,Sigma)
RVar=RobVar(X)
}
\keyword{Mixture Model}


