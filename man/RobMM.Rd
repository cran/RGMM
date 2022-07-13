\name{RobMM}
\alias{RobMM}
\title{RobMM}

\description{Robust Mixture Model
}
\usage{
RobMM(X, nclust=2:5, model="Gaussian", ninit=10,
               nitermax=50, niterEM=50, niterMC=50, df=3,
               epsvp=10^(-4), mc_sample_size=1000, LogLike=-Inf,
               init='genie', epsPi=10^-4, epsout=-20,
               alpha=0.75, c=ncol(X), w=2, epsilon=10^(-8),
               criterion='ICL',methodMC="RobbinsMC", par=TRUE,
               methodMCM="Weiszfeld")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{A matrix giving the data.}
\item{nclust}{A vector of positive integers giving the possible number of clusters.}
\item{model}{The mixture model. Can be \code{'Gaussian'} (by default), \code{'Student'} and \code{'Laplace'}.}
\item{ninit}{The number of random initisalizations. Befault is \code{10}.}
\item{nitermax}{The number of iterations for the Weiszfeld algorithm if \code{MethodMCM= 'Weiszfeld'}.}
\item{niterEM}{The number of iterations for the EM algorithm.}
\item{niterMC}{The number of iterations for estimating robustly the variance of each class if \code{methodMC='FixMC'} or \code{methodMC='GradMC'}.}
\item{df}{The degrees of freedom for the Student law if \code{model='Student'}.}
\item{epsvp}{The minimum values the estimates of the eigenvalues of the Median Covariation Matrix can take. Default is \code{10^-4}.}
\item{mc_sample_size}{The number of data generated for the Monte-Carlo method for estimating robustly the variance.}
\item{LogLike}{The initial loglikelihood to "beat". Defulat is \code{-Inf}.}
\item{init}{Can be \code{F} if no non random initialization of the algorithm is done, \code{'genie'} if the algorithm is initialized with the help of the function \code{'genie'} of the package \code{genieclust} or \code{'Mclust'} if the initialization is done with the function \code{hclass} of the package \code{Mclust}.}
\item{epsPi}{A scalar to ensure the estimates of the probabilities of belonging to a class or uniformly lower bounded by a positive constant.}
\item{epsout}{If the probability of belonging of a data to a class is smaller than \code{exp(epsout)}, this probbility is replaced by \code{exp(epsout)} for calculating the logLikelihood. If the probability is too weak for each class, the data is considered as an outlier. Defautl is \code{-20}.}
\item{alpha}{A scalar between 1/2 and 1 used in the stepsequence for the Robbins-Monro method if \code{methodMC='RobbinsMC'}.}
\item{c}{The constant in the stepsequence if \code{methodMC='RobbinsMC'} or \code{methodMC='GradMC'}.}
\item{w}{The power for the weighted averaged Robbins-Monro algorithm if \code{methodMC='RobbinsMC'}.}
\item{epsilon}{Stoping condition for the Weiszfeld algorithm.}
\item{criterion}{The criterion for selecting the number of cluster. Can be \code{'ICL'} (default) or \code{'BIC'}.}
\item{methodMC}{The method chosen to estimate robustly the variance. Can be \code{'RobbinsMC'}, \code{'GradMC'} or \code{'FixMC'}.}
\item{par}{Is equal to \code{T} if the parallelization of the algorithm is allowed.}
\item{methodMCM}{The method chosen for estimating the Median Covariation Matrix. Can be \code{'Gmedian'} or \code{'Weiszfeld'}}
}

\value{
A list with:
\item{bestresult}{ A list giving all the results fo the best clustering (chosen with respect to the selected criterion.}
\item{allresults}{A list containing all the results.}
\item{ICL}{The ICL criterion for all the number of classes selected.}
\item{BIC}{The ICL criterion for all the number of classes selected.}
\item{data}{The initial data.}
\item{nclust}{A vector of positive integers giving the possible number of clusters.}
\item{Kopt}{The number of clusters chosen by the selected criterion.}
For the lists \code{bestresult} and \code{allresults[[k]]}:
\item{centers}{A matrix whose rows are the centers of the classes.}
\item{Sigma}{A matrix containing all the variance of the classes}
\item{LogLike}{The final LogLikelihood.}
\item{Pi}{A matrix giving the probabilities of each data to belong to each class.}
\item{niter}{The number of iterations of the EM algorithm.}
\item{initEM}{A vector giving the initialized clustering if \code{init='Mclust'} or \code{init='genie'}.}
\item{prop}{A vector giving the proportions of each classes.}
\item{outliers}{A vector giving the detected  outliers.}
}
\references{
Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.

Cardot, H. and Godichon-Baggioni, A. (2017). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Components Analysis.  \emph{Test}, 26(3), 461-480

Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and associated data depth. \emph{Proc. Natl. Acad. Sci. USA}, 97(4):1423-1426.

}
\seealso{
See also \code{\link{Gen_MM}}, \code{\link{RMMplot}} and \code{\link{RobVar}}.
}
\examples{
ech <- Gen_MM()
X <- ech$X
res <- RobMM(X , nclust=3,ninit=0)
RMMplot(res,graph=c('Two_Dim'))
}
\keyword{Robust Mixture Model}
