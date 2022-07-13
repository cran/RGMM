\name{Gen_MM}
\alias{Gen_MM}
\title{Gen_MM}

\description{Generate a sample of a Mixture Model
}
\usage{
Gen_MM(nk=NA, df=3, mu=NA, Sigma=FALSE, delta=0,cont="Student",
                model="Gaussian", dfcont=1, mucont=FALSE, Sigmacont=FALSE,
                minU=-20, maxU=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{nk}{An integer vector containing the desired number of data for each class. The defulat is \code{nk=rep(500,3)}.}
\item{df}{An integer larger (or qual) than \code{3} specifying the degrees of freedom of the Student law, if \code{model='Student'}. Default is \code{3}.}
\item{mu}{A numeric matrix whose raws correspond to the centers of the classes. By default, \code{mu} is generated randomly.}
\item{Sigma}{An array containing the variance of each class. See exemple for more details.}
\item{delta}{A positive scalr between \code{0} and \code{1} giving the proportion of contaminated data. Default is \code{0}.}
\item{cont}{The kind of contamination chosen. Can be equal to \code{'Unif'} or \code{'Student'}.}
\item{model}{A string character specifying the model chosen for the Mixture Model. Can be equal to \code{'Gaussian'} (default) or \code{'Student'}.}
\item{dfcont}{A positive integer specifying the degrees of freedom of the contamination laws if \code{cont='Student'}. Default is \code{1}.}
\item{mucont}{A numeric matrix whose rows correspond to the centers of the contamination laws. By default, \code{mucont} is chosen equal to \code{mu}.}
\item{Sigmacont}{An array containing the variance of each contamination law. By default, \code{sigmacont} is chosen equal to \code{sigma}.}
\item{minU}{A scalar giving the lower bound of the uniform law of the contamination if \code{cont='Unif'}.}
\item{maxU}{A scalar giving the upper bound of the uniform law of the contamination if \code{cont='Unif'}.}
}

\value{
A list with:
\item{Z}{ An integer vector specifying the true classification. If \code{delta} is non nul, the contaminated data are consider as a class.}
\item{C}{A \code{0-1} vector specifying the contaminated data.}
\item{X}{A numerical matrix giving the generated data.}
}

\seealso{
See also \code{\link{RobMM}} and \code{\link{RobVar}}.
}
\keyword{Mixture Model}
\examples{
p <- 3
nk <- rep(50,p)
mu <- c()
for (i in 1:length(nk))
{
  Z <- rnorm(3)
  mu <- rbind(mu,length(nk)*Z/sqrt(sum(Z^2)))
}
Sigma <- array(dim=c(length(nk), p, p))
for (i in 1:length(nk))
{
  Sigma[i, ,] <- diag(p)
}
ech <- Gen_MM(nk=nk,mu=mu,Sigma=Sigma)
}
