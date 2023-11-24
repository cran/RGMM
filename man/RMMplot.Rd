\name{RMMplot}
\alias{RMMplot}
\title{RMMplot}

\description{A plot function for Robust Mixture Model
}
\usage{
RMMplot(a,outliers=TRUE,
    graph=c('Two_Dim','Two_Dim_Uncertainty','ICL','BIC',
    'Profiles','Uncertainty'),bestresult=TRUE,K=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{a}{Output from \code{\link{RobMM}}.}
\item{outliers}{An argument telling if there are outliers or note. In this case, Two dimensional plots and profiles plots will be done without detected outliers. Default is \code{TRUE}.}
\item{graph}{A string specifying the type of graph requested.
Default is \code{c('Two_Dim','Two_Dim_Uncertainty','ICL','BIC',
    'Profiles','Uncertainty')}.}
\item{bestresult}{A logical indicating if the graphs must be done for the result chosen by the selected criterion. Default is \code{TRUE}.}
\item{K}{A logical or positive integer giving the chosen number of clusters for each the graphs should be drawn.}
}


\seealso{
See also \code{\link{RobMM}} and \code{\link{Gen_MM}}.
}
\keyword{Mixture Model}
\examples{
\dontrun{
ech <- Gen_MM(mu = matrix(c(rep(-2,3),rep(2,3),rep(0,3)),byrow = TRUE,nrow=3))
 X <- ech$X
 res <- RobMM(X , nclust=3)
 RMMplot(res,graph=c('Two_Dim'))
 }
}
