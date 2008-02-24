\name{brglm.control}
\alias{brglm.control}
\title{Auxiliary for Controlling BRGLM Fitting}
\description{
Auxiliary function as user interface for \code{\link{brglm}}
fitting. Typically only used when calling \code{brglm} or \code{brglm.fit}.
}
\usage{
brglm.control(br.epsilon = 1e-08, br.maxit = 100, br.trace=FALSE, ...)
}
\arguments{
  \item{br.epsilon}{positive convergence tolerance for the iteration
    described in \code{\link{brglm.fit}}.}
  \item{br.maxit}{integer giving the maximum number of iterations for
  the iteration in \code{\link{brglm.fit}}.}
  \item{br.trace}{logical indicating if output should be prooduced for
    each iteration.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  If \code{br.trace=TRUE} then for each iteration the iteration number
  and the current value of the modified scores is \code{\link{cat}}'ed.
}
\value{
  A list with the arguments as components.
}
\references{
  Kosmidis, I. (2007). Bias reduction in exponential family nonlinear
  models. \emph{PhD Thesis}, Department of Statistics, University of
  Warwick.
}
\author{Ioannis Kosmidis, \email{I.Kosmidis@warwick.ac.uk}}
\seealso{
  \code{\link{brglm.fit}}, the fitting procedure used by
  \code{\link{brglm}}. 
}
\keyword{iteration}
