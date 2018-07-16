\name{vrclust.est}
\alias{vrclust.est}
\title{
Obtain the estimated confusion matrix based on the optimal clustering rule
}
\description{
\code{vrclust.est} calculate the estimated confusion matrix with the clustering object and the optimal clustering rule.
}
\usage{
vrclust.est(clust, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clust}{ clustering object from \code{vrclust.fit}}
  \item{delta}{ A vector containing of cluster label based on \code{FAR.control}}
}

\value{
\item{Vmat}{K*K estimated confusion matrix}
\item{FAR}{Estimated False assignment rate}
}
\references{
Donghwan Lee, Dongseok Choi, and Youngjo Lee. VRclust: An optimal clustering algorithm with varying risks of flase assignments across clusters - 2018 Submitted.
}

\author{
Donghwan Lee, Dongseok Choi, and Youngjo Lee

Maintainer: Donghwan Lee <donghwan.lee@ewha.ac.kr>
}

\examples{
mc3k.clust<-vrclust.fit(x.irop, K=3, family="normal")
## ML rule
vrclust.est(mc3k.clust,mc3k.clust$delta)
}
