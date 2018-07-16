\name{FAR.control}
\alias{FAR.control}
\title{
Construct optimal clustering rule for controling predetermined level of False assignment rate
}
\description{
\code{FAR.control} provide the optimal clustering rule for controlling the false assignment rate at the predetermined level.
}
\usage{
FAR.control(clust, class=1, lambda=NULL, level=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clust}{ clustering object from \code{vrclust.fit}}
  \item{class}{ Label number of targeted cluster which be controlled with the FAR level}
  \item{lambda}{A vector of relative cost lambda. When controlling FAR, the default is NULL.}
  \item{level}{Level of targeted FAR}
}

\value{
\item{delta}{A vector of cluster labels of the optimal clustering rule}
\item{lambda}{relative cost corresponding to FAR level}
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
