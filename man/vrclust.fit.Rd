\name{vrclust.fit}
\alias{vrclust.fit}
\title{
Fit the discrete latent variable model
}
\description{
\code{vrclust.fit} fit the discrete latent variable model. When data is multinomial, \code{poLCA} package is used, and when data is multivariate normal, \code{mclust} package is used for obtaining the maximum likelihood estimate of fixed unknown parameters and the class probabilites.
}
\usage{
vrclust.fit(x,K, family=c("normal","multinomial"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ n*p data matrix}
  \item{K}{ number of clusters}
  \item{family}{type of distributions of \code{x}. }
}

\value{
\item{clust}{The object of clustering}
\item{z}{Class probabilities}
\item{delta}{ML estimate of the clustering rule}
\item{loglik}{Log-likelihood values of the fitted model}
\item{bic}{BIC values}
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
}
