\name{VRclust-package}
\alias{VRclust-package}
\alias{VRclust}
\docType{package}
\title{VRclust: An optimal clustering algorithm with varying risks of false assignments across clusters}
\description{
There are three components to this package. These are: i) Construct and fit the finite mixture model when data are multivariate normal or multinomial (\code{vrclust.fit}). ii) Construct the optimal clustering rule for controlling the false assignment rate (\code{FAR.control}). iii) Obtain the estimated confusion matrix (\code{vrclust.est}).
}

\details{
\tabular{ll}{
Package: \tab VRclust\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2018-06-30\cr
License: \tab GPL version 2.0 or newer\cr
Depend: \tab mclust, poLCA\cr
LazyLoad: \tab yes\cr
}

Index:\cr
\tabular{ll}{
vrclust.fit \tab fit the discret latent variable model\cr
FAR.control \tab construct optimal clustering rule for controling predetermined level of False assignment rate \cr
vrclust.est \tab obtain the estimated confusion matrix based on the optimal clustering rule\cr
}

\author{
Donghwan Lee, Dongseok Choi, and Youngjo Lee
Maintainer: Donghwan Lee <donghwan.lee@ewha.ac.kr>
}

\references{
Donghwan Lee, Dongseok Choi, and Youngjo Lee. VRclust: An optimal clustering algorithm with varying risks of flase assignments across clusters - 2018 Submitted.
}

