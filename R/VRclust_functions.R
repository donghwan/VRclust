library(mclust)
library(poLCA)


vrclust.fit<-function(x,K=2, family=c("normal","multinomial"))    {

  n<-dim(x)[1]; p<-dim(x)[2]

  if (family=="normal"){

                clust<-Mclust(x,G=K)
                prior<-clust$parameters$pro
                z<-clust$z
                delta<-clust$classification
                loglik<-clust$loglik
                bic<-clust$bic

  }

  if (family=="multinomial"){

    model<-do.call(cbind,x)~1
    clust<-poLCA(model,x,nclass=K,nrep=10)

    prior<-clust$probs
    z<-clust$posterior
    delta<-clust$predclass
    loglik<-clust$llik
    bic<-clust$bic
  }



  return(list(clust=clust, prior=prior, z=z, delta=delta, loglik=loglik, bic=bic))
}


Unmap<-function(delta,K){
  n<-length(delta)

  delmat<-matrix(0,n,K)
  for (k in 1:K){
    delmat[delta==k,k]<-1
  }

  delmat
}

reorder<-function(clust, ord=c(1,4,2,3)){
  clust$z<-clust$z[,ord]
  temp<-clust$delta
  for ( i in 1:length(ord)){
    temp[clust$delta==ord[i]]<-i
  }
  clust$delta<-temp
  return(clust)
}


vrclust.est<-function(clust, delta){

  z<-clust$z
  n<-dim(z)[1];K<-dim(z)[2]

  ## Calculation of Error rate

  dmat<-Unmap(delta, K)
  Vmat<-t(z)%*%dmat

  FPR<-FDR<-FAR<-rep(0,K)
  R<-colSums(Vmat)
  N<-rowSums(Vmat)

  for (k in 1:K){
      FAR[k]<-sum(Vmat[-k,k])/R[k]


      FDR[k]<-Vmat[1,k]/R[k]
      if (R[k]==0) {FAR[k]<-0; FDR[k]<-0}

      FPR[k]<-Vmat[1,k]/N[1]
    }

  return(list(Vmat=Vmat, FAR=FAR, FDR=FDR, FPR=FPR))

}

FPR.control<-function(clust, class=1, lambda=NULL, level=0.05){

  z<-clust$z

  n<-dim(z)[1];K<-dim(z)[2]



  if (is.null(lambda)) {lambda<-rep(1,K)}
  if (is.null(level)){level<-0.05}

  conv.code<-FALSE
  niter<-1
  lambda.old<-lambda
  k<-class

  a<-0; b<-9999

  while(conv.code==FALSE ){
    lambda.old<-lambda[k]
    lambda[k]<-(a+b)/2

    delta<-delta.est(clust, lambda)
    err.temp<-vrclust.est(clust, delta)
    FPR<-err.temp$FPR

    dif<-FPR[k]-level

    if (dif>0) {a<-lambda[k]}
    if (dif<0) {b<-lambda[k]}

    lambda.diff<-abs(lambda.old-lambda[k])
    if (abs(dif)<0.001 | niter>100|lambda.diff<1e-05) {conv.code<-TRUE}
    niter<-niter+1
    #print(c(niter, FPR[k], lambda, conv.code))
  }
  print(c(niter, FPR[k], lambda, conv.code))

  return(list(delta=delta,lambda=lambda))


}

FAR.control<-function(clust, class=1, lambda=NULL, level=0.05){

  z<-clust$z

  n<-dim(z)[1];K<-dim(z)[2]



  if (is.null(lambda)) {lambda<-rep(1,K)}
  if (is.null(level)){level<-0.05}

  conv.code<-FALSE
  niter<-1
  lambda.old<-lambda
  k<-class

  a<-0; b<-9999

  while(conv.code==FALSE ){
    lambda.old<-lambda[k]
    lambda[k]<-(a+b)/2

    delta<-delta.est(clust, lambda)
    err.temp<-vrclust.est(clust, delta)
    FAR<-err.temp$FAR

    dif<-FAR[k]-level

    if (dif>0) {a<-lambda[k]}
    if (dif<0) {b<-lambda[k]}

    lambda.diff<-abs(lambda.old-lambda[k])
    if (abs(dif)<0.001 | niter>100|lambda.diff<1e-05) {conv.code<-TRUE}
    niter<-niter+1
    #print(c(niter, FAR[k], lambda, conv.code))
  }
  print(c(niter, FAR[k], lambda, conv.code))

  return(list(delta=delta,lambda=lambda))


}

delta.est<-function(clust, lambda=NULL){

  z<-clust$z

  n<-dim(z)[1];K<-dim(z)[2]

  if (is.null(lambda)) {lambda<-rep(1,K)}

  delta.mat<- (1-z)

  for (i in 1:K){
    delta.mat[,i]<-  lambda[i]*(1-z[,i])
  }

  delta<-c(apply(delta.mat, 1, function(x) which(x == min(x, na.rm = TRUE))) )

  delta
}


delta.est.co<-function(clust, lambda=NULL){

  z<-clust$z

  n<-dim(z)[1];K<-dim(z)[2]

  if (is.null(lambda)) {lambda<-matrix(1,K,K)}

  delta.mat<- matrix(0,n,K)

  for (k in 1:K){
    temp<-0
    for (j in c(1:K)[-k]){
        temp<-temp+lambda[j,k]*z[,j]
    }
    delta.mat[,k]<-  temp
  }

  delta<-c(apply(delta.mat, 1, function(x) which(x == min(x, na.rm = TRUE))) )

  delta
}

## TFAR and WFAR

TFAR<-function(A,lambda.mat){
  n<-sum(A)
  B<-A
  diag(B)<-0
  TFAR<-sum(B)/n
  WFAR<-sum(B*lambda.mat)/n/sum(lambda.mat)
  return(list(raw=A,TFAR=TFAR,WFAR=WFAR))
}
