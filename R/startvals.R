#' Starting values for polytomous logit-normit model
#' @aliases startbetas
#' @description Computes starting values for estimation of polytomous logit-normit model.
#' @usage
#' startalphas(x, ncat, nitem = NULL)
#' startbetas(x, ncat, nitem = NULL)
#' @param x A data matrix. Data can be in one of two formats: 1) raw data
#'   where the number of rows corresponds to the number of raw cases and
#'   each column represents an item, and 2) a matrix of dimensions
#'   \code{nrec}\eqn{\times}{X}(\code{nitem}+1) where each row corresponds to a
#'   response pattern and the last column is the frequency of that response
#'   pattern. A data matrix of the second type requires input for \code{nitem}
#'   and \code{nrec}.
#' @param ncat Number of ordinal categories for each item, coded as
#'    0,\dots,(\code{ncat}-1). Currently supported are items that have the same number of categories.
#' @param nitem Number of items. If omitted, it is assumed that \code{x} contains
#'    a data matrix of the first type (raw data) and the number of
#'    columns in \code{x} will be selected as the number of items.
#' @details \code{startalphas} computes starting values for the (decreasing) cutpoints
#'   for the items based on logit transformed probabilities, assuming independent items.
#'
#'   \code{startbetas} computes starting values for slopes under the polytomous
#'   logit-normit model, using a method based on values that are proportional to the
#'   average correlations of each item with all other items. Starting values are
#'   currently bounded between -.2 and 1.
#' @return
#'  A vector of starting values, depending on which function was called.
#' @author Carl F. Falk \email{cffalk@gmail.com}, Harry Joe
#' @seealso
#'  \code{\link{nrmlepln}}
#'  \code{\link{nrmlerasch}}
#'  \code{\link{nrbcpln}}
#' @examples
#' ### Raw data
#' data(item9cat5)
#' 
#' myAlphas<-startalphas(item9cat5, ncat=5)
#' print(myAlphas)
#' 
#' myBetas<-startbetas(item9cat5, ncat=5)
#' print(myBetas)
#' 
#' nrbcplnout<-nrbcpln(item9cat5, ncat=5, alphas=myAlphas, betas=myBetas, se=FALSE)
#' print(nrbcplnout)
#' 
#' ## Matrix of response patterns and frequencies
#' data(item5fr)
#' 
#' myAlphas<-startalphas(item5fr, ncat=3, nitem=5)
#' print(myAlphas)
#' 
#' myBetas<-startbetas(item5fr, ncat=3, nitem=5)
#' print(myBetas)
#' 
#' nrbcplnout<-nrbcpln(item5fr, ncat=3, nitem=5, alphas=myAlphas, betas=myBetas, se=FALSE)
#' print(nrbcplnout)
#'   
#' @export
startalphas<- function (x, ncat, nitem=NULL) {

  ## input checking & data prep
  myInput<-check.input(x, ncat, nitem)
  
  out <- startpln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX)

  return(out$alphas)

}

startpln.func<-function(nitem, ncat, nrec, myX){

  ## prep return variables
  alphas=rep(0,nitem*(ncat-1))

  ## TO DO: add PACKAGE argument here
  out <- .C("Rstartpln",
            as.integer(nitem), as.integer(ncat), as.integer(nrec), as.double(myX),
            alphas=as.double(alphas))
  return(out)
}

#' @export
startbetas<-function(x, ncat, nitem=NULL){
  myInput<-check.input(x, ncat, nitem)
  x<-myInput$myX
  betas<-startbetas.func(x)
  return(betas)
}

startbetas.func<-function(x){
  nn<-ncol(x)
  nitem<-nn-1
  y<-x[,-nn]
  fr<-x[,nn]
  tot<-sum(fr)
  mnvec<-apply(y*fr,2,sum)/tot ## means
  vvec<-apply(y*y*fr,2,sum)/tot ## variances
  vvec<-vvec-mnvec^2
  cc<-matrix(0,nitem,nitem)
  for(j in 1:(nitem-1))
  { for(k in (j+1):nitem)
    { ss<-sum(y[,j]*y[,k]*fr)/tot
      den<-vvec[j]*vvec[k]
      cc[j,k]<-(ss-mnvec[j]*mnvec[k])/sqrt(den)
   	  cc[k,j]<-cc[j,k]
    }
  }
  avcc<-apply(cc,2,sum)/(nitem-1)
  #print(avcc)
  bvec<-avcc*10; bvec[bvec>=1]=1; bvec[bvec<= -.2]= -.2
  #print(bvec)
  bvec
}