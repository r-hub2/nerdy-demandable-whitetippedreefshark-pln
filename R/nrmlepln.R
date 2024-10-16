#' Full information maximum likelihood and bivariate composite likelihood estimation for
#' polytomous logit-normit (graded logistic) model
#' @aliases nrmlerasch nrbcpln
#' @description Full information maximum likelihood and bivariate composite likelihood estimation for
#' polytomous logit-normit and Rasch models, via Newton Raphson iterations.
#' @usage
#' nrmlepln(x, ncat, nitem=NULL, alphas=NULL, betas=NULL, abound=c(-10,10),
#'          bbound=c(-1,10), nq=48, mxiter=200, m2=TRUE, iprint=FALSE)
#' nrmlerasch(x, ncat, nitem=NULL, alphas=NULL, abound=c(-10,10),
#'          bbound=c(-1,10), nq=48, mxiter=200, m2=TRUE, iprint=FALSE)
#' nrbcpln(x, ncat, nitem=NULL, alphas=NULL, betas=NULL, abound=c(-10,10),
#'          bbound=c(-1,10), nq=48, mxiter=200, se=TRUE, iprint=FALSE)
#' @param x A data matrix. Data can be in one of two formats: 1) raw data
#'  where the number of rows corresponds to an individual's response and
#'  each column represents an item, and 2) a matrix of dimensions
#'  \code{nrec}\eqn{\times}{ X }(\code{nitem}+1) where each row corresponds to a
#'  response pattern and the last column is the frequency of that response
#'  pattern. A data matrix of the second type requires input for \code{nitem}
#'  and \code{nrec}.
#' @param ncat Number of ordinal categories for each item, coded as
#'  0,\dots,(\code{ncat}-1). Currently supported are items that have the same number of categories.
#' @param nitem Number of items. If omitted, it is assumed that \code{x} contains
#'  a data matrix of the first type (raw data) and the number of
#'  columns in \code{x} will be selected as the number of items.
#' @param alphas A vector of length \code{nitem}\eqn{\times}{ X }(\code{ncat}-1)
#'  corresponding to starting values for the (decreasing) cutpoints
#'  for the items. If omitted, these will be computed from the function \code{\link{startalphas}}.
#' @param betas A vector of length \code{nitem} corresponding to starting values for the
#'  beta vectors of slopes. If omitted, these will be computed from the function
#' \code{\link{startbetas}}. For the polytomous logit-normit, there is one slope
#'  for each item; for the Rasch model, there is a common slope beta for
#'  all of the items.
#' @param abound Vector of length 2 that sets upper and lower bounds on parameter estimation for
#'  alphas. Currently experimental; changing defaults it not recommended. Estimation problems
#'  are more likely solved by changing starting values.
#' @param bbound Vector of length 2 that sets upper and lower bounds on parameter estimation for
#'  betas. Currently experimental; changing defaults it not recommended. Estimation problems
#'  are more likely solved by changing starting values.
#' @param nq Number of quadrature points to use during estimation. This argument is currently experimental. It is recommended to use the default of 48.
#' @param mxiter Maximum number of iterations for estimation.
#' @param se Logical. If \code{TRUE}, calculates standard errors for the bivariate composite
#'  likelihood method.
#' @param m2 Logical. If \code{TRUE}, computes goodness-of-fit statistics
#'  from Maydeu-Olivares and Joe (2005, 2006; i.e., \eqn{M_2}).
#' @param iprint Logical. Enables debugging / diagnostic information from C code that conducts
#'  estimation.
#'  
#' @details
#'  Estimation of graded logistic models is performed under the following parameterization:
#'  \deqn{Pr(y_i = k_i| \eta) = \left\{
#'  \begin{array}{ll}
#'  1-\Psi (\alpha_{i,k} + \beta_i \eta) & \mbox{if  } k_i = 0\\
#'  \Psi (\alpha_{i,k} + \beta_i \eta) - \Psi (\alpha_{i,k+1} + \beta_i \eta) & \mbox{if  } 0 < k_i < m-1\\
#'  \Psi (\alpha_{i,k+1} + \beta_i \eta) & \mbox{if  } k_i = m-1
#'  \end{array} \right. }{
#'  Pr(y_i = k_i| \eta) = {
#'  1-\Psi (\alpha_i,k + \beta_i*\eta)  if k_i = 0, 
#'  \Psi (\alpha_i,k + \beta_i*\eta) - \Psi (\alpha_i,k+1 + \beta_i*\eta)  if 0 < k_i < m-1,  
#'  \Psi (\alpha_i,k+1 + \beta_i*\eta)  if k_i = m-1}.
#'  }
#'  Where the items are \eqn{y_i, i = 1, \dots, n}, and response categories are \eqn{k=0, \dots, m-1}. \eqn{\eta} is the latent trait, \eqn{\Psi} is the logistic distribution function, \eqn{\alpha} is an intercept (cutpoint) parameter, and \eqn{\beta} is a slope parameter. When the number of categories for the items is 2, this reduces to the 2PL parameterization:
#'  \deqn{Pr(y_i = 1| \eta) = \Psi (\alpha_1 + \beta_i \eta)}
#'  Both \code{nrmlepln} and \code{nrbcpln} perform estimation under these parameterizations, via Newton Raphson iterations, using full information maximum likelihood (\code{nrmlepln}) and bivariate composite likelihood (\code{nrbcpln}). See Maydeu-Olivares and Joe (2005, 2006) for more information on bivariate composite likelihood estimation (see also Varin, Reid, and Firth, 2011). Under \code{nrmlerasch} a common \eqn{\beta} parameter is estimated for all items.
#'  
#' @return A list containing the following slots.
#' @slot alphas A vector of parameter estimates for alphas. Length is
#'  \code{nitem}\eqn{\times}{ X }(\code{ncat}-1). Estimates are in order by item, e.g., all alphas
#'  for item 1, followed by all alphas for item 2, and so on.
#' @slot betas A vector of parameter estimates for betas. Length is \code{nitem}.
#' @slot nllk Negative (composite) log-likelihood for polytomous 
#'  logit-normit (or Rasch) model.
#' @slot conv Integer indicating whether estimation converged. Currently only returned for
#'  composite likelihood estimation.
#' @slot sealphas A vector of standard errors for the alpha estimates.
#' @slot sebetas A vector of standard errors for the beta estimates.
#' @slot invhes Inverse Hessian matrix for the MLE estimates.
#' @slot vcov Asymptotic covariance matrix for the composite likelihood estimates.
#' @slot teststat Value of \eqn{M_2}.
#' @slot df Degrees of freedom for \eqn{M_2}.
#' @slot pval P-value for \eqn{M_2}.
#'  
#' @references
#'  Bartholomew, D., Knott, M., and Moustaki, I. (2011). \emph{Latent Variable 
#'  Models and Factor Analysis: A Unified Approach}, 3rd Edition. Wiley. 
#'  
#'  Maydeu-Olivares, A., and Joe, H. (2005). Limited and full information estimation
#'  and goodness-of-fit testing in \eqn{2^n} contingency tables: A unified framework.
#'  \emph{Journal of the American Statistical Association, 100}, 1009-1020.
#'  
#'  Maydeu-Olivares, A., and Joe, H. (2006). Limited information and goodness-of-fit
#'  testing in multidimensional contingency tables. \emph{Psychometrika, 71},
#'  713-732.
#'  
#'  Varin, C., Reid, N. and Firth, D. (2011). An overview of composite likelihood
#'  methods. \emph{Statistica Sinica, 21}, 5-42.
#'  
#' @author Carl F. Falk \email{cffalk@gmail.com}, Harry Joe
#' 
#' @examples
#'  ### Matrix of response patterns and frequencies
#'  data(item5fr)
#'  
#'  ## ML estimation
#'  nrmleplnout<-nrmlepln(item5fr, ncat=3, nitem=5)
#'  print(nrmleplnout)
#'  
#'  ## BCL estimation
#'  nrbcplnout<-nrbcpln(item5fr, ncat=3, nitem=5)
#'  print(nrbcplnout)
#'  
#'  ## ML Rasch estimation
#'  nrmleraschout<-nrmlerasch(item5fr, ncat=3, nitem=5)
#'  print(nrmleraschout)
#'  
#'  \donttest{
#'  ### Raw data
#'  data(item9cat5)
#'  
#'  ## ML estimation
#'  nrmleplnout<-nrmlepln(item9cat5, ncat=5)
#'  print(nrmleplnout)
#'  
#'  ## BCL estimation
#'  nrbcplnout<-nrbcpln(item9cat5, ncat=5, se=FALSE)
#'  print(nrbcplnout)
#'  
#'  ## ML Rasch estimation
#'  nrmleraschout<-nrmlerasch(item9cat5, ncat=5)
#'  print(nrmleraschout)
#'  }
#'
#' @seealso  
#'  \code{\link{startalphas}}
#'  \code{\link{startbetas}}
#'
#' @importFrom stats pchisq
#' @useDynLib pln, .registration=TRUE
#' @export
nrmlepln <- function (x, ncat, nitem=NULL, alphas=NULL, betas=NULL, abound=c(-10,10),
    bbound=c(-1,10), nq=48, mxiter=200, m2=TRUE, iprint=FALSE) {

  myInput<-check.input(x, ncat, nitem, nq, mxiter, iprint)

  ## get starting values if not present already
  if(!check.alphas(alphas, myInput$nitem, myInput$ncat)){
    alphas<-startpln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX)$alphas
  }

  ## prep betas
  if(!check.betas(betas, myInput$nitem)){
    betas<-startbetas.func(myInput$myX)
  }
        
  ## check bounds
  abound<-check.bounds(alphas, abound)
  bbound<-check.bounds(betas, bbound)

  nrmleplnout <- nrmlepln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX, alphas,
        betas, abound, bbound, myInput$nq, myInput$mxiter, myInput$iprint)
  alphas<-nrmleplnout$mlePlnOut[1:((myInput$ncat-1)*myInput$nitem)]
  betas<-nrmleplnout$mlePlnOut[((myInput$ncat-1)*myInput$nitem+1):(myInput$ncat*myInput$nitem)]
  
  V<-matrix(nrmleplnout$invHesOut, nrow=myInput$nitem*myInput$ncat,ncol=myInput$nitem*myInput$ncat)
  seVec<-nrmleplnout$seVecOut
    
  sealphas<-seVec[1:((myInput$ncat-1)*myInput$nitem)]
  sebetas<-seVec[((myInput$ncat-1)*myInput$nitem+1):(myInput$ncat*myInput$nitem)]
    
    
  out<-list(alphas=alphas,betas=betas,nllk=nrmleplnout$nllkOut,sealphas=sealphas,
    sebetas=sebetas,invhes=V)

  if(m2){
    m2out<-m2pln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX, alphas, betas,
      myInput$nq, myInput$iprint)
    pval<-pchisq(m2out$m2stat, m2out$df, lower.tail=FALSE)
      ##out<-append(out, list(samplemout=m2out$samplem, m2stat=m2out$m2stat, df=m2out$df, pval=pval))      
    out<-append(out, list(teststat=m2out$m2stat, df=m2out$df, pval=pval))
  }

  return(out)
}

nrmlepln.func <- function(nitem, ncat, nrec, myX, alphas, betas, abound, bbound, nq, mxiter,
    iprint){
  
  ## prep return variables
  nllkOut<-0
  iconv<-0
  np<-ncat*nitem
  mlePlnOut<-rep(0,np)
  seVecOut<-rep(0,np)
  invHesOut<-matrix(0,nrow=np,ncol=np)

  ## TO DO: add PACKAGE argument here
  out <- .C("Rnrmlepln",
            as.integer(nitem), as.integer(ncat), as.integer(nrec), as.double(myX),
            as.double(alphas), as.double(betas), as.double(abound), as.double(bbound),
            nllkOut=as.double(nllkOut), mlePlnOut=as.double(mlePlnOut),
            seVecOut=as.double(seVecOut), invHesOut=as.double(invHesOut),
            as.integer(nq), as.integer(mxiter), iconv=as.integer(iconv),
            as.integer(iprint))
  return(out)
  
}
