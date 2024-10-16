#' Simulate data from polytomous logit-normit (graded logistic) model
#' @description
#'  Simulate data from polytomous logit-normit (graded logistic) model
#' @param n Number of responses to generate.
#' @param nitem Number of items.
#' @param ncat Number of categories for the items.
#' @param alphas A vector of length \code{nitem}\eqn{\times}{ X }(\code{ncat}-1)
#'  corresponding to true values for the (decreasing) cutpoints
#'  for the items.
#' @param betas A vector of length \code{nitem} corresponding to values for the
#'    beta vectors of slopes.
#'
#' @details
#' Data from graded logistic models is generated under the following parameterization:
#'    \deqn{Pr(y_i = k_i| \eta) = \left\{
#'      \begin{array}{ll}
#'      1-\Psi (\alpha_{i,k} + \beta_i \eta) & \mbox{if  } k_i = 0\\
#'      \Psi (\alpha_{i,k} + \beta_i \eta) - \Psi (\alpha_{i,k+1} + \beta_i \eta) & \mbox{if  } 0 < k_i < m-1\\
#'      \Psi (\alpha_{i,k+1} + \beta_i \eta) & \mbox{if  } k_i = m-1
#'      \end{array} \right. }{
#'        Pr(y_i = k_i| \eta) = {
#'          1-\Psi (\alpha_i,k + \beta_i*\eta)  if k_i = 0, 
#'          \Psi (\alpha_i,k + \beta_i*\eta) - \Psi (\alpha_i,k+1 + \beta_i*\eta)  if 0 < k_i < m-1,  
#'          \Psi (\alpha_i,k+1 + \beta_i*\eta)  if k_i = m-1}.
#'      }
#'    Where the items are \eqn{y_i, i = 1, \dots, n}, and response categories are \eqn{k=0, \dots, m-1}. \eqn{\eta} is the latent trait, \eqn{\Psi} is the logistic distribution function, \eqn{\alpha} is an intercept (cutpoint) parameter, and \eqn{\beta} is a slope parameter. When the number of categories for the items is 2, this reduces to the 2PL parameterization:
#'      \deqn{Pr(y_i = 1| \eta) = \Psi (\alpha_1 + \beta_i \eta)}
#'    
#' @return  A data matrix in which each row represents a response pattern and the final column represents the frequency of each response pattern.
#' @author Carl F. Falk \email{cffalk@gmail.com}, Harry Joe
#' @examples
#'    n<-500;
#'    ncat<-3;
#'    nitem<-5
#'    alphas=c(0,-.5,  .2,-1,  .4,-.6,  .3,-.2,  .5,-.5)
#'    betas=c(1,1,1,.5,.5)
#'    
#'    set.seed(1234567)
#'    datfr<-simulpln(n,nitem,ncat,alphas,betas)
#'    nrmleplnout<-nrmlepln(datfr, ncat=ncat, nitem=nitem)
#'    nrmleplnout
#' @seealso
#'    \code{\link{nrmlepln}}
#'    \code{\link{nrmlerasch}}
#'    \code{\link{nrbcpln}}
#' @export
simulpln <- function(n,nitem,ncat,alphas,betas) {
  if(length(betas)!=nitem){stop("Length of betas does not match number of items")}
  if(length(alphas)!=nitem*(ncat-1)){stop("Length of alphas does not match nitem*(ncat-1)")}
    
  ## Size may differ from what C wants to return
  ## But, this is the max size that may be returned
  datvec<-rep(0,n*(nitem+1))
    
  tem <- .C("Rsimulpln", as.integer(nitem), 
     as.integer(ncat),as.integer(n),
     as.double(alphas),as.double(betas),
     datvec=as.double(datvec))
    
  ##nrec<-length(tem$datvec)/(nitems+1)
  ##datfr<-matrix(tem$datvec,nrec,nitem+1)
  datfr<-matrix(tem$datvec, ncol=nitem+1, byrow=TRUE)
    
  ## trim rows with 0 frequencies
  datfr<-datfr[which(datfr[,ncol(datfr)]>0),]    
    
  datfr
}