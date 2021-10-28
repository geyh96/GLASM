
#' Graph Laplacian Augmented screening method 
#'
#' The algorithm based on the effective power iteration method for the partial absorbing random walk algorithm
#'
#' @param amatrix the adjacency matrix with diagonal elements equal 0.
#' @param utility the original screening measure.
#' @param alpha the paramter balances the original ranking measure and the network information
#' @param max_it the maximum of iteration
#' @param threshold make all marignal screening measure that are smaller than the threshold equals 0 to denoise.
#' @return 
#' A vector of size p contains the Graph Laplacian Augmented screening measure
#' @export
PARWaug <- function(amatrix, utility, alpha, max_it = 100, threshold = 0) {
  # uti=utility
  utility <- abs(utility)
  utility[utility < threshold] <- 0
  # utility=utility/sum(utility)
  pdim <- dim(amatrix)[1]
  utility <- matrix(utility, pdim, 1)
  diag(amatrix) <- rep(0, pdim)
  dgc_A <- as(amatrix, "dgCMatrix")
  dgc_D <- as(diag(rowSums(dgc_A), pdim, pdim), "dgCMatrix")
  Mat_1 <- alpha * Matrix::solve(diag(1, pdim) + alpha * dgc_D) %*% dgc_A
  Mat_2 <- Matrix::solve(diag(1, pdim) + alpha * dgc_D)

  rt <- matrix(0, pdim, 1)
  rt_1 <- utility
  error <- sum(abs(rt - rt_1))
  itercount <- 0
  # dgc_Dinv <- as(diag(rowSums(dgc_A), pdim, pdim), "dgCMatrix")
  while (error > 1e-6 && itercount < max_it) {
    rt <- Mat_1 %*% rt_1 + Mat_2 %*% utility
    error <- sum((rt - rt_1)^2)
    rt_1 <- rt
    itercount <- itercount + 1
  }

  return(as.vector(rt))
}


#' PageRank Augmented screening method proposed by Wu., Zhu. and Feng.(2018)
#'
#' The algorithm based on the effective power iteration method.
#'
#' @param A the adjacency matrix with diagonal elements equal 0.
#' @param Orig_RM the original screening measure.
#' @param tol error for  convergence condition.
#' @param d the paramter balances the original ranking measure and the network information
#' @param r the size of the variable set including all active predictors
#' @param max_it the maximum of iteration
#' @param threshold make all marignal screening measure that are smaller than the threshold equals 0 to denoise.
#' @return 
#' A matrix of size p by 1 contains the PageRank Augmented screening measure
#' @export
PageRankaug<-function (A, Orig_RM,d,r,tol=1e-6,max_it=100,threshold=0) 
{
  #####  INPUT #####
  # A: adjacency matrix
  # Orig_RM: original ranking measures
  # d: the parameter balances the original ranking measure and the network information
  # r: the size of the variable set including all active predictors

  ### OUTPUT ####
  # Net_RM: the network-based ranking measure
  Orig_RM=abs(Orig_RM)
  Orig_RM[Orig_RM<threshold]=0
  p=length(Orig_RM)
  diag(A)=rep(0,p)
  degree=colSums(A)
  id=(degree>0)
  # degree[is.infinite(degree)]=0
  D2=1/colSums(A)
  D2[!id]=0
  D2=matrix(D2,1,p)
  A_norm=A*(matrix(1,p,1)%*%D2)

  rank_final=matrix(0,p,1)
  temp_matrix=A_norm
  Orig_RM[sort(Orig_RM,index.return=TRUE)$ix[1:(p-r)]]=1e-20
  Orig_RM_norm= matrix(Orig_RM/sum(Orig_RM),p,1)
  diff=1
  max_iterate=1
  rrank=matrix(1/p,p,1)

  G2=Orig_RM_norm
  temp_matrix2=temp_matrix[id,id]

  while (diff>tol&& max_iterate<max_it){
    G2[id]=matrix(temp_matrix2%*%rrank[id])
    rrank2=(1-d)*Orig_RM_norm+d*G2
    diff=sum(abs(rrank-rrank2))
    rrank=rrank2
    max_iterate=max_iterate+1

  }
  rank_final=rrank
  return(rank_final)
}

