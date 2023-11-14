#' The function implementing the network structured screening method proposed in this manuscripts.
#' @title abrw
#' @description ABsorbed Random Walk Algorithm for network structured screening method.
#' @param list_amatrix list_amatrix is a list contains the matrix information you want. Every matrix is an p by p adjacency matrix here.
#' @param utility utility is the marginal screening meature vector with length p.
#' @param list_alpha is the tuning parameter in the proposed method. It shold have the same length as the number of matrices.
#' @param max_it max_it it the maximum of the iteration
#' @param threshold thresdhold is an feasible choice to deliminate the small value in the utility. Default 0 means use the whole utility vector.
#' @export 
abrw=function(list_amatrix,utility,list_alpha,max_it=100,threshold=0)
{
  utility=abs(utility)
  utility[utility<threshold]=0

  uti_all=utility

  P=dim(list_amatrix[[1]])[1]
  uti_all = matrix(uti_all,P,1)
  
  amatrix_all = matrix(0,P,P)
  for(imat in 1:length(list_amatrix)){
        diag(list_amatrix[[imat]]) = 1
        amatrix_all = amatrix_all + list_alpha[imat] * list_amatrix[[imat]]
  }
  diag(amatrix_all)=rep(0,P,P)
  dgc_A=as(amatrix_all, "dgCMatrix")
  idx_noisolate = which(rowSums(dgc_A)>0)
  


  aamatrix = amatrix_all[idx_noisolate,idx_noisolate]
  utility = uti_all[idx_noisolate,]
  pdim = dim(aamatrix)[1]
  utility = matrix(utility,pdim,1)
  dgc_A=as(aamatrix, "dgCMatrix")
  dgc_D=as(diag(rowSums(dgc_A),pdim,pdim), "dgCMatrix")
  dgc_Dinv=as(diag(1/rowSums(dgc_A),pdim,pdim), "dgCMatrix")

  Mat_B=Matrix::solve(diag(1,pdim)+dgc_D)%*%dgc_D
  Mat_C=Matrix::solve(diag(1,pdim)+dgc_D)

  rt=matrix(0,pdim,1)
  rt_1=utility

  error=sum((rt-rt_1)^2)
  itercount=0
  dim(Mat_B)
  dim(dgc_Dinv)
  dim(dgc_A)
  dim(rt_1)
  while(error>1e-6 && itercount<max_it){
    rt=Mat_B%*%dgc_Dinv%*%dgc_A%*%rt_1 + Mat_C%*%utility
    error=max(abs(rt-rt_1))
    rt_1=rt
    itercount=itercount+1
  }

  result = rep(0,P)
  result = uti_all
  result[idx_noisolate] = rt
  return(as.vector(result))
 }