
#' calculate the conditional correlation among the X's
#'
#' @param X The data matrix for the X
#' @return 
#' A matrix of size p by b contains the conditional correlation among the X.
#' @export
estimate_Corr=function(X){
        P=dim(X)[2]
        N=dim(X)[1]
        BetaMat <- SPMBgraphsqrtKeh(X, lambda = sqrt(log(P) / N), nlambda = as.integer(1), d = (P))
        Corr_cpp = calculate_corrhat(X, BetaMat, N, P)
    return(Corr_cpp)
}


#' calculate the threshold of  Test statistics for the independence of X with control of FDR.
#'
#'
#' @param Corr The estimated conditional correlation among X
#' @param N The sample size N
#' @param P The dimension of X
#' @param alpha_fdr The pregiven FDR level
#' @return 
#' A matrix of size p by b contains the conditional correlation among the X.
#' @export
calculate_Test_threshold=function(Corr,N,P,alpha_fdr=0.01){
        Gt=function(t){
         2 - 2*pnorm(t)
        }
        calculate_FDR=function(thres,P,Tstar){
            (Gt(thres)*(P*P-P)*0.5)/max(1,sum(abs(Tstar)>thres))
        }
        # P=dim(X)[2]
        # N=dim(X)[1]
        # BetaMat <- SPMBgraphsqrtKeh(X, lambda = sqrt(log(P) / N), nlambda = as.integer(1), d = (P))
        # Test_cpp = calculate_Testhat(X, BetaMat, N, P)
        # Corr[upper,tri(Corr)]
        Test_cpp = Corr[upper.tri(Corr)]
        # dim(a)
        # length(a[upper.tri(a)])
        thres=seq(1,10,0.1)
        fdp_for_thres = sapply(thres,calculate_FDR,P=P,Tstar=sqrt(N)*Test_cpp)
        thres_selected=thres[which(fdp_for_thres<alpha_fdr)[1]]/sqrt(N)
    return(thres_selected)
}

#' Estimate the gaussian graphical model with FDR control.
#' @param Corr the calculated conditional correlation of of the data X.
#' @param thres the calculated threshold for the test.
#' @return 
#' A matrix of size p by b contains the 0-1 graph.
#' @export
estimate_network=function(Corr,thres){
    amatrix_hat = ( abs(Corr)>thres) + 0
    return(amatrix_hat)
}
