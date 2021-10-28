#' The KM for censoring time 
#' @param Ti is the vector of the ordered censored Time (ascent)
#' @param censor is the censor indicator for the Ti. 1 means no censor and 0 means censoring.
#' @return weight: The inverse weight of the C time. used in the cMBKRall function.
#' @export
 KMestimate_CTime = function(Ti,censor){
  N = length(Ti)
  risk = N:1
  
  #Ti no same
  risk=(risk-1)/risk
  #risk*censor
  #risk[N]=1
  risk[censor==1]=1
  
  weight=1/cumprod(risk)
  if(weight[N]==Inf){
    weight[N]=0
  }
  return(weight)
}


KMestimate_KM = function(Ti,censor){
  N = length(Ti)
  risk = N:1
  
  #Ti no same
  risk=(risk-1)/risk
  #risk*censor
  #risk[N]=1
  risk[censor==0]=1

  return(cumprod(risk))
  #return(risk)
}