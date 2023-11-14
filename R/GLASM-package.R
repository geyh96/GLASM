#' @useDynLib GLASM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods as
#' @import Rcpp
#' @import RcppEigen
#' @import MASS
#' @importFrom MASS ginv
#' @importFrom stats pnorm
#' @importFrom RSpectra eigs_sym
#' @importFrom Matrix solve
NULL

#' Graphical Laplacian Augmentation feature Screening Method
#' @description A package for Graphical Laplacian Augmentation feature Screening Method

#' \tabular{ll}{
#'   Package: \tab GLASM\cr
#'   Type: \tab Package\cr
#'   Version: \tab 0.0.1\cr
#'   Date: \tab 2020-07-14\cr
#'   License: \tab GPL-2\cr
#'   LazyLoad: \tab yes\cr
#' }
#' The package "GALSMcensor" provides the detailed code for the paper Feng, Ge and Wu(manuscripts):\cr
#' @docType package
#' @aliases GLASMcensor-package
#' @author Yeheng Ge, Xingdong Feng, Tao Li, Hailong  Liu and Mengyun Wu \cr
#' Maintainers: Yeheng Ge<yeheng.ge@polyu.edu.hk>;
#' @references
#' 1.  T. Zhao and H. Liu. The huge Package for High-dimensional Undirected Graph Estimation in R. \emph{Journal of Machine Learning Research}, 2012\cr
#' 4.  Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
#' 6.  H. Liu, J. Lafferty and L. Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research}, 2009 \cr
#' 7.  J. Fan and J. Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, 2008.\cr