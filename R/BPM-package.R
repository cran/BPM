#' BPM software package
#'
#' Bayesian model for purity estimation using DNA
#' methylation data
#'
#' The main function is \code{\link{BayPM}}
#' @docType package
#' @name BPM
#' @rdname BPM-package
#' @references  Jianzhao Gao, Linghao Shen, and Xiaodan Fan,
#' Bayesian model for purity estimation using DNA methylation data.(submitted)
#'
#' @author Jianzhao Gao(gaojz@@nankai.edu.cn), Linghao Shen Xiaodan Fan (xfan@@cuhk.edu.hk)
#'
#'
#' @examples
#' ### need to install package "limma"
#' ### source("https://bioconductor.org/biocLite.R");biocLite("limma");
#' library(BPM);
#' BayPM(simUCEC,20,2);

NULL

if (getRversion()>="2.15.1") utils::globalVariables(c("goodProbes","annotGeneNames"))
