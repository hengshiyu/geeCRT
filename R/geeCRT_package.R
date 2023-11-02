#' geeCRT: a package for implementing the bias-corrected generalized estimating equations in analyzing cluster randomized trials
#' @description geeCRT: a package for implementing the bias-corrected generalized estimating equations in analyzing cluster randomized trials
#' @section geeCRT functions:
#' The simbinPROBIT function performs correlated binary outcome data simulation using the multivariate probit method
#' The simbinCLF function performs correlated binary outcome data simulation using the conditional linear family method
#' The cpgeeSWD function performs cluster-period generalized estimating equations for estimating the marginal mean and correlation parameters in cross-sectional stepped wedge cluster randomized trials
#' The geemaee function performs matrix-adjusted generalized estimating equations on estimating the marginal mean and correlation parameters in cluster randomized trials
#' @docType package
#' @name geeCRT
#' @import MASS
#' @import mvtnorm
#' @import rootSolve
#' @importFrom stats binomial glm rbinom uniroot
#' @importFrom stats qnorm
#' @keywords internal
"_PACKAGE"
