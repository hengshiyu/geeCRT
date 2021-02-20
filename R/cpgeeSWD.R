#' Cluster-Period GEE for Estimating the Mean and Correlation Parameters in Cross-Sectional SW-CRTs
#' @param y a vector specifying the cluster-period means (proportions)
#' @param X design matrix for the marginal mean model, including period indicator and intervention indicator
#' @param id a vector specifying cluster identifier
#' @param n a vector of cluster "sample sizes" after cluster-period aggregation (equal to the number of measurement periods for each cluster)
#' @param m a vector of the cluster-period sizes
#' @param corstr correlation structure specified for the individual-level outcomes, could be \code{'exchangeable'}, \code{'nest_exch'} or \code{'exp_decay'}
#' @param family See corresponding documentation to \code{glm}. The current version only supports \code{family = 'binomial'}
#' @param maxiter maximum number of iterations for Fisher scoring updates
#' @param epsilon tolerance for convergence
#' @param printrange print details of range violations when \code{family = 'binomial'}. The default is \code{TRUE}
#' @param alpadj if \code{TRUE}, performs bias adjustment for the alpha estimating equations. The default is \code{FALSE}
#' @param rho.init user-specified initial value for the decay parameter when \code{corstr = 'exp_decay'}
#' @keywords stepped-wedge-cluster-randomized-trials cluster-period-means generalized-estimating-equations matrix-adjusted-estimating-equations bias-corrected-sandwich-variance
#' @author Hengshi Yu <hengshi@umich.edu>, Fan Li <fan.f.li@yale.edu>, Paul Rathouz <paul.rathouz@austin.utexas.edu>, Elizabeth L. Turner <liz.turner@duke.edu>, John Preisser <jpreisse@bios.unc.edu>
#' @description cpgeeSWD implements the cluster-period GEE developed for cross-sectional stepped wedge cluster randomized trials (SW-CRTs). It provides valid estimation and inference for the
#' treatment effect and intraclass correlation parameters within the GEE framework, and is computationally efficient for SW-CRTs with large cluster sizes. The program currently only allows for
#' a marginal mean model with discrete period effects and the intervention indicator without additional covariates. The program offers bias-corrected ICC estimates
#' as well as bias-corrected sandwich variances for both the treatment effect parameter and the ICC parameters. The technical details of the cluster-period GEE approach are provided in Li et al. (2020+).
#'
#' @references
#'
#' Zhao, L. P., Prentice, R. L. (1990). Correlated binary regression using a quadratic exponential model. Biometrika, 77(3), 642-648.
#'
#' Mancl, L. A., DeRouen, T. A. (2001). A covariance estimator for GEE with improved small sample properties. Biometrics, 57(1), 126-134.
#'
#' Kauermann, G., Carroll, R. J. (2001). A note on the efficiency of sandwich covariance matrix estimation. Journal of the American Statistical Association, 96(456), 1387-1396.
#'
#' Fay, M. P., Graubard, B. I. (2001). Small sample adjustments for Wald type tests using sandwich estimators. Biometrics, 57(4), 1198-1206.
#'
#' Lu, B., Preisser, J. S., Qaqish, B. F., Suchindran, C., Bangdiwala, S. I., Wolfson, M. (2007). A comparison of two bias corrected covariance estimators for generalized estimating equations. Biometrics, 63(3), 935-941.
#'
#' Preisser, J. S., Lu, B., Qaqish, B. F. (2008). Finite sample adjustments in estimating equations and covariance estimators for intracluster correlations. Statistics in Medicine, 27(27), 5764-5785.
#'
#' Li, F., Turner, E. L., Preisser, J. S. (2018). Sample size determination for GEE analyses of stepped wedge cluster randomized trials. Biometrics, 74(4), 1450-1458.
#'
#' Li, F. (2020). Design and analysis considerations for cohort stepped wedge cluster randomized trials with a decay correlation structure. Statistics in Medicine, 39(4), 438-455.
#'
#' Li, F., Yu, H., Rathouz, P., Turner, E. L., Preisser, J. S. (2020+). Marginal modeling of cluster period means and intraclass
#' correlations in stepped wedge designs with binary outcomes. Under Revision at Biostatistics.
#'
#' @export
#' @examples
#'
#' # Simulated SW-CRT example with binary outcome
#' 
#' ########################################################################
#' ### Example 1): simulated SW-CRT with smaller cluster-period sizes (5~10)
#' ########################################################################
#' 
#' sampleSWCRT = sampleSWCRTSmall
#' 
#' #############################################################
#' ### cluster-period id, period, outcome, and design matrix ###
#' #############################################################
#'
#' ### id, period, outcome
#' id = sampleSWCRT$id; period =  sampleSWCRT$period; y =  sampleSWCRT$y_bin
#' X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'treatment')])
#'
#' m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
#' clp_mu<-tapply(y,list(id,period), FUN=mean); y_cp <- c(t(clp_mu))
#'
#' ### design matrix for correlation parameters
#' trt <- tapply(X[, t + 1], list(id, period), FUN=mean); trt <- c(t(trt))
#'
#' time <- tapply(period,list(id, period), FUN = mean); time <- c(t(time)); X_cp <- matrix(0, n * t, t)
#'
#' s = 1
#' for(i in 1:n){for(j in 1:t){X_cp[s, time[s]] <- 1; s = s + 1}}
#' X_cp <- cbind(X_cp, trt); id_cp <- rep(1:n, each= t); n_cp <- rep(t, n); m_cp <-  c(t(m))
#'
#' #####################################################
#' ### cluster-period matrix-adjusted estimating equations (MAEE)
#' ### with exchangeable, nested exchangeable and exponential decay correlation structures ###
#' #####################################################
#'
#' # exchangeable
#' est_maee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                          n = n_cp, m = m_cp,
#'                          corstr = "exchangeable", alpadj = TRUE)
#' print(est_maee_exc)
#' 
#' # nested exchangeable
#' est_maee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                          n = n_cp, m = m_cp,
#'                          corstr = "nest_exch", alpadj = TRUE)
#' print(est_maee_nex)
#' 
#' # exponential decay
#' est_maee_ed <- cpgeeSWD(y  = y_cp, X = X_cp, id = id_cp,
#'                         n = n_cp, m = m_cp,
#'                         corstr = "exp_decay", alpadj = TRUE)
#' print(est_maee_ed)
#' 
#' #####################################################
#' ### cluster-period GEE
#' ### with exchangeable, nested exchangeable and exponential decay correlation structures ###
#' #####################################################
#'
#' # exchangeable
#' est_uee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                         n = n_cp, m = m_cp,
#'                         corstr = "exchangeable",alpadj = FALSE)
#' print(est_uee_exc)
#'
#' # nested exchangeable
#' est_uee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                         n = n_cp, m = m_cp,
#'                         corstr = "nest_exch", alpadj = FALSE)
#' print(est_uee_nex)
#'
#' # exponential decay
#' est_uee_ed <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                        n = n_cp, m = m_cp,
#'                        corstr = 'exp_decay', alpadj = FALSE)
#' print(est_uee_ed)
#' 
#' ########################################################################
#' ### Example 2): simulated SW-CRT with larger cluster-period sizes (20~30)
#' ########################################################################
#' 
#' sampleSWCRT = sampleSWCRTLarge
#' 
#' #############################################################
#' ### cluster-period id, period, outcome, and design matrix ###
#' #############################################################
#'
#' ### id, period, outcome
#' id = sampleSWCRT$id; period =  sampleSWCRT$period; y =  sampleSWCRT$y_bin
#' X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'period5', 'treatment')])
#'
#' m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
#' clp_mu<-tapply(y,list(id,period), FUN=mean); y_cp <- c(t(clp_mu))
#'
#' ### design matrix for correlation parameters
#' trt <- tapply(X[, t + 1], list(id, period), FUN=mean); trt <- c(t(trt))
#'
#' time <- tapply(period,list(id, period), FUN = mean); time <- c(t(time)); X_cp <- matrix(0, n * t, t)
#'
#' s = 1
#' for(i in 1:n){for(j in 1:t){X_cp[s, time[s]] <- 1; s = s + 1}}
#' X_cp <- cbind(X_cp, trt); id_cp <- rep(1:n, each= t); n_cp <- rep(t, n); m_cp <-  c(t(m))
#'
#' #####################################################
#' ### cluster-period matrix-adjusted estimating equations (MAEE)
#' ### with exchangeable, nested exchangeable and exponential decay correlation structures ###
#' #####################################################
#'
#' # exchangeable
#' est_maee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                          n = n_cp, m = m_cp,
#'                          corstr = "exchangeable", alpadj = TRUE)
#' print(est_maee_exc)
#'
#' # nested exchangeable
#' est_maee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                          n = n_cp, m = m_cp,
#'                          corstr = "nest_exch", alpadj = TRUE)
#' print(est_maee_nex)
#'
#' # exponential decay
#' est_maee_ed <- cpgeeSWD(y  = y_cp, X = X_cp, id = id_cp,
#'                         n = n_cp, m = m_cp,
#'                         corstr = "exp_decay", alpadj = TRUE)
#' print(est_maee_ed)
#'
#' #####################################################
#' ### cluster-period GEE
#' ### with exchangeable, nested exchangeable and exponential decay correlation structures ###
#' #####################################################
#'
#' # exchangeable
#' est_uee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                         n = n_cp, m = m_cp,
#'                         corstr = "exchangeable",alpadj = FALSE)
#' print(est_uee_exc)
#'
#' # nested exchangeable
#' est_uee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                         n = n_cp, m = m_cp,
#'                         corstr = "nest_exch", alpadj = FALSE)
#' print(est_uee_nex)
#'
#' # exponential decay
#' est_uee_ed <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp,
#'                        n = n_cp, m = m_cp,
#'                        corstr = 'exp_decay', alpadj = FALSE)
#' print(est_uee_ed)
#' 
#' 
#'
#'
#' @return \code{outbeta} estimates of marginal mean model parameters and standard errors with different finite-sample bias corrections.
#' The current version supports model-based standard error (MB), the sandwich standard error (BC0) extending Zhao and Prentice (2001),
#' the sandwich standard errors (BC1) extending Kauermann and Carroll (2001), the sandwich standard errors (BC2) extending Mancl and DeRouen (2001),
#' and the sandwich standard errors (BC3) extending the Fay and Graubard (2001). A summary of
#' these bias-corrections can also be found in Lu et al. (2007), and Li et al. (2018).
#' @return \code{outalpha} estimates of correlation parameters and standard errors with different finite-sample bias corrections.
#' The current version supports the sandwich standard error (BC0) extending Zhao and Prentice (2001),
#' the sandwich standard errors (BC1) extending Kauermann and Carroll (2001), the sandwich standard errors (BC2) extending
#' Mancl and DeRouen (2001), and the sandwich standard errors (BC3) extending the Fay and Graubard (2001). A summary of
#' these bias-corrections can also be found in Preisser et al. (2008).
#' @return \code{beta} a vector of estimates for marginal mean model parameters
#' @return \code{alpha} a vector of estimates of correlation parameters
#' @return \code{MB} model-based covariance estimate for the marginal mean model parameters
#' @return \code{BC0} robust sandwich covariance estimate of the marginal mean model and correlation parameters
#' @return \code{BC1} robust sandwich covariance estimate of the marginal mean model and correlation parameters with the
#' Kauermann and Carroll (2001) correction
#' @return \code{BC2} robust sandwich covariance estimate of the marginal mean model and correlation parameters with the
#' Mancl and DeRouen (2001) correction
#' @return \code{BC3} robust sandwich covariance estimate of the marginal mean model and correlation parameters with the
#' Fay and Graubard (2001) correction
#' @return \code{niter} number of iterations used in the Fisher scoring updates for model fitting


cpgeeSWD <- function(y, X, id, n, m, corstr, family = "binomial", maxiter = 500, epsilon = 0.001, printrange = TRUE, alpadj = FALSE, rho.init = NULL)
{
  if(corstr == "exchangeable"){
    cpgee_exc(y, X, id, n, m, family, maxiter, epsilon, printrange, alpadj)
  }else if(corstr == "nest_exch"){
    cpgee_nex(y, X, id, n, m, family, maxiter, epsilon, printrange, alpadj)
  }else if(corstr == "exp_decay"){
    cpgee_ed(y, X, id, n, m, family, maxiter, epsilon, printrange, alpadj, rho.init)
  }

}
