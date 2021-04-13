#' GEE and Matrix-adjusted Estimating Equations (MAEE) for Estimating the Marginal Mean and Correlation Parameters in CRTs
#' @param y a vector specifying the outcome variable across all clusters
#' @param X design matrix for the marginal mean model, including the intercept
#' @param id a vector specifying cluster identifier
#' @param Z design matrix for the correlation model, should be all pairs j < k for each cluster
#' @param family See corresponding documentation to \code{glm}. The current version only supports \code{'continuous'} and \code{'binomial'}
#' @param maxiter maximum number of iterations for Fisher scoring updates
#' @param epsilon tolerance for convergence. The default is 0.001
#' @param printrange print details of range violations. The default is \code{TRUE}
#' @param alpadj if \code{TRUE}, performs bias adjustment for the correlation estimating equations. The default is \code{FALSE}
#' @param shrink method to tune step sizes in case of non-convergence including \code{'THETA'} or \code{'ALPHA'}. The default is \code{'ALPHA'}
#' @param makevone if \code{TRUE}, it assumes unit variances for the correlation parameters in the correlation estimating equations. The default is \code{TRUE}
#' @keywords cluster-randomized-trials generalized-estimating-equations matrix-adjusted-estimating-equations bias-corrected-sandwich-variance
#' @author Hengshi Yu <hengshi@umich.edu>, Fan Li <fan.f.li@yale.edu>, Paul Rathouz <paul.rathouz@austin.utexas.edu>, Elizabeth L. Turner <liz.turner@duke.edu>, John Preisser <jpreisse@bios.unc.edu>
#' @description geemaee implements the GEE and matrix-adjusted estimating equations (MAEE) for analyzing cluster randomized trials (CRTs). It supports estimation and inference for the
#' marginal mean and intraclass correlation parameters within the population-averaged modeling framework. With suitable choice of the design matrices, the function can be used to analyze
#' parallel, crossover and stepped wedge cluster randomized trials. The program also offers bias-corrected intraclass correlation estimates, as well as bias-corrected sandwich variances for both the marginal mean and correlation parameters.
#' The technical details of the GEE and MAEE approach are provided in Preisser (2008) and Li et al. (2018, 2019).
#' @references
#' Liang, K. Y., Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1), 13-22.
#'
#' Prentice, R. L. (1988). Correlated binary regression with covariates specific to each binary observation. Biometrics, 1033-1048.
#'
#' Zhao, L. P., Prentice, R. L. (1990). Correlated binary regression using a quadratic exponential model. Biometrika, 77(3), 642-648.
#'
#' Prentice, R. L., Zhao, L. P. (1991). Estimating equations for parameters in means and covariances of multivariate discrete and continuous responses. Biometrics, 825-839.
#'
#' Sharples, K., Breslow, N. (1992). Regression analysis of correlated binary data: some small sample results for the estimating equation approach. Journal of Statistical Computation and Simulation, 42(1-2), 1-20.
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
#' Li, F., Forbes, A. B., Turner, E. L., Preisser, J. S. (2019). Power and sample size requirements for GEE analyses of cluster randomized crossover trials. Statistics in Medicine, 38(4), 636-649.
#'
#' Li, F. (2020). Design and analysis considerations for cohort stepped wedge cluster randomized trials with a decay correlation structure. Statistics in Medicine, 39(4), 438-455.
#'
#' Li, F., Yu, H., Rathouz, P., Turner, E. L., Preisser, J. S. (2021). Marginal modeling of cluster-period means and intraclass correlations in stepped wedge designs with binary outcomes. Biostatistics, kxaa056. 
#' @export
#' @examples
#'
#' # Simulated SW-CRT examples
#'
#' #################################################################
#' ### function to create the design matrix for correlation parameters
#' ### under the nested exchangeable correlation structure
#' #################################################################
#' createzCrossSec = function (m) {
#'       Z = NULL
#'       n = dim(m)[1]
#'       for (i in 1:n) {
#'           alpha_0 = 1; alpha_1 = 2; n_i = c(m[i, ]); n_length = length(n_i)
#'           POS = matrix(alpha_1, sum(n_i), sum(n_i))
#'           loc1 = 0; loc2 = 0
#'           for (s in 1:n_length) {
#'               n_t = n_i[s]; loc1 = loc2 + 1; loc2 = loc1 + n_t - 1
#'               for (k in loc1:loc2) {
#'                   for (j in loc1:loc2) {
#'                       if (k != j) {POS[k, j] = alpha_0} else {POS[k, j] = 0}}}}
#'            zrow = diag(2); z_c = NULL
#'            for (j in 1:(sum(n_i) - 1)) {
#'                for (k in (j + 1):sum(n_i)) {z_c = rbind(z_c, zrow[POS[j,k],])}}
#'            Z = rbind(Z, z_c)}
#'       return(Z)}
#'
#' ########################################################################
#' ### Example 1): simulated SW-CRT with smaller cluster-period sizes (5~10)
#' ########################################################################
#' 
#' sampleSWCRT = sampleSWCRTSmall
#' 
#' ###############################################################
#' ### Individual-level id, period, outcome, and design matrix ###
#' ###############################################################
#' 
#' id = sampleSWCRT$id; period =  sampleSWCRT$period;
#' X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'treatment')])
#'
#' m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]; 
#'
#' ### design matrix for correlation parameters
#' Z = createzCrossSec(m)
#'
#' ################################################################
#' ### (1) Matrix-adjusted estimating equations and GEE
#' ### on continous outcome with nested exchangeable correlation structure
#' ################################################################
#'
#' ### MAEE
#' est_maee_ind_con = geemaee(y = sampleSWCRT$y_con, X = X, id = id,
#'                            Z = Z, family = 'continuous',
#'                            maxiter = 500, epsilon = 0.001,
#'                            printrange = TRUE, alpadj = TRUE,
#'                            shrink = 'ALPHA', makevone = FALSE)
#' print(est_maee_ind_con)
#' 
#' 
#' ### GEE
#' est_uee_ind_con = geemaee(y = sampleSWCRT$y_con, X = X, id = id,
#'                           Z = Z, family = 'continuous',
#'                           maxiter = 500, epsilon = 0.001,
#'                           printrange = TRUE, alpadj = FALSE,
#'                           shrink = 'ALPHA', makevone = FALSE)
#' print(est_uee_ind_con)
#' 
#'
#' ###############################################################
#' ### (2) Matrix-adjusted estimating equations and GEE
#' ### on binary outcome with nested exchangeable correlation structure
#' ###############################################################
#'
#' ### MAEE
#' est_maee_ind_bin = geemaee(y = sampleSWCRT$y_bin, X = X, id = id,
#'                            Z = Z, family = 'binomial',
#'                            maxiter = 500, epsilon = 0.001,
#'                            printrange = TRUE, alpadj = TRUE,
#'                            shrink = 'ALPHA', makevone = FALSE)
#' print(est_maee_ind_bin)
#' 
#' 
#' ### GEE
#' est_uee_ind_bin = geemaee(y = sampleSWCRT$y_bin, X = X, id = id,
#'                           Z = Z, family = 'binomial',
#'                           maxiter = 500, epsilon = 0.001,
#'                           printrange = TRUE, alpadj = FALSE,
#'                           shrink = 'ALPHA', makevone = FALSE)
#' print(est_uee_ind_bin)
#' 
#' 
#' \donttest{
#' ## This will elapse longer. 
#' ########################################################################
#' ### Example 2): simulated SW-CRT with larger cluster-period sizes (20~30)
#' ########################################################################
#' 
#' sampleSWCRT = sampleSWCRTLarge
#' 
#' ###############################################################
#' ### Individual-level id, period, outcome, and design matrix ###
#' ###############################################################
#' 
#' id = sampleSWCRT$id; period =  sampleSWCRT$period;
#' X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'period5', 'treatment')])
#'
#' m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]; 
#'
#' ### design matrix for correlation parameters
#' Z = createzCrossSec(m)
#'
#' ################################################################
#' ### (1) Matrix-adjusted estimating equations and GEE
#' ### on continous outcome with nested exchangeable correlation structure
#' ################################################################
#'
#' ### MAEE
#' est_maee_ind_con = geemaee(y = sampleSWCRT$y_con, X = X, id = id,
#'                            Z = Z, family = 'continuous',
#'                            maxiter = 500, epsilon = 0.001,
#'                            printrange = TRUE, alpadj = TRUE,
#'                            shrink = 'ALPHA', makevone = FALSE)
#' print(est_maee_ind_con)
#' 
#' 
#' ### GEE
#' est_uee_ind_con = geemaee(y = sampleSWCRT$y_con, X = X, id = id,
#'                           Z = Z, family = 'continuous',
#'                           maxiter = 500, epsilon = 0.001,
#'                           printrange = TRUE, alpadj = FALSE,
#'                           shrink = 'ALPHA', makevone = FALSE)
#' print(est_uee_ind_con)
#'
#' ###############################################################
#' ### (2) Matrix-adjusted estimating equations and GEE
#' ### on binary outcome with nested exchangeable correlation structure
#' ###############################################################
#'
#' ### MAEE
#' est_maee_ind_bin = geemaee(y = sampleSWCRT$y_bin, X = X, id = id,
#'                            Z = Z, family = 'binomial',
#'                            maxiter = 500, epsilon = 0.001,
#'                            printrange = TRUE, alpadj = TRUE,
#'                            shrink = 'ALPHA', makevone = FALSE)
#' print(est_maee_ind_bin)
#' 
#' 
#' ### GEE
#' est_uee_ind_bin = geemaee(y = sampleSWCRT$y_bin, X = X, id = id,
#'                           Z = Z, family = 'binomial',
#'                           maxiter = 500, epsilon = 0.001,
#'                           printrange = TRUE, alpadj = FALSE,
#'                           shrink = 'ALPHA', makevone = FALSE)
#' print(est_uee_ind_bin)
#' 
#' }
#' 
#' 
#' 
#' 
#'
#'
#'
#' @return \code{outbeta} estimates of marginal mean model parameters and standard errors with different finite-sample bias corrections.
#' The current version supports model-based standard error (MB), the sandwich standard error (BC0) extending Liang and Zeger (1986),
#' the sandwich standard errors (BC1) extending Kauermann and Carroll (2001), the sandwich standard errors (BC2) extending Mancl and DeRouen (2001),
#' and the sandwich standard errors (BC3) extending the Fay and Graubard (2001). A summary of
#' these bias-corrections can also be found in Lu et al. (2007), and Li et al. (2018).
#' @return \code{outalpha} estimates of intraclass correlation parameters and standard errors with different finite-sample bias corrections.
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


geemaee = function(y, X, id, Z, family, maxiter = 500, epsilon = 0.001, 
    printrange = TRUE, alpadj = FALSE, shrink = "ALPHA", makevone = TRUE) {

    if (family == "continuous") {

        contMAEE(y, X, id, Z, maxiter, epsilon, printrange, alpadj, 
            shrink, makevone)

    } else if (family == "binomial") {

        binMAEE(y, X, id, Z, maxiter, epsilon, printrange, alpadj, 
            shrink, makevone)

    }
}
