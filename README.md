# geeCRT: a package for implementing the bias-corrected generalized estimating equations in analyzing cluster randomized trials
Hengshi Yu, Fan Li, Paul Rathouz, Elizabeth L. Turner, John Preisser 

[[paper]](https://academic.oup.com/biostatistics/advance-article/doi/10.1093/biostatistics/kxaa056/6126172) | [[arXiv]](https://arxiv.org/abs/2101.00484) | [[R package]](https://cran.r-project.org/web/packages/geeCRT/index.html) | [[example code]](https://github.com/lifanfrank/clusterperiod_GEE)

**Maintainer**: Hengshi Yu (<hengshi@umich.edu>)

geeCRT is an R package for implementing the bias-corrected generalized estimating equations in analyzing cluster randomized trials.

Population-averaged models have been increasingly used in the design and analysis of cluster randomized trials (CRTs). To facilitate the applications of population-averaged models in CRTs, we implement the generalized estimating equations (GEE) and matrix-adjusted estimating equations (MAEE) approaches to jointly estimate the marginal mean models correlation models both for general CRTs and stepped wedge CRTs.

Despite the general GEE/MAEE approach, we also implement a fast cluster-period GEE method specifically for stepped wedge CRTs with large and variable cluster-period sizes. The individual-level GEE/MAEE approach becomes computationally infeasible in this setting due to inversion of high-dimensional covariance matrices and the enumeration of a high-dimensional design matrix for the correlation estimation equations. The package gives a simple and efficient estimating equations approach based on the cluster-period means to estimate the intervention effects as well as correlation parameters.

In addition, the package also provides functions for generating correlated binary data with specific mean vector and correlation matrix based on the multivariate probit method (Emrich and Piedmonte, 1991) or the conditional linear family method (Qaqish, 2003). These two functions facilitate generating correlated binary data in future simulation studies. 


## Functions and references
The geeCRT package constains four main functions. In the analysis of individual-level CRT data, users can use the `geemaee()` function to perform joint estimation of the marginal mean model and intraclass correlation parameters. In the analysis of cross-sectional stepped wedge CRTs, users can use the `cpgeeSWD()` function to perform computationally efficient joint estimation of the marginal mean and intraclass correlation parameters simply based on the cluster-period means. For generating binary data with specified mean and correlation structures, user can use the `simbinPROBIT()` function with the multivariate probit method or the `simbinCLF()` function with the conditional linear family method. 


1. geemaee function: GEE and MAEE for estimating the marginal mean and correlation parameters in CRTs
    * Liang, K. Y., Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1), 13-22.
    * Prentice, R. L. (1988). Correlated binary regression with covariates specific to each binary observation. Biometrics, 1033-1048.
    * Zhao, L. P., Prentice, R. L. (1990). Correlated binary regression using a quadratic exponential model. Biometrika, 77(3), 642-648.
    * Prentice, R. L., & Zhao, L. P. (1991). Estimating equations for parameters in means and covariances of multivariate discrete and continuous responses. Biometrics, 825-839.
    * Sharples, K., & Breslow, N. (1992). Regression analysis of correlated binary data: some small sample results for the estimating equation approach. Journal of Statistical Computation and Simulation, 42(1-2), 1-20.
    * Mancl, L. A., DeRouen, T. A. (2001). A covariance estimator for GEE with improved small sample properties. Biometrics, 57(1), 126-134.
    * Kauermann, G., Carroll, R. J. (2001). A note on the efficiency of sandwich covariance matrix estimation. Journal of the American Statistical Association, 96(456), 1387-1396.
    * Fay, M. P., Graubard, B. I. (2001). Small sample adjustments for Wald type tests using sandwich estimators. Biometrics, 57(4), 1198-1206.
    * Lu, B., Preisser, J. S., Qaqish, B. F., Suchindran, C., Bangdiwala, S. I., Wolfson, M. (2007). A comparison of two bias corrected covariance estimators for generalized estimating equations. Biometrics, 63(3), 935-941.
    * Preisser, J. S., Lu, B., Qaqish, B. F. (2008). Finite sample adjustments in estimating equations and covariance estimators for intracluster correlations. Statistics in Medicine, 27(27), 5764-5785.
    * Li, F., Turner, E. L., & Preisser, J. S. (2018). Sample size determination for GEE analyses of stepped wedge cluster randomized trials. Biometrics, 74(4), 1450-1458.
    * Li, F., Forbes, A. B., Turner, E. L., Preisser, J. S. (2019). Power and sample size requirements for GEE analyses of cluster randomized crossover trials. Statistics in Medicine, 38(4), 636-649.
    * Li, F. (2020). Design and analysis considerations for cohort stepped wedge cluster randomized trials with a decay correlation structure. Statistics in medicine, 39(4), 438-455.
    * Li, F., Yu, H., Rathouz, P., Turner, E. L., Preisser, J. S. (2020+). Marginal modeling of cluster period means and intraclass correlations in stepped wedge designs with binary outcomes. Under Revision at Biostatistics.

2. cpgeeSWD function: cluster-period generalized estimating equations for estimating the marginal mean and correlation parameters in cross-sectional stepped wedge CRTs
    * Zhao, L. P., Prentice, R. L. (1990). Correlated binary regression using a quadratic exponential model. Biometrika, 77(3), 642-648.
    * Mancl, L. A., DeRouen, T. A. (2001). A covariance estimator for GEE with improved small sample properties. Biometrics, 57(1), 126-134.
    * Kauermann, G., Carroll, R. J. (2001). A note on the efficiency of sandwich covariance matrix estimation. Journal of the American Statistical Association, 96(456), 1387-1396.
    * Fay, M. P., Graubard, B. I. (2001). Small sample adjustments for Wald type tests using sandwich estimators. Biometrics, 57(4), 1198-1206.
    * Lu, B., Preisser, J. S., Qaqish, B. F., Suchindran, C., Bangdiwala, S. I., Wolfson, M. (2007). A comparison of two bias corrected covariance estimators for generalized estimating equations. Biometrics, 63(3), 935-941.
    * Preisser, J. S., Lu, B., Qaqish, B. F. (2008). Finite sample adjustments in estimating equations and covariance estimators for intracluster correlations. Statistics in Medicine, 27(27), 5764-5785.
    * Li, F., Turner, E. L., & Preisser, J. S. (2018). Sample size determination for GEE analyses of stepped wedge cluster randomized trials. Biometrics, 74(4), 1450-1458.
    * Li, F. (2020). Design and analysis considerations for cohort stepped wedge cluster randomized trials with a decay correlation structure. Statistics in medicine, 39(4), 438-455.
    * Li, F., Yu, H., Rathouz, P., Turner, E. L., Preisser, J. S. (2020+). Marginal modeling of cluster period means and intraclass correlations in stepped wedge designs with binary outcomes. Under Revision at Biostatistics.

3. simbinPROBIT function: generating correlated binary data using the multivariate probit method

    * Emrich, L. J., & Piedmonte, M. R. (1991). A method for generating high-dimensional multivariate binary variates. The American Statistician, 45(4), 302-304.
    * Preisser, J. S., Qaqish, B. F. (2014). A comparison of methods for simulating correlated binary variables with specified marginal means and correlations. Journal of Statistical Computation and Simulation, 84(11), 2441-2452.


4. simbinCLF function: generating correlated binary data using the conditional linear family method
    * Qaqish, B. F. (2003). A family of multivariate binary distributions for simulating correlated binary variables with specified marginal means and correlations. Biometrika, 90(2), 455-463.
    * Preisser, J. S., Qaqish, B. F. (2014). A comparison of methods for simulating correlated binary variables with specified marginal means and correlations. Journal of Statistical Computation and Simulation, 84(11), 2441-2452.

## Installation

The `geeCRT` R package is available on CRAN.
```
install.packages('geeCRT')
```

### `geemaee()` example: matrix-adjusted GEE for estimating the mean and correlation parameters in CRTs

The `geemaee()` function implements the matrix-adjusted GEE or regular GEE developed for analyzing cluster randomized trials (CRTs). It provides valid estimation and inference for the treatment effect and intraclass correlation parameters within the population-averaged modeling framework. The program allows for flexible marginal mean model specifications. The program also offers bias-corrected intraclass correlation coefficient (ICC) estimates as well as bias-corrected sandwich variances for both the treatment effect parameter and the ICC parameters. The technical details of the matrix-adjusted GEE approach are provided in Preisser et al. (2008) and Li et al. (2018).

For the individual-level data, we use the `geemaee()` function to estimate the marginal mean and correlation parameters in CRTs. We use two simulated stepped wedge CRT datasets with true nested exchangeable correlation structure to illustrate the `geemaee()` function examples. We first create an auxiliary function `createzCrossSec()` to help create the design matrix for the estimating equations of the correlation parameters. We then collect design matrix `X` for the mean parameters with five period indicators and the treatment indicator. 

We implement the `geemaee()` function on both the continuous outcome and binary outcome, and consider both matrix-adjusted estimating equations (MAEE) with `alpadj = TRUE` and uncorrected generalized estimating equations (GEE) with `alpadj = FALSE`. For the `shrink` argument, we use the `"ALPHA"` method to tune step sizes and focus on using estimated variances in the correlation estimating equations rather than using unit variances by specifying `makevone = FALSE`. 


```r

### function to create the design matrix for correlation parameters 
### under the nested exchangeable correlation structure of SW-CRTs
createzCrossSec = function (m) {

    Z = NULL
    n = dim(m)[1]
    
    for (i in 1:n) {
        
        alpha_0 = 1; alpha_1 = 2; n_i = c(m[i, ]); n_length = length(n_i)
        POS = matrix(alpha_1, sum(n_i), sum(n_i))
        loc1 = 0; loc2 = 0
        
        for (s in 1:n_length) {
            
            n_t = n_i[s]; loc1 = loc2 + 1; loc2 = loc1 + n_t - 1
            
            for (k in loc1:loc2) {

                for (j in loc1:loc2) {

                    if (k != j) { POS[k, j] = alpha_0 } else { POS[k, j] = 0 }}}}

        zrow = diag(2); z_c = NULL
        
        for (j in 1:(sum(n_i) - 1)) { 

            for (k in (j + 1):sum(n_i)) {z_c = rbind(z_c, zrow[POS[j,k],])}}
        
        Z = rbind(Z, z_c) }

    return(Z)}

########################################################################
### Example 1): simulated SW-CRT with smaller cluster-period sizes (5~10)
########################################################################

sampleSWCRT = sampleSWCRTSmall

### Individual-level id, period, outcome, and design matrix
id = sampleSWCRT$id; period =  sampleSWCRT$period;
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'treatment')])
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]

### design matrix for correlation parameters
Z = createzCrossSec(m) 

### (1) Matrix-adjusted estimating equations and GEE 
### on continous outcome with nested exchangeable correlation structure
 
### MAEE
est_maee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                           X = X, id  = id, Z = Z, 
                           family = "continuous", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_maee_ind_con)

### GEE
est_uee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                          X = X, id = id, Z = Z, 
                          family = "continuous", 
                          maxiter = 500, epsilon = 0.001, 
                          printrange = TRUE, alpadj = FALSE, 
                          shrink = "ALPHA", makevone = FALSE)
print(est_uee_ind_con)

### (2) Matrix-adjusted estimating equations and GEE 
### on binary outcome with nested exchangeable correlation structure

### MAEE
est_maee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                           X = X, id = id, Z = Z, 
                           family = "binomial", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_maee_ind_bin)

### GEE
est_uee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                          X = X, id = id, Z = Z, 
                          family = "binomial", 
                          maxiter = 500, epsilon = 0.001, 
                          printrange = TRUE, alpadj = FALSE, 
                          shrink = "ALPHA", makevone = FALSE)
print(est_uee_ind_bin)


########################################################################
### Example 2): simulated SW-CRT with larger cluster-period sizes (20~30)
########################################################################
## This will elapse longer. 
sampleSWCRT = sampleSWCRTLarge

### Individual-level id, period, outcome, and design matrix
id = sampleSWCRT$id; period =  sampleSWCRT$period;
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'period5', 'treatment')])
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
### design matrix for correlation parameters
Z = createzCrossSec(m) 

### (1) Matrix-adjusted estimating equations and GEE 
### on continous outcome with nested exchangeable correlation structure
 
### MAEE
est_maee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                           X = X, id  = id, Z = Z, 
                           family = "continuous", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_maee_ind_con)

### GEE
est_uee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                          X = X, id = id, Z = Z, 
                          family = "continuous", 
                          maxiter = 500, epsilon = 0.001, 
                          printrange = TRUE, alpadj = FALSE, 
                          shrink = "ALPHA", makevone = FALSE)
print(est_uee_ind_con)

### (2) Matrix-adjusted estimating equations and GEE 
### on binary outcome with nested exchangeable correlation structure

### MAEE
est_maee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                           X = X, id = id, Z = Z, 
                           family = "binomial", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_maee_ind_bin)

### GEE
 est_uee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                           X = X, id = id, Z = Z, 
                           family = "binomial", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = FALSE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_uee_ind_bin)



```

### `cpgeeSWD()` example: cluster-period GEE for estimating the marginal mean and correlation parameters in cross-sectional SW-CRTs

The `cpgeeSWD()` function implements the cluster-period GEE developed for cross-sectional stepped wedge cluster randomized trials (SW-CRTs). It provides valid estimation and inference for the treatment effect and intraclass correlation parameters within the GEE framework, and is computationally efficient for analyzing SW-CRTs with large cluster sizes. The program currently only allows for a marginal mean model with discrete period effects and the intervention indicator without additional covariates. The program offers bias-corrected ICC estimates as well as bias-corrected sandwich variances for both the treatment effect parameter and the ICC parameters. The technical details of the cluster-period GEE approach are provided in Li et al. (2021).

We summarize the individual-level simulated SW-CRT data to cluster-period data and use the `cpgeeSWD()` function to estimate the marginal mean and correlation parameters on cluster-period means of binary outcome. We first transform the variables to get the cluster-period mean outcome `y_cp`, mean parameters' design matrix `X_cp` as well as other arguments.  

We implement the `cpgeeSWD()` function on all the three choices of the correlation structure including `"exchangeable"`, `"nest_exch"` and `"exp_decay"`. We consider both matrix-adjusted estimating equations (MAEE) with `alpadj = TRUE` and uncorrected generalized estimating equations (GEE) with `alpadj = FALSE`. 
```r

# Simulated SW-CRT example with binary outcome

########################################################################
### Example 1): simulated SW-CRT with smaller cluster-period sizes (5~10)
########################################################################

sampleSWCRT = sampleSWCRTSmall

### cluster-period id, period, outcome, and design matrix
### id, period, outcome
id = sampleSWCRT$id; period =  sampleSWCRT$period; y =  sampleSWCRT$y_bin
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'treatment')])
 
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
clp_mu = tapply(y,list(id,period), FUN=mean)
y_cp = c(t(clp_mu))
 
### design matrix for correlation parameters
trt = tapply(X[, t + 1], list(id, period), FUN=mean)
trt = c(t(trt))

time = tapply(period,list(id, period), FUN = mean); time = c(t(time))
X_cp = matrix(0, n * t, t)

s = 1
for (i in 1:n) { for (j in 1:t) { X_cp[s, time[s]] = 1; s = s + 1 }}
X_cp = cbind(X_cp, trt); id_cp = rep(1:n, each= t); m_cp =  c(t(m))

### cluster-period matrix-adjusted estimating equations (MAEE) 
### with exchangeable, nested exchangeable and exponential decay correlation structures 
# exponential
est_maee_exc = cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exchangeable", 
                        alpadj = TRUE)
print(est_maee_exc)

# nested exchangeable
est_maee_nex = cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "nest_exch", 
                        alpadj = TRUE)
print(est_maee_nex)

# exponential decay 
est_maee_ed = cpgeeSWD(y  = y_cp, X = X_cp, id = id_cp, 
                       m = m_cp, corstr = "exp_decay", 
                       alpadj = TRUE)
print(est_maee_ed)
 

### cluster-period GEE 
### with exchangeable, nested exchangeable and exponential decay correlation structures

# exchangeable
est_uee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exchangeable",
                        alpadj = FALSE)
print(est_uee_exc)

# nested exchangeable
est_uee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "nest_exch", 
                        alpadj = FALSE)
print(est_uee_nex)

# exponential decay 
est_uee_ed <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                       m = m_cp, corstr = 'exp_decay', 
                       alpadj = FALSE)
print(est_uee_ed)

########################################################################
### Example 2): simulated SW-CRT with larger cluster-period sizes (20~30)
########################################################################

sampleSWCRT = sampleSWCRTLarge

### cluster-period id, period, outcome, and design matrix
### id, period, outcome
id = sampleSWCRT$id; period =  sampleSWCRT$period; y =  sampleSWCRT$y_bin
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'period5', 'treatment')])
 
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
clp_mu<-tapply(y,list(id,period), FUN=mean)
y_cp <- c(t(clp_mu))
 
### design matrix for correlation parameters
trt <- tapply(X[, t + 1], list(id, period), FUN=mean)
trt <- c(t(trt))

time <- tapply(period,list(id, period), FUN = mean); time <- c(t(time))
X_cp <- matrix(0, n * t, t)

s = 1
for(i in 1:n){for(j in 1:t){X_cp[s, time[s]] <- 1; s = s + 1}}
X_cp <- cbind(X_cp, trt); id_cp <- rep(1:n, each= t); m_cp <-  c(t(m))

### cluster-period matrix-adjusted estimating equations (MAEE) 
### with exchangeable, nested exchangeable and exponential decay correlation structures 
# exponential
est_maee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                         m = m_cp, corstr = "exchangeable", 
                         alpadj = TRUE)
print(est_maee_exc)

# nested exchangeable
est_maee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                         m = m_cp, corstr = "nest_exch", 
                         alpadj = TRUE)
print(est_maee_nex)

# exponential decay 
est_maee_ed <- cpgeeSWD(y  = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exp_decay", 
                        alpadj = TRUE)
print(est_maee_ed)
 

### cluster-period GEE 
### with exchangeable, nested exchangeable and exponential decay correlation structures

# exchangeable
est_uee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exchangeable",
                        alpadj = FALSE)
print(est_uee_exc)

# nested exchangeable
est_uee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "nest_exch", 
                        alpadj = FALSE)
print(est_uee_nex)

# exponential decay 
est_uee_ed <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                       m = m_cp, corstr = 'exp_decay', 
                       alpadj = FALSE)
print(est_uee_ed)

```



### `simbinPROBIT()` example: generating correlated binary data using the multivariate probit method

The `simbinPROBIT()` function generates correlated binary data using the multivariate Probit method (Emrich and Piedmonte, 1991). It simulates a vector of binary outcomes according the specified marginal mean vector and correlation structure. Constraints and compatibility between the marginal mean and correlation matrix are checked.

We use the `simbinPROBIT()` function to generate correlated binary data with different correlation structures. We consider simulating a cross-sectional SW-CRT dataset with 2 clusters, 3 periods with the same cluster-period size of 5. We use two mean vectors for the two clusters and specify the `mu` argument.  

For the exchangeable correlation structure, we specify both the within-period and inter-period correlation parameters to be `0.015`. We use `0.03` and `0.015` for the within-period and inter-period correlations, respectively. The exponential decay correlation structure has an decay parameter `0.8` with the within-period correlation parameter `0.03`. 


```r
#### Simulate 2 clusters, 3 periods and cluster-period size of 5

t = 3; n = 2; m = 5
# means of cluster 1
u_c1 = c(0.4, 0.3, 0.2)
u1 <- rep(u_c1, c(rep(m, t)))
# means of cluster 2
u_c2 = c(0.35, 0.25, 0.2)
u2 <- rep(u_c2, c(rep(m, t)))

# List of mean vectors
mu = list(); mu[[1]] = u1; mu[[2]] = u2;

# List of correlation matrices

## correlation parameters
alpha0 = 0.03; alpha1 = 0.015; rho = 0.8

## (1) exchangeable
Sigma = list()
Sigma[[1]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )

Sigma[[2]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )

y_exc = simbinPROBIT(mu = mu, Sigma = Sigma, n = n)

## (2) nested exchangeable
Sigma = list()
cor_matrix = matrix(alpha1, m * t,  m * t)
loc1 = 0; loc2 = 0
for(t in 1:t){loc1 = loc2 + 1; loc2 = loc1 + m - 1
  for(i in loc1:loc2){for(j in loc1:loc2){
         if(i != j){cor_matrix[i, j] = alpha0}else{cor_matrix[i, j] = 1}}}}
 
Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix

y_nex = simbinPROBIT(mu = mu, Sigma = Sigma, n = n)

## (3) exponential decay
 
Sigma = list()
 
### function to find the period of the ith index
 region_ij<-function(points, i){diff = i - points
     for(h in 1:(length(diff) - 1)){if(diff[h] > 0 & diff[h + 1] <= 0){find <- h}}
  return(find)}

 cor_matrix = matrix(0,  m * t,  m * t)
 useage_m = cumsum(m * t); useage_m = c(0, useage_m)

 for(i in 1:(m * t)){i_reg = region_ij(useage_m, i)
      for(j in 1:(m * t)){j_reg = region_ij(useage_m, j)
          if(i_reg == j_reg & i != j){
              cor_matrix[i, j] = alpha0}else if(i == j){cor_matrix[i, j] = 1
 }else if(i_reg != j_reg){cor_matrix[i,j] = alpha0 * (rho^(abs(i_reg - j_reg)))}}}

Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix

y_ed = simbinPROBIT(mu = mu, Sigma = Sigma, n = n)

```


### `simbinCLF()` example: generating correlated binary data using the conditional linear family method

The `simbinCLF()` function generates correlated binary data using the conditional linear family method (Qaqish, 2003). It simulates a vector of binary outcomes according the specified marginal mean vector and correlation structure. Natural constraints and compatibility between the marginal mean and correlation matrix are checked.

We use the `simbinCLF()` function to generate correlated binary data with different correlation structures. We consider simulating a cross-sectional SW-CRT dataset with 2 clusters, 3 periods with the same cluster-period size of 5. We use two mean vectors for the two clusters and specify the `mu` argument.  

For the exchangeable correlation structure, we specify both the within-period and inter-period correlation parameters to be `0.015`. We use `0.03` and `0.015` for the within-period and inter-period correlations, respectively. The exponential decay correlation structure has an decay parameter `0.8` with the within-period correlation parameter `0.03`. 


```r
##### Simulate 2 clusters, 3 periods and cluster-period size of 5

t = 3; n = 2; m = 5
 
# means of cluster 1
u_c1 = c(0.4, 0.3, 0.2)
u1 <- rep(u_c1, c(rep(m, t)))
# means of cluster 2
u_c2 = c(0.35, 0.25, 0.2)
u2 <- rep(u_c2, c(rep(m, t)))

# List of mean vectors
mu = list()
mu[[1]] = u1; mu[[2]] = u2;

# List of correlation matrices
 
## correlation parameters
alpha0 = 0.03; alpha1 = 0.015; rho = 0.8

## (1) exchangeable
Sigma = list()
Sigma[[1]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )
Sigma[[2]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )
y_exc = simbinCLF(mu = mu, Sigma = Sigma, n = n)

## (2) nested exchangeable
Sigma = list()
cor_matrix = matrix(alpha1, m * t,  m * t)
loc1 = 0; loc2 = 0
for(t in 1:t){loc1 = loc2 + 1; loc2 = loc1 + m - 1
    for(i in loc1:loc2){for(j in loc1:loc2){
         if(i != j){cor_matrix[i, j] = alpha0}else{cor_matrix[i, j] = 1}}}}
 
Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix
y_nex = simbinCLF(mu = mu, Sigma = Sigma, n = n)

## (3) exponential decay
 
Sigma = list()

### function to find the period of the ith index
region_ij<-function(points, i){diff = i - points
     for(h in 1:(length(diff) - 1)){if(diff[h] > 0 & diff[h + 1] <= 0){find <- h}}
  return(find)}

cor_matrix = matrix(0,  m * t,  m * t)
useage_m = cumsum(m * t); useage_m = c(0, useage_m)

for(i in 1:(m * t)){i_reg = region_ij(useage_m, i)
      for(j in 1:(m * t)){j_reg = region_ij(useage_m, j)
          if(i_reg == j_reg & i != j){
              cor_matrix[i, j] = alpha0}else if(i == j){cor_matrix[i, j] = 1
 }else if(i_reg != j_reg){cor_matrix[i,j] = alpha0 * (rho^(abs(i_reg - j_reg)))}}}

Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix

y_ed = simbinCLF(mu = mu, Sigma = Sigma, n = n)

```





