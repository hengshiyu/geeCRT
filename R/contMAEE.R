# Continuous responses: GEE and Matrix-adjusted Estimating Equations (MAEE)
# for Estimating the Marginal Mean and Correlation Parameters in CRTs
# Args:
#  y: a vector specifying the outcome variable across all clusters
#  X: design matrix for the marginal mean model, including the intercept
#  id: a vector specifying cluster identifier
#  Z: design matrix for the correlation model, should be all pairs j < k
#    for each cluster
#  link: a specification for the model link function
#  maxiter: maximum number of iterations for Fisher scoring updates
#  epsilon: tolerance for convergence. The default is 0.001
#  printrange: print details of range violations. The default is TRUE
#  alpadj: if TRUE, performs bias adjustment for the correlation
#    estimating equations. The default is FALSE
#  shrink: method to tune step sizes in case of non-convergence including
#    'THETA' or 'ALPHA'. The default is 'ALPHA'
#  makevone if TRUE, it assumes unit variances for the correlation
#    parameters in the correlation estimating equations. The default is
#    TRUE
# Returns:
#  beta: a vector of estimates for marginal mean model parameters
#  alpha: a vector of estimates of correlation parameters
#  MB: model-based covariance estimate for the marginal mean model parameters
#  BC0: alpha robust sandwich covariance estimate of the marginal mean model
#       and correlation parameters
#  BC1: robust sandwich covariance estimate of the marginal mean model and
#       correlation parameters with the Kauermann and Carroll (2001) correction
#  BC2: robust sandwich covariance estimate of the marginal mean model and
#       correlation parameters with the Mancl and DeRouen (2001) correction
#  BC3: robust sandwich covariance estimate of the marginal mean model and
#       correlation parameters with the Fay and Graubard (2001) correction
#  niter: number of iterations in the Fisher scoring updates for model fitting


continuous_maee <- function(y, X, id, Z, link, maxiter, epsilon, printrange,
                            alpadj, shrink, makevone) {
  
  
  # ginv_scaleup calculate the generalized inverse matrix
  # Args:
  #  input_matrix: matrix for which the Moore-Penrose inverse is required
  #  threshold: threshold limit to decide no scaling
  #  tolerance: tolerance to define positive singular values
  #  maxiter: maximum number of iterations
  # Returns:
  #  Moore-Penrose inverse
  ginv_scaleup <- function(input_matrix, threshold = 0.01, tolerance = 1e-13, maxiter = 20) {
    summary_stat = as.vector(summary(as.vector(input_matrix)))
    low_range_stat = summary_stat[3:1]
    # min, q1, median, mean, q3, max
    scale_factor = 1
    ginv_matrix = tryCatch(ginv(input_matrix  * scale_factor) * scale_factor, error = function(e) e)
    error_message_condition = inherits(tryCatch(ginv(input_matrix * scale_factor) * scale_factor, error = function(e) e),   "error")
    
    niter_ginv = 0
    value_for_scale = low_range_stat[1]
    while (error_message_condition & niter_ginv < maxiter) {
      scale_factor = 1 / (ceiling(value_for_scale) + 1e-8)
      ginv_matrix = tryCatch(ginv(input_matrix  * scale_factor) * scale_factor, error = function(e) e)
      error_message_condition = inherits(tryCatch(ginv(input_matrix * scale_factor) * scale_factor, error = function(e) e),   "error")
      
      niter_ginv = niter_ginv + 1
      if (niter_ginv < 3) {
        value_for_scale = low_range_stat[niter_ginv + 1]
      } else {
        value_for_scale = value_for_scale * 10
      }
      
    }
    if (error_message_condition) {
      return (matrix(0, nrow = ncol(input_matrix), ncol = nrow(input_matrix)))
    } else {
      return(ginv_matrix)
    }
  }

  # begin_end
  # Args:
  #  n: Vector of cluster sample sizes
  # Returns:
  #  matrix with two columns.
  #  first: a vector with starting row for cluster i,
  #  last: a vector with ending row for cluster i
  begin_end <- function(n) {
    last <- cumsum(n)
    first <- last - n + 1
    return(cbind(first, last))
  }

  # is_pos_def
  # Args:
  #  matrix: symmetric matrix
  # Returns:
  #  returns 1 if matrix is positive definite else 0
  is_pos_def <- function(matrix) {
    return(min(eigen(matrix)$values) > 1e-13)
  }

  # get_corr_beta_deriv calculates the derivative of correlation matrix to
  #  marginal mean covariates
  # Args:
  #  mu: marginal means for the cluster
  #  gamma: pairwise correlation for outcome j and k for the cluster
  #  j: indicator for the mean
  #  k: indicator for the mean
  #  X: covariate matrix for the cluster
  #  y: response vector for the cluster
  #  link: a specification for the model link function
  # Returns:
  #  derivative of the empirical correlation over marginal mean covariates
  get_corr_beta_deriv <- function(mu, gamma, j, k, X, y, link) {
    beta_deriv_j <- create_beta_derivative(X[j, ], mu[j], link)
    beta_deriv_k <- create_beta_derivative(X[k, ], mu[k], link)

    row <- gamma * (beta_deriv_j / (y[j] - mu[j]) +
      beta_deriv_k / (y[k] - mu[k]))
    return(row)
  }

  # create_beta_residual creates residual for beta estimating equation
  # Args:
  #  mu: marginal mean vector
  #  y: outcome vector
  # Returns:
  #  residual of marginal outcome
  create_beta_residual <- function(mu, y) {
    return(y - mu)
  }

  # marginal_mean_estimation creates mean estimation
  # Args:
  #  X: covariate matrix for cluster
  #  beta: vector of marginal mean parameters
  #  link: a specification for the model link function
  # Returns:
  #  marginal mean vector
  marginal_mean_estimation <- function(X, beta, link) {
    if (link == "identity") {
      mu <- c(X %*% beta)
    } else if (link == "log") {
      mu <- exp(c(X %*% beta))
    }
    return(mu)
  }


  # create_beta_derivative creates derivative for beta estimating equation
  # Args:
  #  X: covariate matrix for cluster
  #  mu: marginal mean vector
  #  link: a specification for the model link function
  # Returns:
  #  derivative of beta estimating equation
  create_beta_derivative <- function(X, mu, link) {
    if (link == "identity") {
      return(X)
    } else if (link == "log") {
      return(X * mu)
    }
  }

  # create_beta_covariance creates covariance for beta estimating equation
  # Args:
  #  phi: Dispersion parameter
  #  n: sample size for the cluster
  #  gamma: working correlation mean vector
  # Returns:
  #  working covariance matrix for beta estimating equation
  create_beta_covariance <- function(phi, n, gamma = NULL) {
    if (is.null(gamma)) {
      return(as.matrix(phi))
    }

    B <- diag(rep(phi, n))
    l <- 1
    for (j in 1:(n - 1)) {
      for (k in (j + 1):n) {
        B[j, k] <- phi * gamma[l]
        l <- l + 1
      }
    }
    B[lower.tri(B)] <- t(B)[lower.tri(B)]
    return(B)
  }

  # phi_estimate generates moment estimates for dispersion parameter
  # Args:
  #  Ustar: information matrix
  #  beta: vector of marginal mean parameters
  #  alpha: vector of marginal correlation parameters
  #  y: the continuous outcome
  #  X: marginal mean covariates
  #  Z: marginal correlation covariates
  #  n: vector of cluster sample sizes
  #  p: number of marginal mean parameters
  #  q: number of marginal correlation parameters
  #  phi: dispersion parameter
  #  link: a specification for the model link function
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix is
  #                           positive definite
  # Returns:
  #  phi: dispersion parameter for beta estimating equation
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix is
  #                           positive definite
  phi_estimate <- function(Ustar, beta, alpha, y, X, Z, n, p, q,
                           phi, link, not_pos_def_alpadj_flag) {
    resid_sum_square <- 0
    naive_inv_prev <- ginv_scaleup(Ustar[1:p, 1:p])

    loc_x <- begin_end(n)
    loc_z <- begin_end(choose(n, 2))

    for (i in 1:length(n)) {
      X_c <- X[loc_x[i, 1]:loc_x[i, 2], , drop = FALSE]
      y_c <- y[loc_x[i, 1]:loc_x[i, 2]]
      mu_c <- marginal_mean_estimation(X_c, beta, link)

      beta_deriv <- create_beta_derivative(X_c, mu_c, link)
      beta_resid <- create_beta_residual(mu_c, y_c)

      if (loc_x[i, 1] == loc_x[i, 2]) {
        beta_work_cov <- create_beta_covariance(phi, n[i])
      } else {
        Z_c <- Z[loc_z[i, 1]:loc_z[i, 2], , drop = FALSE]
        gamma_c <- c(Z_c %*% alpha)
        beta_work_cov <- create_beta_covariance(phi, n[i], gamma_c)
      }

      square_var <- sqrt(rep(phi, n[i]))
      inv_square_var <- 1 / square_var
      beta_deriv_trans <- t(beta_deriv)

      omega <- beta_deriv %*% naive_inv_prev %*% beta_deriv_trans
      v_min_omega <- beta_work_cov - omega
      psd_vmin <- is_pos_def(v_min_omega)

      if (psd_vmin == 1) {
        Ci <- beta_work_cov %*% ginv_scaleup(v_min_omega)
        Rx <- (y_c - mu_c) * inv_square_var
        Gi <- tcrossprod(Rx)
      } else {
        not_pos_def_alpadj_flag <- 1
        stop("(V - Omega) is not positive definite")
      }

      # moment-based estimation of dispersion (total variance)
      for (j in 1:n[i]) {
        resid_sum_square <- resid_sum_square + phi * (Ci[j, ] %*% Gi[, j])
      }
    }
    phi <- resid_sum_square / (sum(n) - p)
    return(list(phi = phi, not_pos_def_alpadj_flag = not_pos_def_alpadj_flag))
  }

  # score_function generates the score matrix for each cluster and
  # approximate information to be used to estimate parameters and
  # generate standard errors
  # Args:
  #  Ustar_prev: initial values for information matrix
  #  beta: vector of marginal mean parameters
  #  alpha: vector of marginal correlation parameters
  #  y: the continuous outcome
  #  X: marginal mean covariates
  #  Z: marginal correlation covariates
  #  n: vector of cluster sample sizes
  #  p: number of marginal mean parameters
  #  q: number of marginal correlation parameters
  #  phi: dispersion parameter
  #  link: a specification for the model link function
  #  flag: performs an eigen-analysis of covariance to see if positive
  #        definite. only called when computing the variance at the
  #        end (0, 1). Prints warning for each cluster violation.
  #  range_flag: checks to see if correlation is within range based on
  #              marginal means (0 = in range, 1 = out of range).
  #              See Prentice (1988).
  #  alpha_work_cov_flag: checks if all variances for alpha estimating
  #                       equations are positive, terminates if not
  #  not_pos_def_work_cov_flag: checks if working covariance is positive
  #                             definite, terminates if not
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix is
  #                           positive definite, terminates if not
  # Returns:
  #  U: score vector
  #  UUtran: sum of U_i*U_i` across all clusters
  #  Ustar: approximated information matrix
  #  flag: Performs an eigen-analysis of covariance to see if positive
  #        definite. only called when computing the variance at the
  #        end (0, 1). Prints warning for each cluster violation.
  #  range_flag: checks to see if correlation is within range based on
  #              marginal means (0 = in range, 1 = out of range).
  #              See Prentice (1988).
  #  alpha_work_cov_flag: checks if all variances for alpha estimating
  #                       equations are positive, terminates if not
  #  not_pos_def_work_cov_flag: checks if working covariance is positive
  #                             definite, terminates if not
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix is
  #                           positive definite, terminates if not
  score_function <- function(Ustar_prev, beta, alpha, y, X, Z, n, p, q,
                             phi, link, flag, range_flag, alpha_work_cov_flag,
                             not_pos_def_work_cov_flag,
                             not_pos_def_alpadj_flag) {
    # score function
    U <- rep(0, p + q)
    # meat and bread (information) of the sandwich variance estimator
    UUtran <- Ustar <- matrix(0, p + q, p + q)
    naive_inv_prev <- ginv_scaleup(Ustar_prev[1:p, 1:p])

    loc_x <- begin_end(n)
    loc_z <- begin_end(choose(n, 2))

    for (i in 1:length(n)) {
      X_c <- X[loc_x[i, 1]:loc_x[i, 2], , drop = FALSE]
      y_c <- y[loc_x[i, 1]:loc_x[i, 2]]

      U_c <- rep(0, p + q)
      Ustar_c <- matrix(0, p + q, p + q)
      mu_c <- marginal_mean_estimation(X_c, beta, link)

      # case when there was only 1 observation for this cluster
      if (loc_x[i, 1] == loc_x[i, 2]) {
        beta_deriv <- create_beta_derivative(X_c, mu_c, link)
        beta_work_cov <- create_beta_covariance(phi, n[i])
        beta_resid <- create_beta_residual(mu_c, y_c)

        inv_work_cov <- ginv_scaleup(beta_work_cov)
        U_c[1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_resid
        UUtran_c <- tcrossprod(U_c)

        Ustar_c[1:p, 1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_deriv
        U <- U + U_c
        UUtran <- UUtran + UUtran_c
        Ustar <- Ustar + Ustar_c
        next
      }

      Z_c <- Z[loc_z[i, 1]:loc_z[i, 2], , drop = FALSE]
      gamma_c <- c(Z_c %*% alpha)
      # correlation means
      alpha_work_cov <- alpha_resid <- rep(0, choose(n[i], 2))

      beta_deriv <- create_beta_derivative(X_c, mu_c, link)
      beta_work_cov <- create_beta_covariance(phi, n[i], gamma_c)
      beta_resid <- create_beta_residual(mu_c, y_c)

      inv_work_cov <- ginv_scaleup(beta_work_cov)

      if (alpadj) {
        square_var <- sqrt(rep(phi, n[i]))
        inv_square_var <- 1 / square_var
        beta_deriv_trans <- t(beta_deriv)
        omega <- beta_deriv %*% naive_inv_prev %*% beta_deriv_trans
        v_min_omega <- beta_work_cov - omega
        psd_vmin <- is_pos_def(v_min_omega)

        if (psd_vmin == 1) {
          Ci <- beta_work_cov %*% ginv_scaleup(v_min_omega)
          Rx <- (y_c - mu_c) * inv_square_var
          Gi <- tcrossprod(Rx)
        } else {
          not_pos_def_alpadj_flag <- 1
          stop("(V - Omega) is not positive definite")
        }
      } else {
        square_var <- sqrt(rep(phi, n[i]))
        inv_square_var <- 1 / square_var
        Rx <- (y_c - mu_c) * inv_square_var
        Gi <- tcrossprod(Rx)
      }

      # range checks: iterate through the correlations
      # and residual and covariance matrix for the alpha
      # estimating equations
      l <- 1
      for (j in 1:(n[i] - 1)) {
        for (k in (j + 1):n[i]) {
          range_condition <- ((gamma_c[l] > 1) | (gamma_c[l] < -1))
          if (range_condition & (flag == 0)) {
            range_flag <- 1

            if (printrange) {
              warning(cat(
                "Range Violation Detected for Cluster",
                i, "and Pair", j, k, "\n"
              ))
            }
            break
          }

          if (range_condition & (flag == 1)) {
            warning(cat(
              "Last Update Pushes Parameters Out of Range.",
              "\n"
            ))
            warning(cat(
              "Range Violation Detected for Cluster",
              i, "and Pair", j, k, "\n"
            ))
          }
          if (!makevone) {
            alpha_work_cov[l] <- 1 + gamma_c[l]^2
          }

          # insert check that variance is nonnegative
          if (alpha_work_cov[l] <= 0) {
            alpha_work_cov_flag <- 1
            stop("Variance of correlation parameter is negative")
          }

          # Matrix-based multiplicative correction, (I - H_i)^{-1}
          if (alpadj) {
            alpha_resid[l] <- Ci[j, ] %*% Gi[, k] - gamma_c[l]
          } else {
            alpha_resid[l] <- Gi[j, k] - gamma_c[l]
          }

          l <- l + 1
        }
      }

      if (makevone) {
        alpha_work_cov <- rep(1, choose(n[i], 2))
      }

      # Check for positive definite of B;
      if (min(eigen(beta_work_cov)$values) <= 0) {
        not_pos_def_work_cov_flag <- 1
        stop(paste(
          "Var(Y) of Cluster", i, "is not Positive-Definite;",
          "Joint Distribution Does Not Exist and Program terminates"
        ))
      }

      U_c[1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_resid
      U_c[(p + 1):(p + q)] <- t(Z_c) %*% (alpha_resid / alpha_work_cov)

      UUtran_c <- tcrossprod(U_c)
      Ustar_c[1:p, 1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_deriv
      Ustar_c[(p + 1):(p + q), (p + 1):(p + q)] <- t(Z_c) %*%
        (Z_c / alpha_work_cov)

      U <- U + U_c
      UUtran <- UUtran + UUtran_c
      Ustar <- Ustar + Ustar_c
    }
    return(list(
      U = U,
      UUtran = UUtran,
      Ustar = Ustar,
      flag = flag,
      range_flag = range_flag,
      alpha_work_cov_flag = alpha_work_cov_flag,
      not_pos_def_work_cov_flag = not_pos_def_work_cov_flag,
      not_pos_def_alpadj_flag = not_pos_def_alpadj_flag
    ))
  }

  # initial_beta generates initial values for beta
  #  by linear regression with least squares
  # Args:
  #  y: the continuous outcome
  #  X: marginal mean covariates
  #  n: vector of cluster sample sizes
  #  link: a specification for the model link function
  # Returns:
  #  beta: vector of marginal mean parameters
  #  Ustar: approximate information matrix
  #  phi: dispersion parameter
  initial_beta <- function(y, X, n, link) {
    if (link == "identity") {
      beta <- solve(t(X) %*% X, t(X) %*% y)
      mu <- marginal_mean_estimation(X, beta, link)
      beta_deriv <- create_beta_derivative(X, mu, link)
      beta_resid <- create_beta_residual(mu, y)

      phi <- sum(beta_resid^2) / (sum(n) - ncol(X))
      Ustar <- t(X) %*% beta_deriv / phi
    } else {
      beta <- solve(t(X) %*% X, t(X) %*% y)

      for (i in 1:2) {
        mu <- marginal_mean_estimation(X, beta, link)
        beta_resid <- create_beta_residual(mu, y)
        beta_deriv <- create_beta_derivative(X, mu, link)
        phi <- sum(beta_resid^2) / (sum(n) - ncol(X))

        z <- t(X) %*% beta_resid
        Ustar <- t(X) %*% beta_deriv / phi
        d <- solve(Ustar, z)
        beta <- beta + d
      }
    }
    return(list(beta = c(beta), Ustar = Ustar, phi = phi))
  }


  # inv_big_matrix compute (A - mm`)^{-1}c without performing the
  # inverse directly.
  # Args:
  #  ainvc: inverse of matrix A times vector c
  #  ainvm: inverse of matrix A times matrix (with low no. of
  #         columns) M
  #  m: matrix of eigen column vectors m1,m2, ..
  #  c: vector
  #  start: of do loop
  #  end: of do loop, rank of X
  # Returns:
  #  ainvc: inverse of matrix A times vector c
  inv_big_matrix <- function(ainvc, ainvm, m, c, start, end) {
    for (i in start:end) {
      b <- ainvm[, i]
      bt <- t(b)
      btm <- bt %*% m
      btmi <- btm[, i]
      gam <- 1 - btmi
      bg <- b / gam
      ainvc <- ainvc + bg %*% (bt %*% c)
      if (i < end) {
        ainvm <- ainvm + bg %*% btm
      }
    }
    return(ainvc)
  }

  # variance_estimator creates covariance matrix of beta and alpha
  # Args:
  #  Ustar_prev: initial values for information matrix
  #  beta: vector of marginal mean parameters
  #  alpha: vector of marginal correlation parameters
  #  y: the continuous outcome
  #  X: marginal mean covariates
  #  Z: marginal correlation covariates
  #  n: vector of cluster sample sizes
  #  p: number of marginal mean parameters
  #  q: number of marginal correlation parameters
  #  phi: dispersion parameter
  #  link: a specification for the model link function
  #  alpha_work_cov_flag: checks if all variances for alpha estimating
  #                       equations are positive, terminates if not
  #  not_pos_var_est_flag: checks if any robust covariance matrix
  #                        contains non-positive variance
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix is
  #                           positive definite, terminates if not
  # Returns:
  #  naive: Naive (Model-Based) covariance matrix for beta
  #  robust: Robust covariance matrix for beta and alpha
  #  varMD: bias-corrected variance by Mancl and Derouen (2001)
  #  varKC: bias-corrected variance by Kauermann and Carroll (2001)
  #  varFG: bias-corrected variance by Fay and Graubard (2001)
  #  alpha_work_cov_flag: checks if all variances for alpha estimating
  #                       equations are positive, terminates if not
  #  not_pos_def_work_cov_flag: checks if working covariance is positive
  #                             definite, terminates if not
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix is
  #                           positive definite, terminates if not
  variance_estimator <- function(Ustar_prev, beta, alpha, y, X, Z, n, p, q,
                                 phi, link, alpha_work_cov_flag,
                                 not_pos_var_est_flag,
                                 not_pos_def_work_cov_flag,
                                 not_pos_def_alpadj_flag) {
    score_res <- score_function(Ustar_prev, beta, alpha, y, X, Z, n, p,
      q, phi, link,
      flag = 1, range_flag = 0, alpha_work_cov_flag,
      not_pos_def_work_cov_flag, not_pos_def_alpadj_flag
    )

    U <- score_res$U
    UUtran <- score_res$UUtran
    Ustar <- score_res$Ustar
    flag <- score_res$flag
    range_flag <- score_res$range_flag
    alpha_work_cov_flag <- score_res$alpha_work_cov_flag
    not_pos_def_work_cov_flag <- score_res$not_pos_def_work_cov_flag
    not_pos_def_alpadj_flag <- score_res$not_pos_def_alpadj_flag

    naive_beta <- ginv_scaleup(Ustar[1:p, 1:p])
    naive_alpha <- ginv_scaleup(Ustar[(p + 1):(p + q), (p + 1):(p + q)])

    # new commands to compute INV(I - H1)
    eigen_res_1 <- eigen(naive_beta)
    evals_1 <- eigen_res_1$values
    evecs_1 <- eigen_res_1$vectors
    sqrt_evals_1 <- sqrt(evals_1)
    sqrt_e1 <- evecs_1 %*% diag(sqrt_evals_1)

    # new commands to compute INV(I - H2)
    eigen_res_2 <- eigen(naive_alpha)
    evals_2 <- eigen_res_2$values
    evecs_2 <- eigen_res_2$vectors
    sqrt_evals_2 <- sqrt(evals_2)

    ## when there is only 1 alpha parameter
    if (length(sqrt_evals_2) == 1) {
      sqrt_e2 <- evecs_2 %*% sqrt_evals_2
    } else {
      sqrt_e2 <- evecs_2 %*% diag(sqrt_evals_2)
    }

    # Bias-corrected variance
    Ustar_c_array <- UUtran_c_array <- array(0, c(p + q, p + q, length(n)))
    UUtran <- UUbc <- UUbc2 <- UUbc3 <- Ustar <- inv_Ustar <- matrix(
      0, p + q, p + q
    )

    loc_x <- begin_end(n)
    loc_z <- begin_end(choose(n, 2))

    for (i in 1:length(n)) {
      X_c <- X[loc_x[i, 1]:loc_x[i, 2], , drop = FALSE]
      y_c <- y[loc_x[i, 1]:loc_x[i, 2]]
      mu_c <- marginal_mean_estimation(X_c, beta, link)
      U_i <- U_c <- rep(0, p + q)
      Ustar_c <- matrix(0, p + q, p + q)

      # case when there was only 1 observation for this cluster
      if (loc_x[i, 1] == loc_x[i, 2]) {
        beta_deriv <- create_beta_derivative(X_c, mu_c, link)
        beta_work_cov <- create_beta_covariance(mu_c, n[i])
        beta_resid <- create_beta_residual(mu_c, y_c)

        inv_work_cov <- ginv_scaleup(beta_work_cov)
        U_i[1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_resid

        ai1 <- inv_work_cov
        mm1 <- beta_deriv %*% sqrt_e1
        ai1A <- ai1 %*% beta_resid
        ai1m1 <- ai1 %*% mm1
        ai1A <- inv_big_matrix(ai1A, ai1m1, mm1, beta_resid, 1, p)
        U_c[1:p] <- t(beta_deriv) %*% ai1A
        Ustar_c[1:p, 1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_deriv

        UUtran_c <- tcrossprod(U_i)
        UUbc_c <- tcrossprod(U_c)
        UUbc_ic <- tcrossprod(U_c, U_i)

        UUtran <- UUtran + UUtran_c
        UUbc <- UUbc + UUbc_c
        UUbc2 <- UUbc2 + UUbc_ic
        Ustar <- Ustar + Ustar_c


        Ustar_c_array[, , i] <- Ustar_c
        UUtran_c_array[, , i] <- UUtran_c

        next
      }

      Z_c <- Z[loc_z[i, 1]:loc_z[i, 2], , drop = FALSE]
      gamma_c <- c(Z_c %*% alpha)

      # commands for beta
      beta_deriv <- create_beta_derivative(X_c, mu_c, link)
      beta_work_cov <- create_beta_covariance(phi, n[i], gamma_c)
      beta_resid <- create_beta_residual(mu_c, y_c)

      inv_work_cov <- ginv_scaleup(beta_work_cov)
      U_i[1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_resid

      # commands for generalized inverse - beta
      ai1 <- inv_work_cov
      mm1 <- beta_deriv %*% sqrt_e1
      ai1A <- ai1 %*% beta_resid
      ai1m1 <- ai1 %*% mm1
      ai1A <- inv_big_matrix(ai1A, ai1m1, mm1, beta_resid, 1, p)
      U_c[1:p] <- t(beta_deriv) %*% ai1A

      # commands for alpha
      alpha_work_cov <- alpha_resid <- rep(0, choose(n[i], 2))
      corr_beta_deriv <- matrix(0, choose(n[i], 2), p)

      # MAEE
      if (alpadj) {
        square_var <- sqrt(rep(phi, n[i]))
        inv_square_var <- 1 / square_var
        beta_deriv_trans <- t(beta_deriv)
        omega <- beta_deriv %*% naive_beta %*% beta_deriv_trans
        v_min_omega <- beta_work_cov - omega
        psd_vmin <- is_pos_def(v_min_omega)

        if (psd_vmin == 1) {
          Ci <- beta_work_cov %*% ginv_scaleup(v_min_omega)
          Rx <- (y_c - mu_c) * inv_square_var
          Gi <- tcrossprod(Rx)
        } else {
          not_pos_def_alpadj_flag <- 1
          stop("(V - Omega) is not positive definite")
        }
      } else {
        square_var <- sqrt(rep(phi, n[i]))
        inv_square_var <- 1 / square_var
        Rx <- (y_c - mu_c) * inv_square_var
        Gi <- tcrossprod(Rx)
      }

      # residual and covariance matrix for the alpha
      # estimating equations
      l <- 1
      for (j in 1:(n[i] - 1)) {
        for (k in (j + 1):n[i]) {
          if (!makevone) {
            alpha_work_cov[l] <- 1 + gamma_c[l]^2
          }
          corr_beta_deriv[l, ] <- get_corr_beta_deriv(
            mu_c, gamma_c[l], j, k, X_c, y_c, link
          )

          # Matrix-based multiplicative correction, (I - H_i)^{-1}
          if (alpadj) {
            alpha_resid[l] <- Ci[j, ] %*% Gi[, k] - gamma_c[l]
          } else {
            alpha_resid[l] <- Gi[j, k] - gamma_c[l]
          }

          l <- l + 1
        }
      }

      if (makevone) {
        alpha_work_cov <- rep(1, choose(n[i], 2))
      }


      U_i[(p + 1):(p + q)] <- t(Z_c) %*% (alpha_resid / alpha_work_cov)

      mm2 <- Z_c %*% sqrt_e2
      ai2R <- alpha_resid / alpha_work_cov
      ai2m2 <- mm2 / alpha_work_cov
      ai2R <- inv_big_matrix(ai2R, ai2m2, mm2, alpha_resid, 1, q)
      U_c[(p + 1):(p + q)] <- t(Z_c) %*% ai2R

      Ustar_c[1:p, 1:p] <- t(beta_deriv) %*% inv_work_cov %*% beta_deriv
      Ustar_c[(p + 1):(p + q), 1:p] <- t(Z_c) %*% (corr_beta_deriv /
        alpha_work_cov)
      Ustar_c[(p + 1):(p + q), (p + 1):(p + q)] <- t(Z_c) %*% (Z_c /
        alpha_work_cov)
      Ustar <- Ustar + Ustar_c

      UUtran_c <- tcrossprod(U_i)
      UUtran <- UUtran + UUtran_c

      UUbc_c <- tcrossprod(U_c)
      UUbc <- UUbc + UUbc_c

      UUbc_ic <- tcrossprod(U_c, U_i)
      UUbc2 <- UUbc2 + UUbc_ic

      Ustar_c_array[, , i] <- Ustar_c
      UUtran_c_array[, , i] <- UUtran_c
    }

    inv_Ustar[1:p, 1:p] <- ginv_scaleup(Ustar[1:p, 1:p])
    inv_Ustar[(p + 1):(p + q), (p + 1):(p + q)] <- ginv_scaleup(Ustar[(p +
      1):(p + q), (p + 1):(p + q)])
    inv_Ustar[(p + 1):(p + q), 1:p] <- -inv_Ustar[
      (p + 1):(p + q),
      (p + 1):(p + q)
    ] %*% Ustar[(p + 1):(p + q), 1:p] %*%
      inv_Ustar[1:p, 1:p]
    # the minus sign above is crucial, esp. for large correlation;

    inv_Ustar_trans <- t(inv_Ustar)

    # calculating adjustment factor for BC3
    for (i in 1:length(n)) {
      Hi <- diag(1 / sqrt(1 - pmin(
        0.75,
        c(diag(Ustar_c_array[, , i] %*% inv_Ustar))
      )))
      UUbc3 <- UUbc3 + Hi %*% UUtran_c_array[, , i] %*% Hi
    }

    # BC0 or usual Sandwich estimator of Prentice (1988);
    robust <- inv_Ustar %*% UUtran %*% inv_Ustar_trans

    # BC1 or Variance estimator that extends Kauermann and Carroll (2001);
    varKC <- inv_Ustar %*% (UUbc2 + t(UUbc2)) %*% inv_Ustar_trans / 2

    # BC2 or Variance estimator that extends Mancl and DeRouen (2001);
    varMD <- inv_Ustar %*% UUbc %*% inv_Ustar_trans

    # BC3 or Variance estimator that extends Fay and Graubard (2001);
    varFG <- inv_Ustar %*% UUbc3 %*% inv_Ustar_trans

    naive <- inv_Ustar[1:p, 1:p]

    if (min(diag(robust)) <= 0) {
      not_pos_var_est_flag <- 1
    }
    if (min(diag(varMD)) <= 0) {
      not_pos_var_est_flag <- 1
    }
    if (min(diag(varKC)) <= 0) {
      not_pos_var_est_flag <- 1
    }
    if (min(diag(varFG)) <= 0) {
      not_pos_var_est_flag <- 1
    }

    return(list(
      naive = naive,
      robust = robust,
      varMD = varMD,
      varKC = varKC,
      varFG = varFG,
      alpha_work_cov_flag = alpha_work_cov_flag,
      not_pos_var_est_flag = not_pos_var_est_flag,
      not_pos_def_alpadj_flag = not_pos_def_alpadj_flag
    ))
  }


  # model_fit solves generalized estimating equations for marginal mean
  # and correlation parameters
  # Args:
  #  y: the continuous outcome
  #  X: marginal mean covariates
  #  Z: marginal correlation covariates
  #  n: vector of cluster sample sizes
  #  link: a specification for the model link function
  #  maxiter: max number of iterations
  #  epsilon: tolerence for convergence
  #  alpha_work_cov_flag: algorithm terminated due to non-positive variance
  #                       of correlation parameters
  #  singular_naive_var_flag: checks if the naive covariance matrix
  #                           is not positive definite
  #  not_pos_var_est_flag: checks if any robust covariance matrix
  #                        contains non-positive variance
  #  alpha_shrink_flag: checks if the number of alpha shrinkage steps exceeds
  #                     20 in order to make the correlation within bounds
  #  not_pos_def_work_cov_flag: checks if working covariance is positive
  #                             definite, terminates if not
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix
  #                           is positive definite, terminates if not
  # Returns:
  #  beta: a p x 1 vector of marginal mean parameter estimates
  #  alpha: a q x 1 vector of marginal correlation parameter estimates
  #  robust: robust covariance matrix for beta and alpha
  #  naive: naive (model-Based) covariance matrix for beta
  #  varMD: bias-corrected variance due to Mancl and DeRouen (2001)
  #  varKC: bias-corrected variance due to Kauermann and Carroll (2001)
  #  varFG: bias-corrected variance due to Fay and Graubard (2001)
  #  niter: number of iterations required for convergence
  #  converge: whether the algorithm converged (1) or not (0)
  #  alpha_work_cov_flag: algorithm terminated due to non-positive variance
  #                       of correlation parameters
  #  singular_naive_var_flag: checks if the naive covariance matrix
  #                           is not positive definite
  #  not_pos_var_est_flag: checks if any robust covariance matrix
  #                        contains non-positive variance
  #  alpha_shrink_flag: checks if the number of alpha shrinkage steps exceeds
  #                     20 in order to make the correlation within bounds
  #  not_pos_def_work_cov_flag: checks if working covariance is positive
  #                             definite, terminates if not
  #  not_pos_def_alpadj_flag: checks if the alpha adjustment factor matrix
  #                           is positive definite, terminates if not
  model_fit <- function(y, X, Z, n, link, maxiter, epsilon, alpha_work_cov_flag,
                        singular_naive_var_flag, not_pos_var_est_flag,
                        alpha_shrink_flag, not_pos_def_work_cov_flag,
                        not_pos_def_alpadj_flag) {
    p <- ncol(X)
    q <- ncol(Z)
    delta <- rep(2 * epsilon, p + q)
    max_modify <- 20
    converge <- 0

    range_flag <- 0
    alpha <- rep(0.01, q)
    init_res <- initial_beta(y, X, n, link)
    beta <- init_res$beta
    Ustar <- init_res$Ustar
    phi <- init_res$phi

    niter <- 1
    while ((niter <= maxiter) & (max(abs(delta)) > epsilon)) {
      n_modify <- 0
      singular_naive_var_flag <- 0
      alpha_shrink_flag <- 0

      repeat {
        Ustar_prev <- Ustar
        not_pos_def_work_cov_flag <- 0
        not_pos_def_alpadj_flag <- 0

        score_res <- score_function(Ustar_prev, beta, alpha, y, X, Z,
          n, p, q, phi, link,
          flag = 0, range_flag, alpha_work_cov_flag,
          not_pos_def_work_cov_flag,
          not_pos_def_alpadj_flag
        )

        U <- score_res$U
        UUtran <- score_res$UUtran
        Ustar <- score_res$Ustar
        range_flag <- score_res$range_flag
        alpha_work_cov_flag <- score_res$alpha_work_cov_flag

        if (alpha_work_cov_flag == 1) {
          stop("Program terminated due to division by zero in variance")
        }

        if (range_flag == 1) {
          if (shrink == "THETA") {
            if (niter == 1) {
              alpha <- rep(0, q)
            } else {
              theta <- theta - (0.5)^(n_modify + 1) * delta
              beta <- theta[1:p]
              alpha <- theta[(p + 1):(p + q)]
            }
          } else if (shrink == "ALPHA") {
            if (niter == 1) {
              alpha <- rep(0, q)
            } else {
              alpha <- 0.95 * alpha
            }
          }
          n_modify <- n_modify + 1
          if (printrange) {
            warning(cat(
              "Iteration", niter, "and Shrink Number",
              n_modify, "\n"
            ))
          }
        }
        if ((n_modify > max_modify) | (range_flag == 0)) {
          break
        }
      }

      if (n_modify > max_modify) {
        if (printrange) {
          warning(cat("n_modify too large, more than 20 shrinks"))
        }
        alpha_shrink_flag <- 1
      }
      theta <- c(beta, alpha)

      psd_ustar <- is_pos_def(Ustar)
      if (psd_ustar) {
        delta <- solve(Ustar, U)
        theta <- theta + delta
        beta <- theta[1:p]
        alpha <- theta[(p + 1):(p + q)]
        phi_res <- phi_estimate(
          Ustar, beta, alpha, y, X, Z, n, p,
          q, phi, link, not_pos_def_alpadj_flag
        )
        phi <- phi_res$phi
        converge <- (max(abs(delta)) <= epsilon)
      } else {
        singular_naive_var_flag <- 1
      }
      niter <- niter + 1
    }

    Ustar_prev <- Ustar

    # inference
    var_res <- variance_estimator(
      Ustar_prev, beta, alpha, y, X, Z, n,
      p, q, phi, link, alpha_work_cov_flag,
      not_pos_var_est_flag,
      not_pos_def_work_cov_flag,
      not_pos_def_alpadj_flag
    )
    naive <- var_res$naive
    robust <- var_res$robust
    varMD <- var_res$varMD
    varKC <- var_res$varKC
    varFG <- var_res$varFG

    alpha_work_cov_flag <- var_res$alpha_work_cov_flag
    not_pos_var_est_flag <- var_res$not_pos_var_est_flag
    not_pos_def_alpadj_flag <- var_res$not_pos_def_alpadj_flag

    return(list(
      beta = beta,
      alpha = alpha,
      naive = naive,
      robust = robust,
      varMD = varMD,
      varKC = varKC,
      varFG = varFG,
      niter = niter,
      converge = converge,
      alpha_work_cov_flag = alpha_work_cov_flag,
      singular_naive_var_flag = singular_naive_var_flag,
      not_pos_var_est_flag = not_pos_var_est_flag,
      alpha_shrink_flag = alpha_shrink_flag,
      not_pos_def_work_cov_flag = not_pos_def_work_cov_flag,
      not_pos_def_alpadj_flag = not_pos_def_alpadj_flag
    ))
  }

  # beta_alpha_format_results creates printed output to screen of parameters
  #  and other information
  # Args:
  #  beta: vector of marginal mean parameters
  #  alpha: vector of marginal correlation parameters
  #  naive: naive (model-based) covariance matrix for beta
  #  robust: robust covariance matrix for beta and alpha
  #  niter: Number of iterations until convergence
  #  n: vector of cluster sample sizes
  # Returns:
  #  output to screen
  beta_alpha_format_results <- function(beta, alpha, naive, robust, varMD,
                                        varKC, varFG, niter, n) {
    p <- length(beta)
    q <- length(alpha)
    K <- length(n)
    df <- K - p

    beta_numbers <- as.matrix(seq(1:p)) - 1
    bSE <- sqrt(diag(naive))
    bSEBC0 <- sqrt(diag(robust[1:p, 1:p]))
    bSEBC1 <- sqrt(diag(varKC[1:p, 1:p]))
    bSEBC2 <- sqrt(diag(varMD[1:p, 1:p]))
    bSEBC3 <- sqrt(diag(varFG[1:p, 1:p]))

    alpha_numbers <- as.matrix(seq(1:q)) - 1
    if (q == 1) {
      aSEBC0 <- sqrt(robust[(p + 1):(p + q), (p + 1):(p + q)])
      aSEBC1 <- sqrt(varKC[(p + 1):(p + q), (p + 1):(p + q)])
      aSEBC2 <- sqrt(varMD[(p + 1):(p + q), (p + 1):(p + q)])
      aSEBC3 <- sqrt(varFG[(p + 1):(p + q), (p + 1):(p + q)])
    } else {
      # more than 1 correlation parameters
      aSEBC0 <- sqrt(diag(robust[(p + 1):(p + q), (p + 1):(p + q)]))
      aSEBC1 <- sqrt(diag(varKC[(p + 1):(p + q), (p + 1):(p + q)]))
      aSEBC2 <- sqrt(diag(varMD[(p + 1):(p + q), (p + 1):(p + q)]))
      aSEBC3 <- sqrt(diag(varFG[(p + 1):(p + q), (p + 1):(p + q)]))
    }

    outbeta <- cbind(
      beta_numbers, beta, bSE, bSEBC0, bSEBC1,
      bSEBC2, bSEBC3
    )
    outalpha <- cbind(
      alpha_numbers, alpha, aSEBC0, aSEBC1, aSEBC2,
      aSEBC3
    )
    colnames(outbeta) <- c(
      "Beta", "Estimate", "MB-stderr", "BC0-stderr",
      "BC1-stderr", "BC2-stderr", "BC3-stderr"
    )
    colnames(outalpha) <- c(
      "Alpha", "Estimate", "BC0-stderr",
      "BC1-stderr", "BC2-stderr", "BC3-stderr"
    )

    return(list(outbeta = outbeta, outalpha = outalpha))
  }

  # main function operations start
  # reasons for non-results are identified and tallied
  alpha_work_cov_flag <- 0
  singular_naive_var_flag <- 0
  not_pos_var_est_flag <- 0
  alpha_shrink_flag <- 0
  not_pos_def_work_cov_flag <- 0
  not_pos_def_alpadj_flag <- 0

  # sort by id
  id1 <- id[order(id)]
  y <- y[order(id)]
  X <- X[order(id), ]
  id <- id1
  n <- as.vector(table(id))

  # link default
  if (is.null(link)) {
    link <- "identity"
  }

  # fit the GEE2 model
  model_fit_res <- model_fit(
    y, X, Z, n, link, maxiter, epsilon, alpha_work_cov_flag,
    singular_naive_var_flag, not_pos_var_est_flag, alpha_shrink_flag,
    not_pos_def_work_cov_flag, not_pos_def_alpadj_flag
  )

  beta <- model_fit_res$beta
  alpha <- model_fit_res$alpha
  naive <- model_fit_res$naive
  robust <- model_fit_res$robust
  varMD <- model_fit_res$varMD
  varKC <- model_fit_res$varKC
  varFG <- model_fit_res$varFG
  niter <- model_fit_res$niter

  converge <- model_fit_res$converge
  alpha_work_cov_flag <- model_fit_res$alpha_work_cov_flag
  singular_naive_var_flag <- model_fit_res$singular_naive_var_flag
  not_pos_var_est_flag <- model_fit_res$not_pos_var_est_flag
  alpha_shrink_flag <- model_fit_res$alpha_shrink_flag
  not_pos_def_work_cov_flag <- model_fit_res$not_pos_def_work_cov_flag
  not_pos_def_alpadj_flag <- model_fit_res$not_pos_def_alpadj_flag

  # final Results
  if (singular_naive_var_flag == 1) {
    stop("Derivative matrix for beta is singular during updates")
  }
  if (not_pos_var_est_flag == 1) {
    stop("Sandwich variance is not positive definite")
  }
  if (converge == 0 & singular_naive_var_flag == 0) {
    stop("The algorithm did not converge")
  }

  if (converge == 1 & not_pos_var_est_flag == 0) {
    result <- beta_alpha_format_results(
      beta, alpha, naive, robust, varMD, varKC,
      varFG, niter, n
    )
    output_res <- list(
      outbeta = result$outbeta, outalpha = result$outalpha,
      beta = beta, alpha = alpha, MB = naive, BC0 = robust,
      BC1 = varKC, BC2 = varMD, BC3 = varFG, niter = niter
    )
    class(output_res) <- "geemaee"
    return(output_res)
  }
}
