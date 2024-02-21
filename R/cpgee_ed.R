cpgee_ed <- function(
    y, X, id, m, family, maxiter, epsilon, printrange,
    alpadj, rho.init = NULL) {
  ##################################################################################### MODULE: BEGINEND creates two vectors that have the start and
  ##################################################################################### end points for each cluster

  # INPUT n: vector of cluster sizes

  # OUTPUT first: vector with starting row for cluster i last:
  # vector with ending row for cluster i

  BEGINEND <- function(n) {
    last <- cumsum(n)
    first <- last - n + 1
    return(cbind(first, last))
  }

  ##################################################################################### Module: IS_POS_DEF A = symmetric matrix returns 1 if A is
  ##################################################################################### positive definite 0 otherwise

  is_pos_def <- function(A) {
    return(min(eigen(A)$values) > 1e-13)
  }

  ##################################################################################### MODULE: GETROWB generates a row of the -E(d(cov)/dbeta) (ie Q)
  ##################################################################################### matrix

  # INPUT mu: Marginal means for cluster i j: Indicator for mean
  # k: Indicator for mean X: Covariate matrix for cluster i y:
  # Response vector for cluster i

  # OUTPUT row of E(d(cov)/dbeta) (ie Q) matrix

  GETROWB <- function(mu, j, k, X, y) {
    row <- -(y[j] - mu[j]) * X[j, ] * mu[j] * (1 - mu[j]) - (y[k] -
      mu[k]) * X[k, ] * mu[k] * (1 - mu[k])
    return(row)
  }

  ##################################################################################### MODULE: CREATEA creates residual for beta estimating equation,
  ##################################################################################### (Y - mu)

  # INPUT mu: vector of n_i marginal means y: outcome vector for
  # ith cluster

  # OUTPUT residuals for beta estimating equation

  CREATEA <- function(mu, y) {
    return(y - mu)
  }

  ##################################################################################### MODULE: CREATEC creates derivative matrix for beta estimating
  ##################################################################################### equation, dmu/dbeta

  # INPUT X: covariate matrix for cluster i mu: vector of n_i
  # marginal means

  # OUTPUT derivative matrix for beta estimating equation

  CREATEC <- function(X, mu) {
    return(X * (mu * (1 - mu)))
  }

  ##################################################################################### MODULE: CREATEB creates covariance matrix for beta estimating
  ##################################################################################### equation, var(Y)

  # INPUT gamma: vector of s_jk between mean j and mean k for
  # cluster i n: sample size (scalar) for cluster i

  # OUTPUT covariance matrix for beta estimating equation

  CREATEB <- function(gamma, n) {
    B <- matrix(0, n, n)
    l <- 1
    for (j in 1:n) {
      for (k in j:n) {
        B[j, k] <- gamma[l]
        l <- l + 1
      }
    }
    B[lower.tri(B)] <- t(B)[lower.tri(B)]
    return(B)
  }

  ##################################################################################### MODULE: CREATED creates derivative matrix for alpha estimating
  ##################################################################################### equation, deta/dalpha

  # INPUT mu: vector of marginal means alpha: correlation
  # parameters n: sample size (scalar) for cluster i m: vector of
  # cluster-period sizes

  # OUTPUT derivative matrix for alpha estimating equation

  CREATED <- function(mu, alpha, n, m) {
    alpha0 <- alpha[1] # within-period ICC
    rho <- alpha[2] # decay
    v <- mu * (1 - mu)
    D <- NULL
    for (j in 1:(n - 1)) {
      dj <- (v[j] / m[j]) * (m[j] - 1)
      D <- rbind(D, c(dj, 0))
      for (k in (j + 1):n) {
        djk0 <- sqrt(v[j] * v[k]) * rho^(abs(k - j))
        djk1 <- sqrt(v[j] * v[k]) * alpha0 * abs(k - j) *
          rho^(abs(k - j) - 1)
        D <- rbind(D, c(djk0, djk1))
      }
    }
    D <- rbind(D, c((v[n] / m[n]) * (m[n] - 1), 0))
    return(D)
  }

  ##################################################################################### MODULE: gammahat Creates vector of estimated covariances

  # INPUT mu: vector of n_i marginal means alpha: correlation
  # estimates n: cluster size (# of periods) m: vector of
  # cluster-period sizes

  # OUTPUT vector of estimated covariances

  gammahat <- function(mu, alpha, n, m) {
    alpha0 <- alpha[1] # within-period ICC
    rho <- alpha[2] # decay
    v <- mu * (1 - mu)
    gamma_c <- NULL
    for (j in 1:(n - 1)) {
      sj <- (v[j] / m[j]) * (1 + (m[j] - 1) * alpha0)
      gamma_c <- c(gamma_c, sj)
      for (k in (j + 1):n) {
        sjk <- sqrt(v[j] * v[k]) * alpha0 * rho^(abs(k -
          j))
        gamma_c <- c(gamma_c, sjk)
      }
    }
    gamma_c <- c(gamma_c, (v[n] / m[n]) * (1 + (m[n] - 1) * alpha0))
    return(gamma_c)
  }

  ##################################################################################### MODULE: checkalpha Creates compact correlation matrix for
  ##################################################################################### range checks

  # INPUT alpha: correlation estimates n: cluster size (# of
  # periods)

  # OUTPUT A correlation matrix indexed by cluster-periods

  checkalpha <- function(alpha, n) {
    alpha0 <- alpha[1] # within-period ICC
    rho <- alpha[2] # decay
    L <- matrix(0, n, n)
    for (j in 1:(n - 1)) {
      L[j, j] <- alpha0
      for (k in (j + 1):n) {
        L[j, k] <- alpha0 * rho^(abs(k - j))
      }
    }
    L[n, n] <- alpha0
    L[lower.tri(L)] <- t(L)[lower.tri(L)]
    return(L)
  }

  ##################################################################################### MODULE: SCORE generates the score matrix for each cluster and
  ##################################################################################### approximate information to be used to estimate parameters and
  ##################################################################################### generate standard errors

  # INPUT Ustarold: initial values for information matrix beta:
  # vector of marginal mean parameters alpha: vector of marginal
  # correlation parameters y: vector of outcomes (cluster-period
  # means) X: marginal mean covariates m: vector of cluster-period
  # sizes n: vector of cluster sample sizes p: number of marginal
  # mean parameters q: number of marginal correlation parameters
  # NPSDFLAG: checks if B is positive definite, terminates if not
  # NPSDADJFLAG: checks if the alpha adjustment factor matrix is
  # positive definite, terminates if not

  # OUTPUT U: score vector UUtran: sum of U_i*U_i` across all
  # clusters Ustar: approximate information matrix

  SCORE <- function(Ustarold, beta, alpha, y, X, m, n, p, q, NPSDFLAG,
                    NPSDADJFLAG) {
    U <- rep(0, p + q)
    UUtran <- Ustar <- matrix(0, p + q, p + q)
    naiveold <- ginv(Ustarold[1:p, 1:p]) # needed for Hi1 below

    locx <- BEGINEND(n)
    locz <- BEGINEND(choose(n + 1, 2))

    alpdata0 <- NULL
    alpdata1 <- NULL

    for (i in 1:length(n)) {
      X_c <- X[locx[i, 1]:locx[i, 2], , drop = FALSE]
      y_c <- y[locx[i, 1]:locx[i, 2]]
      m_c <- m[locx[i, 1]:locx[i, 2]]

      U_c <- rep(0, p + q)
      Ustar_c <- matrix(0, p + q, p + q)
      mu_c <- 1 / (1 + exp(c(-X_c %*% beta)))
      v_c <- mu_c * (1 - mu_c)
      gamma_c <- gammahat(mu_c, alpha, n[i], m_c)

      R <- rep(0, choose(n[i] + 1, 2))
      DB <- matrix(0, choose(n[i] + 1, 2), p)

      C <- CREATEC(X_c, mu_c)
      B <- CREATEB(gamma_c, n[i])
      A <- CREATEA(mu_c, y_c)
      D <- CREATED(mu_c, alpha, n[i], m_c)

      INVB <- ginv(B)
      CtinvB <- t(C) %*% INVB
      Hi1 <- C %*% naiveold %*% CtinvB
      G_c <- tcrossprod(y_c - mu_c)

      if (alpadj) {
        CT <- t(C)
        omega <- C %*% naiveold %*% CT
        vminomega <- B - omega
        psd_vmin <- is_pos_def(vminomega)


        if (psd_vmin == 1) {
          Ci <- B %*% ginv(vminomega)
          G_c <- Ci %*% G_c
        } else {
          NPSDADJFLAG <- 1
          stop("(V - Omega) is not positive definite")
        }
      }


      l <- 1
      for (j in 1:n[i]) {
        for (k in j:n[i]) {
          R[l] <- G_c[j, k] - gamma_c[l]
          l <- l + 1
        }
      }

      # Check for positive definite of B
      if (min(eigen(B)$values) <= 0) {
        NPSDFLAG <- 1
        stop(paste(
          "Var(Y) of Cluster", i, "is not Positive-Definite;",
          "Joint Distribution Does Not Exist and Program terminates"
        ))
      }

      U_c[1:p] <- t(C) %*% INVB %*% A
      U_c[(p + 1):(p + q)] <- t(D) %*% R
      UUtran_c <- tcrossprod(U_c)
      Ustar_c[1:p, 1:p] <- t(C) %*% INVB %*% C
      Ustar_c[(p + 1):(p + q), (p + 1):(p + q)] <- crossprod(D)
      U <- U + U_c
      UUtran <- UUtran + UUtran_c
      Ustar <- Ustar + Ustar_c

      # build dataset for within-period ICC
      r1 <- cbind(diag(G_c), v_c, m_c)
      alpdata0 <- rbind(alpdata0, r1)

      # build dataset for remaining
      for (j in 1:(n[i] - 1)) {
        for (k in (j + 1):(n[i])) {
          r2 <- c(G_c[j, k], v_c[j] * v_c[k], j, k, i)
          alpdata1 <- rbind(alpdata1, r2)
        }
      }
    }
    return(list(
      U = U, UUtran = UUtran, Ustar = Ustar, alpdata0 = alpdata0,
      alpdata1 = alpdata1, NPSDFLAG = NPSDFLAG, NPSDADJFLAG = NPSDADJFLAG
    ))
  }

  ##################################################################################### MODULE: INITBETA generates initial values for beta
  ##################################################################################### approximates logistic regression using Newton's method

  # INPUT y: vector of cluster-period means m: vector of
  # cluster-period sizes (# of trials) X: marginal mean covariates

  # OUTPUT beta: vector of marginal mean parameters Ustar:
  # approximate score vector

  INITBETA <- function(y, m, X) {
    fit <- glm(cbind(m * y, m * (1 - y)) ~ -1 + X, family = binomial(link = "logit"))
    beta <- as.numeric(fit$coefficients)
    u <- as.numeric(fit$fitted.values)
    v <- u * (1 - u)
    Ustar <- t(X) %*% (X * v * m)
    return(list(beta = c(beta), Ustar = Ustar))
  }

  ##################################################################################### MODULE: getalpha0, getalpha1, getalpha (getalpha0, getalpha1)
  ##################################################################################### updates values for correlation parameters

  # INPUT alpdata0: dataset built for estimating alpha0 alpdata1:
  # dataset built for estimating the rest

  # OUTPUT alpha0: correlation parameter alpha1: correlation decay
  # parameter alpha: c(alpha0, alpha1)
  getalpha0 <- function(alpdata0) {
    s <- alpdata0[, 1]
    v <- alpdata0[, 2]
    m <- alpdata0[, 3]

    d1 <- sum(((m - 1) / m) * (s * v - v^2 / m))
    d0 <- sum(((m - 1) / m)^2 * v^2)
    return(d1 / d0)
  }

  getalpha1 <- function(alpdata1, alpha0) {
    s <- alpdata1[, 1]
    vcross <- alpdata1[, 2]
    d <- abs(alpdata1[, 3] - alpdata1[, 4])

    polyfun <- function(rho) {
      smt1 <- sum(sqrt(vcross) * d * rho^(d - 1) * s)
      smt2 <- sum(vcross * alpha0 * d * rho^(2 * d - 1))
      return(smt1 - smt2)
    }
    opt <- uniroot(polyfun, c(0, 1))
    return(opt$root)
  }

  getalpha <- function(alpdata0, alpdata1, alphaold) {
    s0 <- alpdata0[, 1]
    v <- alpdata0[, 2]
    m <- alpdata0[, 3]

    s1 <- alpdata1[, 1]
    vcross <- alpdata1[, 2]
    d <- abs(alpdata1[, 3] - alpdata1[, 4])

    f <- function(alpha) {
      alpha0 <- alpha[1]
      rho <- alpha[2]
      f0_smt1 <- ((m - 1) / m) * v * (s0 - v / m - alpha0 * v *
        (m - 1) / m)
      f0_smt2 <- sqrt(vcross) * (rho^d) * (s1 - sqrt(vcross) *
        alpha0 * rho^d)
      f0 <- sum(f0_smt1) + sum(f0_smt2)
      f1 <- sum(sqrt(vcross) * d * rho^(d - 1) * (s1 - sqrt(vcross) *
        alpha0 * rho^d))
      return(c(f0, f1))
    }
    res <- multiroot(f, start = alphaold)$root
    return(res)
  }

  ##################################################################################### MODULE: INITALPHA generates initial values for alpha

  # INPUT y: vector of cluster-period means m: vector of
  # cluster-period sizes (# of trials) X: marginal mean covariates
  # n: vector of cluster sizes (# of periods) beta: estimated
  # regression coefficients

  # OUTPUT alpha: estimated correlation values

  INITALPHA <- function(y, X, m, n, beta, maxiter, epsilon, rho.init) {
    alpdata0 <- NULL
    alpdata1 <- NULL
    locx <- BEGINEND(n)

    for (i in 1:length(n)) {
      X_c <- X[locx[i, 1]:locx[i, 2], , drop = FALSE]
      y_c <- y[locx[i, 1]:locx[i, 2]]
      m_c <- m[locx[i, 1]:locx[i, 2]]

      mu_c <- 1 / (1 + exp(c(-X_c %*% beta)))
      v_c <- mu_c * (1 - mu_c)
      G_c <- tcrossprod(y_c - mu_c)

      # build dataset for within-period ICC
      r1 <- cbind(diag(G_c), v_c, m_c)
      alpdata0 <- rbind(alpdata0, r1)

      # build dataset for remaining
      for (j in 1:(n[i] - 1)) {
        for (k in (j + 1):(n[i])) {
          r2 <- c(G_c[j, k], v_c[j] * v_c[k], j, k, i)
          alpdata1 <- rbind(alpdata1, r2)
        }
      }
    }

    # produce initial estimates
    alpha0 <- getalpha0(alpdata0)
    if (is.null(rho.init)) {
      alpha1 <- getalpha1(alpdata1, alpha0)
    } else {
      alpha1 <- rho.init
    }

    return(c(alpha0, alpha1))
  }

  ##################################################################################### MODULE: INVBIG compute (A - mm`)^{-1}c without performing the
  ##################################################################################### inverse directly

  # INPUT ainvc: inverse of matrix A times vector c ainvm: inverse
  # of matrix A times matrix (with low no. of columns) M M: matrix
  # of eigen column vectors m1,m2, ..  c: vector start: of do loop
  # end: of do loop, rank of X

  # OUTPUT ainvc: inverse of matrix A times vector c

  INVBIG <- function(ainvc, ainvm, m, c, start, end) {
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

  ##################################################################################### MODULE: MAKEVAR creates covariance matrix of beta and alpha.

  # INPUT beta: vector of marginal mean parameters alpha: vector
  # of marginal correlation parameters Ustarold: initial values
  # for information matrix y: vector of cluster-period means X:
  # marginal mean covariates m: vector of cluster-period sizes n:
  # vector of cluster sample sizes p: number of marginal mean
  # parameters q: number of marginal correlation parameters
  # NPSDFLAG: checks if B is positive definite, terminates if not
  # NPSDADJFLAG: checks if the alpha adjustment factor matrix is
  # positive definite, terminates if not

  # OUTPUT robust: robust covariance matrix for beta and alpha
  # naive: naive (model-based) covariance matrix for beta varMD:
  # bias-corrected variance by Mancl and Derouen (2001) varKC:
  # bias-corrected variance by Kauermann and Carroll (2001) varFG:
  # bias-corrected variance by Fay and Graubard (2001)

  MAKEVAR <- function(Ustarold, beta, alpha, y, X, m, n, p, q,
                      ROBFLAG, NPSDFLAG, NPSDADJFLAG) {
    SCORE_RES <- SCORE(
      Ustarold, beta, alpha, y, X, m, n, p,
      q, NPSDFLAG, NPSDADJFLAG
    )
    U <- SCORE_RES$U
    UUtran <- SCORE_RES$UUtran
    Ustar <- SCORE_RES$Ustar
    NPSDFLAG <- SCORE_RES$NPSDFLAG
    NPSDADJFLAG <- SCORE_RES$NPSDADJFLAG

    naive <- ginv(Ustar[1:p, 1:p])
    naivealp <- ginv(Ustar[(p + 1):(p + q), (p + 1):(p + q)])

    # new commands to compute INV(I - H1)
    eigenRES1 <- eigen(naive)
    evals1 <- eigenRES1$values
    evecs1 <- eigenRES1$vectors
    sqrevals1 <- sqrt(evals1)
    sqe1 <- evecs1 %*% diag(sqrevals1)

    # new commands to compute INV(I - H2)
    eigenRES2 <- eigen(naivealp)
    evals2 <- eigenRES2$values
    evecs2 <- eigenRES2$vectors
    sqrevals2 <- sqrt(evals2)
    sqe2 <- evecs2 %*% diag(sqrevals2)

    # Bias-corrected variance
    Ustar_c_array <- UUtran_c_array <- array(0, c(
      p + q, p + q,
      length(n)
    ))
    UUtran <- UUbc <- UUbc2 <- UUbc3 <- Ustar <- inustar <- matrix(
      0,
      p + q, p + q
    )

    locx <- BEGINEND(n)
    locz <- BEGINEND(choose(n + 1, 2))

    for (i in 1:length(n)) {
      X_c <- X[locx[i, 1]:locx[i, 2], , drop = FALSE]
      y_c <- y[locx[i, 1]:locx[i, 2]]
      m_c <- m[locx[i, 1]:locx[i, 2]]

      mu_c <- 1 / (1 + exp(c(-X_c %*% beta)))
      v_c <- mu_c * (1 - mu_c)

      U_i <- U_c <- rep(0, p + q)
      Ustar_c <- matrix(0, p + q, p + q)
      gamma_c <- gammahat(mu_c, alpha, n[i], m_c)

      # commands for beta
      C <- CREATEC(X_c, mu_c)
      B <- CREATEB(gamma_c, n[i])
      A <- CREATEA(mu_c, y_c)
      D <- CREATED(mu_c, alpha, n[i], m_c)

      INVB <- ginv(B)
      U_i[1:p] <- t(C) %*% INVB %*% A

      CtinvB <- t(C) %*% INVB
      Hi1 <- C %*% naive %*% CtinvB

      # commands for beta
      ai1 <- INVB
      mm1 <- C %*% sqe1
      ai1A <- ai1 %*% A
      ai1m1 <- ai1 %*% mm1
      ai1A <- INVBIG(ai1A, ai1m1, mm1, A, 1, p)
      U_c[1:p] <- t(C) %*% ai1A

      # commands for alpha
      R <- rep(0, choose(n[i] + 1, 2))
      DB <- matrix(0, choose(n[i] + 1, 2), p)
      G_c <- tcrossprod(y_c - mu_c)

      if (alpadj) {
        # MAEE
        CT <- t(C)
        omega <- C %*% naive %*% CT
        vminomega <- B - omega
        psd_vmin <- is_pos_def(vminomega)


        if (psd_vmin == 1) {
          Ci <- B %*% ginv(vminomega)
          G_c <- Ci %*% G_c
        } else {
          NPSDADJFLAG <- 1
          stop("(V - Omega) is not positive definite")
        }
      }

      # RANGE CHECKS
      rangeflag <- 0
      L_c <- checkalpha(alpha, n[i])
      l <- 1
      for (j in 1:n[i]) {
        for (k in j:n[i]) {
          if ((L_c[j, k] >= min(sqrt((mu_c[j] * (1 - mu_c[k])) / (mu_c[k] *
            (1 - mu_c[j]))), sqrt((mu_c[k] * (1 - mu_c[j])) / (mu_c[j] *
            (1 - mu_c[k]))))) | (L_c[j, k] <= max(-sqrt((mu_c[j] *
            mu_c[k]) / ((1 - mu_c[j]) * (1 - mu_c[k]))), -sqrt(((1 -
            mu_c[j]) * (1 - mu_c[k])) / (mu_c[j] * mu_c[k]))))) {
            rangeflag <- 1
            if (printrange) {
              warning(cat(
                "Range Violation Detected for Cluster",
                i, "and Pair", j, k, "\n"
              ))
            }
            break
          }

          DB[l, ] <- GETROWB(mu_c, j, k, X_c, y_c)

          R[l] <- G_c[j, k] - gamma_c[l]
          l <- l + 1
        }
      }

      # Check for positive definite of B
      if (min(eigen(B)$values) <= 0) {
        NPSDFLAG <- 1
        stop(paste(
          "Var(Y) of Cluster", i, "is not Positive-Definite;",
          "Joint Distribution Does Not Exist and Program terminates"
        ))
      }

      U_i[(p + 1):(p + q)] <- t(D) %*% R
      mm2 <- D %*% sqe2
      ai2R <- R
      ai2m2 <- mm2
      ai2R <- INVBIG(ai2R, ai2m2, mm2, R, 1, q)
      U_c[(p + 1):(p + q)] <- t(D) %*% ai2R

      Ustar_c[1:p, 1:p] <- t(C) %*% INVB %*% C
      Ustar_c[(p + 1):(p + q), 1:p] <- t(D) %*% DB
      Ustar_c[(p + 1):(p + q), (p + 1):(p + q)] <- crossprod(D)
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

    inustar[1:p, 1:p] <- ginv(Ustar[1:p, 1:p])
    inustar[(p + 1):(p + q), (p + 1):(p + q)] <- ginv(Ustar[(p +
      1):(p + q), (p + 1):(p + q)])
    inustar[(p + 1):(p + q), 1:p] <- inustar[
      (p + 1):(p + q),
      (p + 1):(p + q)
    ] %*% Ustar[(p + 1):(p + q), 1:p] %*%
      inustar[1:p, 1:p]

    # the minus sign above is crucial, esp. for large correlation;
    inustartr <- t(inustar)

    # calculating adjustment factor for BC3
    for (i in 1:length(n)) {
      Hi <- diag(1 / sqrt(1 - pmin(0.75, c(diag(Ustar_c_array[
        , ,
        i
      ] %*% inustar)))))
      UUbc3 <- UUbc3 + Hi %*% UUtran_c_array[, , i] %*% Hi
    }

    # BC0 or usual Sandwich estimator;
    robust <- inustar %*% UUtran %*% inustartr

    # BC1 or Variance estimator that extends Kauermann and Carroll
    # (2001);
    varKC <- inustar %*% (UUbc2 + t(UUbc2)) %*% inustartr / 2

    # BC2 or Variance estimator that extends Mancl and DeRouen
    # (2001);
    varMD <- inustar %*% UUbc %*% inustartr

    # BC3 or Variance estimator that extends Fay and Graubard
    # (2001);
    varFG <- inustar %*% UUbc3 %*% inustartr

    # model-based variance
    naive <- inustar[1:p, 1:p]

    if (min(diag(robust)) <= 0) {
      ROBFLAG <- 1
    }
    if (min(diag(varMD)) <= 0) {
      ROBFLAG <- 1
    }
    if (min(diag(varKC)) <= 0) {
      ROBFLAG <- 1
    }
    if (min(diag(varFG)) <= 0) {
      ROBFLAG <- 1
    }

    return(list(
      robust = robust, naive = naive, varMD = varMD,
      varKC = varKC, varFG = varFG, rangeflag = rangeflag,
      ROBFLAG = ROBFLAG, NPSDFLAG = NPSDFLAG, NPSDADJFLAG = NPSDADJFLAG
    ))
  }

  ##################################################################################### MODULE: FITPRENTICE Performs the paired estimating equations
  ##################################################################################### method

  # INPUT y: vector of cluster-period means X: marginal mean
  # covariates n: Vector of cluster sizes m: vector of
  # cluster-period sizes maxiter: maximum number of iterations
  # epsilon: tolerence for convergence SINGFLAG: THE ALGORITHM
  # terminated due to singular MB covariance matrix ROBFLAG: THE
  # ALGORITHM terminated due to singular robust variance matrix
  # NPSDFLAG: checks if B is positive definite, terminates if not
  # NPSDADJFLAG: checks if the alpha adjustment factor matrix is
  # positive definite, terminates if not

  # OUTPUT beta: a p x 1 vector of marginal mean parameter
  # estimates alpha: a q x 1 vector of marginal correlation
  # parameter estimates robust: robust covariance matrix for beta
  # and alpha naive: naive (model-based) covariance matrix for
  # beta varMD: bias-corrected variance due to Mancl and DeRouen
  # (2001) varKC: bias-corrected variance due to Kauermann and
  # Carroll (2001) varFG: bias-corrected variance due to Fay and
  # Graubard (2001) niter: number of iterations required for
  # convergence converge: did the algorithm converge (0 = no, 1 =
  # yes)

  FITPRENTICE <- function(y, X, m, n, maxiter, epsilon, SINGFLAG,
                          ROBFLAG, NPSDFLAG, NPSDADJFLAG) {
    p <- ncol(X)
    converge <- 0
    rangeflag <- 0

    INITRES <- INITBETA(y, m, X)
    beta <- INITRES$beta
    Ustar <- INITRES$Ustar
    alpha <- INITALPHA(y, X, m, n, beta, maxiter, epsilon, rho.init)
    q <- length(alpha)
    delta <- rep(2 * epsilon, p)
    deltaalp <- rep(2 * epsilon, q)

    niter <- 1
    while ((niter <= maxiter) & (max(abs(c(delta, deltaalp))) >
      epsilon)) {
      SINGFLAG <- 0
      NPSDFLAG <- 0
      NPSDADJFLAG <- 0
      Ustarold <- Ustar
      alphaold <- alpha
      SCORE_RES <- SCORE(
        Ustarold, beta, alpha, y, X, m, n,
        p, q, NPSDFLAG, NPSDADJFLAG
      )
      U <- SCORE_RES$U
      UUtran <- SCORE_RES$UUtran
      Ustar <- SCORE_RES$Ustar
      NPSDFLAG <- SCORE_RES$NPSDFLAG
      NPSDADJFLAG <- SCORE_RES$NPSDADJFLAG

      alpdata0 <- SCORE_RES$alpdata0
      alpdata1 <- SCORE_RES$alpdata1
      alpha <- getalpha(alpdata0, alpdata1, alphaold)
      deltaalp <- alpha - alphaold

      psdustar <- is_pos_def(Ustar[1:p, 1:p])
      mineig <- min(eigen(Ustar[1:p, 1:p])$values)
      if (psdustar == TRUE) {
        delta <- solve(Ustar[1:p, 1:p], U[1:p])
        beta <- beta + delta
        converge <- (max(abs(c(delta, deltaalp))) <= epsilon)
      } else {
        SINGFLAG <- 1
      }
      niter <- niter + 1
    }

    Ustarold <- Ustar

    # inference
    MAKEVAR_RES <- MAKEVAR(
      Ustarold, beta, alpha, y, X, m, n,
      p, q, ROBFLAG, NPSDFLAG, NPSDADJFLAG
    )
    robust <- MAKEVAR_RES$robust
    naive <- MAKEVAR_RES$naive
    varMD <- MAKEVAR_RES$varMD
    varKC <- MAKEVAR_RES$varKC
    varFG <- MAKEVAR_RES$varFG
    rangeflag <- MAKEVAR_RES$rangeflag
    ROBFLAG <- MAKEVAR_RES$ROBFLAG
    NPSDFLAG <- MAKEVAR_RES$NPSDFLAG
    NPSDADJFLAG <- MAKEVAR_RES$NPSDADJFLAG

    return(list(
      beta = beta, alpha = alpha, robust = robust,
      naive = naive, varMD = varMD, varKC = varKC, varFG = varFG,
      niter = niter, converge = converge, SINGFLAG = SINGFLAG,
      ROBFLAG = ROBFLAG, NPSDFLAG = NPSDFLAG, NPSDADJFLAG = NPSDADJFLAG
    ))
  }

  ##################################################################################### MODULE: RESULTS creates printed output to screen of parameters
  ##################################################################################### and other information

  # INPUT beta: vector of marginal mean parameters alpha: vector
  # of marginal correlation Parameters robust: robust covariance
  # matrix for beta and alpha naive: naive (model-based)
  # covariance matrix for beta niter: number of iterations until
  # convergence n: vector of cluster sample sizes

  # OUTPUT to screen
  RESULTS <- function(beta, alpha, robust, naive, varMD, varKC,
                      varFG, niter, n) {
    p <- length(beta)
    q <- length(alpha)
    K <- length(n)

    # ** the next message is structure specific **
    corstr <- "Exponential decay"


    beta_numbers <- as.matrix(seq(1:p)) - 1
    bSE <- sqrt(diag(naive))
    bSEBC0 <- sqrt(diag(robust[1:p, 1:p]))
    bSEBC1 <- sqrt(diag(varKC[1:p, 1:p]))
    bSEBC2 <- sqrt(diag(varMD[1:p, 1:p]))
    bSEBC3 <- sqrt(diag(varFG[1:p, 1:p]))

    alpha_numbers <- as.matrix(seq(1:q)) - 1
    aSEBC0 <- sqrt(diag(robust[(p + 1):(p + q), (p + 1):(p +
      q)]))
    aSEBC1 <- sqrt(diag(varKC[(p + 1):(p + q), (p + 1):(p + q)]))
    aSEBC2 <- sqrt(diag(varMD[(p + 1):(p + q), (p + 1):(p + q)]))
    aSEBC3 <- sqrt(diag(varFG[(p + 1):(p + q), (p + 1):(p + q)]))

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

  # reasons for non-results are identified and tallied
  SINGFLAG <- 0
  ROBFLAG <- 0
  NPSDFLAG <- 0
  NPSDADJFLAG <- 0

  id1 <- id[order(id)]
  y <- y[order(id)]
  X <- X[order(id), ]
  m <- m[order(id)]
  id <- id1

  n <- as.vector(table(id))
  # Fit the GEE/MAEE algorithm
  PRENTICE_RES <- FITPRENTICE(
    y, X, m, n, maxiter, epsilon, SINGFLAG,
    ROBFLAG, NPSDFLAG, NPSDADJFLAG
  )
  beta <- PRENTICE_RES$beta
  alpha <- PRENTICE_RES$alpha
  robust <- PRENTICE_RES$robust
  naive <- PRENTICE_RES$naive
  varMD <- PRENTICE_RES$varMD
  varKC <- PRENTICE_RES$varKC
  varFG <- PRENTICE_RES$varFG
  niter <- PRENTICE_RES$niter
  converge <- PRENTICE_RES$converge
  SINGFLAG <- PRENTICE_RES$SINGFLAG
  ROBFLAG <- PRENTICE_RES$ROBFLAG
  NPSDFLAG <- PRENTICE_RES$NPSDFLAG
  NPSDADJFLAG <- PRENTICE_RES$NPSDADJFLAG

  # Final Results
  if (SINGFLAG == 1) {
    stop("Derivative matrix for beta is singular during updates")
  }
  if (ROBFLAG == 1) {
    stop("Sandwich variance is not positive definite")
  }
  if (converge == 0 & SINGFLAG == 0) {
    stop("The algorithm did not converge")
  }
  if (converge == 1 & ROBFLAG == 0) {
    result <- RESULTS(
      beta, alpha, robust, naive, varMD, varKC,
      varFG, niter, n
    )
    outList <- list(
      outbeta = result$outbeta, outalpha = result$outalpha,
      beta = beta, alpha = alpha, MB = naive, BC0 = robust,
      BC1 = varKC, BC2 = varMD, BC3 = varFG, niter = niter
    )
    class(outList) <- "cpgeeSWD"
    return(outList)
  }
}
