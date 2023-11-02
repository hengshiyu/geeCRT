#' Generating Correlated Binary Data using the Conditional Linear Family Method.
#' @param mu a mean vector when \code{n = 1} or is \code{NULL}, otherwise a list of mean vectors for the \code{n} clusters
#' @param Sigma a correlation matrix when \code{n = 1} or is \code{NULL}, otherwise a list of correlation matrices for the \code{n} clusters
#' @param n number of clusters. The default is \code{1}
#' @keywords cluster-randomized-trials correlated-binary-data conditional-linear-family
#' @author Hengshi Yu <hengshi@umich.edu>, Fan Li <fan.f.li@yale.edu>, Paul Rathouz <paul.rathouz@austin.utexas.edu>, Elizabeth L. Turner <liz.turner@duke.edu>, John Preisser <jpreisse@bios.unc.edu>
#' @description simbinCLF generates correlated binary data using the conditional linear family method (Qaqish, 2003).
#' It simulates a vector of binary outcomes according the specified marginal mean vector and correlation structure.
#' Natural constraints and compatibility between the marginal mean and correlation matrix are checked.
#'
#'
#' @references
#'
#' Qaqish, B. F. (2003). A family of multivariate binary distributions for simulating correlated binary variables with
#' specified marginal means and correlations. Biometrika, 90(2), 455-463.
#'
#' Preisser, J. S., Qaqish, B. F. (2014). A comparison of methods for simulating correlated binary variables with specified marginal means
#' and correlations. Journal of Statistical Computation and Simulation, 84(11), 2441-2452.
#'
#' @export
#' @examples
#'
#' #####################################################################
#' # Simulate 2 clusters, 3 periods and cluster-period size of 5 #######
#' #####################################################################
#'
#' t <- 3
#' n <- 2
#' m <- 5
#'
#' # means of cluster 1
#' u_c1 <- c(0.4, 0.3, 0.2)
#' u1 <- rep(u_c1, c(rep(m, t)))
#' # means of cluster 2
#' u_c2 <- c(0.35, 0.25, 0.2)
#' u2 <- rep(u_c2, c(rep(m, t)))
#'
#' # List of mean vectors
#' mu <- list()
#' mu[[1]] <- u1
#' mu[[2]] <- u2
#' # List of correlation matrices
#'
#' ## correlation parameters
#' alpha0 <- 0.03
#' alpha1 <- 0.015
#' rho <- 0.8
#'
#' ## (1) exchangeable
#' Sigma <- list()
#' Sigma[[1]] <- diag(m * t) * (1 - alpha1) + matrix(alpha1, m * t, m * t)
#' Sigma[[2]] <- diag(m * t) * (1 - alpha1) + matrix(alpha1, m * t, m * t)
#'
#' y_exc <- simbinCLF(mu = mu, Sigma = Sigma, n = n)
#'
#' ## (2) nested exchangeable
#' Sigma <- list()
#' cor_matrix <- matrix(alpha1, m * t, m * t)
#' loc1 <- 0
#' loc2 <- 0
#' for (t in 1:t) {
#'     loc1 <- loc2 + 1
#'     loc2 <- loc1 + m - 1
#'     for (i in loc1:loc2) {
#'         for (j in loc1:loc2) {
#'             if (i != j) {
#'                 cor_matrix[i, j] <- alpha0
#'             } else {
#'                 cor_matrix[i, j] <- 1
#'             }
#'         }
#'     }
#' }
#'
#' Sigma[[1]] <- cor_matrix
#' Sigma[[2]] <- cor_matrix
#'
#' y_nex <- simbinCLF(mu = mu, Sigma = Sigma, n = n)
#'
#' ## (3) exponential decay
#'
#' Sigma <- list()
#'
#' ### function to find the period of the ith index
#' region_ij <- function(points, i) {
#'     diff <- i - points
#'     for (h in 1:(length(diff) - 1)) {
#'         if (diff[h] > 0 & diff[h + 1] <= 0) {
#'             find <- h
#'         }
#'     }
#'     return(find)
#' }
#'
#' cor_matrix <- matrix(0, m * t, m * t)
#' useage_m <- cumsum(m * t)
#' useage_m <- c(0, useage_m)
#'
#' for (i in 1:(m * t)) {
#'     i_reg <- region_ij(useage_m, i)
#'     for (j in 1:(m * t)) {
#'         j_reg <- region_ij(useage_m, j)
#'         if (i_reg == j_reg & i != j) {
#'             cor_matrix[i, j] <- alpha0
#'         } else if (i == j) {
#'             cor_matrix[i, j] <- 1
#'         } else if (i_reg != j_reg) {
#'             cor_matrix[i, j] <- alpha0 * (rho^(abs(i_reg - j_reg)))
#'         }
#'     }
#' }
#'
#' Sigma[[1]] <- cor_matrix
#' Sigma[[2]] <- cor_matrix
#'
#' y_ed <- simbinCLF(mu = mu, Sigma = Sigma, n = n)
#'
#' @return \code{y} a vector of simulated binary outcomes for \code{n} clusters.




simbinCLF <- function(mu, Sigma, n = 1) {
    # a[1:n, 1:n] is the input covariance matrix of Y[1:n].  Returns
    # b[1:n,1:n] such that b[1:t-1, t] are the slopes for regression
    # of y[t] on y[1:t-1], for t=2:n.  Diagonals and lower half of
    # b[,] are copied from a[,].  a[,] is assumed +ve definite
    # symmetric, not checked.

    allreg <- function(a) {
        n <- nrow(a)
        b <- a
        for (t in 2:n) {
            t1 <- t - 1
            gt <- a[1:t1, 1:t1]
            st <- a[1:t1, t]
            bt <- solve(gt, st)
            b[1:t1, t] <- bt
        }
        return(b)
    }

    # returns variance matrix of binary variables with mean vector
    # u[] and corr matrix r[,].

    cor2var <- function(r, u) {
        s <- diag(sqrt(u * (1 - u)))
        return(s %*% r %*% s)
    }

    # r[1:n, 1:n] is the corr mtx u[1:n] is the mean of a binary
    # vector checks that pairwise corrs are in-range for the given
    # u[] only upper half of r[,] is checked.  return 0 if ok return
    # 1 if out of range

    chkbinc <- function(r, u) {
        n <- length(u)
        s <- sqrt(u * (1 - u))
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                uij <- u[i] * u[j] + r[i, j] * s[i] * s[j]
                ok <- ((uij <= min(u[i], u[j])) & (uij >= max(
                    0,
                    u[i] + u[j] - 1
                )))
                if (!ok) {
                    return(1)
                }
            }
        }
        return(0)
    }

    # Multivariate Binary Simulation by Linear Regression.  Simulate
    # a single vector.  Returns a simulated binary random vector
    # y[1:n] with mean u[1:n] and regression coefs matrix b[1:n,1:n]
    # (obtained by calling allreg() above).  y[] and u[] are column
    # vectors.  Returns -1 if the cond. linear family not
    # reproducible

    mbslr1 <- function(b, u) {
        n <- nrow(b)
        y <- rep(-1, n)
        y[1] <- rbinom(1, 1, u[1])
        for (i in 2:n) {
            i1 <- i - 1
            r <- y[1:i1] - u[1:i1] # residuals
            ci <- u[i] + sum(r * b[1:i1, i]) # cond.mean
            if (ci < 0 | ci > 1) {
                stop(paste("mbslr1: ERROR:", ci))
                return(-1)
            }
            y[i] <- rbinom(1, 1, ci)
        }
        return(y)
    }

    if (is.null(n)) {
        n <- 1

        if (!is.list(mu)) {
            meanList <- list()
            meanList[[1]] <- mu
        } else {
            meanList <- mu
        }

        if (!is.list(Sigma)) {
            corList <- list()
            corList[[1]] <- Sigma
        } else {
            corList <- Sigma
        }
    } else if (n == 1) {
        n <- 1

        if (!is.list(mu)) {
            meanList <- list()
            meanList[[1]] <- mu
        } else {
            meanList <- mu
        }

        if (!is.list(Sigma)) {
            corList <- list()
            corList[[1]] <- Sigma
        } else {
            corList <- Sigma
        }
    } else {
        meanList <- mu
        corList <- Sigma
    }



    # Simulate correlated binary outcomes
    y <- NULL
    for (i in 1:n) {
        r <- corList[[i]]
        u <- meanList[[i]]

        v <- cor2var(r, u)
        oor <- chkbinc(r, u)
        if (oor) {
            stop("ERROR: Correlation out of range for given mean")
        }
        b <- allreg(v)
        y <- c(y, mbslr1(b, u))
    }
    return(y)
}
