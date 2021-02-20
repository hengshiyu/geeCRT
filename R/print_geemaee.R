#' The print format for geemaee output
#' 
#' @description The print format for geemaee output
#' @param x The object of geemaee output
#' @param ... further arguments passed to or from other methods
#' @return The output from \code{\link{print}}
#' @rdname print.geemaee
#' @export
#'
#'
print.geemaee <- function(x, ...) {
  
    p = length(x$beta)
    q = length(x$alpha)
    niter = x$niter
    outbeta = x$outbeta
    outalpha = x$outalpha

    cat("GEE for correlated Gaussian data","\n",
        "Number of Iterations:",niter,"\n")


    cat("Results for marginal mean parameters \n")
    print(outbeta)
    cat("\n")

    cat("Results for correlation parameters \n")
    print(outalpha)
    cat("\n")
}

