#' The print format for cpgeeSWD output
#' @description The print format for cpgeeSWD output
#' @param x The object of cpgeeSWD output
#' @param ... further arguments passed to or from other methods
#' @return The output from \code{\link{print}}
#' @rdname print.cpgeeSWD
#' @export
#'
#'
print.cpgeeSWD <- function(x, ...) {
    p <- length(x$beta)
    q <- length(x$alpha)

    niter <- x$niter
    outbeta <- x$outbeta
    outalpha <- x$outalpha

    cat(
        "GEE and MAEE for Cluster-Period Summaries", "\n", "Number of Iterations:",
        niter, "\n"
    )


    cat("Results for marginal mean parameters \n")
    print(outbeta)
    cat("\n")

    cat("Results for correlation parameters \n")
    print(outalpha)
    cat("\n")
}
