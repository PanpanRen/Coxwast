#' The score testing method for for change-plane Cox models similar in Kang et al.(2017).
#'
#' @param data A list, including \eqn{y} (time response), \eqn{x} (predictors), \eqn{z} (predictors), \eqn{u} (change-plane variables), status (censoring indicator).
#' @param Gamma A matrix for initial gamma values. If not given then it will be automatically generated based on the data.
#' @param B A constant. Number of bootstrap samples. Default is 1000.
#' @param K A constant. The number of the initial gamma values. Default is 1000.
#' @param qlb A constant. The lower quantile specified for Z%*%gamma.initials. Default is 0.1.
#' @param tol A constant. The precision of the Newton method. Default is \code{1e-8}.
#'
#' @return A list.
#' \itemize{
#'   \item TestR - The value of test statistic.
#'   \item TestB - B values of test statistic obtained from the bootstrap.
#'   \item Pvalue - The p-value of the test.
#'   \item time - Running time.
#' }
#' @export
#'
#' @examples
#' n = 100
#' p1 = 2
#' p2 = 1
#' p3 = 3
#' alpha = rep(1, p1)
#' beta  = rep(1, p2)/2
#' gamma = c(1, seq(-1,1,length.out = p3-1)) 
#' rho = 0.3
#' cenRate = 0.2
#' set.seed(100)
#' data = generate_cox_data(n, alpha, beta, gamma, rho, cenRate = cenRate)
#' fit <- ST(data)
#' fit$Pvalue

ST <- function(data, Gamma, B = 1000, K = 1000, qlb = 0.1, tol = 1e-8) {
    tic = proc.time()
    
    X = data$x
    Y = data$y
    Status = data$status
    Z = data$z
    U = data$u
    
    if (missing(Gamma)) {
        cols = apply(U, 2, var) != 0
        Gamma = gam.init(K, U[,cols], lb.quantile=qlb, ub.quantile=1-qlb, ss=1)
        rm(cols)
    }
    
    n  = length(Y)
    px = ncol(X)
    Z  = matrix(Z, n)
    pz = ncol(Z)
    K  = nrow(Gamma)
    p = px + pz
    
    yo 	    = order(Y)
    Y       = Y[yo]
    X       = X[yo,,drop=F]
    Z       = Z[yo,,drop=F]
    U       = U[yo,,drop=F]
    Status	= Status[yo]
    
    yrmax = rank(Y, ties.method="max") - 1
    yrmin = rank(Y, ties.method="min") - 1
    rvmatrix = matrix(rnorm(n*B), B, n)
    index = U%*%t(Gamma) > 0 # index of sample in subgroup
    
    fit = .Call(
        "_calculate_ST00",
        as.integer(Status),
        as.integer(yrmin),
        as.integer(yrmax),
        as.numeric(t(X)),
        as.numeric(t(Z)),
        as.integer(index),
        as.numeric(t(rvmatrix)),
        as.integer(n),
        as.integer(pz),
        as.integer(px),
        as.integer(B),
        as.integer(K),
        as.numeric(tol))
    
    toc = proc.time()
    fit$time = toc[3] - tic[3]
    
    return(fit)
}