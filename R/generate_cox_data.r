#' Function for generating data from Cox proportional hazard model with a change plane.
#'
#' @param n A constant. The sample size.
#' @param alpha A vector. The true parameter for baseline covariates.
#' @param beta A vector. The true parameter denoting the heterogeneous effect of the subgroup.
#' @param gamma A vector. The true parameter for grouping variables.
#' @param rho A constant. The strength of correlation among covariates.
#' @param cenRate A constant. Censoring rate. Default is 0.1.
#' @param censortype Censroing type, including "RightCensor" and "RandomCensor".
#'
#' @return A list
#' \itemize{
#'   \item y - A length \eqn{n} vector. The survival time.
#'   \item x - A matrix. The baseline covariates.
#'   \item z - A matrix. The baseline covariates.
#'   \item u - A matrix. The grouping variables.
#'   \item status - A length \eqn{n} vector. Censoring indicator.
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

generate_cox_data <- function(
    n, alpha, beta, gamma, rho, cenRate = 0.1, censortype = c("RightCensor", "RandomCensor")
) {
    censortype = match.arg(censortype)
    
    p1 = length(alpha)
    p2 = length(beta)
    p3 = length(gamma) - 1

    p = max(p1, p2) + p3
    sigma = diag(rep(1-rho,p)) + rho
    RV = mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=sigma, method="chol")

    x = matrix(as.numeric(RV[,1:p1,drop=F]>0), n, p1)
    z = matrix(as.numeric(RV[,1:p2,drop=F]>0), n, p2)
    u = cbind(1, RV[,-c(1:max(p1,p2)),drop=F])
    
    etaX = x %*% alpha
    etaZ = z %*% beta
    
    # generate U
    qz = 0.6
    etaU = u[,-1,drop=FALSE]%*%gamma[-1]
    threshold = quantile(etaU, probs = qz)
    etaU = (etaU>threshold)

    eta = as.numeric(etaX + etaZ*etaU)
    
    # error   ==================================================
    # the natural log of a Weibull random time is an extreme value random observation.   
    e = log(rweibull(n, shape = 1, scale = 1)) - digamma(1)
    y0 = exp(-eta + e)
    
    # censor  ==================================================
    if (is.null(cenRate)){
        y = y0
        status = rep(1, n)
    }else{
        if (censortype == "RightCensor") {
            cens = quantile(y0, 1-cenRate)
            y = pmin(y0, cens)
            status = 1 * (y0<=cens)
        } else if (censortype == "RandomCensor") {
            status     = rbinom(n = n, prob = 1-cenRate, size = 1)
            y         = ys
        }
    }
    
    # In the transformation model, y may be infinity.
    y[which(y=="NaN"|y=="Inf")] =  max(y[which(y!="NaN"&y!="Inf")]) 

    return(list(
        y=y, 
        x=x, 
        z=z,
        u=u,
        status=status
    ))
}