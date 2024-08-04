#' Function to select initial gamma values spanning its space
#'
#' @param n.initials Number of gamma values.
#' @param Z The grouping variables.
#' @param lb.quantile The lower quantile specified for Z%*%gamma.initials.
#' @param ub.quantile The upper quantile specified for Z%*%gamma.initials.
#' @param ss A positive integer with n.initials/ss indicating how many sets of gamma are chosen, default to 1.
#'
#' @return A matrix. A set of gamma values selected for defining subgroup.
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
#' K = 1000
#' qlb = 0.1
#' cols = apply(data$u, 2, var) != 0
#' Gamma = gam.init(K, data$u[,cols], lb.quantile=qlb, ub.quantile=1-qlb, ss=1)

gam.init.sub = function(q,n.initials){
    out = matrix(rnorm(n.initials*q), n.initials, q)
    dd 	= apply(out^2,1,sum)
    out	= out/sqrt(dd)
    return(out)
}

gam.init = function(n.initials, Z, lb.quantile, ub.quantile, ss=1){
    q = ifelse(is.matrix(Z), ncol(Z), 1)
    if(q==1){
        gamma.initials = matrix(1,n.initials,q+1)
        gamma.initials[,1] = -quantile(Z,seq(lb.quantile,ub.quantile,length=n.initials))
    }else{
        gamma.initials = gam.init.sub(q, n.initials/ss)
        Z.gamma.initials = Z %*% t(gamma.initials)
        ll=round(n.initials/ss)
        qtile = sample(seq(lb.quantile,ub.quantile,length=n.initials),n.initials)
        gamma.initials.1 = sapply(1:n.initials,function(x)return(
                            -quantile(Z.gamma.initials[,x-floor((x-0.1)/ll)*ll],qtile[x])
                        ))

        gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,min),gamma.initials.1-0.001,gamma.initials.1)
        gamma.initials.1=ifelse(gamma.initials.1==(-1)*apply(Z.gamma.initials,2,max),gamma.initials.1+0.001,gamma.initials.1)
        gamma.initials.aug=do.call("rbind", rep(list(gamma.initials), ss))
        gamma.initials = cbind(gamma.initials.1, gamma.initials.aug)
    }
    return(gamma.initials)
}