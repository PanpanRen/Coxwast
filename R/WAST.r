#' The weighted-sum-based subgrouping testing method for change-plane Cox models
#'
#' @param data A list, including \eqn{y} (time response), \eqn{x} (predictors), \eqn{z} (predictors), \eqn{u} (change-plane variables), status (censoring indicator).
#' @param B A constant. Number of bootstrap samples. Default is 1000.
#' @param CoreNum A constant. Number of processes to be used in parallel operation. Default is 4.
#' @param par A logical value. Whether to conduct the parallel operation. Default is TRUE.
#'
#' @return A list.
#' \itemize{
#'   \item alpha_hat - The estimator of alpha under the null hypothesis.
#'   \item test_statistic - The value of test statistic.
#'   \item test_statistic_boot - B values of test statistic obtained from the bootstrap.
#'   \item pval - The p-value of the test.
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
#' fit <- WAST(data)
#' fit$pval

WAST <- function(data, B = 1000, CoreNum = 4, par = TRUE){
    y 	= data$y
    x 	= data$x
    z   = data$z
    u   = data$u
    status = data$status
    n = length(y)
    tic = proc.time()
    
    if (is.vector(u)) {
        u = as.matrix(u, ncol = 1)
    }
    if (is.vector(z)) {
        z = as.matrix(z, ncol = 1)
    }
    if (is.vector(x)) {
        x = as.matrix(x, ncol = 1)
    }
    
    y.order = order(y)
    x       = x[y.order, ,drop=FALSE]
    z       = z[y.order, ,drop=FALSE]
    u       = u[y.order, ,drop=FALSE]
    y       = y[y.order]
    status  = status[y.order]
    
    data = list(
        y=y, 
        x=x, 
        z=z,
        u=u,
        status=status
    )
    
    param = c(dim(x), ncol(z), ncol(u))
    
    datasur <- list( time=y, status=status, tx=x )
    fit_cox = coxph(Surv(time, status) ~ tx, datasur)
    alpha_hat <- coef(fit_cox)
    
    yrmax = rank(y, ties.method="max") - 1
    yrmin = rank(y, ties.method="min") - 1
    
    fit <- .Call(
        "_COX_Subgroup_Test",
        as.integer(yrmax),
        as.integer(yrmin),
        as.numeric(t(x)),
        as.numeric(t(z)),
        as.numeric(t(u)),
        as.integer(status),
        as.numeric(alpha_hat),
        as.integer(param))
    
    fit$omega = matrix(fit$omega, nrow=n, ncol=n, byrow = TRUE)
    fit$alpha_hat = alpha_hat
    fit$data = data


    n = length(fit$data$y)
    
    df <- data.frame(
        time = fit$data$y,
        status = fit$data$status,
        X = fit$data$x,
        Z = fit$data$z,
        U = fit$data$u
    )
    
    get_ts <- function(df) {
        y   = df$time
        x 	= as.matrix(df[,grep("X", colnames(df))])
        z   = as.matrix(df[,grep("Z", colnames(df))])
        u   = as.matrix(df[,grep("U", colnames(df))])
        status = df$status
        
        y.order = order(y)
        x       = x[y.order, ,drop=FALSE]
        z       = z[y.order, ,drop=FALSE]
        u       = u[y.order, ,drop=FALSE]
        y       = y[y.order]
        status  = status[y.order]
        
        param = c(dim(x), ncol(z), ncol(u))
        
        datasur <- list( time=y, status=status, tx=x )
        fit_cox = coxph(Surv(time, status) ~ tx, datasur)
        alpha_hat <- coef(fit_cox)
        
        yrmax = rank(y, ties.method="max") - 1
        yrmin = rank(y, ties.method="min") - 1
        
        result <- .Call(
            "_COX_Subgroup_Test",
            as.integer(yrmax),
            as.integer(yrmin),
            as.numeric(t(x)),
            as.numeric(t(z)),
            as.numeric(t(u)),
            as.integer(status),
            as.numeric(alpha_hat),
            as.integer(param))
        
        result$alpha_hat = alpha_hat
        
        return(
            c(result$test_statistic, result$alpha_hat)
        )
    }
    
    # related input
    fit_cox = coxph(Surv(y, status) ~ x, fit$data)
    fit_surv = survfit(fit_cox)
    fit_cens = survfit(Surv(y-0.001*status, 1-status)~1, data = fit$data)

    if (par) {
        clnum <- min(detectCores(), CoreNum) 

        boot5 = censboot(
            data = df, statistic = get_ts, R = B, F.surv = fit_surv, G.surv = fit_cens, 
            cox = fit_cox, sim = "model", parallel = "multicore", ncpus = clnum)
        boot5 = boot5$t
    } else {
        boot5 = censboot(
            data = df, statistic = get_ts, R = B, F.surv = fit_surv, G.surv = fit_cens, 
            cox = fit_cox, sim = "model")
        boot5 = boot5$t
    }
    
    tsVec4 = boot5[,1]
    alphaMat4 = boot5[,-1]
    pval4 = mean(tsVec4 > fit$test_statistic)
    
    boot = list(
        alpha = alphaMat4,
        ts = tsVec4,
        pval = pval4
    )

    toc = proc.time()
    output = list(
        alpha_hat = fit$alpha_hat,
        test_statistic = fit$test_statistic,
        test_statistic_boot = boot$ts,
        pval = boot$pval,
        time = toc[3] - tic[3]
    )
    return(output)
}