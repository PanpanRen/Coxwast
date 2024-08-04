# Coxwast
Subgroup testing method for change-plane Cox models.

Weighted-sum-based subgroup testing method is proposed for change-plane Cox models. Different testing methods are also provided in this package. 

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("PanpanRen/Coxwast")

# Usage

   - [x] [Coxwast-manual.pdf](https://github.com/PanpanRen/Coxwast/inst/Coxwast-manual.pdf) ---------- Details of the usage of the package.
# Example
    library(Coxwast)
    library(survival)
    library(simsurv)
    library(boot)
    library(mvtnorm)
    library(parallel)

    n = 100
    p1 = 2
    p2 = 1
    p3 = 3
    alpha = rep(1, p1)
    beta  = rep(1, p2)/2
    gamma = c(1, seq(-1,1,length.out = p3-1)) 
    rho = 0.3
    cenRate = 0.2
    set.seed(100)
    data = generate_cox_data(n, alpha, beta, gamma, rho, cenRate = cenRate)
    fit <- WAST(data)
    fit$pval


# References
Deng, Y., Cai, J., and Zeng, D. (2022). Maximum Likelihood Estimation for Cox Proportional Hazards Model with a Change Hyperplane. Statistica Sinica, 32(2), 983.

Kang, S., Lu, W., and Song, R. (2017). Subgroup detection and sample size calculation with proportional hazards regression for survival data. Statistics in medicine, 36(29), 4646-4659.

Ren, P., Zhang, X., Li, Y., and Liu, X. (2024). Weighted-sum-based subgroup testing procedure for change-plane Cox models. Manuscript.
