# R functions to compute O/E values, fit simple log-linear models 
# by maximum likelihood (ML) or maximum a-posteriori (MAP) estimation,
# and compare models with the likelihood-ratio test (LRT), AIC, BIC, or the Laplace approximation.
# A short example appears below the function definitions.
# Colin Wilson, Johns Hopkins University, colin.wilson@jhu.edu
# February 2022

require(tidyverse)
require(matrixStats)

# # # # # Functions # # # # #

# Observed counts divided by counts that are expected,  
# given the assumption of complete independence of 
# rows and columns, in a two-dimensional table dat.
OE = function(dat) {
    n = sum(dat)
    O = dat / n
    E = outer(rowSums(O), colSums(O))
    return (O/E)
}

# Probability distribution defined by the log-linear 
# model with design matrix X and parameters lambda.
# Optionally returns log probabilities.
loglinear = function(lambda, X, logprob=TRUE) {
    log_pstar = X %*% lambda
    log_Z = logSumExp(log_pstar)
    log_p = log_pstar - log_Z
    if (logprob)
        return (log_p)
    return (exp(log_p))
}

# Negative log-likelihood of the log-linear model with 
# design matrix X and parameters lambda given data 
# (vector of counts) dat.
negLogLik = function(lambda, X, dat) {
    log_p = loglinear(lambda, X)
    L = -sum(dat * log_p)
    return (L)
}

# Negative log posterior of the log-linear model 
# with design matrix X and parameters lambda given 
# data (vector of counts) dat and isotropic Gaussian 
# prior with precision (inverse variance) l2; the prior 
# is equivalent to L2 regularization with scale l2/2.
negLogPost = function(lambda, X, dat, l2) {
    L = negLogLik(lambda, X, dat)
    P = l2/2 * sum(lambda^2)
    return (L + P)
}

# Matrix of second derivatives (Hessian) of the 
# negative log posterior of the log-linear model 
# with design matrix X, evaluated at parameters lambda, 
# given data (vector of counts) dat and Gaussian prior 
# with precision (inverse variance) l2.
# Compare with output of optim(..., hessian=TRUE).
dd_negLogPost = function(lambda, X, dat, l2) {
    # Probability distribution give lambda, X
    p = loglinear(lambda, X, logprob=FALSE)
    # Double derivatives of negative log-likelihood
    k = ncol(X)
    ddL = matrix(0, nrow=k, ncol=k)
    for (i in 1:k) {
        for (j in 1:k) {
            # Cov(Xi,Xj) = E[Xi * Xj] - E[Xi] * E[Xj]
            ddL[i,j] = sum(p * X[,i] * X[,j]) - sum(p * X[,i]) * sum(p * X[,j])
        }
    }
    ddL = sum(dat) * ddL
    # Double derivatives of neg Gaussian prior
    ddG = l2 * diag(1, k, k)
    # Hessian
    A = ddL + ddG
    return (A)
}

# Numerical fit of log-linear model with design 
# matrix X given data (count vector) dat.
# Return maximum-likelihood estimate if l2 = 0. 
# Otherwise (l2 > 0), return maximum a-posteriori 
# esimate given Gaussian prior that has precision 
# (inverse variance) l2.
loglinear_fit = function(X, dat, l2=0, ...) {
    k = ncol(X)
    lambda0 = rnorm(k, 0, .1)
    fit = optim(
        lambda0,
        negLogPost,
        X=X,
        dat=dat,
        l2=l2,
        method='L-BFGS-B',
        hessian=TRUE,
        ...)
    return (fit)
}

# Laplace approximation for the log-linear model with 
# design matrix X given data (vector of counts) dat 
# and Gaussian prior with precision (inverse variance) l2.
laplace_approx = function(X, dat, l2) {
    # MAP estimation of parameters
    fit = loglinear_fit(X, dat, l2)
    lambda = fit$par

    # Laplace approximation
    k = ncol(X)
    I = diag(nrow=k, ncol=k)

    L = negLogLik(lambda, X, dat)
    A = dd_negLogPost(lambda, X, dat, l2)
    val = -L - l2/2 * sum(lambda^2) -
          log(sqrt(det(I * 1/l2))) - 
          1/2 * log(det(A))
    return (val)
}

# AIC value of model with k free parameters and 
# minimum negative log-likelihood L.
aic = function(L, k) {
    val = 2 * (k + L)
    return (val)
}

# BIC value of model with k free parameters and 
# minimum negative log-likelihood L, fit to 
# data with n observations (i.e., total count).
bic = function(L, k, n) {
    val = k * log(n) + 2 * L
    return (val)
}

# Likelihood ratio test (LRT) for nested models 
# given mimimum negative log-likelihood values 
# (L1, L2) and number of free parameters (k1, k2). 
# Model M1 is assumed to be nested within model M2.
lrt = function(L1, k1, L2, k2) {
    q = 2 * (L1 - L2) # Equiv. -2 * ((-L1) - (-L2))
    df = (k2 - k1)
    return (pchisq(q, df, lower.tail=FALSE))
}

# # # # # Example # # # # #

# Domain of consonant place pairs.
D = expand.grid(
    C1 = c('p', 't', 'k'),
    C2 = c('p', 't', 'k'))
D = D[order(D$C1),]

# Constraints on consonant pairs.
D %>% mutate(
    'p1' = (C1 == 'p'),
    't1' = (C1 == 't'),
    'k1' = (C1 == 'k'),
    'p2' = (C2 == 'p'),
    't2' = (C2 == 't'),
    'k2' = (C2 == 'k'),
    'OCP' = (C1 == C2),
    'OCPlab' = OCP * p1,
    'OCPcor' = OCP * t1,
    'OCPdor' = OCP * k1) %>%
    select(-C1, -C2) -> 
    X
X = 1 * as.matrix(X)

# Hypothetical probability distribution (see Wilson 
# & Obdeyn 2009, p. 14 and Wilson 2022, Tables 3, 5).
theta = c('p1' = 2/3, 't1' = 1, 'k1' = 1/3,
          'p2' = 2/3, 't2' = 1, 'k2' = 1/3,
          'OCP' = 1/2, 'OCPlab' = 1, 'OCPcor' = 1, 'OCPdor' = 1)
lambda = log(theta)

p = loglinear(lambda, X, logprob=FALSE)

# Most probable sample of consonant cooccurrence counts 
# from the hypothetical distribution, with total sample 
# size n and rounding fractional counts.
n = 5000
dat = floor(n * p)

# Observed/Expected values for the sample (see Wilson
# & Obdeyn 2009, p. 14 and Wilson 2022, Table 4).
dat_oe = OE(matrix(dat, nrow=3))
round(dat_oe, 3)

# Fit log-linear model M1 with a single OCP constraint 
# for all same-place consonant combinations.
X1 = X[,c('p1', 'k1', 'p2', 'k2', 'OCP')]
fit1 = loglinear_fit(X=X1, dat=dat)
data.frame(
    constraint = colnames(X1),
    lambda = round(fit1$par, 3),
    theta = round(exp(fit1$par), 3))

# Laplace approximation for the log-linear model with 
# design matrix X given data (vector of counts) dat 
# and Gaussian prior with precision (inverse variance) l2.
laplace_approx = function(X, dat, l2) {
    # MAP estimation of parameters
    fit = loglinear_fit(X, dat, l2)
    lambda = fit$par

    # Laplace approximation
    k = ncol(X)
    I = diag(nrow=k, ncol=k)

    L = negLogLik(lambda, X, dat)
    A = dd_negLogPost(lambda, X, dat, l2)
    val = -L - l2/2 * sum(lambda^2) -
          log(sqrt(det(I * 1/l2))) - 
          1/2 * log(det(A))
    return (val)
}

# AIC value of model with k free parameters and 
# minimum negative log-likelihood L.
aic = function(L, k) {
    val = 2 * (k + L)
    return (val)
}

# BIC value of model with k free parameters and 
# minimum negative log-likelihood L, fit to 
# data with n observations (i.e., total count).
bic = function(L, k, n) {
    val = k * log(n) + 2 * L
    return (val)
}

# Likelihood ratio test (LRT) for nested models 
# given mimimum negative log-likelihood values 
# (L1, L2) and number of free parameters (k1, k2). 
# Model M1 is assumed to be nested within model M2.
lrt = function(L1, k1, L2, k2) {
    q = 2 * (L1 - L2) # Equiv. -2 * ((-L1) - (-L2))
    df = (k2 - k1)
    return (pchisq(q, df, lower.tail=FALSE))
}

# # # # # Example # # # # #

# Domain of consonant place pairs.
D = expand.grid(
    C1 = c('p', 't', 'k'),
    C2 = c('p', 't', 'k'))
D = D[order(D$C1),]

# Constraints on consonant pairs.
D %>% mutate(
    'p1' = (C1 == 'p'),
    't1' = (C1 == 't'),
    'k1' = (C1 == 'k'),
    'p2' = (C2 == 'p'),
    't2' = (C2 == 't'),
    'k2' = (C2 == 'k'),
    'OCP' = (C1 == C2),
    'OCPlab' = OCP * p1,
    'OCPcor' = OCP * t1,
    'OCPdor' = OCP * k1) %>%
    select(-C1, -C2) -> 
    X
X = 1 * as.matrix(X)

# Hypothetical probability distribution (see Wilson 
# & Obdeyn 2009, p. 14 and Wilson 2022, Tables 3, 5).
theta = c('p1' = 2/3, 't1' = 1, 'k1' = 1/3,
          'p2' = 2/3, 't2' = 1, 'k2' = 1/3,
          'OCP' = 1/2, 'OCPlab' = 1, 'OCPcor' = 1, 'OCPdor' = 1)
lambda = log(theta)

p = loglinear(lambda, X, logprob=FALSE)

# Most probable sample of consonant cooccurrence counts 
# from the hypothetical distribution, with total sample 
# size n and rounding fractional counts.
n = 5000
dat = floor(n * p)

# Observed/Expected values for the sample (see Wilson
# & Obdeyn 2009, p. 14 and Wilson 2022, Table 4).
dat_oe = OE(matrix(dat, nrow=3))
round(dat_oe, 3)

# Fit log-linear model M1 with a single OCP constraint 
# for all same-place consonant combinations.
X1 = X[,c('p1', 'k1', 'p2', 'k2', 'OCP')]
fit1 = loglinear_fit(X=X1, dat=dat)
data.frame(
    constraint = colnames(X1),
    lambda = round(fit1$par, 3),
    theta = round(exp(fit1$par), 3))

# Model comparison values for M1.
n = sum(dat)
k1 = length(fit1$par) # Equiv. number of columns in X1
L1 = round(fit1$value, 2)
aic1 = aic(L1, k1)
bic1 = bic(L1, k1, n)
laplace1 = laplace_approx(X=X1, dat=dat, l2=1/10)

# Fit log-linear model M2 with a different OCP constraint 
# for each same-place combination. (Rewrite as weights 
# of place-specific constraints with OCPcor' = OCP, 
# OCPlab' = OCP + OCPlab, OCPdor' = OCP + OCPdor.)
X2 = X[,c('p1', 'k1', 'p2', 'k2', 'OCP', 'OCPlab', 'OCPdor')]
fit2 = loglinear_fit(X=X2, dat=dat)
data.frame(
    constraint = colnames(X2),
    lambda = round(fit2$par, 3),
    theta = round(exp(fit2$par), 3))

# Model comparison values for M2.
n = sum(dat)
k2 = length(fit2$par) # Equiv. number of columns in X2
L2 = round(fit2$value, 2) # Same as L1 in this case
aic2 = aic(L2, k2)
bic2 = bic(L2, k2, n)
laplace2 = laplace_approx(X=X2, dat=dat, l2=1/10)

# Model comparison summary table and likelihood ratio test. 
# The more restrictive/parsimonious model M1, which is true 
# of the underlying population, is preferred by all values. 
data.frame(
    model = c('M1', 'M2'),
    AIC = c(aic1, aic2),    # Lower is better
    BIC = c(bic1, bic2),    # Lower is better
    Laplace = c(laplace1, laplace2) # Higher is better
)
lrt(L1, k1, L2, k2) # Fail to reject M1 in favor of M2

bic1 = bic(L1, k1, n)
laplace1 = laplace_approx(X=X1, dat=dat, l2=1/10)

# Fit log-linear model M2 with separate OCP constraints 
# for each same-place combination. (Recover the weights 
# of place-specific constraints with OCPcor' = OCP, 
# OCPlab' = OCP + OCPlab, OCPdor' = OCP + OCPdor.)
X2 = X[,c('p1', 'k1', 'p2', 'k2', 'OCP', 'OCPlab', 'OCPdor')]
fit2 = loglinear_fit(X=X2, dat=dat)
data.frame(
    constraint = colnames(X2),
    lambda = round(fit2$par, 3),
    theta = round(exp(fit2$par), 3))

# Model comparison values for M2.
n = sum(dat)
k2 = length(fit2$par) # Equiv. number of columns in X2
L2 = round(fit2$value, 2) # Same as L1 in this case
aic2 = aic(L2, k2)
bic2 = bic(L2, k2, n)
laplace2 = laplace_approx(X=X2, dat=dat, l2=1/10)

# Model comparison summary table and likelihood ratio test. 
# The more restrictive/parsimonious model M1, which is true 
# of the underlying population, is preferred by all values. 
data.frame(
    model = c('M1', 'M2'),
    AIC = c(aic1, aic2),    # Lower is better
    BIC = c(bic1, bic2),    # Lower is better
    Laplace = c(laplace1, laplace2) # Higher is better
)
lrt(L1, k1, L2, k2) # Fail to reject M1 in favor of M2
