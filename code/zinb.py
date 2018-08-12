# zero-inflated negative binomial model
# based on Ahbishek's implementation
# https://users.rcc.uchicago.edu/~aksarkar/singlecell-qtl/zinb.html

import numpy as np
import scipy.optimize as so
import scipy.special as sp
import pandas as pd

def log(x):
    return np.log(x + 1e-8)

def sigmoid(x):
    lim = np.log(np.finfo(np.float64).resolution)
    return np.clip(sp.expit(x), lim, -lim)

def nb(theta, x, size):
    mean, inv_disp = np.exp(theta)
    mean *= size
    assert mean.shape == x.shape
    return (x * log(mean / inv_disp) -
            x * log(1 + mean / inv_disp) -
            inv_disp * log(1 + mean / inv_disp) +
            sp.gammaln(x + inv_disp) -
            sp.gammaln(inv_disp) -
            sp.gammaln(x + 1))

def _nb(theta, x, size):
    return -nb(theta, x, size).mean()

def zinb(theta, x, size):
    theta, logodds = theta[:2], theta[2]
    case_non_zero = -np.log1p(np.exp(logodds)) + nb(theta, x, size)
    case_zero = np.logaddexp(logodds - np.log1p(np.exp(logodds)), case_non_zero)
    return -np.where(x < 1e-8, case_zero, case_non_zero).mean()

def _fit_gene(chunk):
    if chunk[:,0].sum() == 0:
        # Without stringent QC, we need to handle the case of all zero observations
        return [-np.inf, -np.inf, np.inf, 0, 0]
    res0 = so.minimize(_nb, x0=[0, 0], args=(chunk[:,0], chunk[:,1]))
    pi0 = (chunk[:,0] == 0).sum() / chunk.shape[0]
    res = so.minimize(zinb, x0=list(res0.x) + [sp.logit(pi0 + 1e-8)], args=(chunk[:,0], chunk[:,1]))
    if res0.fun < res.fun:
    # This isn't a likelihood ratio test. Numerically, our implementation of
    # ZINB can't represent pi = 0, so we need to use a separate implementation
    # for it
        log_mu, neg_log_phi = res0.x
        logit_pi = -np.inf
    else:
        log_mu, neg_log_phi, logit_pi = res.x
    mean_by_sample = chunk[:,1].ravel() * np.exp(log_mu)
    var_by_sample = mean_by_sample + np.square(mean_by_sample) * np.exp(-neg_log_phi)
    mean_by_ind = mean_by_sample.mean()
    var_by_ind = (np.square(mean_by_sample - mean_by_ind) + var_by_sample).mean()
    return [log_mu, -neg_log_phi, logit_pi, mean_by_ind, var_by_ind]

# for each gene
# y_i ~ Pois(scale_i * lambda_i)
# lambda_i ~ pi_0*delta_0 + (1-pi_0)*Gamma(alpha, beta) or Gamma(1/phi, 1/(phi*mu))
# output mean, activemean, cv(coef of variation), loglike
def _fit_gene1(y, scale):
    if y.sum() == 0:
        # Without stringent QC, we need to handle the case of all zero observations
        return [1, np.nan, np.nan, 0, 0, np.nan, np.nan]
    res0 = so.minimize(_nb, x0=[0, 0], args=(y, scale))
    pi0 = (y == 0).sum() / len(y)
    res = so.minimize(zinb, x0=list(res0.x) + [sp.logit(pi0 + 1e-8)], args=(y, scale))
    if res0.fun < res.fun:
    # This isn't a likelihood ratio test. Numerically, our implementation of
    # ZINB can't represent pi = 0, so we need to use a separate implementation
    # for it
        log_mu, neg_log_phi = res0.x
        logit_pi = -np.inf
        loglike = -res0.fun*len(y)
    else:
        log_mu, neg_log_phi, logit_pi = res.x
        loglike = -res.fun*len(y)
    mu = np.exp(log_mu)
    pi0 = np.exp(logit_pi)/(np.exp(logit_pi)+1)
    alpha = 1/np.exp(-neg_log_phi)
    beta = alpha/mu
    mean = (1-pi0)*mu
    activemean = mu
    cv = np.sqrt(alpha/beta**2)/activemean
    return pd.Series([pi0, alpha, beta, mean, activemean, cv, loglike])
