library(ashr)
library(descend)

##############
lik_pois = function(y, scale=1, link=c("identity","log")){
  link = match.arg(link)
  if (link=="identity"){
    list(name = "pois",
         const = TRUE,
         lcdfFUN = function(x){stats::pgamma(abs(x),shape=y+1,rate=scale,log.p=TRUE)-log(scale)},
         lpdfFUN = function(x){stats::dgamma(abs(x),shape=y+1,rate=scale,log=TRUE)-log(scale)},
         etruncFUN = function(a,b){-my_etruncgamma(-b,-a,y+1,scale)},
         e2truncFUN = function(a,b){my_e2truncgamma(-b,-a,y+1,scale)},
         data=list(y=y,link=link))
  }else if (link=="log"){
    y1 = y+1e-5 # add pseudocount
    list(name = "pois",
         const = TRUE,
         lcdfFUN = function(x){stats::pgamma(exp(-x),shape=y1,rate=scale,log.p=TRUE)-log(y1)},
         lpdfFUN = function(x){stats::dgamma(exp(-x),shape=y1,rate=scale,log=TRUE)-log(y1)},
         etruncFUN = function(a,b){-my_etruncgamma(exp(-b),exp(-a),y1,scale)},
         e2truncFUN = function(a,b){my_e2truncgamma(exp(-b),exp(-a),y1,scale)},
         data=list(y=y,link=link))
  }
}

# For unimix object (uniform mixture distribution) g,
# compute its density on grid x
dens_unimix = function(g, x){
  sapply(x, dens_unimix_sing, pi=g$pi, a=g$a, b=g$b)
}
dens_unimix_sing = function(x,pi,a,b){
  sum((x>=a & x<b)/(b-a)*pi,na.rm=TRUE)
}

# Compute mean for unimix object g
mean_unimix = function(g, include_g0=TRUE){
  if (include_g0==FALSE){
    g$pi = g$pi[-1]
    g$pi = g$pi/sum(g$pi)
    g$a = g$a[-1]
    g$b = g$b[-1]
  }
  sum(g$pi*(g$a+g$b)/2)
}

# Compute mean for a unimix object g
sd_unimix = function(g){
  sqrt(sum(g$pi*(g$b-g$a)^2/12)+sum(g$pi*(g$b+g$a)^2/4)-mean_unimix(g)^2)
}

# Compute coefficient of variation (CV) for unimix object g
cv_unimix = function(g){
  sd_unimix(g)/mean_unimix(g)
}


# Number of interior modes for vector x
n_modes = function(x){
  l = length(x)
  sum((x[2:(l-1)]-x[1:(l-2)]>0) & (x[2:(l-1)]-x[3:l]>0))
}

# Deconvolution for single vector y
# Suppose y~Pois(scale*lambda), lambda~g,
# use DESCEND, Poisson ASH and nonparametric deconvolution
# to estimate the unknown prior g
# Input: 
#   data: gene expression matrix, rows are genes and columns are samples
#   scale: scaling factors for samples (eg. library size)
#   output.fit: flag, if output the DESCEND/ASH fitted object
#   plot.dens: flag, if plot the fitted prior g's density (exclude delta_0)
#   plot.cdf: flag, if plot the fitted prior g's cdf
# Output: properties of fitted g
#   pi0: mixture proportion of the pointmass at zero (delta_0)
#   mean: mean of fitted g
#   activemean: mean of the non-zero components, so mean=(1-pi0)*activemean
#   cv: coefficient of variation of fitted g
#   loglike: log-likelihood log(P(y;g))
deconv_sing = function(y, scale, output.fit=FALSE, plot.dens=FALSE, plot.cdf=FALSE){
  y = as.numeric(y)
  
  ##### Nonparametric deconvolution (fitted by ASH)
  
  # set uniform mixture components for g:
  # [0,c],[c,2c],[2c,3c]..., where c=gridmax/ncomp
  gridmax = (max(y/scale)+0.01)*3  # maximum for the grid
  ncomp = 100  # number of mixture components
  grida = c(0,seq(0, gridmax-gridmax/ncomp, by=gridmax/ncomp)) # left limits
  gridb = c(0,seq(gridmax/ncomp, gridmax, by=gridmax/ncomp)) # right limits
  g.nonp = unimix(pi=rep(1,ncomp+1)/(ncomp+1), a=grida, b=gridb)
  
  # fit ASH with given prior grid 
  fit_ashnonp = try(ash.workhorse(rep(0, length(y)), 1,
                              lik = lik_pois(y, scale=scale, link="identity"), 
                              prior="uniform", pointmass=FALSE,
                              g = g.nonp, control=list(maxiter=20000)))
  if (class(fit_ashnonp)!="try-error"){
    g_ashnonp = fit_ashnonp$fitted_g
    pi0_ashnonp = g_ashnonp$pi[1]
    mean_ashnonp = mean_unimix(g_ashnonp)
    activemean_ashnonp = mean_unimix(g_ashnonp, include_g0 = FALSE)
    cv_ashnonp = cv_unimix(g_ashnonp)
    loglike_ashnonp = fit_ashnonp$loglik
  }else{
    pi0_ashnonp = NA
    mean_ashnonp = NA
    activemean_ashnonp = NA
    cv_ashnonp = NA
    loglike_ashnonp = NA
  }
  
  ##### Poisson ASH (unimodal prior with point mass at 0)
  
  # first fit the model with only unimodal prior (non-zero mode estimated by ASH)
  fit_ashuni = try(ash.workhorse(rep(0, length(y)), 1,
                              lik = lik_pois(y, scale=scale, link="identity"),
                              mixcompdist = "halfuniform", prior="uniform", pointmass=FALSE))
  # fit the model with unimodal prior (mode at 0)
  fit_ashuni0 = try(ash.workhorse(rep(0, length(y)), 1,
                               lik = lik_pois(y, scale=scale, link="identity"),
                               mixcompdist = "+uniform", prior="uniform", pointmass=TRUE))
  
  if (class(fit_ashuni0)!="try-error"){
    g_ashuni0 = fit_ashuni0$fitted_g
    pi0_ashuni0 = g_ashuni0$pi[1]
    mean_ashuni0 = mean_unimix(g_ashuni0)
    activemean_ashuni0 = mean_unimix(g_ashuni0, include_g0 = FALSE)
    cv_ashuni0 = cv_unimix(g_ashuni0)
    loglike_ashuni0 = fit_ashuni0$loglik
  }else{
    pi0_ashuni0 = NA
    mean_ashuni0 = NA
    activemean_ashuni0 = NA
    cv_ashuni0 = NA
    loglike_ashuni0 = NA
  }
  
  if (class(fit_ashuni)!="try-error"){
    g_ashuni = fit_ashuni$fitted_g
    pi0_ashuni = NA 
    mean_ashuni = mean_unimix(g_ashuni)
    activemean_ashuni = mean_unimix(g_ashuni, include_g0 = FALSE)
    cv_ashuni = cv_unimix(g_ashuni)
    loglike_ashuni = fit_ashuni$loglik
    
    # then fit the model with unimodal prior and point mass at 0
    g = unimix(rep(1/(length(g_ashuni$a)+1), length(g_ashuni$a)+1), 
               c(0, g_ashuni$a), c(0, g_ashuni$b)) # add delta_0
    
    fit_ash = ash.workhorse(rep(0, length(y)), 1,
                            lik = lik_pois(y, scale=scale, link="identity"), 
                            g=g, prior="uniform")
    g_ash = fit_ash$fitted_g
    
    pi0_ash = g_ash$pi[1]
    mean_ash = mean_unimix(g_ash)
    activemean_ash = mean_unimix(g_ash, include_g0 = FALSE)
    cv_ash = cv_unimix(g_ash)
    loglike_ash = fit_ash$loglik
  }else{
    pi0_ashuni = NA
    mean_ashuni = NA
    activemean_ashuni = NA
    cv_ashuni = NA
    loglike_ashuni = NA
    
    pi0_ash = NA
    mean_ash = NA
    activemean_ash = NA
    cv_ash = NA
    loglike_ash = NA
  }
  
  ##### DESCEND
  
  fit_descend = try(deconvSingle(y, scaling.consts = scale, plot.density = FALSE,
                                 control=list(max.sparse=c(1,1))))
  
  if (class(fit_descend)=="DESCEND"){
    pi0_descend = 1-fit_descend@estimates[1,1]
    mean_descend = fit_descend@estimates[3,1]
    activemean_descend = fit_descend@estimates[2,1]
    cv_descend = fit_descend@estimates[4,1]
    nmodes_descend = n_modes(fit_descend@density.points[,2])

    # DESCEND outputs the unpenalized log-likelihood for the full model
    # and non-zero model (assuming pi0=0)
    msg = try(capture.output(deconvSingle(y, scaling.consts = scale, 
                                          plot.density = FALSE,control=list(max.sparse=c(1,1)))))
    msg = try(msg[grep("neg",msg)])
    loglike_descfull = try(-as.numeric(strsplit(msg," ")[[1]][6]))
    loglike_descnon0 = try(-as.numeric(strsplit(msg," ")[[2]][6]))
    if (class(loglike_descfull)=="try-error"){
      loglike_descfull = NA
    }
    if (class(loglike_descnon0)=="try-error"){
      loglike_descnon0 = NA
    }
    
    # DESCEND also outputs the likelihood ratio test 
    # p-value comparing the full model and non-zero model
    LRTpval_descend = fit_descend@pval
  }else{
    pi0_descend = NA
    mean_descend = NA
    activemean_descend = NA
    cv_descend = NA
    nmodes_descend = NA
    LRTpval_descend = NA
    loglike_descfull = NA
    loglike_descnon0 = NA
  }
  
  # plot density of fitted prior (exclude delta_0)
  if (plot.dens==TRUE){
    x = fit_descend@density.points[-1,1]
    dens_ash = dens_unimix(g_ash, x)
    plot(x, dens_ash, type='l',
         main="Fitted H", xlab="x", ylab="density",
         ylim=c(0,max(dens_ash)*1.5))
    lines(x, fit_descend@density.points[-1,2], col=3)
    legend("topright", legend=c("Poisson ash", "DESCEND") ,
           col=c(1,3), lty=1)
  }
  
  # plot cdf of fitted prior
  if (plot.cdf==TRUE){
    x = fit_descend@density.points[-1,1]
    dens_descend = fit_descend@density.points[-1,2]
    middens = (dens_descend+c(dens_descend[1], dens_descend[1:(length(x)-1)]))/2
    cdf_descend = cumsum((x[1:length(x)]-c(0,x[1:(length(x)-1)]))*middens)+pi0_descend
    
    plot(x, cdf.ash(fit_ash, x)$y, type='l', ylim=c(0,1),
         main="Fitted G", xlab="x", ylab="cdf")
    lines(x, cdf_descend, col=3)
    legend("bottomright", legend=c("Poisson ash", "DESCEND") ,col=c(1,3), lty=1)
  }
  
  res = cbind(pi0_ashnonp, mean_ashnonp, activemean_ashnonp, cv_ashnonp, loglike_ashnonp,
              pi0_ash, mean_ash, activemean_ash, cv_ash, loglike_ash,
              pi0_ashuni, mean_ashuni, activemean_ashuni, cv_ashuni, loglike_ashuni,
              pi0_ashuni0, mean_ashuni0, activemean_ashuni0, cv_ashuni0, loglike_ashuni0,
              pi0_descend, mean_descend, activemean_descend, cv_descend, nmodes_descend, LRTpval_descend,
              loglike_descfull, loglike_descnon0
              )
  
  if (output.fit==TRUE){
    return(list(res=res, fit_descend=fit_descend, 
                fit_ash=fit_ash, fit_ashuni=fit_ashuni, fit_ashuni0=fit_ashuni0,
                fit_ashnonp=fit_ashnonp))
  }else{
    return(res)
  }
}

# Deconvolution for gene expression matrix
# Input: 
#   data: gene expression matrix, rows are genes and columns are samples
#   scale: scaling factors for samples (eg. library size)
#   plot.dens: flag, if plot the fitted prior g's density (exclude delta_0)
#   plot.cdf: flag, if plot the fitted prior g's cdf
# Output: properties of fitted g
#   pi0: mixture proportion of the pointmass at zero (delta_0)
#   mean: mean of fitted g
#   activemean: mean of the non-zero components, so mean=(1-pi0)*activemean
#   cv: coefficient of variation of fitted g
#   loglike: log-likelihood log(P(y;g))
deconv_data = function(data, scale=NULL, plot.dens=FALSE, plot.cdf=FALSE){
  # normalized library sizes as scaling factors
  if (missing(scale)){
    libsize = colSums(data)
    scale = libsize/mean(libsize)
  }
  res = apply(data, 1, deconv_sing, scale=scale, 
              only_nonparam=only_nonparam, 
              plot.dens=plot.dens, plot.cdf=plot.cdf)
  
  row.names(res) = c("pi0_ashnonp", "mean_ashnonp", "activemean_ashnonp", "cv_ashnonp", "loglike_ashnonp",
                     "pi0_ash", "mean_ash", "activemean_ash", "cv_ash", "loglike_ash",
                     "pi0_ashuni", "mean_ashuni", "activemean_ashuni", "cv_ashuni", "loglike_ashuni",
                     "pi0_ashuni0", "mean_ashuni0", "activemean_ashuni0", "cv_ashuni0", "loglike_ashuni0",
                     "pi0_descend", "mean_descend", "activemean_descend", "cv_descend", 
                     "nmodes_descend", "LRTpval_descend", "loglike_descfull", "loglike_descnon0")
  return(t(res))
}
