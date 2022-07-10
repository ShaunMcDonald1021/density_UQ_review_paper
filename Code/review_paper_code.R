# CODE FOR EXAMPLES IN REVIEW PAPER
# (McDonald/Campbell: A Review of Uncertainty Quantification in Density Estimation)
#TODO: Markdown on Github


# Simulate data from a mixture of normals
set.seed(24601)
par(mfrow = c(2,2))
n = 1000 # Sample size
eval_grid = seq(0, 1, by = 0.002)
# coverage_mat = matrix(0, nrow = 100, ncol = 1)
mu = c(0.5, 5/7)
sd = c(1, sqrt(0.1))/7
f0 = function(x) 0.5*(dnorm(x, mu[1], sd[1]) + dnorm(x, mu[2], sd[2]))

mix_comp = sample(1:2, n, replace = TRUE)
X = rnorm(n, mu[mix_comp], sd[mix_comp])

#### KDE methods ####
# First try pointwise intervals based on robust bias correction (Calonico, Catteneo, Farrell 2018)
library(nprobust)
# Unlike the paper, we'll use IMSE-optimal bandwidth instead of pointwise MSE-optimal bandwidth
# We'll also keep the bandwidth constant throughout (i.e. not adjusting it to improve ESS near boundaries)
# Makes it more straightforward to compare to the stuff in the paper imo
ccf_est = kdrobust(X, eval = eval_grid, bwselect = 'imse-dpi', bwcheck = 0)

# Variable-width bootstrap confidence bands based on RBC(Cheng and Chen 2019)
h = ccf_est$Estimate[1,'h']
B = 1000 # Number of bootstrap iterations
piv_quants = numeric(B)
for(b in 1:B){
  X_star = sample(X, n, replace = TRUE) 
  ccf_est_boot = kdrobust(X_star, eval = eval_grid, h = h, bwcheck = 0)
  piv_quants[b] = max(abs(ccf_est_boot$Estimate[,'tau.bc']-ccf_est$Estimate[,"tau.bc"]))
}


plot(eval_grid, ccf_est$Estimate[,"tau.us"], type = 'l', ylim = c(-0.5, 6.65), lty = 2, lwd = 2.5,
     xlab = expression(italic(x)), main = 'KDE methods', ylab = '')
lines(eval_grid, ccf_est$Estimate[,"tau.bc"] + qnorm(0.975)*ccf_est$Estimate[,"se.rb"], lty = 3)
lines(eval_grid, ccf_est$Estimate[,"tau.bc"] - qnorm(0.975)*ccf_est$Estimate[,"se.rb"], lty = 3)
lines(eval_grid, f0(eval_grid), lwd = 1.75)
polygon(c(eval_grid, rev(eval_grid)), 
        c(ccf_est$Estimate[,'tau.bc'] + quantile(piv_quants, 0.95), 
          rev(ccf_est$Estimate[,'tau.bc'] - quantile(piv_quants, 0.95))), 
        border = NA, col = rgb(0.85, 0.85, 0.85))
legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))), 
                             "Pointwise confidence interval", "Confidence band"), 
       lty = c(1, 2, 3, 0), lwd = c(1.75, 2.5, 1, 0), pch = c(NA,NA,NA,22), 
       col = c('black','black','black',rgb(0.85, 0.85, 0.85)), 
       pt.bg = c(NA,NA,NA,rgb(0.85, 0.85, 0.85)), pt.cex = 2)


#### Bernstein polynomial methods ####
# Running the algorithm of Petrone (1999) for a variable-dimension basis expansion density estimate
# Inference is Bayesian and simply comes from MCMC
source('unknown_dim_basis_dens_functions.R')
petrone_est = run_mcmc(list(X), K0 = 4, rate = 1, K_log_prior = numeric(150), 
                       basis_type = 'Bernstein', M = 10, S = 5000, x_min = 0, x_max = 1)
# Note: the first 1000 samples will be discarded as burn-in

# Variable-width credible bands based on MAD's (Edwards, Meyer, and Christensen 2019)
post_med = apply(petrone_est$mcmc_chains$f[1001:5000,,], 2, median)
abs_devs = abs(sweep(petrone_est$mcmc_chains$f[1001:5000,6:496,], 2, post_med[6:496]))
mad = apply(abs_devs, 2, median)
width_quant = quantile(apply(sweep(abs_devs, 2, mad, "/"), 1, max), 0.95)

eval_grid2 = seq(0, 1, by = 0.002)
plot(eval_grid2, post_med, type = 'l', ylim = c(-0.5, 6.65), lty = 2, lwd = 2.5,
     xlab = expression(italic(x)), main = 'Basis expansion methods', ylab = '')
lines(eval_grid2, apply(petrone_est$mcmc_chains$f[1001:5000,,], 2, quantile, 0.025), lty = 3)
lines(eval_grid2, apply(petrone_est$mcmc_chains$f[1001:5000,,], 2, quantile, 0.975), lty = 3)
lines(eval_grid2, f0(eval_grid2), lwd = 1.75)
polygon(c(eval_grid2[6:496], rev(eval_grid2[6:496])), 
        c(post_med[6:496] + width_quant*mad, 
          rev(post_med[6:496] - width_quant*mad)), 
        border = NA, col = rgb(0.85, 0.85, 0.85))
legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))), 
                             "Pointwise credible interval", "Credible band"), 
       lty = c(1, 2, 3, 0), lwd = c(1.75, 2.5, 1, 0), pch = c(NA,NA,NA,22), 
       col = c('black','black','black',rgb(0.85, 0.85, 0.85)), 
       pt.bg = c(NA,NA,NA,rgb(0.85, 0.85, 0.85)), pt.cex = 2)


#### Logspline methods ####
# Using a bootstrap standard error estimate for approximate Gaussian pointwise intervals (Kooperberg and Stone 2003/2004)
library(logspline)
logspline_est = logspline(X, 0, 1, maxknots = 20)

logspline_est_boot_mat = matrix(0, 25, length(eval_grid))

for(b in 1:25){
  X_star = sample(X, n, replace = TRUE)
  logspline_est_boot = logspline(X_star, 0, 1, maxknots = 20)
  logspline_est_boot_mat[b,] = log(dlogspline(eval_grid2, logspline_est_boot))
  #lines(eval_grid, exp(logspline_est_boot_mat[b,]), col = b)
}

st_err = apply(logspline_est_boot_mat, 2, sd)

plot(eval_grid2, dlogspline(eval_grid2, logspline_est), type = 'l', ylim = c(-0.5, 6.65), lty = 2, lwd = 2.5,
     xlab = expression(italic(x)), main = 'Logspline methods', ylab = '')
lines(eval_grid2, exp(log(dlogspline(eval_grid2, logspline_est)) + qnorm(0.975)*st_err), lty = 3)
lines(eval_grid2, exp(log(dlogspline(eval_grid2, logspline_est)) - qnorm(0.975)*st_err), lty = 3)
lines(eval_grid, f0(eval_grid), lwd = 1.75)
legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))), 
                             "Pointwise confidence interval"), 
       lty = c(1, 2, 3), lwd = c(1.75, 2.5, 1), pch = c(NA,NA,NA), 
       col = c('black','black','black'), 
       pt.bg = c(NA,NA,NA), pt.cex = 2)


#### Dirichlet process mixture ####
library(dirichletprocess)
dp = DirichletProcessGaussian((X-mean(X))/sd(X))
dp = Fit(dp, its=5000)

posteriorFit <- sapply(1001:5000, function(i) PosteriorFunction(dp, i)((eval_grid-mean(X))/sd(X)))
plot(eval_grid, apply(posteriorFit, 1, mean)/sd(X), type = 'l', ylim = c(-0.5, 6.65), lty = 2, lwd = 2.5,
     xlab = expression(italic(x)), main = 'DPM methods', ylab = '')
lines(eval_grid, apply(posteriorFit, 1, quantile, 0.975)/sd(X), lty = 3)
lines(eval_grid, apply(posteriorFit, 1, quantile, 0.025)/sd(X), lty = 3)
lines(eval_grid, f0(eval_grid), lwd = 1.75)
legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))), 
                             "Pointwise credible interval"), 
       lty = c(1, 2, 3), lwd = c(1.75, 2.5, 1), pch = c(NA,NA,NA), 
       col = c('black','black','black'), 
       pt.bg = c(NA,NA,NA), pt.cex = 2)


#Actual plotting code
setEPS()
postscript(file = '../Documentation/mcdonald_density_uq_fig.eps', width = 8, height = 8*11/14, family = 'LM Roman 10')
#loadfonts(device = 'win')

par(mfrow = c(2,2), mar = c(2.6, 3.1, 2.1, 0.6), cex.axis = 1.25, cex.lab = 1.25,
    cex.main = 1.25, mgp = c(2.1, 1, 0), family = 'LM Roman 10')

## KDE stuff
plot(eval_grid, ccf_est$Estimate[,"tau.us"], type = 'l', ylim = c(-0.5, 7.125), col = NA,
     ylab = '', xaxt = 'n', xlab = '')
title('KDE methods', line = 0.5)
# Plot bias-corrected simultaneous band first so lines can be overlaid on top
polygon(c(eval_grid, rev(eval_grid)),
        c(ccf_est$Estimate[,'tau.bc'] + quantile(piv_quants, 0.95),
          rev(ccf_est$Estimate[,'tau.bc'] - quantile(piv_quants, 0.95))),
        border = NA, col = rgb(0.85, 0.85, 0.85))

# Plot the KDE itself
lines(eval_grid, ccf_est$Estimate[,"tau.us"], lty = 2, lwd = 2.5)

# Bias-corrected pointwise confidence intervals
lines(eval_grid, ccf_est$Estimate[,"tau.bc"] + qnorm(0.975)*ccf_est$Estimate[,"se.rb"],
      lty = 3, lwd = 1.5)
lines(eval_grid, ccf_est$Estimate[,"tau.bc"] - qnorm(0.975)*ccf_est$Estimate[,"se.rb"],
      lty = 3, lwd = 1.5)

# The true density
lines(eval_grid, f0(eval_grid), lwd = 1.75)

legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))),
                             "Pointwise conf. int.", "Conf. band"),
       lty = c(1, 2, 3, 0), lwd = c(1.75, 2.5, 1.5, 0), pch = c(NA, NA, NA, 22),
       col = c('black', 'black', 'black', NA),
       pt.bg = c(NA, NA, NA, rgb(0.85, 0.85, 0.85)), pt.cex = 2, y.intersp = 1.1)

## Bernstein polynomial stuff
par(mar = c(2.6, 2.1, 2.1, 1.6))
plot(eval_grid, post_med, type = 'l', ylim = c(-0.5, 7.125), col = NA,
     ylab = '', yaxt = 'n', xaxt = 'n', xlab = '')
title('Basis expansion methods', line = 0.5)

# Simultaneous bands based on MAD's, but only over the interval [0.01, 0.99]
# Again, plot this first to overlay lines
polygon(c(eval_grid[6:496], rev(eval_grid[6:496])),
        c(post_med[6:496] + width_quant*mad,
          rev(post_med[6:496] - width_quant*mad)),
        border = NA, col = rgb(0.85, 0.85, 0.85))

# Plot the posterior median
lines(eval_grid, post_med, lty = 2, lwd = 2.5)

# Pointwise credible intervals
lines(eval_grid, apply(petrone_est$mcmc_chains$f[1001:5000,,], 2, quantile, 0.025),
      lty = 3, lwd = 1.5)
lines(eval_grid, apply(petrone_est$mcmc_chains$f[1001:5000,,], 2, quantile, 0.975),
      lty = 3, lwd = 1.5)

# The true density
lines(eval_grid, f0(eval_grid), lwd = 1.75)

legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))),
                             "Pointwise cred. int.", "Cred. band"),
       lty = c(1, 2, 3, 0), lwd = c(1.75, 2.5, 1.5, 0), pch = c(NA,NA,NA,22),
       col = c('black', 'black', 'black', NA),
       pt.bg = c(NA, NA, NA, rgb(0.85, 0.85, 0.85)), pt.cex = 2, y.intersp = 1.1)

## Logspline stuff
par(mar = c(3.6, 3.1, 1.1, 0.6))
# Plot the logspline estimate
plot(eval_grid, dlogspline(eval_grid, logspline_est), type = 'l', ylim = c(-0.5, 7.125),
     lty = 2, lwd = 2.5, xlab = expression(italic(x)), ylab = '')
title('Logspline methods', line = 0.5, xpd = NA)

# Exponentiated pointwise intervals from bootstrap std. err. estimates
lines(eval_grid, exp(log(dlogspline(eval_grid, logspline_est)) + qnorm(0.975)*st_err),
      lty = 3, lwd = 1.5)
lines(eval_grid, exp(log(dlogspline(eval_grid, logspline_est)) - qnorm(0.975)*st_err),
      lty = 3, lwd = 1.5)

# The true density
lines(eval_grid, f0(eval_grid), lwd = 1.75)

legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))),
                             "Pointwise conf. int."),
       lty = c(1, 2, 3), lwd = c(1.75, 2.5, 1.5), col = c('black','black','black'),
       pt.cex = 2, y.intersp = 1.1)

## Dirichlet process mixture stuff
par(mar = c(3.6, 2.1, 1.1, 1.6))
# Plot posterior mean. Note that it must be scaled since data was linearly transformed
plot(eval_grid, apply(posteriorFit, 1, mean)/sd(X), type = 'l', ylim = c(-0.5, 7.125),
     lty = 2, lwd = 2.5, xlab = expression(italic(x)), ylab = '', yaxt = 'n')
title('DPM methods', line = 0.5, xpd = NA)

# Pointwise credible intervals
lines(eval_grid, apply(posteriorFit, 1, quantile, 0.975)/sd(X), lty = 3, lwd = 1.5)
lines(eval_grid, apply(posteriorFit, 1, quantile, 0.025)/sd(X), lty = 3, lwd = 1.5)

# The true density
lines(eval_grid, f0(eval_grid), lwd = 1.75)
legend('topleft', legend = c(expression(italic(f[0])), expression(italic(hat(f))),
                             "Pointwise cred. int."),
       lty = c(1, 2, 3), lwd = c(1.75, 2.5, 1.5), col = c('black','black','black'),
       pt.cex = 2, y.intersp = 1.1)
dev.off()




