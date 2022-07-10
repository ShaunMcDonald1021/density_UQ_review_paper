##### R functions for basis expansion-based densities with shared unknown dimensionality K.
##### Using the sampler (and most of the ideas) from the following two papers:
##### Petrone, S. (1999). Random Bernstein Polynomials. Scandinavian Journal of Statistics 26 373-393.
##### Petrone, S. (1999). Bayesian density estimation using Bernstein polynomials. Canadian Journal of Statistics 27 105-126.
##### While these papers used Bernstein polynomial bases, we also allow for normalized cubic B-spline and
##### constant (i.e. histogram) bases with equally-spaced knots.
##### We also allow for hierarchical fitting of multiple densities.


# TODO: add diagnostic functions using some rstan machinery (PSIS-LOO, IMSE, n_eff, R-hat/split R-hat)
# Also add plotting functions, etc. into here

#### Install these libraries...unless we already did ####
tryCatch(library(MCMCpack), error=function(c) install.packages('MCMCpack'),
         finally=library(MCMCpack))
tryCatch(library(evd), error=function(c) install.packages('evd'),
         finally=library(evd))
tryCatch(library(fda), error=function(c) install.packages('fda'),
         finally=library(fda))
tryCatch(library(Matrix), error=function(c) install.packages('Matrix'),
         finally=library(Matrix))

#### Functions to evaluate basis functions ####
# X = a vector of samples (MUST BE SCALED TO [0, 1] IN ADVANCE)
# Z = vector of continuous latent variables in [0, 1]
# k = dimensionality of the basis
# labs = a vector of discrete mixture component labels, calculated in advance given Z, k

.max_label_one = function(Z, k){
  # Helper function. Returns discrete mixture component labels corresponding to continuous latent var's
  return(pmax(ceiling(k * Z), 1))
}

bernstein_eval = function(X, labs, k){
  return(dbeta(X, labs, k - labs + 1))
}

bspline_eval = function(X, labs, k, inv_ints){
  # inv_ints = a vector containing reciprocals of B-spline integrals
  # (which are required to normalize)
  if(k < 4){ # Must treat this separately since we can't define a cubic basis with < 4 functions
    basisvals = rep(0, length(X))
  }
  else{ # Does the same thing as eval.fd but MUCH faster
    basisvals = t(((t(splineDesign(X, knots = c(rep(0, 3), seq(0, 1, length.out = k - 2),
                                                rep(1, 3)))) *
                      inv_ints))[cbind(labs, 1:length(X))])
  }
  
  return(basisvals)
}

histogram_eval = function(X, labs, k){
  # For each sample in X, nonzero iff the corresponding Z-value is in the same bin
  # Then scale to normalize
  return(as.numeric(.max_label_one(X, k) == labs) * k)
}

#### MCMC sampling functions for each parameter (original Petrone algo) ####
# Using the Monte Carlo/Gibbs hybrid sampling techniques by Petrone (1999)
# In this version, the DP concentration parameter M is fixed a priori and the base distribution
# F0 is assumed to be uniform for simplicity
K_sampler = function(X, Z, K_log_prior, basis_fn, inv_ints = NULL){
  # Sample from "full" conditional of dimensionality K.
  # K_log_prior = logs of prior probabilities for K
  # basis_fn = one of the eval functions defined above, depending on basis choice
  # inv_ints = list containing reciprocals of B-spline integrals for every K-value if we're using those
  Xu = unlist(X)
  Zu = unlist(Z)
  #print(range(Zu - Xu))
  log_vals = log(sapply(1:length(K_log_prior),
                        FUN = function(k) do.call(basis_fn,
                                                  c(list(inv_ints = inv_ints[[k]])[!is.null(inv_ints)],
                                                    list(X = Xu, labs = .max_label_one(Zu, k), k = k)))))
  #print(log_vals)
  #print(.max_label_one(Zu, 150))
  #print(colSums(log_vals))
  K_log_conditional = K_log_prior + colSums(log_vals) # Log values of "full" conditional

  # Use "Gumbel max trick" to get a sample from the correct discrete distribution
  gumbel_samps = rgumbel(n = length(K_log_conditional))
  #print(K_log_conditional+gumbel_samps)
  K = which.max(K_log_conditional + gumbel_samps)
  
  return(K)
}

latent_sampler = function(X_i, Z0_i, M, K, basis_fn, inv_ints = NULL){
  # Sample Z_i (latent variables for i^th density) from "full" conditional.
  # X_i = samples from i^th density (SCALED TO [0, 1])
  # Z0_i = previous value for Z_i
  # TODO: speed up by calculating all necessary basis function evaluations at the beginning,
  # then simply picking out the ones corresponding to the Z-labels you need
  Z_i = Z0_i
  labels = .max_label_one(Z_i, K)
  
  all_probs = outer(X_i, 1:K, function(x, y) do.call(basis_fn,
                                                     c(list(X = x, labs = y, k = K), list(inv_ints = inv_ints)[!is.null(inv_ints)])))
  # Chinese restaurant process-style sampling
  for(l in 1:length(Z_i)){
    assign_probs = all_probs[l, labels]
    assign_probs[l] = M
    
    assignment = sample(1:length(Z_i), size = 1, prob = assign_probs)
    
    # Either sample a new Z_i[l] from a piecewise constant density...
    if(assignment == l){
      subint_probs = all_probs[l,]
      # Must multiply piecewise values by integrals if we're using a B-spline basis
      if(identical(basis_fn, bspline_eval)){
        subint_probs = subint_probs / inv_ints
      }
      
      subinterval = sample(1:K, size = 1, prob = subint_probs)
      Z_i[l] = runif(1, (subinterval - 1) / K, subinterval / K)
      
    }
    else{
      # ...or set it equal to one of the other components of Z_i
      Z_i[l] = Z_i[assignment]
    }
    
    
    # Update the l^th label
    labels[l] = .max_label_one(Z_i[l], K)
  }
  
  return(Z_i)
}

density_sampler = function(Z_i, M, K, inv_ints = NULL){
  # Sample density coefficients for i^th density
  labels = .max_label_one(Z_i, K)
  N = numeric(K)
  
  for(k in 1:K){
    N[k] = sum(labels == k) # N[j] = number of mixture component labels equal to k
  }
  
  # Different prior parameters for Dirichlet iff we're using B-spline bases
  # iff we passed a non-null inv_ints
  if(!is.null(inv_ints)){
    prior_probs = M / inv_ints
  }
  else{
    prior_probs = M / K
  }
  
  phi_i = rdirichlet(1, N + prior_probs)[1,]
  
  return(phi_i)
}

#### Setup function ####
prelim_setup = function(X, x_min, x_max, K_log_prior, basis_type){
  # X = list of vectors containing UNSCALED/UNSTANDARDIZED samples
  # x_max and x_min = optional scalars giving known boundaries for the densities
  # (e.g. if we know they're all exponential, we could set x_min = 0)
  # K_log_prior = Prior log probabilities for K. We must have K_log_prior[k] = p(k) for k = 1, 2, ... 
  # basis_type = One of "B-spline", "Bernstein", or "Histogram"
  
  # Scale data to (0, 1) and define density supports
  X_scaled = vector("list", length(X))
  # Allowing densities to extend slightly beyond sample range (based on handwavey
  # order statistic theory) if endpoints not pre-specified
  if(is.null(x_min)){
    x_min = sapply(1:length(X), 
                   function(n) min(X[[n]]) - (max(X[[n]]) - min(X[[n]])) / (length(X[[n]]) - 1))
  }
  else{
    x_min = rep(x_min, length(X))
  }

  if(is.null(x_max)){
    x_max = sapply(1:length(X), 
                   function(n) max(X[[n]]) + (max(X[[n]]) - min(X[[n]])) / (length(X[[n]]) - 1))
  }
  else{
    x_max = rep(x_max, length(X))
  }

  for(i in 1:length(X)){
    X_scaled[[i]] = (X[[i]] - x_min[i]) / (x_max[i] - x_min[i])
  }
  
  # Set up the basis evaluation function
  basis_fn = switch(basis_type, 'B-spline' = bspline_eval,
                    'Bernstein' = bernstein_eval, 'Histogram' = histogram_eval)
  if(is.null(basis_fn)){
    stop("basis_type must be 'B-spline', 'Bernstein', or 'Histogram'.")
  }
  
  K_max = length(K_log_prior)
  
  # For the spline basis type, we need normalizing constants for all K
  if(basis_type == 'B-spline'){
    inv_ints = vector('list', K_max)
    for(k in 1:3){ # Placeholders since we can't have k < 4 for cubic splines
      inv_ints[[k]] = 0
    }
    
    for(k in 4:K_max){
      temp_basis = create.bspline.basis(c(0, 1), norder = 4, nbasis = k)
      inv_ints[[k]] = 1/inprod(temp_basis)[,1]
    }
  }
  else{
    inv_ints = NULL
  }
  
  # Also need a list of basis values at a fine grid for each K
  basis_mat_list = vector('list', K_max)
  for(k in 1:K_max){
    basis_mat_list[[k]] = outer(seq(0, 1, 0.002), 1:k,
                                function(x, y) do.call(basis_fn,
                                                       c(list(X = x, labs = y, k = k),
                                                         list(inv_ints =
                                                                inv_ints[[k]])[!is.null(inv_ints)])))
  }
  
  return(list(X_scaled = X_scaled, x_min = x_min, x_max = x_max, inv_ints = inv_ints,
              basis_mat_list = basis_mat_list, basis_fn = basis_fn))
}

#### MCMC functions (original Petrone algo) ####
samp_step = function(X, Z0, K0, K_log_prior, M, basis_fn, inv_ints = NULL, x_min, x_max,
                     basis_mat_list){
  # Function to perform a single MCMC sampling step, combining all parameter-sampling functions.
  # x_min and x_max = vectors containing endpoints for X-ranges on unstandardized scale
  # basis_mat_list = Premade list of basis function evaluations at a fine grid for every k
  # TODO: Add log-likelihood for future PSIS-LOO CV
  
  # Sample density dimensionality
  K = K_sampler(X, Z0, K_log_prior, basis_fn, inv_ints)
  
  # Sample latent variables
  Z = mapply(X, Z0, FUN = latent_sampler, SIMPLIFY = FALSE,
             MoreArgs = list(M = M, K = K, basis_fn = basis_fn, inv_ints = inv_ints[[K]]))
  
  # Sample density coefficients
  phi = sapply(Z, FUN = density_sampler, M = M, K = K, inv_ints = inv_ints[[K]])
  
  # Evaluate densities on unstandardized scale
  f = t(t(basis_mat_list[[K]] %*% phi) / (x_max - x_min))
  
  # Return a single sample of everything
  return(list(phi = phi, Z = Z, K = K, f = f))
}

single_chain = function(X, Z0, K0, K_log_prior, M, basis_fn, inv_ints = NULL, x_min, x_max,
                        basis_mat_list, S = 5000, r){
  # A function to run a single MCMC chain, given PRE-SCALED data and all necessary prior stuff.
  # S = number of iterations
  # Z0 = inital values for latent variables
  # r = chain ID
  
  # Time it
  start  = proc.time()
  
  # Set up objects to hold samples
  K = numeric(S)
  phi = list()
  Z = list()
  f = array(0, dim = c(S, 501, length(X)))
  
  # First step with initial values for latents
  first_step = samp_step(X, Z0, K0, K_log_prior, M, basis_fn, inv_ints, x_min, x_max,
                         basis_mat_list)
  K[1] = first_step$K
  phi[[1]] = first_step$phi
  Z[[1]] = first_step$Z
  f[1,,] = first_step$f
  
  # And awaaaay we go
  for(s in 2:S){
    next_step = samp_step(X, Z[[s-1]], K[s-1], K_log_prior, M, basis_fn, inv_ints, x_min, x_max,
                          basis_mat_list)
    K[s] = next_step$K
    phi[[s]] = next_step$phi
    Z[[s]] = next_step$Z
    f[s,,] = next_step$f
    
    print(paste('Chain ', r, ', iteration ', s, sep = ''))
  }
  
  tot_time = proc.time() - start
  
  return(list(K = K, phi = phi, Z = Z, f = f, tot_time = tot_time, inits = Z0))
}

run_mcmc = function(X, K0, rate, K_log_prior, basis_type,
                    M, S = 5000, n_chains = 1, x_max = NULL, x_min = NULL){
  # Sets everything up and runs the MCMC top to bottom. This is the "main" function.
  # X = list containing vectors of unstandardized samples
  # x_max and x_min = optional scalars giving known boundaries for the densities
  # (e.g. if we know they're all exponential, we could set x_min = 0)
  # K_log_prior = Prior log probabilities for K. We must have K_log_prior[k] = p(k) for k = 1, 2, ... 
  # basis_type = One of "B-spline", "Bernstein", or "Histogram"
  # M = concentration parameter for Dirichlet Process
  # S = number of MCMC iterations
  # n_chains = Number of parallel chains to run.
  # TODO: allow for > 1 chains
  
  # Get the data all set up (including scaling)
  start = proc.time()
  print('Setting up data...')
  setup_stuff = prelim_setup(X, x_min, x_max, K_log_prior, basis_type)
  print(paste('Data setup complete in', format((proc.time() - start)[3], 4), 'sec'))
  
  # Extract everything from setup_stuff for convenience
  mapply(assign, names(setup_stuff), setup_stuff, MoreArgs = list(envir = as.environment(-1)))
  
  mcmc_chains = vector('list', n_chains)
  
  for(r in 1:n_chains){ # Will parallelize this later
    # Initialize latent variables Z by adding random noise with random variance to each X-vector
    Z0 = vector('list', length(X))
    for(i in 1:length(X)){
      Z0[[i]] = pmin(pmax(X_scaled[[i]] + rnorm(length(X[[i]]), 0, sd(X_scaled[[i]]) / rate),
                          0), 1)
    }
    
    # Run the chain
    mcmc_chains[[r]] = single_chain(X_scaled, Z0, K0, K_log_prior, M, basis_fn, inv_ints, x_min, x_max,
                                    basis_mat_list, S, r)
  }
  
  # Extract the stuff in setup_stuff in the output again for convenience
  return(c(mapply(assign, names(setup_stuff), setup_stuff, MoreArgs = list(envir = as.environment(-1))),
           mcmc_chains = mcmc_chains))
}

#### EM algorithm functions ####
# Functions to implement the mostly-frequentist "hybrid estimator" approach using the EM
# algorithm and IC weighting, as in Petrone/Wasserman (2002)

EM_step = function(X, K, phi_init, basis_vals, inv_ints = NULL){
  # Combined E- and M-step given data X (a list of samples), dimensionality K, and phi_init
  # (a list containing initial weights/coefficients for each entry in X).
  # basis_vals = a list of length length(X), each entry of which contains all basis function evaluations
  # at all points in the X-sample
  # inv_ints = Normalizing constants, only needed for B-spline basis type
  
  R = vector('list', length(X))
  for(i in 1:length(X)){
    # dens_vals = outer(X[[i]], 1:K,
    #                   function(x, z) do.call(basis_fn, c(list(X = x, labs = z, k = K),
    #                                                      list(inv_ints = inv_ints)[!is.null(inv_ints)])))
    resp_mat = t(t(basis_vals[[i]]) * phi_init[[i]])
    R[[i]] = resp_mat / rowSums(resp_mat)
  }
  
  phi = t(sapply(R, function(x) colSums(x) / nrow(x)))
  phi = split(phi, row(phi))
  
  return(list(R = R, phi = phi))
}

EM_alg = function(X, K, basis_fn, rel_tol, inv_ints = NULL, phi_init = NULL, max_iter = 2000){
  # Wrapper function to iteratre the E- and M-steps until convergence, as measured by relative tolerance rel_tol.
  # Need SCALED list of samples X. Can optionally supply an initial MATRIX of mixture coef's phi_init
  # (coercion to list is done within this function); otherwise the rows are all sampled from a uniform distribution
  # on the K-simplex
  
  if(is.null(phi_init)){
    phi_init = rdirichlet(length(X), rep(1, K))
  }
  
  phi_init = split(phi_init, row(phi_init))
  
  basis_vals = vector('list', length(X))
  for(i in 1:length(X)){
    basis_vals[[i]] =  outer(X[[i]], 1:K,
                             function(x, y) do.call(basis_fn,
                                                    c(list(X = x, labs = y, k = K),
                                                      list(inv_ints = inv_ints[[K]])[!is.null(inv_ints)])))
  }
  
  loglik = sum(log(unlist(mapply("%*%", basis_vals, phi_init))))
  
  step_one = EM_step(X, K, phi_init, basis_vals, inv_ints)
  phi = step_one$phi
  loglik = c(loglik, sum(log(unlist(mapply("%*%", basis_vals, phi)))))
  improvement = abs((loglik[2] - loglik[1]) / loglik[1])
  iter_count = 1
  
  while(improvement > rel_tol & iter_count < max_iter){
    
    step = EM_step(X, K, phi, basis_vals, inv_ints)
    phi = step$phi
    loglik = c(loglik, sum(log(unlist(mapply("%*%", basis_vals, phi)))))
    improvement = abs((loglik[length(loglik)] - loglik[length(loglik) - 1]) / loglik[length(loglik) - 1])
    
    iter_count = iter_count + 1
    # print(iter_count)
  }
  
  if(improvement > rel_tol & iter_count >= max_iter){
    print('Convergence not achieved after maximum no. of iterations')
  }
  
  return(list(phi = phi, iter_count = iter_count, loglik = loglik, final_resp = step$R))
}
