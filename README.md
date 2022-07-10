# density_UQ_review_paper
This repo contains materials for [A review of uncertainty quantification for density estimation](https://projecteuclid.org/journals/statistics-surveys/volume-15/issue-none/A-review-of-uncertainty-quantification-for-density-estimation/10.1214/21-SS130.full) (McDonald and Campbell, 2021).

**DOCUMENTATION**

Contains the full text of the review paper, as well as the supplement and the figure.

**CODE**

`review_paper_code.R`: the code used to generate the figure and supplementary material.

`unknown_dim_basis_dens_functions.R`: a custom implementation of the sampler described in [Bayesian Density Estimation Using Bernstein Polynomials](https://www.jstor.org/stable/3315494) (Petrone 1999), for estimating densities as mixtures of Bernstein polynomials (see also [Random Bernstein Polynomials](https://www.jstor.org/stable/4616563), Petrone 1999).
Our version is extended to allow for other basis types (including B-splines and step functions), as well as multiple densities (all sharing the same dimensionality).
