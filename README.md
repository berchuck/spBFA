# spBFA
Implements a spatial Bayesian factor analysis model with inference in a Bayesian non-parametric setting using Markov chain Monte Carlo (MCMC). Spatial correlation is introduced in the columns of the factor loadings matrix using a Bayesian non-parametric prior, the probit stick-breaking process. Areal spatial data is modeled using a conditional autoregressive (CAR) prior and point-referenced spatial data is treated using a Gaussian process. The response variable can be modeled as Gaussian, probit, Tobit, or Binomial (using Polya-Gamma augmentation). Temporal correlation is introduced for the latent factors through a hierarchical structure and can be specified as exponential or first-order autoregressive. See `vignette('spBFA-example')` for usage. A manuscript by Berchuck et al. 2019 that details the theory corresponding to the functions in the `spBFA` package is forthcoming.