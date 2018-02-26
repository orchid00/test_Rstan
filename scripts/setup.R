
# To implement Bayesian models, like linear models we use rstan packages
# Statistical Thinking, A bayesian 
# Statistical modeling to deal with uncertainty

# Statistical modeling: Try to estimate the likelihood of the data given 
# the model parameters
# in bayesian modeling we add a prior, and choosing th prior is somehow
# arbitrary
# In bayes we try to identify the posterior distribution
# based on the prior distribution and the data
# Likelihood is one number
# prior (distribution) is a range of parameters (usually normally distributed)

# In Bayesian statistics, the parameters (e.i. Alpha and betha)
# have their own distributions, unlike lm
# parameter pooling can help regularise the model
# The use of Cauchy distribution is common iyesian statistn Baics, but we can use
# uniform priors and gamma distribution
# if the priors are loose and flexible it doesno't mandate the data 

# Stan is a programming langage for probabilistic analysis
# open source mc-stan.org
# first choose parameters
# then fit the model
# stan is cross-programming (R/python/julia)

# library load
library(rstan)
library(rstanarm)
# there is no package called ‘DT’ before loading rstanarm
# install.packages("DT")
library(rstantools)
library(microbiome)

# https://microbiome.github.io/
# https://microbiome.github.io/microbiome/rstanarm.html

# Exercise
# In that example, the use of Bayesian prior helps to avoid overfitting and 
# finds a more acurate solution when the sample size is small compared to the
# model complexity.

# Install the example data set used in rstanarm tutorials
#install.packages("HSAUR3")
library(HSAUR3)

