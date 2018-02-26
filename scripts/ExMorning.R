# Load example data
library(microbiome)
#data(dietswap)

# Some necessary data preparations
#sample_data(dietswap)$time <- sample_data(dietswap)$timepoint

# Pick jist the baseline data points (first time point 1)
#x <- baseline(dietswap)

# Pick metadata from the phyloseq object
#df <- meta(x) # return sample metadata as a data frame

# Add microbiome diversity to the metadata
#df$diversity <- diversities(x, "shannon")$shannon

#head(df)
#summary(df)

# Exercise session 1: introduction to rstanarm -------------------------------
# 1.1 Check the tutorial, and try out rstanarm examples with ANOVA in rstanarm.

# Using the Weightgain data
data("weightgain", package = "HSAUR3")
head(weightgain)
summary(weightgain)

# try anova with weightgain
standard_aov <- coef(aov(weightgain ~ source * type, data = weightgain))
standard_aov

# try anova with rstanarm
library(rstanarm)

# adjusted location from 0.5 to 0.2
# adapt_delta = 0.999 to decrease the stepsize and largely prevent 
# divergent transitions. 
post1 <- stan_aov(weightgain ~ source * type, data = weightgain, 
                  prior = R2(location = 0.2), adapt_delta = 0.999)
post1
summary(post1)
plot(post1)
# convert as data frame
df_post1 <- as.data.frame(post1)
summary(df_post1)

# Interpretation
plot(post1)
# select sourceCereal parameter
hist(df_post1$sourceCereal, breaks = 10)

d_post_Cereal <- density(df_post1$sourceCereal) # returns the density data 
plot(d_post_Cereal)

# 1.2 Explore the use of shinystan R package. This provides tools to diagnose 
# stan models visually and numerically, including model parameters and convergence
# diagnostics for MCMC simulations.
library(shinystan)
launch_shinystan_demo()
launch_shinystan(object = post1)

# 
# 1.3 Compare the results between rstanarm and standard alternatives 
# (lm vs. stan_lm or aov vs. stan_aov) etc.

standard_aov
summary(df_post1)


# Test lm
fit <- lm(weightgain ~ source * type,
          data = weightgain)
fit
standard_lm <- coef(fit) # point estimates

par(mfrow = c(2,2)) # Change the panel layout to 2 x 2
plot(fit)
par(mfrow = c(1,1)) # Change back to 1 x 1
summary(fit)

# lm with stan
fit_stan <- stan_lm(weightgain ~ source * type,
                    data = weightgain, 
                    prior = R2(location = 0.2, what = "mean"), 
                    adapt_delta = 0.999)

fit_stan
summary(fit_stan)
df_fitstan <- as.data.frame(fit_stan)
head(df_fitstan)

summary(df_fitstan)
standard_lm

hist(df_fitstan$sourceCereal, breaks = 15)

# try lm with a different parameter
fit_stan_round2 <- stan_lm(weightgain ~ source * type,
                           data = weightgain, 
                           prior = R2(location = 0.3, 
                                      what = "mean"), 
                           adapt_delta = 0.999)

fit_stan_round2
pairs(fit_stan_round2)

df_fitstan_round2 <- as.data.frame(fit_stan_round2)
head(df_fitstan_round2)
summary(df_fitstan_round2)
summary(df_fitstan)

# compare models use the loo package
library(loo)
looR2_0_2 <- loo(fit_stan, cores = 2)
looR2_0_3 <- loo(fit_stan_round2, cores = 2)
compare_models(looR2_0_2, looR2_0_3 )

kf1_R2_0_2 <- kfold(fit_stan, K = 10)
kf2_R2_0_3 <- kfold(fit_stan_round2, K = 10)
compare_models(kf1_R2_0_2, kf2_R2_0_3 )

# when little data the bayesian model performs better, 
# because of the prior

# R2 is the square of the correlation 
# (location and scale) instead as the prior
# is the betha distribution between 0 and 1
# location the mode or mean 
# what = 



# intercept: 100 to 90.06
# souceCereal: -14.1 to -12.027
# typeLow: -20.8 to -17.64
# sourceCereal:typeLow: 18.8 to 16.10

# 1.4 If the time allows, you can also make rstanarm examples with linear models 
# in rstanarm.

# Try with df
# head(df)
# summary(df)
# coef(aov(diversity ~ bmi_group * nationality, data = df ))
# post1 <- stan_aov(diversity ~ bmi_group * nationality, data = df, 
# prior = R2(location = 0.3), adapt_delta = 0.999)
# post1

launch_shinystan(object = post1)
