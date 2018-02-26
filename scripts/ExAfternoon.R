# extending the models
# https://microbiome.github.io/microbiome/rstanarm.html

# Exercise session 2: extend and vary the linear model
# Try out alternative formulations of the linear model 
# (stan_lm; stan_glm; stan_lmer)

my_stan_lm <- stan_lm(weightgain ~ source * type,
                           data = weightgain, 
                           prior = R2(location = 0.3, 
                                      what = "mean"), 
                           adapt_delta = 0.999)

my_stan_lm

my_stan_glm <- stan_glm(weightgain ~ source * type,
                      data = weightgain, 
                      prior = cauchy(), 
                      adapt_delta = 0.999)

my_stan_lm

my_stan_lmer <- stan_lmer(weightgain ~ 1 + (1|source:type),
                   data = weightgain, prior_intercept = cauchy(),
                   prior_covariance = decov(shape = 1, scale = 3),
                   adapt_delta = 0.999, chains = 4, cores = 4)

my_stan_lmer
launch_shinystan(object = my_stan_lmer)

my_stan_lmer1 <- stan_lmer(weightgain ~ 1 + (1|source:type),
                          data = weightgain, prior_intercept = cauchy(),
                          prior_covariance = decov(shape = 2, scale = 2),
                          adapt_delta = 0.999, chains = 4, cores = 4)

my_stan_lmer1

# compare the two lmer
summary(as.data.frame(my_stan_lmer))
summary(as.data.frame(my_stan_lmer1))
launch_shinystan(object = my_stan_lmer1)
# PPcheck
# posterior predictive check

plot(my_stan_lmer)
plot(my_stan_lmer1)

# loo
# leave one out cross validation
(l_my_stan_lmer <- loo(my_stan_lmer))
(l1_my_stan_lmer <- loo(my_stan_lmer1))

compare_models(l_my_stan_lmer, l1_my_stan_lmer)

# WAIC 
?waic
(waic1 <- waic(my_stan_lmer))
(waic2 <- waic(my_stan_lmer1))
print(compare(waic1, waic2), digits = 3)

par(mfrow = c(2,1))
plot(l_my_stan_lmer, label_points = TRUE)
plot(l1_my_stan_lmer, label_points = TRUE)

# the number of parameters in a model gives a mesuarement of 
# complexity
# if some parameters are correlated then it means, there are
# less number of efficient parameters

# LPPD = log posterior density LPPD
# P = parameters
# waic = -2 (LPPD - P)

# MCMC
# Marcov Models
# samples are correlated between consecutive points.
# that's why we take every 10th or every 100th 
# 10 thousand 10^4 points is common to select
# we burn the first 1000 
# if we have 3 chains of values, then you realize the
# tend to converge
# tipycally at least 3 chains
# Rhat measures how similar the chains are 
# (if it is similar to 1 they are convereged )

# SECOND EXERCISE ---------------------------------------------------------
# Increase the number of predictors in the linear model
# (instead of gender, use gender, nationality, and BMI,
#   for instance) and investigate how this affects the inference

data(dietswap)

# Some necessary data preparations
sample_data(dietswap)$time <- sample_data(dietswap)$timepoint

# Pick jist the baseline data points (first time point 1)
x <- baseline(dietswap)

# Pick metadata from the phyloseq object
df <- meta(x) # return sample metadata as a data frame

rm(x)
# Add microbiome diversity to the metadata
df$diversity <- diversities(x, "shannon")$shannon

head(df)
summary(df)

# Increase the number of predictors in the linear model
# (instead of gender, use gender, nationality, and BMI,
#   for instance) and investigate how this affects the inference
my_gnb_stan_glm <- stan_glm(diversity ~  nationality * sex  ,
                           data = df, 
                           prior = normal(),
                           adapt_delta = 0.999, chains = 4, cores = 4)

my_gnb_stan_glm
launch_shinystan(object = my_gnb_stan_glm)
summary(as.data.frame(my_gnb_stan_glm))
 
# -----------------------------------------------------------------------
# Read about priors in linear models. Experiment with 
# different priors and see how this affects the results 
# (gaussian, cauchy, uniform, different parameter values). 
# You can also read more on priors in rstanarm.

my_dns_stan_glm <- stan_glm(diversity ~  nationality * sex  ,
                            data = df, 
                            prior = student_t(),
                            adapt_delta = 0.999, chains = 4, cores = 4)

my_dns_stan_glm
summary(as.data.frame(my_dns_stan_glm))


# Example
data("clouds", package = "HSAUR3")
head(clouds)
ols <- lm(rainfall ~ seeding * (sne + cloudcover + prewetness + echomotion) +
            time, data = clouds)
round(coef(ols), 3)
post_clouds <- 
  stan_lm(rainfall ~ seeding * (sne + cloudcover + prewetness + echomotion) + 
            time, 
          data = clouds,
          prior = R2(location = 0.2), 
          chains = 4, cores = 4)
post_clouds


head(df)
str(df)
summary(df)
table(df$sex, df$nationality)
table(df$sex, df$nationality, df$bmi_group)
# Increase the number of predictors in the linear model
# (instead of gender, use gender, nationality, and BMI,
#   for instance) and investigate how this affects the inference
library(ggplot2)
ggplot(data = df, aes( x = bmi_group, fill = nationality)) +
  geom_bar() + 
  facet_grid(~ sex) +
  theme_bw()

#
summary(df$diversity)
ggplot(data = df, aes( x = diversity)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sex)


my_stan_lm_small <- stan_lm(diversity ~ sex,
                      data = df, 
                      prior = R2(location = 0.5, what = "mean"), 
                      adapt_delta = 0.999)

my_stan_lm_small
summary(as.data.frame(my_stan_lm_small))
plot(my_stan_lm_small)

head(df)
my_stan_lm_small1 <- stan_lm(diversity ~ bmi_group * ( nationality + bmi_group ),
                            data = df, 
                            prior = R2(location = 0.5, what = "mean"), 
                            adapt_delta = 0.999)

my_stan_lm_small1
summary(as.data.frame(my_stan_lm_small1))
plot(my_stan_lm_small1)

# We can estimate or visualize the average treatment effect (ATE) 
# using rstanarm’s posterior_predict function.
clouds
clouds_cf <- clouds
(clouds_cf$seeding[] <- "yes")
y1_rep <- posterior_predict(post_clouds, newdata = clouds_cf)
clouds_cf$seeding[] <- "no"
y0_rep <- posterior_predict(post_clouds, newdata = clouds_cf)
qplot(x = c(y1_rep - y0_rep), geom = "histogram", xlab = "Estimated ATE")
# test data clouds_df in this example is a copy of the clouds
# tries to see the relationship of the new point in the model

post_clouds_glm <- stan_glm(rainfall ~ seeding * (sne + cloudcover + prewetness + 
                                           echomotion) + time,
                   data = clouds, family = gaussian(), 
                   prior = cauchy(), prior_intercept = cauchy(),
                   chains = 4, cores = 4)
post_clouds_glm

(loo_post <- loo(post_clouds))
# .Warning message:
#   Found 2 observation(s) with a pareto_k > 0.7. We recommend calling 'loo' 
# again with argument 'k_threshold = 0.7' in order to calculate the ELPD without
# the assumption that these observations are negligible. This will refit the model 
# 2 times to compute the ELPDs for the problematic observations directly.
(loo_post <- loo(post_clouds, pareto_k = 0.8, k_threshold = 0.7))
(loo_post_glm <-  loo(post_clouds_glm))

compare_models(loo_post,loo_post_glm )


