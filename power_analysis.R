install.packages("pwrss")
library(pwrss)

# do bees choose the salty nectar significantly more than chance? (repeated for each species)
# => is p > .5 (p = probability bees choose salty nectar)
# => this is a one-sample binomial test
# key decisions: 
  # this test will be one-sided because the literature shows no evidence of bees avoiding salty nectar, so I do not think it is plausible for bees to choose the control MORE than the treatment
  # i'm leaving power at 0.80 as a default number. while salting out may be subtle chemically, the effect on behavior should not be subtle, so the larger sample size a higher power value would require is not worth it
  # I'm using 0.65 as the expected probability, but this number is not founded in anything other than 1) a desire to see more subtle (closer to average) choice differences and 2) setting a number that allows me to collect a manageable amount of bees
expected_prob <- 0.60
desired_power <- 0.8 
alpha <- 0.05
binom_test_power <- power.exact.oneprop(prob = expected_prob,
                    null.prob = 0.5,
                    power = desired_power,
                    alpha = alpha,
                    alternative = "one.sided")
num_bees_per_species <- binom_test_power$size
actual_power <- binom_test_power$power
# ========================================

# is species a predictor of the probability of choosing salty nectar?
# there's no clean function to do the power analysis here like there is for binom test above
# so, i am going to run a bunch of simulated trials and count what proportion got statistically significant results 
# in order to simulate this, i must guess at the proportion of each species that will choose salty nectar
est_honey_p = 0.7
est_bumble_p = 0.6
est_sweat_p = 0.5
# then run a bunch of simulations and accumulate the significant results
num_sims <- 1000
significant_results <- logical(num_sims)

for (i in 1:num_sims) {
  # simulate data
  species <- rep(c("honey", "bumble", "sweat"), each = num_bees_per_species)
  probs <- rep(c(est_honey_p, est_bumble_p, est_sweat_p), each = num_bees_per_species)
  choice <- rbinom(length(probs), size = 1, prob = probs) 
  
  # fit the model
  model <- glm(choice ~ species, family = binomial)
  p_val <- anova(model, test = "Chisq")[2, "Pr(>Chi)"] # p-value for species predictor
  
  significant_results[i] <- (p_val < alpha)
}

# estimate power using the simulation
power_estimate <- mean(significant_results)
power_estimate

