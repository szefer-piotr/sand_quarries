
source("data/data_processing.R")

library(brms)
# library(lme4)
# library(lmerTest)

# Estimating probability of obtaining vulnerable individual
# vulsp0 <- glm(vuln~1, family="bernoulli", 
#       data = rdf)


raredf$fStage <- "Stage I"
raredf[raredf$stage == 2, ]$fStage <- "Stage II"
raredf[raredf$stage == 3, ]$fStage <- "Stage III"
raredf$vuln <- raredf$vuln == 1


# vulsp0c <- brm(vuln ~ lifeh*fStage, 
#                data = raredf, 
#                family = 'bernoulli', 
#                prior = set_prior("uniform(2,4)", lb = 2, ub = 4),
#                iter = 1000, 
#                chains = 4)

# Ordinary glm
library(emmeans)
library(multcomp)
glm1 <- (glm(vuln ~ lifeh*fStage, 
    family = binomial(link = "logit"),
    data =raredf))

# Predict the data and prepare plots

it1 <- emmeans(glm1, pairwise~ fStage|lifeh )
emmip(glm1, lifeh~fStage)

plot(it1)
# 
# preds <- predict(vulsp0c, robust=F) # 95% credible intervals
# hist(preds[,1])
# conditional_effects(vulsp0c)

# https://mjskay.github.io/tidybayes/articles/tidy-brms.html

# pp = brms::pp_check(m)
# pp
# plot(vulsp0b)
# conditional_effects(vulsp1)
# 
# library(dplyr)
# d <- tibble(z=rbinom(100, 1, 0.6))
# m <- brm(z ~ 1, d, family=bernoulli())
# summary(m)
# predict(m, robust=F)
# 
