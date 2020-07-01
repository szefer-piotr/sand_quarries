
source("data/data_processing.R")

# library(brms)
# library(lme4)
# library(lmerTest)

# Estimating probability of obtaining vulnerable individual
raredf$fStage <- "Stage I"
raredf[raredf$stage == 2, ]$fStage <- "Stage II"
raredf[raredf$stage == 3, ]$fStage <- "Stage III"
raredf$vuln <- raredf$vuln == 1

# Bayesian model
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
as.data.frame(it1$contrasts)
# emmip(glm1, lifeh~fStage)

rarespdf <- as.data.frame(it1$emmeans)

invlog <- function(logitp){exp(logitp)/(1 + exp(logitp))}
preds <- predict(glm1, robust=T) # 95% credible intervals
invlog(preds)

newdat <- expand.grid(fStage = c("Stage I",
                                 "Stage II",
                                 "Stage III"),
                      lifeh = c("herbivore",
                                "predator",
                                "kleptoparasite"))
newdata.pred = predict(glm1, type = "response",
                       newdata = newdat,se.fit = T,
                       interval = "confidence")
newdata.pred
plotdat <- data.frame(Stage = newdat$fStage,
                      Group = newdat$lifeh,
                      Proportion = newdata.pred$fit,
                      SE = newdata.pred$se.fit,
                      CLR = rep(colvec[1:3], each=3))

p<- ggplot(plotdat, aes(x=Stage, y=Proportion, 
                        group = Group,
                        color = alpha(CLR, 0.5))) + 
  geom_point()+
  geom_errorbar(aes(ymin=Proportion-SE, 
                    ymax=Proportion+SE), width=.2,
                position=position_dodge(0.05))+
  facet_wrap(~Group)
p+ theme(legend.position = "none")

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
