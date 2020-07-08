# Rare species analysis

source("data/data_processing.R")

# library(brms)
# library(lme4)
# library(lmerTest)

# Estimating probability of obtaining vulnerable individual
raredf$fStage <- "Early"
raredf[raredf$stage == 2, ]$fStage <- "Mid"
raredf[raredf$stage == 3, ]$fStage <- "Late"
raredf$vuln <- raredf$vuln == 1

raredf$lifeh <- factor(raredf$lifeh, levels = c('herbivore',
                                                 'kleptoparasite',
                                                 'predator'),
                        labels = c('Herbivore',
                                   'Kleptoparasite',
                                   'Predator'))

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

raredf$lifeStage <- paste(raredf$lifeh, raredf$fStage, sep="_")

# Same model but formatted speifically for pairwise comparisons
glm2 <- glm(vuln ~ lifeStage,
            family = binomial(link = "logit"),
            data =raredf)
it2 <- emmeans(glm2, "lifeStage")
it2let <- cld(it2, Letter="abcdefghijklm")
ltrs <- c(it2let$.group[3],
          it2let$.group[6],
          it2let$.group[2],
          it2let$.group[4],
          it2let$.group[5],
          it2let$.group[7],
          it2let$.group[1],
          it2let$.group[8],
          it2let$.group[9])
ltrs <- gsub(" ", "", ltrs)

# Predict the data and prepare plots
it1 <- emmeans(glm1, pairwise~ fStage|lifeh )
as.data.frame(it1$contrasts)

rarespdf <- as.data.frame(it1$emmeans)

newdat <- expand.grid(fStage = c("Early",
                                 "Mid",
                                 "Late"),
                      lifeh = c("Herbivore",
                                "Kleptoparasite",
                                "Predator"
                                ))
newdata.pred = predict(glm1, type = "response",
                       newdata = newdat,se.fit = T,
                       interval = "confidence")
newdata.pred
table(raredf$lifeh, raredf$fStage)
plotdat <- data.frame(Stage = newdat$fStage,
                      Group = newdat$lifeh,
                      Proportion = newdata.pred$fit,
                      N = as.vector(t(table(raredf$lifeh, raredf$fStage))),
                      SE = newdata.pred$se.fit,
                      CLR = rep(colvec[3:5], each=3),
                      sig = ltrs)
#Lower limit = p - (z) (Estimated σp) - 0.5/N
plotdat$ucl <- plotdat$Proportion + qnorm(0.975)*plotdat$SE - 0.5/plotdat$N
plotdat$lcl <- plotdat$Proportion - qnorm(0.975)*plotdat$SE - 0.5/plotdat$N

p<- ggplot(plotdat, aes(x=Stage, y=Proportion, 
                        group = Group,
                        color = Group,
                        label = sig)) + 
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=.2,lwd=1,
                position=position_dodge(0.05))+
  stat_summary(fun=mean, geom="text", 
               col = rgb(10,10,10,180,maxColorValue = 255),
               hjust = 1.7)+
  facet_wrap(~Group)+
  ylim(0, 0.21)

# jpeg("fig_4_rare_species.jpg", width = 800, height = 300)
p+ theme_bw()+ theme(legend.position = "none")+
  scale_color_manual(values=c(alpha(colvec[2],1), 
                              alpha(colvec[1],1), 
                              alpha(colvec[3],1)))
# dev.off()

