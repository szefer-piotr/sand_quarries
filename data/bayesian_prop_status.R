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
                        labels = c('Herbivores',
                                   'Kleptoparasites',
                                   'Predators'))

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
glmNULL <- (glm(vuln ~ 1, 
               family = binomial(link = "logit"),
               data =raredf)) 
glm0 <- (glm(vuln ~ fStage, 
             family = binomial(link = "logit"),
             data =raredf))
summary(glm0)
anova(glmNULL, glm0, test = "Chi")

glm1 <- (glm(vuln ~ lifeh*fStage, 
    family = binomial(link = "logit"),
    data =raredf))

raredf$lifeStage <- paste(raredf$lifeh, raredf$fStage, sep="_")

# Same model but formatted specifically for pairwise comparisons
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
                      lifeh = c("Herbivores",
                                "Kleptoparasites",
                                "Predators"
                                ))
newdata.pred = predict(glm1, type = "response",
                       newdata = newdat,se.fit = T,
                       interval = "confidence")

newdata.pred0 = predict(glm0, type = "response",
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

library(dplyr)
plotdat$Group <- recode(plotdat$Group, Kleptoparasites = "Parasitic species")

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

pdf("revision_1/figures/Fig6.pdf", width = 10, height = 4, onefile = FALSE)
p+ theme_bw()+ theme(legend.position = "none")+
  scale_color_manual(values=c(alpha(colvec[1],1), 
                              alpha(colvec[2],1), 
                              alpha(colvec[3],1)))
dev.off()

# General proportion for stages ----
it0 <- emmeans(glm0, "fStage")
it0let <- cld(it0, Letter="abcdefghijklm")
newdat0 <- expand.grid(fStage = c("Early",
                                 "Mid",
                                 "Late"))
newdata.pred0 = predict(glm0, type = "response",
                        newdata = newdat0,se.fit = T,
                        interval = "confidence")

plotdat0 <- data.frame(Stage = newdat0$fStage,
                      Proportion = newdata.pred0$fit,
                      N = as.vector(t(table(raredf$fStage))),
                      SE = newdata.pred0$se.fit,
                      CLR = colvec[3:5],
                      sig = it0let$.group)
#Lower limit = p - (z) (Estimated σp) - 0.5/N
plotdat0$ucl <- plotdat0$Proportion + qnorm(0.975)*plotdat0$SE - 0.5/plotdat0$N
plotdat0$lcl <- plotdat0$Proportion - qnorm(0.975)*plotdat0$SE - 0.5/plotdat0$N

p<- ggplot(plotdat0, aes(x=Stage, y=Proportion,
                        label = sig)) + 
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lcl, 
                    ymax=ucl), width=.2,lwd=1,
                position=position_dodge(0.05))+
  stat_summary(fun=mean, geom="text", 
               col = rgb(10,10,10,180,maxColorValue = 255),
               hjust = 1.7)+
  ylim(0, 0.21)

# jpeg("fig_4_rare_species.jpg", width = 800, height = 300)
p+ theme_bw()+ theme(legend.position = "none")

