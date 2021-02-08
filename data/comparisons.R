# Comparisons

# Source data processing script ----
source("data/data_processing.R")
# library(nlme)
library(MASS)
library(emmeans)
library(multcomp)
# 1. Abundance ----
# >>> Analysis ----
# This should probably benegative binomial for abundance

# General model
abund_gen <- glm.nb(abun~1+succession, data = desasc)
anova(abund_gen, test = "Chisq")
anova(abund_gen, test = "F")
inter.test_gen <- emmeans(abund_gen, "succession")
phabu_gen <- cld(inter.test_gen, Letter="abcdefghijklm")



abund_int1 <- glm(abun~1+GS, 
             family = "poisson", data = desasc)

abund_int2 <- glm.nb(abun~1+GS, data = desasc)

abund_int3 <- glm(abun~1+GS, data = desasc)

# Select the best model
AIC(abund_int1, abund_int2,abund_int3)
# negative binomial model fits best

# >>> Post-hoc ----
inter.test1 <- emmeans(abund_int2, "GS")
phabu <- cld(inter.test1, Letter="abcdefghijklm")
contrast(inter.test1, type = "response")

# 2. Diversity ----
# * SIMPSON ----
# >>> Analysis ----
simpson <- glm(simp~1+GS, data = desasc)

# >>> Post-hoc ----
# inter.test1 <- emmeans(simpson, "GS")
# gps <- cld(inter.test1, Letter="abcdefghijklm")

# * SHANNON ----
# >>> Analysis ----
sw <- glm(sw~1+GS, data = desasc)
sw_gen <- glm(sw~succession, data = desasc)
anova(sw_gen)

# >>> Post-hoc ----
inter.test2 <- emmeans(sw, "GS")
phsw <- cld(inter.test2, Letter="abcdefghijklm")

inter.test2_gen <- emmeans(sw_gen, "succession")
phsw_gen <- cld(inter.test2_gen, Letter="abcdefghijklm")

# 3. Richness ----
bg <- glm(rich~1+GS, data = desasc, family = "poisson")
bg_gen <- glm(rich~1+succession, data = desasc, family = "poisson")
summary(bg_gen)
anova(bg_gen, test="Chisq")

# >>> Post-hoc ----
inter.test3 <- emmeans(bg, "GS")
phrich <- cld(inter.test3, Letter="abcdefghijklm")

inter.test3_gen <- emmeans(bg_gen, "succession")
phrich_gen <- cld(inter.test3_gen, Letter="abcdefghijklm")


# Facet plot ----
facet_dat <-data.frame( 
  val = c(log(desasc$abun+1),
          desasc$rich,
          desasc$sw),
  group = rep(desasc$group,3),
  stage = rep(desasc$succession, 3),
  type = rep(c("No. of individuals",
               "Richness",
               "Diversity"), each=96)
)

facet_dat$fsucc <- 0
facet_dat[facet_dat$stage == 1,]$fsucc <- "Stage I"
facet_dat[facet_dat$stage == 2,]$fsucc <- "Stage II"
facet_dat[facet_dat$stage == 3,]$fsucc <- "Stage III"
facet_dat$int <- as.character(desasc$GS)
facet_dat$type <- as.character(facet_dat$type)
facet_dat$post_hoc <- "x"

# Generate labels
# For "No. of individuals"
abugroups <- data.frame(int = phabu$GS, 
                        gr = gsub(" ", "", phabu$.group))
abugroups$gr <- as.character(abugroups$gr)
for (itr in abugroups$int){
  indices <- (facet_dat$type == "No. of individuals" & facet_dat$int == itr)
  facet_dat[indices,]$post_hoc <- abugroups[abugroups$int == itr,]$gr
}

# For "Diversity"
divgroups <- data.frame(int = phsw$GS, 
                        gr = gsub(" ", "", phsw$.group))
divgroups$gr <- as.character(divgroups$gr)
for (itr in divgroups$int){
  indices <- (facet_dat$type == "Diversity" & facet_dat$int == itr)
  facet_dat[indices,]$post_hoc <- divgroups[divgroups$int == itr,]$gr
}

# For "Richness"
richgroups <- data.frame(int = phrich$GS, 
                        gr = gsub(" ", "", phrich$.group))
richgroups$gr <- as.character(richgroups$gr)
for (itr in richgroups$int){
  indices <- (facet_dat$type == "Richness" & facet_dat$int == itr)
  facet_dat[indices,]$post_hoc <- richgroups[richgroups$int == itr,]$gr
}

fpl <- ggplot(facet_dat, aes(x = fsucc, y = val, 
                          col = group,
                          group = group,
                          label = post_hoc))
fullflp <- fpl + geom_jitter(width=0.1, alpha=0.3, cex=2) +
  facet_wrap(~type, scales = "free") +
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange", lwd=0.8) +
  stat_summary(fun=mean, geom="point",cex = 2) +
  stat_summary(fun=mean, geom="line",lwd=1, lty=2)+
  stat_summary(fun=mean, geom="text", 
               col = rgb(10,10,10,180,maxColorValue = 255),
               hjust = 1.4,
               vjust = -1.5) +
  theme() +
  xlab("") + ylab("") +
  theme_bw() +
  scale_color_manual(values=c(colvec[1], 
                              colvec[2], 
                              colvec[3]))
          

# pdf("fig1_glm.pdf", width = 12, height=4)
# fullflp
# dev.off()

# Table 1 - results from mixed effect models
summary(abund_int2)
summary(sw)
summary(bg)

phrich
phsw
phabu

phrich_gen
phsw_gen
phabu_gen

# Predict abundance
abund_int2
py1 <- data.frame(GS = unique(desasc$GS))
p1  <- predict(abund_int2, newdata=py1, se.fit=TRUE, type='response')
py1 <- data.frame(py1, p1)

py2 <- data.frame(GS = unique(desasc$GS))
p2  <- predict(sw, newdata=py2, se.fit=TRUE, type='response')
py2 <- data.frame(py2, p2)

py3 <- data.frame(GS = unique(desasc$GS))
p3  <- predict(bg, newdata=py3, se.fit=TRUE, type='response')
py3 <- data.frame(py3, p3)
