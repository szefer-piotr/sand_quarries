# Comparisons

# Source data processing script ----
source("data/data_processing.R")
library(nlme)
library(MASS)
library(emmeans)
library(multcomp)
# 1. Abundance ----
# >>> Analysis ----
# This should probably benegative binomial for abundance

abund_int1 <- glm(abun~1+GS, 
             family = "poisson", data = desasc)

abund_int2 <- glm.nb(abun~1+GS, data = desasc)

abund_int3 <- glm(abun~1+GS, data = desasc)

# Select the best model
AIC(abund_int1, abund_int2,abund_int3)
# negative binomial model fits best

# See the summary
summary(abund_int2)

# >>> Post-hoc ----
inter.test1 <- emmeans(abund_int2, "GS")
phabu <- cld(inter.test1, Letter="abcdefghijklm")

# >>> Plot ----
# generate_letters
# ap <- ggplot(desasc, aes(x = succession, y = abun, 
#                          col = group,
#                          group = group))
# ap <- ap + geom_jitter(width=0.1, alpha=0.3) + 
#   stat_summary(fun.data=mean_cl_boot, 
#              geom="pointrange", lwd=1) +
#   stat_summary(fun.y=mean, geom="point",cex = 4) +
#   stat_summary(fun.y=mean, geom="line",cex = 0.5, lty=2)+
#   theme_bw()
# ap

# 2. Diversity ----
# * SIMPSON ----
# >>> Analysis ----
simpson <- glm(simp~1+GS, data = desasc)
summary(simpson)
# >>> Plot ----
# si <- ggplot(desasc, aes(x = succession, y = simp, 
#                          col = group,
#                          group = group))
# si <- si + geom_jitter(width=0.1, alpha=0.3) + 
#   stat_summary(fun.data=mean_cl_boot, 
#                geom="pointrange", lwd=0.8) +
#   stat_summary(fun.y=mean, geom="point",cex = 4) +
#   stat_summary(fun.y=mean, geom="line",cex = 1)
# si

# >>> Post-hoc ----
# inter.test1 <- emmeans(simpson, "GS")
# gps <- cld(inter.test1, Letter="abcdefghijklm")

# * SHANNON ----
# >>> Analysis ----
sw <- glm(sw~1+GS, data = desasc)
summary(sw)
# # >>> Plot ----
# sp <- ggplot(desasc, aes(x = succession, y = sw, 
#                          col = group,
#                          group = group))
# sp <- sp + geom_jitter(width=0.1, alpha=0.3) + 
#   stat_summary(fun.data=mean_cl_boot, 
#                geom="pointrange", lwd=0.8) +
#   stat_summary(fun.y=mean, geom="point",cex = 4) +
#   stat_summary(fun.y=mean, geom="line",cex = 1)
# 
# sp

# >>> Post-hoc ----
inter.test2 <- emmeans(sw, "GS")
phsw <- cld(inter.test2, Letter="abcdefghijklm")


# 3. Richness ----
bg <- glm(rich~1+GS, data = desasc, family = "poisson")

# # >>> Plot ----
# bgp <- ggplot(desasc, aes(x = succession, y = rich, 
#                          col = group,
#                          group = group))
# bgp <- bgp + geom_jitter(width=0.1, alpha=0.3) + 
#   stat_summary(fun.data=mean_cl_boot, 
#                geom="pointrange", lwd=0.8) +
#   stat_summary(fun.y=mean, geom="point",cex = 4) +
#   stat_summary(fun.y=mean, geom="line",cex = 1)
# 
# bgp


# >>> Post-hoc ----
inter.test3 <- emmeans(bg, "GS")
phrich <- cld(inter.test3, Letter="abcdefghijklm")

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
  stat_summary(fun.y=mean, geom="point",cex = 2) +
  stat_summary(fun.y=mean, geom="line",lwd=1, lty=2)+
  stat_summary(fun.y=mean, geom="text", 
               col = rgb(10,10,10,180,maxColorValue = 255),
               hjust = 1.2,
               vjust = -1.5) +
  theme() +
  xlab("") + ylab("") +
  theme_bw() +
  scale_color_manual(values=c(colvec[1], 
                              colvec[2], 
                              colvec[3]))
          
