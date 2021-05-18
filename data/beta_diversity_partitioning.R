# Beta-diveristy partitioning ----

# Load data ----
source("data/data_processing.R")

# packages
library(betapart)
require(vegan)
library(ggplot2)
library(codyn)
library(emmeans)
library(multcomp)
library(brms)
library(betareg)
library(glmmTMB)
library(lmtest)
library(ggeffects)
library(dplyr)

# Convert matrix into incidence matrix (0-1)

# Make pairs and calcualte beta.pair.abund and beta.pair
api1 <- decostand(api[stages$succession == 1, ], "hel") 
api2 <- decostand(api[stages$succession == 2, ], "hel")
api3 <- decostand(api[stages$succession == 3, ], "hel")

sph1 <- decostand(sph[stages$succession == 1, ], "hel")
sph2 <- decostand(sph[stages$succession == 2, ], "hel")
sph3 <- decostand(sph[stages$succession == 3, ], "hel")

chr1 <- decostand(chr[stages$succession == 1, ], "hel")
chr2 <- decostand(chr[stages$succession == 2, ], "hel")
chr3 <- decostand(chr[stages$succession == 3, ], "hel")

randompairs <- function(reps, ind1, ind2){
  chain1 <- sample(1:ind1, reps, replace = TRUE)
  chain2 <- sample(1:ind2, reps, replace = TRUE)
  return(data.frame(chain1,chain2))
}

# Returns indices of unique site combination
get_indices <- function(mat1, mat2){
  mat <- matrix(1, nrow = dim(mat1)[1], ncol=dim(mat2)[1])
  mat <- lower.tri(mat)
  return(which(mat, arr.ind = T))
}

# mat <- matrix(1, nrow = dim(api2)[1], ncol=dim(api3)[1])
# which(mat==1, arr.ind = T)

get_indices2 <- function(mat1, mat2){
  # Because I am comparing two different stages of succession, any comparison of sites will be valid, (e.g. even 1 vs 1) because these are from different successional stages.
  mat <- matrix(1, nrow = dim(mat1)[1], ncol=dim(mat2)[1])
  return(which(mat==1, arr.ind = T))
}

# Main function
compute_pairs <- function(api1, api2, 
                          ...){
  
  # Get only unique set of pairs
  indices <- get_indices2(api1,api2)
  
  result <- data.frame()
  for(pair in 1:dim(indices)[1]){
    bpa <- beta.pair.abund(rbind(api1[indices[pair,1],],
                                 api2[indices[pair,2],]))
    result <- rbind(result, 
                    c(bpa$beta.bray.bal,
                      bpa$beta.bray.gra,
                      bpa$beta.bray,
                      indices[pair,1]))
  }
  colnames(result) <- c("Balanced", "Gradient", "Bray", "InitSite")
  return(result)
}

a12 <- compute_pairs(api1,api2)
a13 <- compute_pairs(api1,api3)
a23 <- compute_pairs(api2,api3)

s12 <- compute_pairs(sph1, sph2)
s13 <- compute_pairs(sph1, sph3)
s23 <- compute_pairs(sph2, sph3)

c12 <- compute_pairs(chr1, chr2)
c13 <- compute_pairs(chr1, chr3)
c23 <- compute_pairs(chr2, chr3)

dset <- data.frame(value = c(a12$Balanced,
                             a13$Balanced,
                             a23$Balanced,
                             a12$Gradient,
                             a13$Gradient,
                             a23$Gradient,
                             s12$Balanced,
                             s13$Balanced,
                             s23$Balanced,
                             s12$Gradient,
                             s13$Gradient,
                             s23$Gradient,
                             c12$Balanced,
                             c13$Balanced,
                             c23$Balanced,
                             c12$Gradient,
                             c13$Gradient,
                             c23$Gradient),
                   initSite = c(a12$InitSite,
                             a13$InitSite,
                             a23$InitSite,
                             a12$InitSite,
                             a13$InitSite,
                             a23$InitSite,
                             s12$InitSite,
                             s13$InitSite,
                             s23$InitSite,
                             s12$InitSite,
                             s13$InitSite,
                             s23$InitSite,
                             c12$InitSite,
                             c13$InitSite,
                             c23$InitSite,
                             c12$InitSite,
                             c13$InitSite,
                             c23$InitSite),
                   type = rep(c("Balanced I  vs II", "Balanced I vs III",
                                "Balanced II vs III",
                                "Gradient I vs II", "Gradient I vs III",
                                "Gradient II vs III",
                                "Balanced I  vs II", "Balanced I vs III",
                                "Balanced II vs III",
                                "Gradient I vs II", "Gradient I vs III",
                                "Gradient II vs III",
                                "Balanced I  vs II", "Balanced I vs III",
                                "Balanced II vs III",
                                "Gradient I vs II", "Gradient I vs III",
                                "Gradient II vs III"),
                              c(dim(a12)[1],dim(a13)[1],dim(a23)[1],
                                dim(a12)[1],dim(a13)[1],dim(a23)[1],
                                dim(s12)[1],dim(s13)[1],dim(s23)[1],
                                dim(s12)[1],dim(s13)[1],dim(s23)[1],
                                dim(c12)[1],dim(c13)[1],dim(c23)[1],
                                dim(c12)[1],dim(c13)[1],dim(c23)[1])),
                   comparison = rep(c("I vs II", "I vs III","II vs III",
                                      "I vs II", "I vs III","II vs III",
                                      "I vs II", "I vs III","II vs III",
                                      "I vs II", "I vs III","II vs III",
                                      "I vs II", "I vs III","II vs III",
                                      "I vs II", "I vs III","II vs III"),
                                    c(dim(a12)[1],dim(a13)[1],dim(a23)[1],
                                      dim(a12)[1],dim(a13)[1],dim(a23)[1],
                                      dim(s12)[1],dim(s13)[1],dim(s23)[1],
                                      dim(s12)[1],dim(s13)[1],dim(s23)[1],
                                      dim(c12)[1],dim(c13)[1],dim(c23)[1],
                                      dim(c12)[1],dim(c13)[1],dim(c23)[1])),
                   group = rep(c("Herbivores", "Predators", "Kleptoparasites"),
                               c(2*dim(a12)[1]+2*dim(a13)[1]+2*dim(a23)[1],
                                 2*dim(s12)[1]+2*dim(s13)[1]+2*dim(s23)[1],
                                 2*dim(c12)[1]+2*dim(c13)[1]+2*dim(c23)[1])))


dset$diff <- substr(dset$type,1,8)

# TESTS ----
# logit <-function(x){log(x/(1-x))}

# I am making here an assumption that no comparison can be compeletly dissimilar 1 or similar 0
dset[dset$value == 0, ]$value <- 0.001
dset[dset$value == 1, ]$value <- 0.999

# lk <- logit(dset[dset$group == "Kleptoparasites", ]$value)
# lp <- logit(dset[dset$group == "Predators", ]$value)
# lh <- logit(dset[dset$group == "Herbivores", ]$value)

# # Logit transformation didn't do a good job
# qqnorm(lk)
# qqline(lk)
# qqnorm(lp)
# qqline(lp)
# qqnorm(lh)
# qqline(lh)

# Make an interaction term for pairtwise comparisons
dset$fInt<- paste(substr(dset$type, 1,1),
                     gsub(" ", "", dset$comparison),
                     substr(dset$group, 1, 1),
                     sep="")

# Comparison of II-nd and III-rd stage should always be lower
# see what i can get if I only analyze IvsII and IvsIII.

# Use beta regression for pairwise comparisons
# Balanced component comparison

# Data for genaral Bray-Curtis index comparison ----
bs1 <- dset[dset$diff == "Balanced", c("value","initSite",
                                       "comparison",
                                       "group",
                                       "diff","fInt")]
gs2 <- dset[dset$diff == "Gradient", c("value","initSite",
                                       "comparison",
                                       "group",
                                       "diff", "fInt")]

names(bs1) <- paste("b", names(bs1), sep="_")
names(gs2) <- paste("g", names(gs2), sep="_")

gset <- cbind(bs1, gs2)
gset$bcvals <- gset$b_value+gset$g_value

gset[gset$bcvals == 0, ]$bcvals <- 0.001
gset[gset$bcvals == 1, ]$bcvals <- 0.999

names(gset)[5] <- "fInt"

# Add a random effect
gset$reff <- paste(gset$b_initSite, 
                   substr(gset$b_comparison, 1, 2), sep = "")

# Get rid of the I vs III comparison
gset <- gset[gset$b_comparison != "I vs III", ]

# BRAY_CURTIS DISSIMILARITY ----
# br <- betareg(bcvals~b_group*b_comparison, data=gset, link="logit")

br <- glmmTMB(bcvals ~ b_group*b_comparison, data = gset, family= beta_family(link = "logit"))
brrand <- glmmTMB(bcvals ~ b_group*b_comparison + (1|reff), data = gset, family= beta_family(link = "logit"))

brrandtest <- emmeans(brrand, pairwise ~ b_comparison | b_group)
# brtest <- emmeans(br, pairwise ~ b_comparison | b_group)

# Compare model withouth and with a random effect
AIC(br, brrand)
lrtest(br, brrand)
# Seems like trandom effect model fits data better

# Predict values for the plot
bcpred <- ggpredict(brrand, terms = c("b_group", "b_comparison"),
                   type = "re")

cfs <- as.data.frame(bcpred)
names(cfs) <- c("b_group",
                "bcvals",
                "SE",
                "lcl",
                "ucl",
                "b_comparison")

gset$b_group <- recode(gset$b_group, Kleptoparasites = "Parasitic species")
cfs$b_group <- recode(cfs$b_group, Kleptoparasites = "Parasitic species")

bgcb <- ggplot(gset, aes(b_comparison, bcvals)) +
  geom_jitter(width = 0.1,aes(color=b_group), alpha = 0.2)+
  geom_line(aes(group = 1, color = b_group), 
            lty = c(2,2,1,1,1,1), data = cfs)+
  geom_errorbar(aes(ymin = lcl, ymax = ucl, color = b_group),
                data = cfs, width = 0, lwd=1) +
  geom_point(data = cfs, size = 3, aes(color = b_group))+
  facet_wrap(~b_group)+
  scale_color_manual(values=c(colvec[1], 
                              colvec[2], 
                              colvec[3]))+
  ylab("Bray-Curtis dissimilarity")+xlab("")
 
bgcb  + theme_bw()+theme(legend.position = "none")

# pdf("fig3b_bc.pdf", width=12, height=4)
# jpeg("fig3b_bc.jpg", width = 1200, height = 400, res = 150)
# pdf("revision_1/figures/Fig5.pdf", width=12, height=4,
#     onefile = FALSE)
# bgcb  + theme_bw()+theme(legend.position = "none")
# dev.off()

# GRADIENT COMPONENT----
# condition <- dset$diff == "Balanced"
# condition <- dset$diff == "Gradient"
# condition <- (dset$diff == "Balanced" & dset$comparison != "I vs III")
condition <- (dset$diff == "Gradient" & dset$comparison != "I vs III")
dset$reff <- paste(dset$initSite, 
                   gsub(" ","",substr(dset$comparison,1,2)), 
                   sep="_")

null_gra_rand <- glmmTMB(value ~ 1 + (1|reff), 
                    data = dset[condition,], 
                    family= beta_family(link = "logit"))
gra_rand <- glmmTMB(value ~ comparison*group + (1|reff), 
                    data = dset[condition,], 
                    family= beta_family(link = "logit"))

summary(gra_rand)
gra_rand_test <- emmeans(gra_rand, pairwise ~ comparison | group)
summary(gra_rand_test)

grapred <- ggpredict(gra_rand, terms = c("group", "comparison"),
                    type = "re")

cfs <- as.data.frame(grapred)
names(cfs) <- c("group",
                "value",
                "SE",
                "lcl",
                "ucl",
                "comparison")


library(dplyr)

data_used <- dset[condition,]
gset$b_group <- recode(gset$b_group, Kleptoparasites = "Parasitic species")
cfs$group <- recode(cfs$group, Kleptoparasites = "Parasitic species")

data_used$group <- recode(data_used$group,
                                 Kleptoparasites = "Parasitic species")

bgcb <- ggplot(data_used, aes(comparison, value)) +
  geom_jitter(width = 0.1,aes(color=group), alpha = 0.2)+
  geom_line(aes(group = 1, color = group), lty = c(2,2,2,
                                                   2,2,2), data = cfs)+
  geom_errorbar(aes(ymin = lcl, ymax = ucl, color = group),
                data = cfs, width = 0, lwd=1) +
  geom_point(data = cfs, size = 3, aes(color = group))+
  facet_wrap(~group)+
  scale_color_manual(values=c(alpha(colvec[1],1), 
                              alpha(colvec[2],1), 
                              alpha(colvec[3],1)))+
  # ggtitle("Balanced component")
  ylab("")+xlab("")+
  ggtitle("Gradient component")

# pdf("fig3c_gradient.pdf", width=12, height=4)
# jpeg("fig3c_gradient.jpg", width = 1200, height = 400, res = 150)
pdf("revision_1/figures/FigA3.pdf", width=12, height=4,
    onefile = FALSE)
bgcb  + theme_bw()+theme(legend.position = "none")
dev.off()

# br1 <- betareg(value~comparison*group, 
               # data = dset[condition,], 
               # link = "logit")

# BALANCED ----

condition <- (dset$diff == "Balanced" & dset$comparison != "I vs III")

null_bal_rand <- glmmTMB(value ~ 1 + (1|reff), 
                    data = dset[condition,], 
                    family= beta_family(link = "logit"))
bal_rand <- glmmTMB(value ~ comparison*group + (1|reff), 
                    data = dset[condition,], 
                    family= beta_family(link = "logit"))

summary(bal_rand)
bal_rand_test <- emmeans(bal_rand, pairwise ~ comparison | group)
summary(bal_rand_test)

balpred <- ggpredict(bal_rand, terms = c("group", "comparison"),
                     type = "re")

cfs <- as.data.frame(balpred)
names(cfs) <- c("group",
                "value",
                "SE",
                "lcl",
                "ucl",
                "comparison")

bgcb <- ggplot(dset[condition,], aes(comparison, value)) +
  geom_jitter(width = 0.1,aes(color=group), alpha = 0.2)+
  geom_line(aes(group = 1, color = group), lwd = 1, lty = c(2,2,1,
                                                   1,1,1), data = cfs)+
  geom_errorbar(aes(ymin = lcl, ymax = ucl, color = group),
                data = cfs, width = 0, lwd=1) +
  geom_point(data = cfs, size = 3, aes(color = group))+
  facet_wrap(~group)+
  scale_color_manual(values=c(alpha(colvec[1],1), 
                              alpha(colvec[2],1), 
                              alpha(colvec[3],1)))+
  ylab("")+xlab("")+
  ggtitle("Balanced component")

# pdf("fig3c_gradient.pdf", width=12, height=4)
jpeg("fig3c_balanced.jpg", width = 1200, height = 500, res = 150)
bgcb + theme_bw()+theme(legend.position = "none")
dev.off()

# RESULTS TABLE ----
br <- summary(bal_rand)
anova(null_bal_rand, bal_rand)
summary(gra_rand)
anova(null_bal_rand, gra_rand)
summary(brrand)

write.table(rbind(as.data.frame(balpred),
      as.data.frame(grapred),
      as.data.frame(bcpred)), "beta_random_model.txt")
library(insight)
get_variance(bal_rand)
get_variance(gra_rand)
get_variance(brrand)

# fixed effects variance (manually - not sure how to get
# predict.glmmTMB to do this)
linear.predictor <- model.matrix(bal_rand) %*% fixef(bal_rand)$cond
fixed.var <- var(linear.predictor)
# sum variance of all random effects ***excluding OLRE***
all.ranef.var <- unlist(VarCorr(bal_rand)$cond)
ranef.var <- all.ranef.var[!names(all.ranef.var) %in% "obs"]

linear.predictor <- model.matrix(gra_rand) %*% fixef(gra_rand)$cond
fixed.var <- var(linear.predictor)
# sum variance of all random effects ***excluding OLRE***
all.ranef.var <- unlist(VarCorr(gra_rand)$cond)
ranef.var <- all.ranef.var[!names(all.ranef.var) %in% "obs"]

# Pairwise cumulative ----
ba12 <- beta.pair.abund(rbind(colSums(api1),colSums(api2)))
ba13 <- beta.pair.abund(rbind(colSums(api1),colSums(api3)))
ba23 <- beta.pair.abund(rbind(colSums(api2),colSums(api3)))

bs12 <- beta.pair.abund(rbind(colSums(sph1),colSums(sph2)))
bs13 <- beta.pair.abund(rbind(colSums(sph1),colSums(sph3)))
bs23 <- beta.pair.abund(rbind(colSums(sph2),colSums(sph3)))

bc12 <- beta.pair.abund(rbind(colSums(chr1),colSums(chr2)))
bc13 <- beta.pair.abund(rbind(colSums(chr1),colSums(chr3)))
bc23 <- beta.pair.abund(rbind(colSums(chr2),colSums(chr3)))

pwcomp <- data.frame(val = as.numeric(c(ba12,ba13,ba23,
                           bs12,bs13,bs23,
                           bc12,bc13,bc23)),
                     comp = rep(c("Stage I vs II",
                                  "Stage I vs III",
                                  "Stage II vs III",
                                  "Stage I vs II",
                                  "Stage I vs III",
                                  "Stage II vs III",
                                  "Stage I vs II",
                                  "Stage I vs III",
                                  "Stage II vs III"),each=3),
                     group = rep(c("Herbivores", 
                                   "Predators", 
                                   "Kleptoparasites"),
                                 c(9,9,9)),
                     beta = rep(c("Balanced",
                                  "Gradient",
                                  "Bray-Curtis"),9))

pwcomp$comp <- as.character(pwcomp$comp)
pwcomp$beta <- factor(pwcomp$beta,levels=c("Balanced",
                                           "Gradient",
                                           "Bray-Curtis"))

pwplot <- ggplot(pwcomp, aes(fill=beta, y=val, x=comp)) + 
  geom_bar(position="dodge", stat="identity") + 
  ylim(c(0,0.65)) +
  xlab("")+ylab("")+
  scale_fill_manual(values=c(colvec[1], 
                              colvec[2], 
                              colvec[3])) +
  facet_wrap(~group,
             drop = T,
             scales = "free") + theme_bw()

pdf("fig3_beta.pdf", width = 12, height = 4)
pwplot
dev.off()

# Incidence matrix ----

# Convert matrix into incidence matrix (0-1)
apiI <- api
apiI[apiI > 0] <- 1
sphI <- sph
sphI[sphI > 0] <- 1
chrI <- chr
chrI[chrI > 0] <- 1

# Make pairs and calcualte beta.pair.abund and beta.pair
api1 <- apiI[stages$succession == 1, ]
api2 <- apiI[stages$succession == 2, ]
api3 <- apiI[stages$succession == 3, ]

sph1 <- sphI[stages$succession == 1, ]
sph2 <- sphI[stages$succession == 2, ]
sph3 <- sphI[stages$succession == 3, ]

chr1 <- chrI[stages$succession == 1, ]
chr2 <- chrI[stages$succession == 2, ]
chr3 <- chrI[stages$succession == 3, ]

a12 <- compute_pairs(api1,api2)
a13 <- compute_pairs(api1,api3)

s12 <- compute_pairs(sph1, sph2)
s13 <- compute_pairs(sph1, sph3)

c12 <- compute_pairs(chr1, chr2)
c13 <- compute_pairs(chr1, chr3)

dset <- data.frame(value = c(a12$Balanced,
                             a13$Balanced,
                             a12$Gradient,
                             a13$Gradient,
                             s12$Balanced,
                             s13$Balanced,
                             s12$Gradient,
                             s13$Gradient,
                             c12$Balanced,
                             c13$Balanced,
                             c12$Gradient,
                             c13$Gradient),
                   type = rep(c("Balanced I  vs II", "Balanced I vs III",
                                "Gradient I vs II", "Gradient I vs III",
                                "Balanced I  vs II", "Balanced I vs III",
                                "Gradient I vs II", "Gradient I vs III",
                                "Balanced I  vs II", "Balanced I vs III",
                                "Gradient I vs II", "Gradient I vs III"),
                              c(dim(a12)[1],dim(a13)[1],
                                dim(a12)[1],dim(a13)[1],
                                dim(s12)[1],dim(s13)[1],
                                dim(s12)[1],dim(s13)[1],
                                dim(c12)[1],dim(c13)[1],
                                dim(c12)[1],dim(c13)[1])),
                   comparison = rep(c("I vs II", "I vs III",
                                      "I vs II", "I vs III",
                                      "I vs II", "I vs III",
                                      "I vs II", "I vs III",
                                      "I vs II", "I vs III",
                                      "I vs II", "I vs III"),
                                    c(dim(a12)[1],dim(a13)[1],
                                      dim(a12)[1],dim(a13)[1],
                                      dim(s12)[1],dim(s13)[1],
                                      dim(s12)[1],dim(s13)[1],
                                      dim(c12)[1],dim(c13)[1],
                                      dim(c12)[1],dim(c13)[1])),
                   group = rep(c("Apoidae", "Spheciformes", "Chrysididae"),
                               c(dim(a12)[1]+dim(a13)[1]+dim(a12)[1]+dim(a13)[1],
                                 dim(s12)[1]+dim(s13)[1]+dim(s12)[1]+dim(s13)[1],
                                 dim(c12)[1]+dim(c13)[1]+dim(c12)[1]+dim(c13)[1])))

plotapi <- ggplot(dset, aes(y=log(value/(1-value)), x = type,color = type)) +
  geom_jitter(width = 0.2) +
  facet_wrap(~group, scales = "free")
plotapi

# individual group
apibal <- dset[dset$group == "Apoidae",]
apibal <- apibal[grep("Gradient", dset$type),]
apibal <- apibal[complete.cases(apibal),]
dist1 <- apibal[apibal$comparison == "I vs II",]
dist2 <- apibal[apibal$comparison == "I vs III",]

apibal <- dset[dset$group == "Apoidae",]
apibal <- apibal[grep("Balanced", dset$type),]
apibal <- apibal[complete.cases(apibal),]
dist1 <- apibal[apibal$comparison == "I vs II",]
dist2 <- apibal[apibal$comparison == "I vs III",]

t.test(dist1$value, dist2$value)

plotapi <- ggplot(dset, aes(value, color = type, 
                            linetype=comparison)) +
  # geom_density(aes(y = ..scaled..), adjust=5) +
  geom_density(adjust=5) +
  facet_wrap(~group, scales = "free") + 
  xlab("") + ylab("") +
  xlim(c(0,1))+
  theme_classic() +
  scale_color_manual(values=c(colvec[1], 
                              colvec[1], 
                              colvec[2],
                              colvec[2]))
# plotapi

# Pairwise cumulative
ba12 <- beta.pair.abund(rbind(colSums(api1),colSums(api2)))
ba13 <- beta.pair.abund(rbind(colSums(api1),colSums(api3)))
ba23 <- beta.pair.abund(rbind(colSums(api2),colSums(api3)))

bs12 <- beta.pair.abund(rbind(colSums(sph1),colSums(sph2)))
bs13 <- beta.pair.abund(rbind(colSums(sph1),colSums(sph3)))
bs23 <- beta.pair.abund(rbind(colSums(sph2),colSums(sph3)))

bc12 <- beta.pair.abund(rbind(colSums(chr1),colSums(chr2)))
bc13 <- beta.pair.abund(rbind(colSums(chr1),colSums(chr3)))
bc23 <- beta.pair.abund(rbind(colSums(chr2),colSums(chr3)))

pwcomp <- data.frame(val = as.numeric(c(ba12,ba13,ba23,
                                        bs12,bs13,bs23,
                                        bc12,bc13,bc23)),
                     comp = rep(c("Stage I vs II",
                                  "Stage I vs III",
                                  "Stage II vs III",
                                  "Stage I vs II",
                                  "Stage I vs III",
                                  "Stage II vs III",
                                  "Stage I vs II",
                                  "Stage I vs III",
                                  "Stage II vs III"),each=3),
                     group = rep(c("Apoidae", 
                                   "Spheciformes", 
                                   "Chrysididae"),
                                 c(9,9,9)),
                     beta = rep(c("Balanced",
                                  "Gradient",
                                  "Bray-Curtis"),9))

pwcomp$comp <- as.character(pwcomp$comp)
pwcomp$beta <- factor(pwcomp$beta,levels=c("Balanced",
                                           "Gradient",
                                           "Bray-Curtis"))

pwplot <- ggplot(pwcomp, aes(fill=beta, y=val, x=comp)) + 
  geom_bar(position="dodge", stat="identity") + 
  ylim(c(0,0.65)) +
  xlab("")+ylab("")+
  scale_fill_manual(values=c(colvec[1], 
                             colvec[2], 
                             colvec[3])) +
  facet_wrap(~group,
             drop = T,
             scales = "free") + theme_bw()


# pwplot



# Grouped

# # Partition of diversity for each group separately

# >>> Apiformes
partapi1 <- beta.multi.abund(api[stages$succession == 1, ])
partapi2 <- beta.multi.abund(api[stages$succession == 2, ])
partapi3 <- beta.multi.abund(api[stages$succession == 3, ])

# >>> Spheciformes
partsph1 <- beta.multi.abund(sph[stages$succession == 1, ])
partsph2 <- beta.multi.abund(sph[stages$succession == 2, ])
partsph3 <- beta.multi.abund(sph[stages$succession == 3, ])

# >>> Chrysididae
partchr1 <- beta.multi.abund(chr[stages$succession == 1, ])
partchr2 <- beta.multi.abund(chr[stages$succession == 2, ])
partchr3 <- beta.multi.abund(chr[stages$succession == 3, ])

# Combine the results into a data.frame

# <<<needs to be fixed... columns are lists >>>
table <- data.frame(rbind(partapi1,
                             partapi2,
                             partapi3,
                             partsph1,
                             partsph2,
                             partsph3,
                             partchr1,
                             partchr2,
                             partchr3
                             ),
                    group = rep(c("Herbivores",
                                     "Predators",
                                     "Kleptoparasites"),
                                   each=3),
                    stage = rep(c(1,2,3),3))