# Beta-diveristy partitioning ----

# Load data ----
source("data/data_processing.R")

# bray.part function
# Baselga, A. (2013). Separating the two components of abundance-based dissimilarity: balanced changes in abundance vs. abundance gradients. Methods in Ecology and Evolution, 4(6), 552–557. doi: 10.1111/2041-210X.12029

# Abundance‐based dissimilarity can be derived from two antithetic patterns: (i) balanced variation in abundance, whereby the individuals of some species in one site are substituted by the same number of individuals of different species in another site or, in other words, species abundances change from site to site with different signs for different species and changes balance each other; and (ii) abundance gradients, whereby some individuals are lost from one site to the other or, in other words, all the species that change their abundance from one site to the other make it with the same sign. The widely used Bray‐Curtis index of dissimilarity is the result of summing these two sources of dissimilarity, and therefore might consider equivalent patterns that are markedly different. The partition of dBC dissimilarity into two components of balanced abundance variation and abundance gradient, respectively, may thus be useful to assess biodiversity patterns and to explore their causes, as substitution and loss of individuals are patterns that can derive from completely different processes. 

library(betapart)
require(vegan)
library(ggplot2)
library(codyn)
library(emmeans)
library(multcomp)
library(brms)
library(betareg)



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

# Main function
compute_pairs <- function(api1, api2, reps=999, 
                          plt = 2, ...){
  
  # Get only unique set of pairs
  
  indices <- get_indices(api1,api2)
  
  result <- data.frame()
  for(pair in 1:dim(indices)[1]){
    bpa <- beta.pair.abund(rbind(api1[indices[pair,1],],
                                 api2[indices[pair,2],]))
    result <- rbind(result, 
                    c(bpa$beta.bray.bal,
                      bpa$beta.bray.gra,
                      bpa$beta.bray))
  }
  colnames(result) <- c("Balanced", "Gradient", "Bray")
  return(result)
  
  # plot(density(result[,plt], adjust = 1, kernel = "gaussian"),
  #      xlim = c(0,  0.9), ylim = c(0,5),
  #      xlab = "",ylab = "",main="",...)
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


#PLOT ----

# plotapi <- ggplot(dset, aes(y=log(value/(1-value)), x = type,color = type)) +
#   geom_jitter(width = 0.2) +
#   facet_wrap(~group, scales = "free")
# plotapi

dset$diff <- substr(dset$type,1,8)

# plotapi <- ggplot(dset, aes(value, color = type,
#                             fill = type,
#                             linetype=comparison)) +
#   # geom_density(aes(y = ..scaled..), adjust=5) +
#   geom_density(adjust=5) +
#   facet_wrap(~group*diff, scales = "free") + 
#   xlab("") + ylab("") +
#   xlim(c(0,1))+
#   theme_classic() +
#   scale_color_manual(values=c(colvec[1], 
#                               colvec[2], 
#                               colvec[5],
#                               colvec[1],
#                               colvec[2],
#                               colvec[5])) +
#   scale_fill_manual(values = alpha(c(colvec[1], 
#                                      colvec[2], 
#                                      colvec[5],
#                                      colvec[1],
#                                      colvec[2],
#                                      colvec[5]),0.2))
# plotapi


# TESTS ----
logit <-function(x){log(x/(1-x))}

# I am making here an assumption that no comparison can be compeletly dissimilar 1 or similar 0
dset[dset$value == 0, ]$value <- 0.001
dset[dset$value == 1, ]$value <- 0.999

# lk <- logit(dset[dset$group == "Kleptoparasites", ]$value)
# lp <- logit(dset[dset$group == "Predators", ]$value)
# lh <- logit(dset[dset$group == "Herbivores", ]$value)
# 
# # Logit transformation doesn't do a good job ----
# qqnorm(lk)
# qqline(lk)
# qqnorm(lp)
# qqline(lp)
# qqnorm(lh)
# qqline(lh)

# Make an interaction term for pairtwise comparisons ----
dset$fInt<- paste(substr(dset$type, 1,1),
                     gsub(" ", "", dset$comparison),
                     substr(dset$group, 1, 1),
                     sep="")
# Comparison of II-nd and III-rd stage should always be lower
# see what i can get if I only analyze IvsII and IvsIII.

# Use beta regression for pairwise comparisons
# Balanced component comparison

# Data for genarl comparison ----
bs1 <- dset[dset$diff == "Balanced", c("value",
                                       "comparison",
                                       "group",
                                       "diff","fInt")]
gs2 <- dset[dset$diff == "Gradient", c("value",
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

# Get rid of the I vs III comparison
gset <- gset[gset$b_comparison != "I vs III", ]

# General bray-curtis ----
br <- betareg(bcvals~b_group*b_comparison, data=gset, link="logit")

summary(br)
brtest <- emmeans(br, pairwise ~ b_comparison | b_group)
# emmip(br, b_group ~ b_comparison)
summary(brtest)
names(brtest)
cfs <-rbind(brtest$emmeans[1], brtest$emmeans[2],
      brtest$emmeans[3], brtest$emmeans[4],
      brtest$emmeans[5], brtest$emmeans[6])
cfs <- as.data.frame(cfs)
names(cfs) <- c("b_comparison",
                "b_group",
                "bcvals",
                "SE",
                "df",
                "lcl",
                "ucl") 

bgcb <- ggplot(gset, aes(b_comparison, bcvals)) +
  geom_jitter(width = 0.1,aes(color=b_group), alpha = 0.2)+
  geom_line(aes(group = 1, color = b_group), lwd = c(0,0,1,1,1,1), data = cfs)+
  geom_errorbar(aes(ymin = lcl, ymax = ucl, color = b_group),
                data = cfs, width = 0, lwd=1) +
  geom_point(data = cfs, size = 3, aes(color = b_group))+
  facet_wrap(~b_group)+
  scale_color_manual(values=c(colvec[1], 
                              colvec[2], 
                              colvec[3]))
 
# bgcb

# pdf("fig3b_bc.pdf", width=12, height=4)
# bgcb
# dev.off()

# CONDITIONS ----
# condition <- dset$diff == "Balanced"
# condition <- dset$diff == "Gradient"
# condition <- (dset$diff == "Balanced" & dset$comparison != "I vs III")
condition <- (dset$diff == "Gradient" & dset$comparison != "I vs III")

# Components ----
br1 <- betareg(value~comparison*group, 
                        data = dset[condition,], 
                        link = "logit")
summary(br1)
brtest <- emmeans(br1, pairwise ~ comparison | group)
# emmip(br, b_group ~ b_comparison)
summary(brtest)
names(brtest)
cfs <-rbind(brtest$emmeans[1], brtest$emmeans[2],
            brtest$emmeans[3], brtest$emmeans[4],
            brtest$emmeans[5], brtest$emmeans[6])
cfs <- as.data.frame(cfs)
names(cfs) <- c("comparison",
                "group",
                "value",
                "SE",
                "df",
                "lcl",
                "ucl") 

bgcb <- ggplot(dset[condition,], aes(comparison, value)) +
  geom_jitter(width = 0.1,aes(color=group), alpha = 0.2)+
  geom_line(aes(group = 1, color = group), lty = c(2,2,2,
                                                   2,2,2), data = cfs)+
  geom_errorbar(aes(ymin = lcl, ymax = ucl, color = group),
                data = cfs, width = 0, lwd=1) +
  geom_point(data = cfs, size = 3, aes(color = group))+
  facet_wrap(~group)+
  scale_color_manual(values=c(colvec[1], 
                              colvec[2], 
                              colvec[3]))+
  # ggtitle("Balanced component")
  ggtitle("Gradient component")

# pdf("fig3c_balanced.pdf", width=12, height=4)
# bgcb
# dev.off()

# pdf("fig3c_gradient.pdf", width=12, height=4)
# bgcb
# dev.off()



# Mean beta-comparisons ----

condition <- dset$diff == "Balanced"
# condition <- dset$diff == "Gradient"

br2 <- betareg(value~fInt, 
                      data = dset[condition, ], 
                      link = "logit")
brtest <- emmeans(br2, "fInt")
brtest <- cld(brtest, Letter="abcdefghijklm")
brsum <- as.data.frame(brtest)
grlist <- strsplit(as.character(brsum$fInt), "I")
for(i in 1:length(grlist)){
  print(grlist[[i]])
  grlist[[i]] <- grlist[[i]][-(which(grlist[[i]] == ""))]
}
groups <- c()
for(i in 1:length(grlist)){
  gi <- grlist[[i]][length(grlist[[i]])]
  groups <- c(groups, gi)
}

brsum$group <- groups
dset$fComp <- as.numeric(dset$comparison)

mean_plot <- ggplot(brsum, aes(x = fInt, y = emmean, 
                  color = group,
                  label = .group)) + 
  geom_point()+
  geom_errorbar(aes(ymin=asymp.LCL, 
                    ymax=asymp.UCL), 
                width=.2,
                position=position_dodge(0.05))+
  stat_summary(fun=mean, geom="text", 
               col = rgb(10,10,10,180,maxColorValue = 255),
               hjust = 1.4,
               vjust = -1.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# pdf("fig3d_mean_comparison.pdf", width=8, height=8)
# mean_plot
# dev.off()

# pdf("fig3d_mean_comparison_gradient.pdf", width=8, height=8)
# mean_plot
# dev.off()


# Bayesian estimation ----

# That would have to be broken into two groups
get_prior(value~comparison*group, family="beta", 
          data = dset[condition, ])
brm1 <- brm(value~comparison*group, family="beta", 
            data = dset[condition, ])

pp = brms::pp_check(brm1)
pp
plot(brm1)
conditional_effects(brm1)

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



# Notes ----

# Grouped

# # Partition of diversity for each group separately
# 
# # >>> Apiformes

partapi1 <- beta.multi.abund(api[stages$succession == 1, ])
partapi2 <- beta.multi.abund(api[stages$succession == 2, ])
partapi3 <- beta.multi.abund(api[stages$succession == 3, ])
# 
# # >>> Spheciformes
partsph1 <- beta.multi.abund(sph[stages$succession == 1, ])
partsph2 <- beta.multi.abund(sph[stages$succession == 2, ])
partsph3 <- beta.multi.abund(sph[stages$succession == 3, ])
# 
# # >>> Chrysididae
partchr1 <- beta.multi.abund(chr[stages$succession == 1, ])
partchr2 <- beta.multi.abund(chr[stages$succession == 2, ])
partchr3 <- beta.multi.abund(chr[stages$succession == 3, ])
# 
# # Combine the results into a data.frame
# 
# # <<<needs to be fixed... columns are lists >>>
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

table
# 
# # Species turnover between stages
# 
# # >>> using beta.multi.abund
# 
# # Baselga, A. 2017. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution 8: 799-808
# 
# bapi12 <- beta.multi.abund(rac_mx[c(1,2), groups == "Apiformes"])
# bapi23 <- beta.multi.abund(rac_mx[c(2,3), groups == "Apiformes"])
# bapi13 <- beta.multi.abund(rac_mx[c(1,3), groups == "Apiformes"])
# 
# datapi <- rbind(BetaB = c(bapi12$beta.BRAY.BAL, bapi23$beta.BRAY.BAL, bapi13$beta.BRAY.BAL),
#                 BetaG = c(bapi12$beta.BRAY.GRA, bapi23$beta.BRAY.GRA, bapi13$beta.BRAY.GRA),
#                 BRAY = c(bapi12$beta.BRAY, bapi23$beta.BRAY, bapi13$beta.BRAY))
# colnames(datapi) <- c("I vs. II", "II vs. III", "I vs. III")
# 
# 
# # ggplot(datapi) + geom_bar()
# # beta.multi.abund(rac_mx[, groups == "Spheciformes"])
# # beta.multi.abund(rac_mx[, groups == "Chrysididae"])
# 


# # Plot
# api1 <- api[stages$succession == 1, ]
# api2 <- api[stages$succession == 2, ]
# api3 <- api[stages$succession == 3, ]
# 
# st <- 5
# rep <- 999
# bsa1 <- beta.sample.abund(api1, index.family="bray", sites=st, samples=rep)
# bsa2 <- beta.sample.abund(api2, index.family="bray", sites=st, samples=rep)
# bsa3 <- beta.sample.abund(api3, index.family="bray", sites=st, samples=rep)
# # 
# dat1 <- bsa1$sampled.values$beta.BRAY.BAL
# dat2 <- bsa2$sampled.values$beta.BRAY.BAL
# dat3 <- bsa3$sampled.values$beta.BRAY.BAL
# 
# par(mfrow = c(1,2))
# 
# plot(density(dat1, adjust = 6, kernel = "gaussian"),
#      xlim = c(0.5,  0.95), ylim = c(0,20),
#      lwd=2, col="navyblue")
# par(new=TRUE)
# plot(density(dat2,adjust = 6, kernel = "gaussian"),
#      xlim = c(0.5,  0.95), ylim = c(0,20),
#      lwd=2, col="orange")
# par(new=TRUE)
# plot(density(dat3,adjust = 6,kernel = "gaussian"),
#      xlim = c(0.5,  0.95), ylim = c(0,20),
#      lwd=2, col="red")
# 
# dat1 <- bsa1$sampled.values$beta.BRAY.GRA
# dat2 <- bsa2$sampled.values$beta.BRAY.GRA
# dat3 <- bsa3$sampled.values$beta.BRAY.GRA
# 
# plot(density(dat1, adjust = 6, kernel = "gaussian"),
#      lty=2,xlim = c(0,  0.2), ylim = c(0,20),
#      lwd=2, col="navyblue")
# par(new=TRUE)
# plot(density(dat2, adjust = 6, kernel = "gaussian"),
#      lty=2,xlim = c(0,  0.2), ylim = c(0,20),
#      lwd=2, col="orange")
# par(new=TRUE)
# plot(density(dat3, adjust = 6, kernel = "gaussian"),
#      lty=2,xlim = c(0,  0.2), ylim = c(0,20),
#      lwd=2, col="red")

# <<<<
# Monte Carlo comparison
# set.seed(11)                              # Today's date in the US - no cherry-picking!
# r = 1:6                                   # The possible ranks of our non-zero differences
# nsim = 1e5 
# V = 0 
# for (i in 1:nsim){ 
#   rank = sample(r)                        # Sampling the ranks...
#   sign = sample(c(1, -1), 6, replace = T) #... and the signs for each rank.
#   V[i] = sum(rank[sign > 0])              # Doing the sum to get the V.
# } 
# (p_value <- sum(V <= 13) / nsim)

# runs <- 1000000
# s1 <- sample(d1d2[,1], runs, replace = T)
# s2 <- sample(d1d2[,2], runs, replace = T)
# mc.p.value <- sum((s2-s1) < 0)/runs
# mc.p.value
# hist(s1/s2)