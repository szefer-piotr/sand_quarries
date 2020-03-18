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

# Make pairs and calcualte beta.pair.abund and beta.pair

api1 <- api[stages$succession == 1, ]
api2 <- api[stages$succession == 2, ]
api3 <- api[stages$succession == 3, ]

sph1 <- sph[stages$succession == 1, ]
sph2 <- sph[stages$succession == 2, ]
sph3 <- sph[stages$succession == 3, ]

chr1 <- chr[stages$succession == 1, ]
chr2 <- chr[stages$succession == 2, ]
chr3 <- chr[stages$succession == 3, ]

randompairs <- function(reps, ind1, ind2){
  chain1 <- sample(1:ind1, reps, replace = TRUE)
  chain2 <- sample(1:ind2, reps, replace = TRUE)
  return(data.frame(chain1,chain2))
}

get_indices <- function(mat1, mat2){
  mat <- matrix(1, nrow = dim(mat1)[1], ncol=dim(mat2)[1])
  mat <- lower.tri(mat)
  return(which(mat, arr.ind = T))
}

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
# plotapi

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

rank(sort(dist1$value) - sort(dist2$value))

# KS test
ks.test(dist1$value,dist2$value)
ks.test(dist1$value,dist2$value)

d1d2 <- cbind(sort(dist1$value), sort(dist2$value))
ks.test(d1d2[,1],d1d2[,2])
# diff <- d1d2[,1]-d1d2[,2]
# rank(diff[-7])
# duplicated(diff)
# d1d2rem <- d1d2[-7,]
# ks.test(d1d2rem[,1],d1d2rem[,2])
# ks.test(diff[-7])
# duplicated(d1d2rem[,1]-d1d2rem[,2])
# 
# wilcox.test(diff[-7])

###################################
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

#################

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

# # Partition of diversity for each group separately ----
# 
# # >>> Apiformes ----
# partapi1 <- beta.multi.abund(api[stages$succession == 1, ])
# partapi2 <- beta.multi.abund(api[stages$succession == 2, ])
# partapi3 <- beta.multi.abund(api[stages$succession == 3, ])
# 
# # >>> Spheciformes ----
# partsph1 <- beta.multi.abund(sph[stages$succession == 1, ])
# partsph2 <- beta.multi.abund(sph[stages$succession == 2, ])
# partsph3 <- beta.multi.abund(sph[stages$succession == 3, ])
# 
# # >>> Chrysididae ----
# partchr1 <- beta.multi.abund(chr[stages$succession == 1, ])
# partchr2 <- beta.multi.abund(chr[stages$succession == 2, ])
# partchr3 <- beta.multi.abund(chr[stages$succession == 3, ])
# 
# # Combine the results into a data.frame ----
# 
# # <<<needs to be fixed... columns are lists >>>
# table <- data.frame(rbind(partapi1,
#                              partapi2,
#                              partapi3,
#                              partsph1,
#                              partsph2,
#                              partsph3,
#                              partchr1,
#                              partchr2,
#                              partchr3
#                              ),
#                     group = rep(c("Apiformes", 
#                                      "Spheciformes",
#                                      "Chrysididae"),
#                                    each=3),
#                     stage = rep(c(1,2,3),3))
# 
# # Species turnover between stages ----
# 
# # >>> using beta.multi.abund ----
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



# # Plot ----
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
