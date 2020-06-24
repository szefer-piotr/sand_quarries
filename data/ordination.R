# Source data processing script ----
source("data/data_processing.R")
library(vegan)

# Transformation
ascnt <- asc
asc <- decostand(asc, "hel")


# 1. RDA ordination ----
rda1 <- rda(asc~as.factor(succession), data=stages)
rdaa <- rda(asc[,groups=="Herbivores"]~as.factor(succession), data=stages)
rdas <- rda(asc[,groups=="Predators"]~as.factor(succession), data=stages)
rdac <- rda(asc[,groups=="Kleptoparasites"]~as.factor(succession), data=stages)
# cca1 <- cca(asc~succession, data=stages)

# m <- cca1
m <- rda1

#>>> Statistical test ----
anova(m, by="terms", permutations = 999)
anova(m, by="axis", permutations = 999)
# both axes are significant
# Succession is a statistically significant factor

# >>> R2 -----
coef(m)
R2 <- RsquareAdj(m)$r.squared
R2

# RDA: Only 12.45 % of the variance explained for the whole dataset
# CCA: Only 10.96 % of the variance explained for the whole dataset


# 2. Variance explained -----
sumeig <- summary(eigenvals(m))
round(sumeig[2,1],3)
round(sumeig[2,2],3)
# >>> goodness of fit ----

# GOF - total inertia explained for given number of axes
# gof <- "goodness"(m, display = "species", 
#                   statistic = "explained",
#                   model="CCA",
#                   summarize = FALSE)

# Fit species onto ordination
spc <- envfit(m, asc)
# write.table(cbind(rownames(spc$vectors),
#                   spc$vectors$arrows,
#                   spc$vectors$r,
#                   spc$vectors$pvals), "envfit.txt")
# selected_gof <- names(spc$vectors$pvals)[which(spc$vectors$pvals <= 0.05)]

# pick both axes
# specs <-sort(round(gof[,1]+gof[,2], 3), decreasing = T)
# Firrst 10 speecies
# perc <- 0.25
# selected_gof <- names(specs[specs>=perc])

# decomposes the inertia into partial,  constrained and unconstrained components for each site or species. Legendre & De Cáceres (2012) called these inertia components as local contributions to beta-diversity (LCBD) and species contributions to beta-diversity (SCBD),and they give these as relative contributions summing up to unity 

# gof2sp <- inertcomp(rda1, display = "species", proportional = F)
# gof2si <- inertcomp(rda1, display = "sites", proportional = T)

# plot(colSums(log(asc+1))~ gof2[,2])
# par(new=T)
# plot(colSums(log(asc[1:10,]+1))~ gof2[,2], col="red")

# names(colSums(log(asc[1:10,]+1))) == names(gof2[,2])
# Check out how species selected in goodness fit here
# gof2[rownames(gof2) %in% selected_gof,]


#>>>Plot ----
# par(mfrow=c(1,1))
# x11(10,10)

# png("figs/fig3.png", width=7, height=7, units = "in", res = 800)

# Main plot - species
# scl = 3
# plot(m, type = "n", scaling = scl, 
#      xlab = paste("RDA1 [", 
#                   round(sumeig[2,1]*100,2),
#                   "%]", 
#                   sep=""),
#      ylab = paste("RDA2 [",
#                   round(sumeig[2,2]*100,2),
#                   "%]", 
#                   sep=""), 
#      cex.lab = 1.5)

# plot(m, type = "n", scaling = scl, xlab = "CCA1 [xxx%]",
#      ylab = "CCA2 [xxx%]", cex.lab = 1.5)

# Species

# withouth selection
as.factor(groups)
names(summary(rda1))
which(colnames(asc) == rownames(summary(rda1)$species))

make_rda_plot <- function(m, scl=3, type = "Apiformes", color){
  plot(m, type = "n", scaling = scl, 
       xlab = paste("RDA1 [", 
                    round(sumeig[2,1]*100,2),
                    "%]", 
                    sep=""),
       ylab = paste("RDA2 [",
                    round(sumeig[2,2]*100,2),
                    "%]", 
                    sep=""), 
       cex.lab = 1.5)
  
  points(m, display = "species",
         col = colvec[color],
         scaling = scl, pch=19, cex = 1.5,
         select = colnames(asc[groups == type]))
  
  points(m, display="cn", scaling = scl,
         col = "red",
         pch = 25)
  
  text(m, display="cn", scaling = scl,
       labels = c("Early",
                  "Mid",
                  "Late"),
       col = "red", lwd=1.5,
       pos = 3)
  
  grnames <- colnames(asc[,groups == type])
  
  lgn <- length(selected_gof[selected_gof %in% grnames])
  
  if(lgn != 0){
    text(m, display="species", scaling=scl,
         select = selected_gof[selected_gof %in% grnames],
         col=cgray)
  }
}


pdf("fig2_rda.pdf", width = 12, height = 4)
par(mfrow = c(1,3))
make_rda_plot(m, 3, "Herbivores", 1)
make_rda_plot(m, 3, "Predators", 3)
make_rda_plot(m, 3, "Kleptoparasites", 2)
dev.off()

# points(m, display = "species",
#        col = colvec[1],
#        scaling = scl, pch=19, cex = 1.5,
#        select = colnames(asc[groups == "Apiformes"]))
# 
# points(m, display = "species",
#        col = colvec[3],
#        scaling = scl, pch=19, cex = 1.5,
#        select = colnames(asc[groups == "Spheciformes"]))
# 
# points(m, display = "species",
#        col = colvec[2],
#        scaling = scl, pch=19, cex = 1.5,
#        select = colnames(asc[groups == "Chrysididae"]))
# 
# points(m, display="cn", scaling = scl,
#        col = "red",
#        pch = 25)
# 
# text(m, display="cn", scaling = scl,
#      labels = c("Early",
#                 "Mid",
#                 "Late"),
#      col = "red", lwd=1.5,
#      pos = 3)
# # 
# text(m, display="species", scaling=scl,
#      select = selected_gof, col=cgray)

# # Legend
# with(test, legend("topright",
#                   legend = c("C","F","I","P","H2","H1"), bty = "n",
#                   col = colvec, pch = pch, pt.bg = colvec))
# 
# barplot(m$CA$eig/m$tot.chi, names.arg = 1:m$CA$rank, cex.names = 0.5, 
#         ylab="Proportion of variance explained", xlab="CA axis")
# m$CCA$tot.chi/m$tot.chi


# Are some groups responding stronger to succession?

# What are the species associated with different stages?

# 3. Separate ordination for each group ----

# makerda <- function(mtrx, stages, selcol, perc=0.25){
#   
#   # Model specification
#   m <- rda(log(mtrx+1)~as.factor(succession), data=stages)
#   
#   # Test
#   print(anova(m, by="terms", permutations = 999))
#   print(anova(m, by="axis", permutations = 999) )
#   
#   # R square
#   print(RsquareAdj(m)$r.squared)
#   
#   # Variance explained
#   sev <- summary(eigenvals(m))
#   varx1 <- sev[2,1]*100
#   varx2 <- sev[2,2]*100
#   
#   # GOF - total inertia explained for given number of axes
#   gof <- goodness(m, display = "species", statistic = "explained",
#                   summarize = FALSE)
#   
#   # Sum up two axes
#   specs <-sort(round(gof[,1]+gof[,2], 3), decreasing = T)
#   # Firrst 10 speecies
#   selected_gof <- names(specs[specs>=perc])
#   
#   # gof2 <- inertcomp(m, display = "species", proportional = T)
#   # Check out how species selected in goodness fit here
#   # gof2[rownames(gof2) %in% selected_gof,]
#   
#   # The above decomposes the inertia into partial,constrained and 
#   # unconstrained components for each site or species. 
#   # Legendre & De Cáceres (2012) called these inertia components 
#   # as local contributions to beta-diversity (LCBD) and species 
#   # contributions to beta-diversity (SCBD),and they give these as 
#   # relative contributions summing up to unity 
# 
#   # Main plot - species
#   scl = 3
#   plot(m, type = "n", scaling = scl, 
#        xlab = paste("RDA1","[",round(varx1,3),"%]", sep=""),
#        ylab = paste("RDA2","[",round(varx2,3),"%]", sep=""), 
#        cex.lab = 1)
#   
#   points(m, display = "species",
#          col = selcol,
#          scaling = scl, pch=19, cex = 1.5)
#   
#   points(m, display="cn", scaling = scl, 
#          col = cgray,
#          pch = 25)
#   
#   text(m, display="cn", scaling = scl,
#        labels = c("Early",
#                   "Mid",
#                   "Late"),
#        col = "red", lwd=1.5)
#   
#   text(m, display="species", scaling=scl,
#        select = selected_gof, col=cgray)
# }
# 
# par(mfrow=c(1,3))
# makerdaplot(api, stages, colvec[1])
# makerdaplot(sph, stages, colvec[3])
# makerdaplot(chr, stages, colvec[2]) 
# par(mfrow=c(1,1))

# 4. NMDS ----
# library(ggplot2)
# library(RColorBrewer)
# cols <- brewer.pal(3, "Set1")
# fGroups <- as.numeric(as.factor(groups))
# par(mfrow=c(1,1))
# nmds <-metaMDS(asc)
# png("nmds.png", width = 500, height=500)
# plot(nmds, display = "species", type = "n")
# points(nmds, display = "species", 
#        col = alpha(cols[fGroups], 0.9), 
#        pch = stages$succession+14, cex = 1.5)
# with(nmds, legend("topleft",
#                   legend = c("Apoidea Stage I",
#                              "Apoidea Stage II",
#                              "Apoidea Stage III",
#                              "Spheciformes Stage I",
#                              "Spheciformes Stage II",
#                              "Spheciformes Stage III",
#                              "Chrysididae Stage I",
#                              "Chrysididae Stage II",
#                              "Chrysididae Stage III"
#                              ), bty = "n",
#                   col = rep(c(cols[1], cols[3], cols[2]), each=3), 
#                   pch = rep(unique(stages$succession)+14,3), 
#                   pt.bg = colvec, cex=0.5))
# dev.off()
# site_rich <- rowSums(asc>0)
# gof2si
# 
# dat <- as.data.frame(cbind(site_rich, gof2si))
# colnames(dat) <- c("rich", "constr", "unconstr")
# 
# lmcon <- lm(rich~constr, data=dat)
# summary(lmcon)
# plot(rich~constr, data=dat)
# abline(lmcon)
# lmuncon <- lm(rich~unconstr, data=dat)
# summary(lmuncon)
# plot(rich~unconstr, data=dat)
# abline(lmuncon)

# LCBD

# gof2spa <- inertcomp(rdaa, display = "species", unity = T)
# gof2sia <- inertcomp(rdaa, display = "sites", unity = T)
# 
# dat <- cbind(gof2sia, stages)
# ggplot(dat, aes(x=succession, y=CCA)) +
#   geom_point()
# ggplot(dat, aes(x=succession, y=CA)) +
#   geom_point() +
#   stat_summary(fun.data=mean_cl_boot,
#                geom="pointrange", lwd=0.8) +
#   stat_summary(fun.y=mean, geom="point",cex = 2) +
#   stat_summary(fun.y=mean, geom="line",lwd=1, lty=2)
# 
# gof2sp[rownames(gof2sp) == nam]
# gof2si[1:10,]
# 
# 
# 
# NOTES
# nam <- selected_gof[2]
# nam <- "tet_sal"
# dsts <- cbind(stages, ascnt[,colnames(ascnt) == nam])
# colnames(dsts) <- c("site","succession",nam)
# pt <- ggplot(dsts, aes(y = log(dsts[,3]+1),
#                  x = as.factor(succession),
#                  group = 1)) +
#   geom_jitter(width=0.1)+
#   stat_summary(fun.data=mean_cl_boot,
#                geom="pointrange", lwd=0.8)# +
# 
# pt + stat_summary(fun.y=mean, geom="point",cex = 2)
# pt +  stat_summary(fun.y=mean, geom="line",lwd=1, lty=2)+
#   theme_bw()
# 
# 
# stat_summary(fun.data=mean_cl_boot, 
#              geom="pointrange", lwd=0.8) +
#   stat_summary(fun.y=mean, geom="point",cex = 2) +
#   stat_summary(fun.y=mean, geom="line",lwd=1, lty=2)+
#   stat_summary(fun.y=mean, geom="text", 
#                col = rgb(10,10,10,180,maxColorValue = 255),
#                hjust = 1.2,
#                vjust = -1.5)


# # Which group had stronger relationships
# sum(gof2sp[groups == "Spheciformes","CCA"])
# sum(gof2sp[groups == "Apiformes","CCA"])
# sum(gof2sp[groups == "Chrysididae","CCA"])
