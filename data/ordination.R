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


# 2. Variance explained -----
sumeig <- summary(eigenvals(m))
round(sumeig[2,1],3)
round(sumeig[2,2],3)
# >>> goodness of fit ----

# Fit species onto ordination
set.seed(1234)
gps <- data.frame(Species = colnames(asc), Group = groups)
spc <- envfit(m, asc, permutations = 9999)

spcdf <- data.frame(Spec = names(spc$vectors$pvals),
                    R = spc$vectors$r,
                    Pvals = spc$vectors$pvals,
                    gps[names(spc$vectors$pvals), ]$Group)

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


# pdf("fig2_rda.pdf", width = 12, height = 4)
par(mfrow = c(1,3))
make_rda_plot(m, 3, "Herbivores", 1)
make_rda_plot(m, 3, "Predators", 3)
make_rda_plot(m, 3, "Kleptoparasites", 2)
# dev.off()


# 5. Contingency table for envfit ----
# https://www.datacamp.com/community/tutorials/contingency-tables-r

efs <- names(spc$vectors$pvals[spc$vectors$pvals <= 0.05])
rownames(gps) <- gps[,1]
gps$Significance <- "NS"
# gps[!(rownames(gps) %in% efs), ]$ns <- "NS"
gps[efs,]$Significance <- "SI" 

Ag.3.5raw <- table(gps$Group,  gps$Significance)
Ag.3.5 <- as.table(matrix(Ag.3.5raw, nrow = 3, byrow = FALSE, dimnames = list(Group = c('Herb', 'Klep', 'Pred'),Significance = c('NS', 'SI'))))
Ag.3.5
addmargins(Ag.3.5)

# Independence between Stage and Group
library(MASS)
Ag.3.5.loglm <- loglm( ~ Significance + Group, data = Ag.3.5) 
Ag.3.5.loglm

# Two columns Early and Late that have similar proportions
Ag.3.5.prop.mar.1.cols.12.table <- prop.table(Ag.3.5[,1:2], margin = 1)
# transform to percentages
Ag.3.5.percent.mar.1.cols.12.table <- 100*Ag.3.5.prop.mar.1.cols.12.table

# Herbivores and kleptoparasites have similar proportions of significant to nonsignificant species.

Ag.3.6.A <- Ag.3.5[1:2, ]
Ag.3.6.A
Ag.3.6.A.loglm <- loglm( ~ Significance + Group, data = Ag.3.6.A)
Ag.3.6.A.loglm
# Homogenous matrix

margin.table(Ag.3.6.A, margin = 2)
Ag.3.6.B <- as.table(matrix(c(margin.table(Ag.3.6.A, margin = 2),
                              Ag.3.5[3,]),
                            nrow = 2,
                            byrow = TRUE,
                            dimnames = list(Group = c('Herb+Klep','Pred'),
                                            Significance = c('NS', 'SI')
                                            )))
Ag.3.6.B
Ag.3.6.B.loglm <- loglm( ~  Significance + Group, data = Ag.3.6.B)
Ag.3.6.B.loglm
# H0 rejected - nonhomogenous

# Test for 12 decimal places
round(Ag.3.5.loglm$lrt,12) == round((Ag.3.6.A.loglm$lrt+Ag.3.6.B.loglm$lrt),12)

Ag.3.6.A
Ag.3.6.B

# Two columns Early and Late that have similar proportions
Ag.3.6.B.prop.mar.1.cols.12.table <- prop.table(Ag.3.6.B[,1:2], margin = 1)
# transform to percentages
Ag.3.6.B.percent.mar.1.cols.12.table <- 100*Ag.3.6.B.prop.mar.1.cols.12.table
