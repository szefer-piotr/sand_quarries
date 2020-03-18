# Data processiong script

# Load packages ----
library(xlsx)
library(ggplot2)
library(vegan)

# Color palette
alpha = 200
alpha2 = 100
colvec <- c(rgb(230,159,0,alpha, max=255),
            rgb(86,180,233,alpha, max=255),
            rgb(0,158,115,alpha, max=255),
            rgb(0,114,178,alpha2, max=255),
            rgb(213,94,0,alpha2, max=255),
            rgb(240,228,66,alpha2, max=255))
cgray <- rgb(10,10,10,alpha2, maxColorValue = 255)

# Read the insect table ----
api <- read.xlsx("data/abundance_table.xls", 1)
sph <- read.xlsx("data/abundance_table.xls", 2)
chr <- read.xlsx("data/abundance_table.xls", 3)

# Extract successional stages ----
stages <- data.frame(site = 1:32, succession = api[,2])

rownames(api) <- api[,1]
api <- api[, 3:192]
colnames(api)

# Check for duplicates
colnames(api)[duplicated(colnames(api))]

rownames(sph) <- sph[,1]
sph <- sph[, 3:74]
colnames(sph)

# Check for duplicates
colnames(sph)[duplicated(colnames(sph))]

# Deal with duplicate names
# ‘ect_con’ #real error, 
# ‘oxy_qua’ # real duplicate concatenate these rows

# convert names
for(i in 1:dim(sph)[2]){
  teststr <- as.character(colnames(sph[i]))
  teststr <- tolower(teststr)
  nm <- strsplit(teststr,"_")
  new_name <- paste(substr(nm[[1]][1], 1,3),
                    substr(nm[[1]][2], 1,3), sep="_")
  if(new_name %in% colnames(sph)){
    new_name <- paste(substr(nm[[1]][1], 1,3),
                      substr(nm[[1]][2], 2,4), sep="_")
  }
  colnames(sph)[i] <- new_name
}

rownames(chr) <- chr[,1]
chr <- chr[, 3:13]
colnames(chr)
# convert names
for(i in 1:dim(chr)[2]){
  teststr <- as.character(colnames(chr[i]))
  teststr <- tolower(teststr)
  nm <- strsplit(teststr,"_")
  new_name <- paste(substr(nm[[1]][1], 1,3),
                    substr(nm[[1]][2], 1,3), sep="_") 
  colnames(chr)[i] <- new_name
}


# Create a combined table ----
asc <- cbind(api, sph, chr)
# Indicate group name
groups <- rep(c("Apiformes", "Spheciformes","Chrysididae"),
              times=c(dim(api)[2], dim(sph)[2], dim(chr)[2]))
group_col <- rep(c("red", "orange","black"),
                 times=c(dim(api)[2], dim(sph)[2], dim(chr)[2]))

# Insect community descriptors - data frame ----
# Full data set.

desasc <- rbind(stages,stages,stages)
desasc$group <- rep(c("Apiformes", "Spheciformes","Chrysididae"), 
                     each = dim(stages)[1])
ind <-c("shannon")
desasc$sw <- c(diversity(api, index = ind, MARGIN = 1, base = exp(1)),
               diversity(sph, index = ind, MARGIN = 1, base = exp(1)),
               diversity(chr, index = ind, MARGIN = 1, base = exp(1)))
ind <-c("simpson")
desasc$simp <- c(diversity(api, index = ind, MARGIN = 1, base = exp(1)),
               diversity(sph, index = ind, MARGIN = 1, base = exp(1)),
               diversity(chr, index = ind, MARGIN = 1, base = exp(1)))
desasc$rich <- c(rowSums(api>0),
                 rowSums(sph>0),
                 rowSums(chr>0))
desasc$abun <- c(rowSums(api),
                 rowSums(sph),
                 rowSums(chr))

desasc$succession <- as.factor(desasc$succession)

desasc$s1 <- as.numeric(desasc$succession == 1)
desasc$s2 <- as.numeric(desasc$succession == 2)
desasc$s3 <- as.numeric(desasc$succession == 3)

desasc$GS <- interaction(desasc$succession, desasc$group)

s1 <- stages$succession == 1
s2 <- stages$succession == 2
s3 <- stages$succession == 3

stage1 <- list("Apiformes" = as.numeric(sort(colSums(api[s1,]),
                                                   decreasing = T)),
               "Spheciformes" = as.numeric(sort(colSums(sph[s1,]),
                                                      decreasing = T)),
               "Chrysididae" = as.numeric(sort(colSums(chr[s1,]),
                                                     decreasing = T)))
# stage2 <-
# stage3 <-

rac_df <- expand.grid(stage = c("Stage I", "Stage II", "Stage III"),
            group = unique(groups))
rac_mxa <- rbind(
  
  # Apiformes
  colSums(api[stages$succession == 1, ]),
  colSums(api[stages$succession == 2, ]),
  colSums(api[stages$succession == 3, ])
)  

rac_mxs <- rbind(
  
  # Spheciformes
  colSums(sph[stages$succession == 1, ]),
  colSums(sph[stages$succession == 2, ]),
  colSums(sph[stages$succession == 3, ])
)  

rac_mxc <- rbind(

  # Chrysididae
  colSums(chr[stages$succession == 1, ]),
  colSums(chr[stages$succession == 2, ]),
  colSums(chr[stages$succession == 3, ])
)  
rac_mx <- cbind(rac_mxa,
                rac_mxs,
                rac_mxc)

# df for multivariate_difference codyn, 
# this won't work as we do not have an experimental design

dfa <- data.frame(group = rep("Apiformes", dim(api)[2]),
                  stage = rep(c("SI","SII","SIII"), 
                              each=dim(api)[2]))
t1 <- rac_mxa[1,]
t2 <- rac_mxa[2,]
t3 <- rac_mxa[3,]

dfa <- cbind(dfa, abund = c(t1,t2,t3))
dfa$spec <- c(names(t1),names(t2),names(t3))
