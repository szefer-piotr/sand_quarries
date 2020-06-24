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
# api <- read.xlsx("data/abundance_table.xls", 1)
# sph <- read.xlsx("data/abundance_table.xls", 2)
# chr <- read.xlsx("data/abundance_table.xls", 3)

clean_dat <- function(dat, clean.row = TRUE){
  ind <- grep("NA", colnames(dat))
  if(length(ind) != 0){ccdat <- dat[, -ind]} else {ccdat <- dat}
  if(clean.row){
    ccdat <- ccdat[complete.cases(ccdat), ]
  }
  return(ccdat)
}

# Read the new dataset
apir <- read.xlsx("data/life_history.xls", "herbivore", startRow = 3)
sphr <- read.xlsx("data/life_history.xls", "predator", startRow = 3)
chrr <- read.xlsx("data/life_history.xls", "cleptoparasite", startRow = 3)

api <- clean_dat(apir)
sph <- clean_dat(sphr)
chr <- clean_dat(chrr)

# Extract successional stages ----
stages <- data.frame(site = 1:32, succession = api[,2])

rownames(api) <- api[,1]
# api <- api[, 3:192]
# colnames(api)

# Check for duplicates
colnames(api)[duplicated(colnames(api))]

rownames(sph) <- sph[,1]
# sph <- sph[, 3:74]
colnames(sph)

# Check for duplicates
colnames(sph)[duplicated(colnames(sph))]

# Deal with duplicate names
# ‘ect_con’ #real error, 
# ‘oxy_qua’ # real duplicate concatenate these rows

# convert names
# for(i in 1:dim(sph)[2]){
#   teststr <- as.character(colnames(sph[i]))
#   teststr <- tolower(teststr)
#   nm <- strsplit(teststr,"_")
#   new_name <- paste(substr(nm[[1]][1], 1,3),
#                     substr(nm[[1]][2], 1,3), sep="_")
#   if(new_name %in% colnames(sph)){
#     new_name <- paste(substr(nm[[1]][1], 1,3),
#                       substr(nm[[1]][2], 2,4), sep="_")
#   }
#   colnames(sph)[i] <- new_name
# }

rownames(chr) <- chr[,1]
# chr <- chr[, 3:13]
colnames(chr)

# # convert names
# for(i in 1:dim(chr)[2]){
#   teststr <- as.character(colnames(chr[i]))
#   teststr <- tolower(teststr)
#   nm <- strsplit(teststr,"_")
#   new_name <- paste(substr(nm[[1]][1], 1,3),
#                     substr(nm[[1]][2], 1,3), sep="_") 
#   colnames(chr)[i] <- new_name
# }


api <- api[, -grep("site|succession", colnames(api))]
sph <- sph[, -grep("site|succession", colnames(sph))]
chr <- chr[, -grep("site|succession", colnames(chr))]

# Create a combined table ----
asc <- cbind(api, sph, chr)
# Indicate group name
groups <- rep(c("Herbivores", "Predators","Kleptoparasites"),
              times=c(dim(api)[2], dim(sph)[2], dim(chr)[2]))
group_col <- rep(c("red", "orange","black"),
                 times=c(dim(api)[2], dim(sph)[2], dim(chr)[2]))

# Insect community descriptors - data frame ----
# Full data set.

desasc <- rbind(stages,stages,stages)
desasc$group <- rep(c("Herbivores", "Predators","Kleptoparasites"), 
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

stage1 <- list("Herbivores" = as.numeric(sort(colSums(api[s1,]),
                                                   decreasing = T)),
               "Predators" = as.numeric(sort(colSums(sph[s1,]),
                                                      decreasing = T)),
               "Kleptoparasites" = as.numeric(sort(colSums(chr[s1,]),
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

dfa <- data.frame(group = rep("Herbivores", dim(api)[2]),
                  stage = rep(c("SI","SII","SIII"), 
                              each=dim(api)[2]))
t1 <- rac_mxa[1,]
t2 <- rac_mxa[2,]
t3 <- rac_mxa[3,]

dfa <- cbind(dfa, abund = c(t1,t2,t3))
dfa$spec <- c(names(t1),names(t2),names(t3))

# Rare species data preparation
raredatr <- read.xlsx("data/life_history.xls", "rare_species", startRow = 2)
raredat <- clean_dat(raredatr, clean.row = FALSE)
names(raredat)[1:3] <- c("status", "lifeh", "spec")

matrix.to.data2 <- function(mat, col.shift=0){
  # This function takes matrix and converts each entry into a 
  # data frame with rows and column naes and entry values.
  drow <- 1
  data <- data.frame(row=0,col=0,val=0)
  for (irow in rownames(mat)){
    for (icol in colnames(mat)){
      if (mat[irow,icol] != 0) {
        #print(paste(irow,icol,mat[irow,icol]))
        data[drow, ]$row = irow
        data[drow, ]$col = icol
        data[drow, ]$val = as.numeric(mat[irow,icol])
        drow <- drow+1
      }
    }
  }
  return(data)
}


raredf = data.frame()
lgth = dim(raredat)[2]
rw = 1
for(rw in 1:dim(raredat)[1]){
  rowdes = raredat[rw,1:3]
  rowsp=raredat[rw,4:lgth]
  specs <- colnames(rowsp)[which(rowsp != 0)]
  specno = rowsp[which(rowsp != 0)]
  rowdes
  rowsp
  specs
  specno
  dfunit = data.frame(site=rep(specs, specno),
                      status=rowdes$status,
                      lifeh= rowdes$lifeh,
                      spec = rowdes$spec)
  raredf <- rbind(raredf,dfunit)
}
raredf
raredf$vuln <- as.numeric(!is.na(raredf$status))
raredf$cr <- 0
raredf$cr[grep("CR", raredf$status)] <- 1

raredf$rm <- 0
raredf$rm[grep("rm|vrm", raredf$status)] <- 1

raredf$stage <- 1
raredf$stage[grep("X2", raredf$site)] <- 2
raredf$stage[grep("X3", raredf$site)] <- 3
raredf$sno <- 1
raredf$sno[grep(".1", raredf$site, fixed = T)] <- 2
raredf$sno[grep(".2", raredf$site, fixed = T)] <- 3
raredf$sno[grep(".3", raredf$site, fixed = T)] <- 4
raredf$sno[grep(".4", raredf$site, fixed = T)] <- 5
raredf$sno[grep(".5", raredf$site, fixed = T)] <- 6
raredf$sno[grep(".6", raredf$site, fixed = T)] <- 7
raredf$sno[grep(".7", raredf$site, fixed = T)] <- 8
raredf$sno[grep(".8", raredf$site, fixed = T)] <- 9
raredf$sno[grep(".9", raredf$site, fixed = T)] <- 10

raredf[raredf$lifeh == "cleptoparasite",]$lifeh <- "kleptoparasite"
raredf$lifeh <- as.factor(as.character(raredf$lifeh))
