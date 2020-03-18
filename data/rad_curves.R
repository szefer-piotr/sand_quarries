# Rank abundance curves

library(codyn) # Avolio et al. ECOSPHERE

# Detailed species-level community changes. Studying the changes in the shape of RAC curves. RAC and multivariate measures are not sensitive to species richness and evenness and all measures detail uniqe aspects of temporal and spatial differences. *Species reordering* is the strongest correlate of a multivariate measure of compositional change. 

# The problem here is taht we dont have grouped sites. We have only randomly taken sites from stage1, stage II and stage III of succession. WE could make comparisons between all sites but this would result in dependencies between mesured changes (Avolio et al 2019, page 4)

# We can only compare how RAC change in different stages as a single value. Withouth statistical test.

# Vignette
data(pplots)
# Without replicates
df <- subset(pplots, plot == 25)
RAC_change(df = df,
           species.var = "species",
           abundance.var = "relative_cover",
           time.var = "year")
# With replicates
df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
RAC_change(df = df,
           species.var = "species",
           abundance.var = "relative_cover",
           replicate.var = "plot",
           time.var = "year")