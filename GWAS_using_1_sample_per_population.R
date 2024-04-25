### Script for categorizing lon/lat data into factors and sampling 1 plant from each factor ###
#Author Patrick Woods #
# Date: April 19, 2024 #


# Read in the raw metadata file supplied by UMN and subset for only genotyeps that are considered 'Wild' #
cords <- read.csv('meta_for_patrick.csv',header=T)
wild <- subset(cords, cords$Status == 'Wild')

# Make a new column that combines every longitude + latitude coordinate into a single column #
wild$cat <- paste(wild$lon, wild$lat, sep="_")

wild$cat <- as.factor(wild$cat)

str(wild) # 249 unique lon/lat coordinate combinations.

# Randomly select 1 sample from each unique longitude_latitude combination #

library(dplyr)

wild_sub <- wild %>%
  group_by(cat) %>%
  slice_sample(n = 1)

# Subset the extracted bioclimate variable data for the subset of 249 genotypes #

phe <- read.csv('final_climate_data_resolution_10_trim.csv',header=T)

wild_sub_phe <- left_join(wild_sub, phe, by = 'Genotype')
wild_sub_phe <- as.data.frame(wild_sub_phe)
write.csv(wild_sub_phe, file = 'wild_sub_phe.csv', row.names = F, quote = F)

# Test each of the 19 bioclimate variables for normality amongst the 249 gneotype subset #
p_vals <- data.frame()

traits <- 10:ncol(wild_sub_phe)

for (i in traits){
  shap <- shapiro.test(wild_sub_phe[,i])
  p_vals <- rbind(p_vals, shap$p.value)
}

max(p_vals) # greatest p value is still less than 0.05. Will quantile normalize the traits.

# all traits show significant departure from normality, will use quantile normalization to transform #

trans <- data.frame()

library(bestNormalize)
library(rowr)

for (i in traits) {
  trait_qn <- orderNorm(wild_sub_phe[,i])
  df <- as.data.frame(trait_qn$x.t)
  colnames(df) <- colnames(wild_sub_phe[i])
  trans<- cbind.fill(trans,df,fill = NA)
}

p_vals_trans <- data.frame()

traits <- 2:ncol(trans)

for (i in traits){
  shap <- shapiro.test(trans[,i])
  p_vals_trans <- rbind(p_vals_trans, shap$p.value)
}

max(p_vals_trans)

write.csv(trans, file = 'bioclim_qn_unique_locations.csv', row.names = F, quote=F)

# Convert the .csv file to a .txt file using excel #

### load GAPIT functions ###
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

myGD=read.delim(file="wild_sub_maf05.hmp.txt",header=FALSE,sep="\t") # VCF converted to HapMap format, subsetted to the 249 wild genotypes and filtered for 
a minor allele frequency of 5%
myY <- read.table('bioclim_qn_unique_locations_final_GAPIT.txt',header=T)
#install.packages("bigmemory",repos = "http://cran.us.r-project.org")
#install.packages("biganalytics",repos = "http://cran.us.r-project.org")
#install.packages("devtools",repos = "http://cran.us.r-project.org")
#install.packages("scatterplot3d",repos = "http://cran.us.r-project.org")

### GWAS using BLINK model and 5 PCs ###
myGAPIT <- GAPIT(
  Y=myY[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)],
  G=myGD,
  PCA.total=5,
  model=c("Blink"), SNP.MAF = 0.05)
