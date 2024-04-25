###Script for downloading Climate Data from World Clim for Use in GWAS###
##Author: Patrick Woods##
##Date: April 25, 2024##

setwd("Desktop/Projects/eGWAS/")
install.packages("raster")

#install.packages(file.choose(), repos=NULL, type='source') #Needed to download and manually install terra_1.7-71.tar
library(raster)

#need a file with coordinates
coords <- read.csv("meta_for_patrick.csv",header=T) #Need a .csv file with longitude and latitude coordinates for individual samples

#Subset file for 'Wild' status accessions only
coords_wild <- subset(coords, coords$Status == 'Wild')

#import all climate data
r <- getData('worldclim',var='bio', res=10) #resolution set to 10 arcminutes to account for uncertainty in geographic locations

#Extract climate data for coordinates of wild accessions
clim1 <- coords_wild$lon
clim2 <- coords_wild$lat
clim1x <- c(clim1)
clim2y <- c(clim2)
df <- data.frame(x = clim1x, y = clim2y)
points <- SpatialPoints(df, proj4string = r@crs)
extrct <- extract(r, points)
cods <- cbind.data.frame(coordinates(points), extrct)

cods <- cbind(coords_wild$Genotype, cods)

write.csv(final_clim_data, file = "final_climate_data_resolution_10_trim.csv")
