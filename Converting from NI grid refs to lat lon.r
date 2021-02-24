## Converting from NI grid refs to lat lon
## ==========================================================================================
## Only need this the very first time to get the packages installed:
install.packages("sp","maptools","rgeos","rgdal") 

###Open libraries

library(sp)
library(maptools)
library(rgeos)
library(rgdal)

## Set the working directory so the computer knows where to look for your data
setwd("Desktop/R Stuff/Irish badgers")   ## (Change to indicate the folder where your data are)

## This code is useful to search for specific map projections: 
EPSG <- make_EPSG()      ## Creates object containing all possible projections
EPSG[grep("Irish", EPSG$note), ]   ## Search for projections associated with specific word, here "Irish"

## Set the relevant map projections
latlon_CRS <- CRS("+proj=longlat +ellps=WGS84")
Irish_CRS  <- CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1.000035 +x_0=200000 +y_0=250000 
                   +datum=ire65 +units=m +no_defs")

## Import your data on the locations: change this to the name of your datafile (make sure it's in .csv format)
data <- read.csv("Badger coords.csv")
## Change this to a "spatial points" object 
data_sp <- data
coordinates(data_sp) <- cbind(data$X, data$Y)  ## (assumes location columns are labelled "X" and "Y")
## Tell it that the locations are as Irish grid references
proj4string(data_sp) <- Irish_CRS

## Transform Irish grid refs to lat lon
data_latlon <- spTransform(data_sp, latlon_CRS)
## Extract the transformed coordinates
llcoords <- coordinates(data_latlon)
## Name the columns of the lat lon coords
colnames(llcoords) <- c("X_transformed","Y_transformed") 

## Add the new coords as two extra columns in your original data table (cbind=bind columns to a dataframe)
data <- cbind(data, llcoords)
## This will save the results as a .csv file in the folder you specified at the top:
write.csv(data, "data_plus_latlon.csv", row.names=F, quote=F)

