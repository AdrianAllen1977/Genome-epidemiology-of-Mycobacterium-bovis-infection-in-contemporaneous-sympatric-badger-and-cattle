#### Generic Script for Building Distance matrix for times
setwd("Desktop/R stuff/TVR genomes/Red dog outputs/Badger genetic data")
library(MASS)

### Import table with sample names and dates (This works with just years - precise dates require different approach)

x<-read.table("263_Pos_Badger_year.txt", header=T)

### Make a vector of the column of dates

y<-x$Year

### Create a distance object that compares all pairs of dates - include the abs function as this records differences between 2014-2015 and 2015-2014

z<-dist(abs(y))

## Make the distance output a matrix and write to file
matrix<-as.matrix(z)
write.matrix(matrix, "263_Pos_Badger_time_matrix.txt")
