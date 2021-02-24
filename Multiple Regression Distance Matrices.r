#### Multiple Regression of Distance Matrices (MRM)

### Packages
library(ecodist)

## You want to compare multiple sets of distance Matrices
## To see which explanatory variables correlate with the outcome variables
## First import all your matrices

gen<-read.table("TVR pos badger 263 genomes SNP dist matrix.txt", header=T)

### Check dimensions to see it is as expected
dim(gen)

geo<-read.table("TVR_pos_badger_map_distance_matrix.txt", header=T)
dim(geo)

micro<-read.table("TVR_263_pos_badger_microsat_dist_matrix.txt", header=T)
dim(micro)

time<-read.table("263_Pos_Badger_time_matrix.txt", header=T)
dim(time)

### Convert all the imported tables to distance matrices
gend<-as.dist(gen, upper=FALSE)
geod<-as.dist(geo, upper=FALSE)
microd<-as.dist(micro, upper=FALSE)
timed<-as.dist(time, upper=FALSE)

### Now, run MRM.
## The mrank argument must be set to F for the Pearson coefficient for continuous variables.
## Pearson (mrank=F) is the default setting.

x<-MRM(gend~geod+microd+timed, nperm=10000, mrank=F)


## Check the summary data by typing 'x'
