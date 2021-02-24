#### Determining pairwise geographic distance between pairs of GPS points.
### First, ensure all GPS coords are converted to decimal lat and lon using the appropriate R script.
### Create two files - one for the initial location for each sample, the second for the follow up location for the sample.

###set working drive
setwd("Desktop/R stuff/TVR genomes/recapture")

### Open geosphere package
library("geosphere")

###Import two data files
p1<-read.table("2014 TVR coords.txt")
p2<-read.table("2015 TVR coords.txt")

###Run analysis
x<-distGeo(p1, p2, a=6378137, f=1/298.257223563)

###Write distances between points to table
write.table(x, file="TVR badgers 2014-2015 recapture distances.txt")

###Plot histogram of results
hist(x, xlim=range(0:7000), ylim=range(0:100), main="Histogram of TVR  recaptured badger movements 2014-2015 ", xlab="Distance (metres)")
