###Making a pairwise distance matrix between multiple GPS points.

### First, make a 2 column matrix / table in Excel of the GPS co-ordinates of interest.  These MUST be Lon Lat decimal data.  Use the other R script to convert XY co-ords to this format.

###Set working drive
setwd("Desktop/R stuff/TVR genomes")

### Open geosphere package.
library("geosphere")

###Import a 2 column table (lon lat)
x<-read.table("TVR_pos_badger_locations.txt")

### Run distm function.
a<-distm(x, fun=distGeo)

###Write matrix to table
write.table(a, file="TVR_pos_badger_map_distance_matrix.txt")
