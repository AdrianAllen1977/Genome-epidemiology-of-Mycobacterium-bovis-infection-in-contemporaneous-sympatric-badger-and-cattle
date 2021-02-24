### GIS in R

##Load libraries

library(rgdal)
library(sp)
library(raster)

## First, import your shape files for the landscape of interest.
### A layer is defined by more than one file in ARC GiS, so to drag them all into R, in one consolidated file, you need to use the ‘readOGR’ function defining that you want all files that share the same prefix to be bundled together.

TVRbuff<-readOGR(dsn=".", layer="2km_Buffer")

TVRzone<-readOGR(dsn=".", layer="TVR")

Irelandshp<-readOGR(dsn=".", layer="Ireland")

plot(TVRbuff)

plot(TVRzone)

plot(Irelandshp)

##Check what projection the shape file is in.  All layers MUST be in the same projection, and so must any points you want to plot on the layers.

proj4string(Irelandshp)

##The Ireland projection is different than the two other layers, so you can alter it.
## Simply determine the projection of the other two layers using the ‘proj4string’ function, copy it, and paste it into the function below in the ‘CRS’ command.

Irelandshp2<-spTransform(Irelandshp, CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1.000035 +x_0=200000 +y_0=250000 +datum=ire65 +units=m +no_defs +ellps=mod_airy +towgs84=482.530,-130.596,564.557,-1.042,-0.214,-0.631,8.15"))

### You can now layer multiple layers on top of one another, provided they are in the same projection. Use the command below:

plot(Irelandshp2 + TVRzone + TVRbuff)

### to add colours to layers

plot(Irelandshp2 + TVRzone + TVRbuff, col=c("green", "yellow", "red"))

### Importing data points - create a text file with the X and Y coords in two columns. Then import them thus.

locs<-read.table("Point locations/All TVR locs.txt", header=T)

## Specifically read in the two columns of interest.

xy<-locs[,c(2,3)]

### Now, make a spatial data frame of these point locations.

xy2<-SpatialPointsDataFrame(coords = xy, data=locs, proj4string=CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1.000035 +x_0=200000 +y_0=250000 +datum=ire65 +units=m +no_defs +ellps=mod_airy +towgs84=482.530,-130.596,564.557,-1.042,-0.214,-0.631,8.15"))

### Now, add the points to the layer

points(xy2, pch=16, col="red")

## You can add in multiple layers of points by just creating a separate points data file for each variable.

## NOTE - the spTransform function only works on layers, not on point data.  You'll have to move between Irish Grid and WGS84 data using the other script.

### add a scale bar
scalebar(100, xy=click(), type='bar', divs=2, below="kilometres")
