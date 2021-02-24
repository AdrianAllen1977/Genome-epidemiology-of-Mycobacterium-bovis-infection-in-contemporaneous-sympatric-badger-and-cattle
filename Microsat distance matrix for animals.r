#### Creating a pairwise genetic distance matrix between microsat genotyped animals.

## First, install packages PopGenReport, adegenet and all dependencies. Then open all packages
library(adegenet)
library(PopGenReport)

### Set working drive
setwd("Desktop/R stuff/TVR genomes")

### Import genepop file (importing converts it to a genind object) ncode = no of numerals in each genotype of an allele.
a<-read.genepop("TVR pos.gen", ncode=3)

## Alternatively, you may already have imported Excel data to genind format via another script.

### Now perform the pairwise difference analysis using the Smouse and Peakall method.
b<-gd.smouse(a, verbose=TRUE)

### Now write the resulting file to a matrix.
library(MASS)
write.matrix(b, file="TVR_263_pos_badger_dist_matrix.txt")
