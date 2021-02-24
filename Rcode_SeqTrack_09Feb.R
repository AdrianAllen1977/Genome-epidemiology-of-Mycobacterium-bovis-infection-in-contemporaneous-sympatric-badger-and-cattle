## Reconstruction of transmission trees in an endemic phase of 
## Mycobacterium bovis infection in cattle and badgers in NI

## 1.) PREPARATIONS

## 1.1) load packages
library(ape)
library(adegenet)
library(igraph)
library(tcltk)

## 1.2) Set your working directory on your computer 
## Import note: Check that you use slashes / (and not backslash \) 
setwd("C:/Users/Assel/Desktop/ASEKE/wgs TVR Aug 2020")

## 1.3) LOAD DATA
## First look at data
cases <- read.csv("TVR302_input_SeqTrack_v1.csv")

#Add "-01-01" (date January 1st) to all of the isolates
cases$Date <- paste(cases$Year,"-01-01", sep = "")

#Order by date
cases <- cases[order(cases$Date),]

## NOte: the samples ("cases") are sorted by sampling date. '1 is the first/oldest case,
## 1 is the last/most recent case. 

## Define the columms with dates as GMT dates in POSIXct format.
SampleDate <- as.POSIXct(cases$Date, tz = "GMT")       # time zone (tz) in Universal Time, Coordinated (UTC), which is GMT


## In this work we will select the sampling date as the infection date.
Date <- SampleDate
head(Date)
range(Date)

#days=how many days of the infection
days <- as.integer(difftime(Date, min(Date), unit="days"))
days


## Load the SNP sequence data (ordered from oldest to newest)
dna <- fasta2DNAbin("TVR_263clade_n=302_subset_Idonly_06102020.fasta")

#######################################################################################

## 2) EXPLORE SNP DATA
## To have an idea of the existing diversity in these sequences, we compute the
## simple pair-wise Hamming distances and plot their distribution:
D <- dist.dna(dna, model="N")
hist(D, col="grey", nclass=30,
     main="Distribution of pairwise genetic distances",
     xlab="Number of differing nucleotides")

## check if this remarkable polymorphism is distributed randomly across the genome
snps <- seg.sites(dna)
head(snps)
length(snps)

## There are 323 polymorphic sites in the sample. We can visualize their position, and
## try to detect hotspots of polymorphism by computing the density of SNPs as we
## move along the genome
snpposi.plot(dna, codon=FALSE)

#########################################################################################################################

## 3.) SEQTRACK ANALYSIS

## 3.1) BAsic analysis based on genomic data (without including epidemiological information or disease biology) 
## Here, we use SeqTrack on the matrix of pairwise distances (distmat), indicating the labels of the cases
## (x.names=cases$label), the collection dates (x.dates=dates), the mutation rate per site per day (mu_MAP) and 
## the number of nucleotides = length of analysed sequences (nbNucl):
distmat <- as.matrix(D)
nbNucl <- ncol(as.matrix(dna))

#Substitution rate for M.bovis
mu_mbovis <- 8.054e-8/365    ## 0.35 substitutions per genome per year (from paper section 3.5)
                              ## 4,345,492-bp M bovis genome

#seqtrack algorithm
res <- seqTrack(distmat, x.names=cases$Id, x.dates=Date, mu=mu_mbovis, haplo.le=nbNucl)
class(res)
res

## The result res is a data.frame with the special class seqTrack, containing the
## following information:
## res$id: the indices of the cases.
## res$ances: the indices of the putative ancestors of the cases.
## res$weight: the number of mutations corresponding to the ancestries.
## res$date: the collection dates of the cases.
## res$ances.date: the collection dates of the putative ancestors.

#Assign hosts to the SeqTrack output data
res$Hosts <- cases$Host

#Create badgers dataset
res_badgers <- subset(res, Hosts == "B")

#Identify ancestor hosts
res$ances.hosts <- ifelse(res$ances %in% res_badgers$id, 'B', 'C')

## seqTrack objects can be plotted simply using the "plot" command.
## Each sequence/case is a node of the graph, and arrows model putative ances-
## tries/transmissions. The number of mutations between ancestors and descendents
## are indicated by the color of the arrows (red = no/few mutations; light grey= many
## mutations) and the numbers in blue. Time is represented on the y axis (up: ancient;
## down: recent). 

res_network <- plot(res, main="SeqTrack reconstruction of the transmission tree")
mtext(side=3, text="dark red: no/few mutations; light grey: many mutations")
res_network

#Rename res_network Ids for Hosts
V(res_network)$name <- res$Hosts

#Create interactive plot for SeqTrack tree
tkplot(res_network, vertex.color="yellow", layout=layout_with_lgl, 
       vertex.size=7, vertex.label.cex = 0.8,
       edge.label.cex = 0.7)

## Store result files ("res") on your computer  
getwd()                                ## this command will give you the path of the default working directory
write.csv(res, "SeqTrack_outputfile_Mbovis network_10Feb_final.csv")  ## this will save the "res" object in your working directory

#Examine SeqTrack output
#Create barplot with distribution of mutations between ancestors and descendants (Supplementary Figure S10)
range(res$weight, na.rm=TRUE)
barplot(table(res$weight)/sum(!is.na(res$weight)), ylab="Proportions",
        xlab="Mutations between inferred ancestor and descendent", col="grey",
        cex.lab = 1.3)

#Create transmissions plot -> change hosts names for Badger, Cattle
res$Hosts <- ifelse(res$Hosts == "B", 'Badger', 'Cattle')
res$ances.hosts <- ifelse(res$ances.hosts == "B", 'Badger', 'Cattle')
res$transmissions <- paste(res$ances.hosts, res$Hosts, sep = "_")
#Create barplot with transmission events (Supplementary Figure S11)
barplot(table(res$transmissions)/length(res$transmissions)*100, col="lightblue", 
        xlab = "Transmission from/to", ylab = "Number (%) from total",
        cex.names = 0.75,
        cex.lab = 1.3)
