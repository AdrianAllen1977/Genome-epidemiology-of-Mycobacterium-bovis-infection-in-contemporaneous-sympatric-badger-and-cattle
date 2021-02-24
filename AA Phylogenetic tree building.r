###Adapted from www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/

####Running phylogenies in R with the ape package.
##First open the libraries for ape, phangorn and seqinr

library(ape)
library(phangorn)
library(seqinr)

###set the working drive
setwd("Desktop/R stuff/TVR genomes")

### Next, import your DNA alignment file of SNPs in the fasta format.

genomes <- read.dna("29Jun17bestsnp.fasta", format="fasta")

###Now convert this alignment data to a “phyDat” object for use in phangorn

genomes_phyDat<- phyDat(genomes, type = "DNA", levels= NULL)

###Now, make a distance metric between all samples - you need to test which model of nucleotide substitution / evolution best fits the data. So, do a model test.

mt<-modelTest(genomes_phyDat)

### Print out the results of the test
print(mt)

###Taking the model that fits best, now create the distance matrix. IN this case, it’s the F81 model - lowest AIC, BIC and log likelihood and useable for nucleotides. And better than Jukes Cantor.

dna_dist <- dist.ml(genomes_phyDat, model="F81")

###Now, plot your trees using UPGMA or neighbour joining methods.

genomes_UPGMA <-upgma(dna_dist)
genomes_NJ<-NJ(dna_dist)

plot(genomes_UPGMA, main="UPGMA", cex=0.3)
plot(genomes_NJ, main="Neighbour Joining", cex=0.3)

###Which of these two tree methods is best? Use the parsimony function to check. The lowest score is the best tree.

parsimony(tree=genomes_UPGMA, data=genomes_phyDat)
parsimony(tree=genomes_NJ, data=genomes_phyDat)

###Next, export the tree so you can view it in other programmes.  You can export it in the Newick format which works in figtree.  Newick format is like the Nexus format - to do a conversion between textsummarjes of trees, go to this website http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi

write.tree(phy=genomes_NJ, file="TVR genomes_phylo")