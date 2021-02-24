#### Simple Phylogeny Script for Figures.

#### PACKAGES ####
library(ape)
library(phangorn)
library(phytools)

# Read in the FASTA file
sequences<- read.dna("/yourfilename",
                            format = "fasta")

# Note the number of sites in the FASTA file
nSitesInFasta <- length(as.list(sequences)[[1]])

#### CHOOSE SUBSTITUTION MODEL ####

# Convert the sequences from a DNAbin format into a PhyDat format
sequencesPhyDat <- phyDat(sequences, type = "DNA", levels = NULL)

# Run model testing to select appropriate substitution model
modelTestResults <- modelTest(sequencesPhyDat, model = c("JC", "HKY", "GTR"))

# Get best model
cat(paste("Best substitution model:", 
          modelTestResults$Model[which.min(modelTestResults$AIC)], "\n"))


#### BUILD PHYLOGENY ####

# Build the distance matrix
distanceMatrix <- dist.dna(sequences, model="JC69")

# Build a neighbour joining tree - an initial tree
njTree <- nj(distanceMatrix)

# Compute likelihood of the initial Neighbour Joining tree given sequences
likelihoodObject <- pml(njTree, sequencesPhyDat)

# Set the controls for the Maximum Likelihood algorithm
controls <- pml.control(maxit=100000, trace=0)

# Run maximum likelihood
fittingOutput <- optim.pml(likelihoodObject,
                           optNni = TRUE, # Optimise topology
                           optInv = TRUE, # Optimise proportion of variable sites
                           model = "GTR", # Substitution model
                           rearrangement="NNI", # Nearest Neighbour Interchanges
                           control=controls)

# Bootstrap the result of maximum likelihood
bootstrapResults <- bootstrap.pml(
  fittingOutput, # Use Maximium Likelihood settings on bootstrapped sequences
  bs = 100, # Number times to bootstrap sequences
  optNni = TRUE, # Use Nearest Neighbour Interchanges in tree building
  jumble=TRUE) # Jumble bootstrapped sequences before building trees

# Get phylogenetic tree with bootstrap values - plotBS will return a tree that can be used
tree <- plotBS(
  fittingOutput$tree,
  bootstrapResults,
  p = 50, # Plot bootstrap values if node in >=50 bootstrap trees
  type="phylogram", cex=0.1) # Type of phylogenetic tree shape to plot

### To ensure you can get a SNP difference scale bar on subsequent phylogenies.
### You need to alter the branch lengths to SNPs and not proportions.
tree$edge.length<-tree$edge.length*nSitesInFasta

###If you want to write the tree to the Newick format that can be read by Figtree, do the following.
write.tree(phy=tree, file="TVR phylogeny Nov19 - n=303 single isolates all years")

### Annotating tree with tip shapes and colours

# Remove any non-dichotomous branching nodes (nodes with more than two daughters)
rootedTree <- multi2di(tree, random=TRUE)

rootedTree <- root(rootedTree, 
                   outgroup="TVR014_Badger_2015", # Set the isolate to root the tree on - can also use node
                   resolve.root=TRUE) # Ensures new root will be bifurcating
                   
# Make branches with lengths of 0 very very small
rootedTree$edge.length[rootedTree$edge.length <= 0] <- 0.000001



### Set plotting margins wide for phyogeny resolution
currentMar <- par()$mar
par(mar=c(0, 0, 0, 0))


### Plot rooted phylogenetic tree
plot.phylo(rootedTree, cex=0.1)


# Get tip states from their names
tipStates <- c()
for(index in 1:length(rootedTree$tip.label)){
  tipStates[index] <- strsplit(rootedTree$tip.label[index], split="_")[[1]][2]
}

### Give the refernce tip a state
### The lead state is the "Reference" genome - it has been entered as an NA because the scripts above were only pulling out the label badger or bovine.
### So, relabel this using the following script.
tipStates[1]<-"Bovine"

rootedTreeWithoutRef <- drop.tip(rootedTree, tip="Reference")
tipStates <- tipStates[-1]



# Assign colours to the states
stateColours <- setNames(c("red", "blue"), c("Badger", "Bovine"))

### Plot phylogenetic tree in fan format
plot.phylo(rootedTreeWithoutRef, type="fan", # Set the shape of the tree to be a fan
           show.tip.label=FALSE) # Don't plot the tip labels

# Add circles to indicate the states of the tips
tiplabels(pie=to.matrix(tipStates,c("Badger", "Bovine")), # Provide the state of each tip
          piecol=stateColours, # Provide the colour assigned to each state
          cex=0.1) # Set the size of the tip circle

### Add scale bar
add.scale.bar(-20, 0, length=5, lcol="black", lwd=1, cex=0.5)

# Add a legend
legend("bottomright",
       legend=c("Badger", "Bovine"), # Name the labels used in the legend
       pch=19, # Set the shape of the points used in legend
       col=c("red", "blue"), # Set the colour of each label point
       bty="n") # Remove box surrounding legend
                   
