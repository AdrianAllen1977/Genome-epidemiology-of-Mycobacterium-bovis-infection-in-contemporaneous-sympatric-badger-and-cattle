### ggtree for phylogeny annotation
### Installing ggtree

BiocManager::install("ggtree")

## Open libraries

library(ggtree)
library(ape)
library(Biostrings)
library(ggplot2)

### Read in tree
### ggtree accepts multiple format trees newick, RaxML, Beast etc etc

tree<-read.tree(yourfilename)

### Check branch length totals and correct if necessary to actual SNP lengths.
## L is the number of sites from the fasta that made the tree:

sum(tree$edge.length)
tree$edge.length<-tree$edge.length*L


### Determine node numbers so you can highlight clades and root the tree etc

ggtree(tree) + geom_tippoint() + geom_text(aes(label=node))

### Root tree to specific out-group sequence
# Set the isolate to root the tree on - can also use node

rootedTree <- root(tree, outgroup="tip_label_ID", resolve.root=TRUE)
                   
## OR

rootedTree <- root(tree, node="node_number",resolve.root=TRUE)

### Print off a copy of the rooted phylogeny and check which nodes correspond to MRCA of clades for highlighting.

ggtree(rootedTree) + geom_tippoint() + geom_text(aes(label=node))
_________________________________________________________________________________
#######################################
### Several ways to highlight clades###
#######################################

## METHOD 1 - Adding a coloured line opposite all clade members with a label.

x<-ggtree(rootedTree, size=0.5)

clade_tree<- x + geom_cladelabel(node=739, label="6.263", color=("green"), barsize=2, offset=0.8, align = TRUE)

## Add in more geom_cladelabels as needed.

### METHOD 2 - Adding a coloured block to each clade - defined by node number

x<-ggtree(rootedTree, size=0.5)

clade_tree<- x + geom_hilight(node=739, fill="green")

## Add in more geom_hilight fields as needed.

_________________________________________________________________________________

### ADD A SCALE BAR

tree1<-clade_tree + geom_treescale(x=100, y=100, width=100, linesize=0.5, fontsize=0, color="black")

## x and y help place the bar in exact coordinates in the plot
## width is no of SNPs you want the scale bar to be for
## Linesize is how heavy the scale line will be
## Fontsize - scale bar font size.

_________________________________________________________________________________

### Add tip shapes and colour them by meta data entry.

## First, prepare and read in a tab separated txt file of metadata values beside tip names

## Extract tip names

labels<-tree$tip.label

### Make an empty vector called tipdates
metadata_labels<-c() 

### String split the name of the label - all samples MUST have same string format!
## Make sure fasta used to build tree has sample names in same format.
for(index in 1:length(labels)){metadata_labels[index] <-strsplit(labels[index], split="_")}

### Make the vector into a matrix with columns - make sure no of columns is right!
metadata_labels2<-matrix(unlist(metadata_labels), ncol=3, byrow=TRUE)

## Now extract the column which has the metadata you want in it.
host<-metadata_labels2[,2]

tip_metadata <- read.table("meta_data.txt", sep="\t", header=TRUE,check.names=FALSE, stringsAsFactor=F)

### Now add metadata tips to phylogeny - note meta data will be taken alphabetically, so list preferred colours in order of preference below

tree2<- tree1 %<+% tip_metadata + geom_tippoint(aes(color=host), size=0.05) + scale_color_manual(values=c("red", "blue","grey"))
