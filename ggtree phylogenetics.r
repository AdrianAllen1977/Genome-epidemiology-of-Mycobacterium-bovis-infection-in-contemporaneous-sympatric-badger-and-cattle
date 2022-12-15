### ggtree for phylogeny annotation
### Installing ggtree

BiocManager::install("ggtree")

## Open libraries

library(ggtree)
library(ape)
library(Biostrings)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)
library(tidytree)
library(dplyr)
library(ggstar)

### Read in tree
### ggtree accepts multiple format trees newick, RaxML, Beast etc etc

tree<-read.tree(yourfilename)

### Check branch length totals and correct if necessary to actual SNP lengths.
## L is the number of sites from the fasta that made the tree:

sum(tree$edge.length)
tree$edge.length<-tree$edge.length*L


### Determine node numbers so you can highlight clades and root the tree etc

ggtree(tree) + geom_tippoint() + geom_text(aes(label=node), size=1)

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

### Add labels to colour blocked clades

clade_tree + geom_cladelabel(node=48, label="2.142", align=TRUE, color="blue", offset.text=0.1)

_________________________________________________________________________________

### ADD A SCALE BAR

tree1<-clade_tree + geom_treescale(x=100, y=100, width=100, linesize=0.5, fontsize=0, color="black")

## x and y help place the bar in exact coordinates in the plot
## width is no of SNPs you want the scale bar to be for
## Linesize is how heavy the scale line will be
## Fontsize - scale bar font size.

_________________________________________________________________________________

### Add tip shapes and colour them by meta data entry.

### Add a matrix of meta data for samples that will sit alongside the tree tips

## First, prepare and read in a tab separated txt file of metadata values beside tip names
### This depends on you having sample file names / tip file names which contain metadata fields separated by '_'

## Import your tree

tree<-read.tree(yourfilename)

### Check branch length totals and correct if necessary to actual SNP lengths.
## L is the number of sites from the fasta that made the tree:

sum(tree$edge.length)
tree$edge.length<-tree$edge.length*L

### Determine node numbers so you can highlight clades and root the tree etc

ggtree(tree) + geom_tippoint() + geom_text(size=0.5, aes(label=node), col="blue")

### Root tree to specific out-group sequence
# Set the isolate to root the tree on - can also use node

rootedTree <- root(tree, outgroup="tip_label_ID", resolve.root=TRUE)

## OR

rootedTree <- root(tree, node="node_number",resolve.root=TRUE)


### You may wish to drop the reference tip as it's name isn't in the same format as sample tips.
## And this format messes with the tip name extraction detailed below.

rootedTree_noRef <- drop.tip(rootedTree, tip="Ref")

## Plot the tree in ggtree
p<-ggtree(rootedTree_noRef)

### you can also do a circular tree if you want to save space
#p<-ggtree(rootedTree_noRef, layout='circular')

## Extract tip names

labels<-rootedTree_noRef$tip.label

### Make an empty vector called metadata_labels
metadata_labels<-c() 

### String split the name of the label - all samples MUST have same string format!
## Make sure fasta used to build tree has sample names in same format.
for(index in 1:length(labels)){metadata_labels[index] <-strsplit(labels[index], split="_")}

### Make the vector into a matrix with columns - make sure no of columns is right!
metadata_labels2<-matrix(unlist(metadata_labels), ncol=7, byrow=TRUE)

### Write the data extracted to a tab delimited text file
write.table(metadata_labels2, "meta_data.txt", quote=F, sep="\t", col.names = c("organism", "survey", "ID", "host", "year", "geno", "country"), row.names = F)

## Import the metadata as a dataframe

meta_data<-read.table("meta_data.txt", header=T)

### MAKE SURE COLUMN 1 HAS THE FULL TIP NAMES OR ELSE IT WON'T WORK!

meta_data2<-cbind(labels, meta_data)

## Make whatever columns you need to into cateogrical variables using as.factor.
## Categories based on numeric classifiers will be seen as continuous not discrete variables if you don't specific now.
meta_data2$spol<-as.factor(meta_data2$spol)

### Check structure of meta data file
str(meta_data2)

### Attach the meta data file to the tree using the %<+% operator from tidyverse

p2<-p %<+% meta_data2

### Now that the data frame is attached, we can use it to add tip colours and labels.
## First, the tip labels

p3<-p2+geom_tippoint(aes(color=host), size=1.0) + scale_color_manual(name="Host", values=c("red", "blue", "green")) 

## Next, labels beside the tree that sit next to each tip

p4<-p3 +  new_scale_fill() + geom_fruit(geom=geom_tile, mapping=aes(fill=country), width=5) + scale_fill_manual(name="Country", values=c("steelblue", "green"), guide=guide_legend(keywidth=0.3, keyheight=0.3))

### And we can add additional ones as required.
### For a classifer that has a numeric value, but you wish it to be used categorically, you need to specify in the dataframe that it is a factor - see above.

p5<-p4 + new_scale_fill() +geom_fruit(geom=geom_tile, mapping=aes(fill=spol), width=5) + scale_fill_manual(name="Spoligotype", values=c("red", "orange", "yellow", "green", "cyan", "darkgreen", "blue", "darkblue", "purple", "violet", "black"), guide=guide_legend(keywidth=0.3, keyheight=0.3))

## Add a scale bar

p5 + geom_treescale(x=10, y=150, width=100, offset=5, color="red", linesize=0.5)

### If tree is outside of plotting limits, use ggsave to adjust and save a PDF
ggsave("infantisplot.pdf", width=50, height=20, units="cm", limitsize =FALSE)