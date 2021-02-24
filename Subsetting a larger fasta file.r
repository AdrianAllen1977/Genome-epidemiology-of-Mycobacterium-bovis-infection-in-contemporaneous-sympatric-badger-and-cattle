##### Sorting a Fasta File and Subsetting a larger fasta file.

###  Create a two column text file of the names of the sequences you want to subset.  Make column 1 "index" and column 2 "ids".  Don't include the > symbol in the sequence name

###First import the large fasta you want to subset from
### THE LARGE FASTA MUST HAVE THE SEQUENCE NAMES IN ALPHANUMERIC ORDER
### IF NOT ORDERED, THE SUBSETTING WILL NOT WORK - YOU'LL GET SUBSETTED GENOMES BUT WITH INCORRECT DATA

## To reorder a fasta by alphanumeric title of isolates, use the script below:
library(ape)
seq<-read.FASTA("yourfasta.fasta")

sorted_seq<-seq[sort(labels(seq))]

write.FASTA(sorted_seq, file="sortedfastaname.fasta")

## Now you can subset the FASTA further.

library(seqinr)
genomes<-read.fasta("TVR_Red_Dog_Trim1.fasta", seqtype="DNA", as.string=TRUE)

###Now import the subset list file with sequence names
## The subset list must also be in alphanumeric order

subsetlist<-read.table("Badger_Subset_IDs.txt", header=T)

### Check the lists of names to makes sure they match between files
check<-names(genomes) %in% subsetlist$ids
summary(check)

### If the number of trues (matches of ids between lists) is less or more than you expected
### There is a difference in the names of seqs.
### Check for suffixes / prefixes - like 'b' in the TVR set.
### Otherwise, you extract / subset fewer seqs but then mislabel them with the final write.fasta script here
## When you try to then do a SNP dist matrix, it will be wrong.
### Now, if you're content, run the subsetting script

a<-genomes[c(which(names(genomes) %in% subsetlist$ids))]

### Write the new fasta file to your working drive

write.fasta(a, names=subsetlist$ids, file.out="263_all_Pos_Badger.fasta")
