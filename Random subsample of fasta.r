###############################################################
# FastaUtils - subsetting a fasta and random sampling from one#
###############################################################

##########
#Packages#
##########

# install FastaUtils via bioconductor

devtools::install_github("GuillemSalazar/FastaUtils")

# require package

library(FastaUtils)

#########
#Scripts#
#########

# Random sampling from a larger fasta file without replacement

fasta.sample(infile="yourname.fasta", nseq=x, file.out="yourname.fasta",replacement=FALSE)

