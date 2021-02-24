### Extracting details from fasta name strings

## Packages

library("phytools")

## Read fasta
sequences<-read.dna("yourfasta", format="fasta")

## Extract the sequence names from the fasta
names<-labels(sequences)

## Make sure the number of labels is what you expect
length(names)

### The labels may be of this format "TVR850_Bovine_21/6/2012"
## You want to extract the date section after the second underscore
## Use a loop script to do this

## First create and empty vector you can add the date sections to
dates<-c()

## Then populate the vector with the dates using the loop script.

for(index in 1:length(names)){      
dates[index] <- strsplit(names[index], split="_")[[1]][3]
}

## LINE 25 - Creates an index from position 1 in the 'names' object to the final position
## LINE 26 - The loop will create an index in the 'dates' vector which is populated from a string split of the index of the names object you just made
### LINE 26 - end - [[1]][2] is the recursive part of the loop - the [[1]] tells the loop to start at name index position 1, whilst the [3] tells it to only collect the 3rd group split by an underscore - the date section.
## The loop recursively then goes through each line of the names index, and populates the dates vector with the appropriate string split.

## Write the dates vector to a txt file
write.table(dates, "Dates.txt", row.names=F, col.names=F, quote=F)

## Changing column names in a dataframe

colnames(dataframe)<-c("Name1", "Name2", "Name3")
