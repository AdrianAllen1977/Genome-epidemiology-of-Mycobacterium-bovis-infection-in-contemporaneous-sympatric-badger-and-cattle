##################################
## TransPhylo for phylodynamics ##
##################################

# Refs;  Didelot et al (2017) Mol Biol Evol 34(4): 997-1007; 
# Based on tutorial - Didelot et al. 2021 Current Protocols e60 Volume 1

###############
### PACKAGES ##
###############

library(TransPhylo)
library(ape)
library(coda)
library(lattice)

########################################################
### Simulating a simple dataset to introduce outputs ###
########################################################

set.seed(0) ## sets random number generator so results are reproducible
s<-simulateOutbreak(nSampled=10) ## Simulates a ten case outbreak - creates 'ctree' class of object
pt<-extractPTree(s) ## Extracts the phylogeny for the outbreak 's'
plot(pt) # plots the phylogeny
tt<-extractTTree(s) # Extracts the transmission tree of the outbreak 's'
plot(tt) # plots the transmission tree
plot(s) ## plots the simulated colour tree (ctree) of the outbreak

## Transmission tree shows unfilled circles - hypothesised un-sampled cases.
## Also gives hypothesised who infected whom data

## ctree object contains all info of the phylogeny AND transmission tree
## Transmission events are marked as asterisks / stars - 
## NODES ARE NOT TRANSMISSION - they are just ancestors

# So, that is the introductory bit to show the different outputs.
# In simulated data sets - the transmission tree and ctree are known
# In the real world, all we have is the phylogeny, and need to generate the other two outputs


#############################################
## Further options in simulating outbreaks ##
#############################################

## Multiple processes affect epidemic modelling.  TransPhylo can simulate these.
## Simulation is useful to appreciate how tweaking parameters of each process affects epidemic

# 1 - The Transmission Process
# dateStartOutbreak parameter is self explanatory.
# An index case goes on to infect others in a branching process with the reproduction number (R0) distribution informing this process.
# Neg binomial dist is used to represent R0
# Parameters of R0 are off.r and off.p
# R0 = off.r*off.p/(1-off.p)
# Interval between infection and transmission is based on a gamma distribution
# Gamma parameters are w.shape and w.scale
# You can also use w.mean and w.std instead of w.scale and w.shape to set intervals.

# 2 - Within host evolution process
# Within host evolution is assumed to follow coalescent process with constant Ne
# Single parameter is neg - which is average time of coalsecence between two lineages
# Bottleneck at time of transmission is assumed - only a single lucky clone escapes to infect new host

# 3 - Sampling of the epidemic process
# Sampling of epidemics is imperfect - not all positives are found by diag methods.
# Sampling of cases assumed to have probability pi for each case
# Assumed to happen at a time, post infection, that is drawn from distribution.
# This distribution is gamma dist with parameters ws.shape and ws.scale
# Limit to sampling times can be introduced by specifying parameter dateT
# Any case after dateT specification is not sampled.
# You can set dateT to dateT=Inf to leave the outbreak ongoing - but this means the simulation continues (options to handle this are presented later)

# You can set all of these different parameters in the simulateOutbreak function.
# Defaults are assumed for all unless specified.
# dateStartOutbreak = 2000; off.r=1, off.p=0.5, w.shape=2, w.scale=1, neg=0.25, pi=0.5, ws.shape=w.shape, wss.scale=w.scale, dateT=Inf
# Optional parameter nSampled can be included to force a simualtion with only x sampled cases (as above).

# Simulating an outbreak with R0>1 can lead to a lot of computational time and effort.
# R0>1 can go exponential in growth and TransPhylo tries to simulate the whole thing!
# Simulating an R0>1 outbreak therefore needs some careful choice in parameters to limit number of infectees

{setTimeLimit(elapsed=5, transient=T) 
s=simulateOutbreak(off.r=3)}

# elapsed is number of secs command runs, transient=T limits to this command only

# OR you can just use specific exact parameters to limit no. of infectees
# Hypothetical outbreak with pathogen of generation time 1 year
# In gamma dsitrbution w.shape=10, w.scale=0.1 are equal to w.mean=1 and w.std=0.316

set.seed(0)
s=simulateOutbreak(off.r=3, dateStartOutbreak=2020, dateT=2023, w.shape=10, w.scale=0.1)

## plot a more detailed transmission tree
plot(extractTTree(s), type='detailed', w.shape=10, w.scale=0.1)

## Each sampled and inferred case appears on separate line
## Vertical arrows between lines infer who infected whom.
## There were 4 unsampled cases. Index case was unsampled

# plot the colour tree
plot(s)

# Unsampled cases are the pink, light blue and dark blue lines plus an extra host on the dark blu line leading to sampled case 5.  Limits in colour pallet affect visualisation of all hosts on a ctree.

# So use the ctree in conjunction with the transmission tree to understand the outbreak more.

## In the example above you don't get a feel for the exponential increase - but by increasing the number of years, you can.

set.seed(0)
s<-simulateOutbreak(off.r=3, dateStartOutbreak=2020, pi=0.1, dateT=2027, w.shape=10, w.scale=0.1)
plot(extractTTree(s), showLabels=F)

## NOTE - TransPhylo only retains unsampled cases which lead directly to sampled ones.
## So, tree pruning occurs that leads to exclusion of some sub-lineages that don't fit the above criteria.
## So, it is important to have a systematic sampling of the epidemic you are investigating.
## This pruning is why you can have seemingly a higher proportion of sampled cases to unsampled than pi infers.

## Simulating outbreaks with realistic parameters can help you determine if an actual dataset is producing sensible results.

## You can export a simulated outbreak's phylogeny for use in other programmes.
# You can export as either a Newick or Nexus format.

phy<-phyloFromPTree(extractPTree(s))
write.tree(phy, file="tree.nwk")
write.nexus(phy, file="tree.nex")
print(dateLastSample(s))

# Newick and Nexus trees store dates of tips in relative time.
# You might need exact dates - so you can extract 1 date and then work out the rest from the tree.


#############################################
## Inferring transmission from actual data ##
#############################################

## A dated phylogeny is needed to infer transmission in TransPhylo.
## The output from BEAST analyses - an MCC tree is most useful here.
## However, for this tutorial, we'll simulate an outbreak and export it's phylogeny to interrogate.

set.seed(0)
s<-simulateOutbreak(off.r=3, dateStartOutbreak=2020, dateT=2025, pi=0.3, w.shape=5, w.scale=0.3, neg=1)
print(s) ## Prints summary of simulated outbreak - no cases, sampled and unsampled.
phy<-phyloFromPTree(extractPTree(s))
dls<-round(dateLastSample(s), digits=4) # returns 2024.9646
phy$edge.length=round(phy$edge.length, digits=4) # converts tree edge lengths to 4 sig figs to date tips from last date
write.tree(phy, file="input.nwk")

## Reimport the tree you just made to be the input for the inference process.
## Your tree must be in the Newick format - export from FigTree in this format

tree<-read.tree("input.nwk")

## If importing an MCC tree from BEAST - be aware, treeannotator can introduce negative branch lengths.
## So a descendant node can appear to be older than its ancestor
## This occurs because an MCC tree is the average of multiple trees - and some clades don't appear in all trees.
## As a result, the averaging process can place these clades with negative branches from ancestors who have a better established position on the tree.
## You can manually correct branch lengths in R - before you do this makes sure of several things
## Your tree has good support, MCMC converged etc.
## Also keep account of how many negative branch lengths you manually had to set to zero.
## Use script below to set negative branch lengths to a small non zero number (zero isn't allowed!)
## Branch lengths if zero can make polytomies, which again aren't allowed in TransPhylo.
## Beast MCC trees don't include polytomies, always seeking bifurcation.
## BUT - add in branch lengths of zero can introduce them to your MCC tree.

tree$edge.length<-pmax(tree$edge.length, 1/365) ## day used as branch length
## tree2<-multi2di(tree) ## removes polytomies - shouldn't be needed for a beast MCC tree.
plot(p)

# Output your new edge adjusted tree

write.tree(tree, "TVR_Skyline_edge_MCC.nwk")

## Import new tree and set latest date of sampling

tree<-read.tree("input.nwk")
p<-ptreeFromPhylo(tree, dateLastSample=2024.9646)

#############################
### Running the inference ###
#############################

## Assumptions

# You need to define the last sample date dateT - if outbreak is ongoing do not use dateT=Inf
# You also need to provide w.shape and w.scale or w.mean and w.std parameters for the gamma dist
# For a real dataset, you need to estimate these using prior publications etc
# Dist of time between infection and sampling (ws.scale and ws.shape) is by default set to be same as dist for the time between being infected and transmitting.
# Latter assumes diagnosis is typically associated with onset of symptoms - mmm maybe not true for endemic livestock disease.
# MCMC needs to be run until convergence with sufficient sampling / thinning to store good estimates.

set.seed(0)
r<-inferTTree(p, w.shape=5, w.scale=0.3, dateT=2025, mcmcIterations=1e5, thinning =10)

##TVR code ##r<-inferTTree(p, w.shape=1.3, w.scale=3.33, ws.shape=1.05, ws.scale=2.85, dateT=2017, mcmcIterations = 1e7, thinning=10)

## A progress bar for the MCMC will appear in R terminal.

# Before you do anything else, check the MCMC has converged.

plot(r)

# Look for the characteristic hairy caterpillar trace to determine convergence.

# BUT, it's worth digging a little deeper to check the outputs of the MCMC for quality
# Check the ESS for each trace to ensure your sampling of the MCMC is representative
# Increase the mcmcIteration to increase ESS - it will take longer to run though!

mcmc=convertToCoda(r)
effectiveSize(mcmc)

## Aim for ESS above 200 to ensure representativeness of parameter estimates

######################################################
## Checking & understanding the inferred parameters ##
######################################################

print(r) # mean and 95% credible interval for all parameters inferred from outbreak.

## Check the inferred transmission tree.
## TransPhylo's MCMC explores multiple transmission trees - quite like BEAST
## You can summarise all these trees using the script below:

med<-medTTree(r)
plot(med) # plots summary ctree
plot(extractTTree(med), type='detailed', w.shape=5, w.scale=0.3) # plots detailed who infected whom tree

## The trouble with condensing and summarising all trees is that some nuance is lost.
## Specifically, some of the uncertainty around transmission dynamics gets lost.
## In real life, there can be multiple transmission trees with decent statistical support.
## Computing the probability of transmission between pairs of cases can be done.
## The TTree can do this by summarising the average number of transmission links between close phylogenetic pairs.

matWIW<-computeMatWIW(r) # infers who infected whom matrix for all sampled cases
levelplot(matWIW) # plots the who infected whom matrix as a heat map.
matTDist<-computeMatTDist(r) # infers average no. of transmission links between pairs of sampled cases
levelplot(matTDist) ## plots the inferred transmission links matrix as a heat map

## Inspect cases where posterior probability of transmission between pairs is >50%

### Further breakdown of individual sampled cases

# You can see the distribution across multiple trees for when an individual case got infected and how many daughter infections it cassed.
# Can potentially identify superspreaders.

time<-getInfectionTimeDist(r, k='1', show.plot=T) # for sampled case no. 1
daughterinf<-getOffspringDist(r, k='1', show.plot=T )

## Summary stats for wider dataset ##

# Plot distribution of time between getting infected and transmitting
getGenerationTimeDist(r, show.plot=T, maxi=5) # maxi sets upper limit on times considered.

# Plot distribution of time between getting infected, and being sampled / detected
getSamplingTimeDist(r, show.plot=T, maxi=4)

