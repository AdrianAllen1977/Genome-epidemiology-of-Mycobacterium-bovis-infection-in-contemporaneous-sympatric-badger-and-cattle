## Randomising Tip Dates for BEAST analyses to asses temporal signal.## Rieux, A. and Khatchikian, C.E. (2017), TIPDATINGBEAST: an R package to assist the implementation of phylogenetic tip-dating tests using BEAST. Molecular Ecology Resources, 17: 608-613. ## To test for significance of temporal signal, randomise tip dates in an XML and run multiple times in BEAST.## Packagesinstall.packages("TipDatingBeast")library(TipDatingBeast)

## Set your working drive to the location of the XML created in BEAUTi

setwd("your location")

## You do not need to import the xml, just run the script below with XML file's name minus the xml suffix.
## The reps argument indicates how many date randomised xmas you want to produce

RandomDates("your_xml_name", reps=10)

## The script now outputs the desired number of date randomised xml files.