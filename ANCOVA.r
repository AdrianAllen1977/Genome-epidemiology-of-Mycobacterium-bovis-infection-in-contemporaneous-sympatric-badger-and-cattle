### ANCOVA analysis

### Having run several linear regressions in datasets stratified by a categorical factor (eg. year).
## You may want to compare the separate regressions to each-other to search for differences in  slope / intercept
## To do this, you need to run an ANCOVA (Analysis of covariance).

## First, visualise your regression lines all in one plot.
setwd("working_drive_path")


## require packages
library(lattice)

## Import data - data should have the response variable, the explanatory variable and the categorical variable for each regression.

data<-read.table("your_data_file_name", header=T)

## Now plot the curves - note type=c("r) tells the plot function to only plot the regression lines.
## type=c("r", "p") would plot both the lines and the points.
## auto.key just makes each line a different color and adds a key.
## Add axes limits if needed example - ylim=c(10,20)

xyplot(STR_Dist ~ Map_Dist, data=data, groups=Year, type=c("r"), auto.key=T, ylab="Pairwise STR distance between badgers", xlab="Pairwise badger map dist (m)", main="TVR 2014-2018 Individual Years IBD Plots")

## Now, to perform the ANCOVA run a linear model of the "Map_Dist" response variable vs the "STR_Dist" explanatory variable interacting with the categorical "Year" variable.
## The treatment group (in this case Year) must be a factor!
## Otherwise downstream analyses are impossible
## So in your data table, if the treatment is numeric, enter a letter code before the number to make it a factor.

IBD.lm<-lm(STR_Dist~Map_Dist*Year, data=data)

## Then, run an ANOVA of the just run linear model

anova(IBD.lm)

### The resulting table will tell you about the comparison of the curves.  The interaction term is the covariate in this ANCOVA - how STR DIST and Year covary with eachother - it specifically tells you whether there is a difference in model slopes.
## Example table below

Analysis of Variance Table
Response: STR_Dist
Df  Sum Sq Mean Sq   F value    Pr(>F)    
Map_Dist                  1   16512 16512.4 1002.6553 < 2.2e-16 ***
factor(Year)              4   10726  2681.5  162.8249 < 2.2e-16 ***
Map_Dist:factor(Year)     4     491   122.8    7.4589 5.302e-06 ***
Residuals             68489 1127922    16.5                        
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### Very low p value for 'year' also tells us STR_dists vary across the years
### Low p value for Map_Dist shows us STR_Dist increases with Map distance.
## So the interaction Map_Dist:Year does exhibit significant covariance.
### So some regression lines have differing slopes.
## Which ones though?  Need to do multiple comparison post hoc testing.


library(multcomp)
test<-glht(model=IBD.lm,linfct=mcp(Year="Tukey"))

Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts


Fit: lm(formula = STR_Dist ~ Map_Dist * Year, data = data)

Linear Hypotheses:
  Estimate Std. Error t value Pr(>|t|)    
B2015 - A2014 == 0 -1.14463    0.09706 -11.794   <0.001 ***
C2016 - A2014 == 0 -0.71178    0.14063  -5.061   <0.001 ***
D2017 - A2014 == 0 -0.33598    0.12294  -2.733   0.0461 *  
E2018 - A2014 == 0 -0.95680    0.10104  -9.469   <0.001 ***
C2016 - B2015 == 0  0.43285    0.15786   2.742   0.0450 *  
D2017 - B2015 == 0  0.80865    0.14233   5.682   <0.001 ***
E2018 - B2015 == 0  0.18783    0.12390   1.516   0.5367    
D2017 - C2016 == 0  0.37580    0.17497   2.148   0.1898    
E2018 - C2016 == 0 -0.24502    0.16034  -1.528   0.5287    
E2018 - D2017 == 0 -0.62082    0.14507  -4.279   <0.001 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- single-step method)
