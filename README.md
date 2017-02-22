# carrot-diallel

## Synopsis

This repository contains data and corresponding code to perform diallel analysis of shoot growth in carrots. Code outlines analyses using
Griffing's Method I Model I (Griffing 1956) and BayesDiallel (Lenarcic et al. 2012, Crowley et al. 2014). 

## Data sets (.csv)
1. **parentalData.csv**  
phenotypic data for parents used to construct the diallel  
2. **diallelRawData.csv**  
phenotypic observations for diallel progenies, includes variables for year, replication, male parent, and female parent  
3. **Diallel tau prior.csv**  
prior information for BayesDiallel

## Scripts

1. **diallel parental means.R**  
computes and plots means + 95% confidence intervals for parent phenotypes  
2. **diallel missing data heatmap.R**  
plots phenotypic values and location of missing data by cross, year, and replication
3. **diallel correlations.R**  
calculates and plots Pearson's correlations for phenotypes
4. **diallel imputation.R**  
uses mice package to impute missing data (see van Buuren and Groothuis-Oudshoorn 2011); generates imputed datasets
5. **pooled diallel analysis.R**  
computes Griffing's ANOVA and effect estimates (i.e. GCA, SCA, reciprocal) for each imputed dataset and combines results. 
Calls functions from three other scripts (listed below)  
      5A. **Diallel_analysis_functions_ST.R**   
          modified diallele1 function from plantbreeding package (Rosyara 2014); 
          allows calculation of Griffing's ANOVA with multiple environments; includes additional edits to 
          SCA calculation from G. Ramstein  
      5B. **calcSCA.R**  
          function to calculate SCA estimates for multiple environments and pooling  
      5C. **calcRecip.R**  
          function to calculate reciprocal effect estimates for multiple environments and pooling  
6. **BayesDiallel full diallel analyze.R**  
applies and stores AFD objects from DiallelAnalyzer function in BayesDiallel  
7. **read in AFD objects.R**  
reads in AFD objects for subsequent analyses in BayesDiallel  
8. **BayesDiallel plots.R**  
plots observed vs. expected phenotypes, highest posterior density (HPD) intervals, and strawplots 
(see Lenarcic et al. 2012 and BayesDiallel documentation)   
9. **create psq data frame.R**  
exports table of posterior PSq values (diallel variance projection [VarP] in Crowley et al. 2014) 
for all traits  
10. **p2 diallel variance projection plot.R**  
plots relative contribution of inhertiance classes to VarP 

## References
Crowley JJ, Kim Y, Lenarcic AB, Quackenbush CR, Barrick C, Adkins DE, Shaw GS, Miller DR, Pardo Manuel de Villena F, Sullivan PF, 
Valdar W (2014) Genetics of adverse reactions to haloperidol in a mouse diallel: A drug-placebo experiment and Bayesian causal analysis. 
*Genetics* 196(1):321-47.  

Griffing B (1956) Concept of general and specific combining ability in relation to diallel crossing systems. *Aust. J. Biol. Sci.* 9:463-493.   

Lenarcic AB, Svenson KL, Churchill GA, Valdar W (2012) A general Bayesian approach to analyzing diallel crosses of inbred strains. 
*Genetics* 190:413-435. doi: 10.1534/genetics.111.132563   

Rosyara U (2014) plantbreeding: Analysis and visualization of data from plant breeding and genetics experiments. 
http://R-Forge.R-project.org/projects/plantbreeding/   

van Buuren S and Groothuis-Outshoorn K (2011) mice: multivariate imputation by chained equations in R. *J. Stat. Softw.* 45:1-67  

BayesDiallel website: http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html

