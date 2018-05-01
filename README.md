# carrot-diallel

## Synopsis

This repository contains data and corresponding code to perform diallel analysis of shoot growth in carrots. Scripts include analyses using
Griffing's Method I Model I (Griffing 1956) and BayesDiallel (Lenarcic et al. 2012, Crowley et al. 2014). 

## Data sets (.csv)
1. **parentalData.csv**  
phenotypic data for parents used to construct the diallel  
2. **diallelRawData.csv**  
phenotypic observations for diallel progenies, includes variables for year, replication, male parent, and female parent  
3. **Diallel tau prior.csv**  
prior information for BayesDiallel

## Scripts

**01_diallel_parental_means.R**  
computes and plots means + 95% confidence intervals for parent phenotypes  
**02_missing_data_heatmap.R**  
plots phenotypic values and location of missing data by cross, year, and replication
**03_phenotypic_correlations.R**  
calculates and plots Pearson's correlations for phenotypes
**04_impute_missing_data.R**  
uses mice package to impute missing data (see van Buuren and Groothuis-Oudshoorn 2011); generates imputed datasets for pooled diallel analysis
**05_pooled_diallel_analysis.R**  
computes Griffing's ANOVA and effect estimates (i.e. GCA, SCA, reciprocal) for each imputed dataset and combines results. 
Calls functions from three other scripts (listed below)  
      **05A_diallel_analysis_functions.R**   
      modified diallele1 function from plantbreeding package (Rosyara 2014); 
      allows calculation of Griffing's ANOVA with multiple environments; includes additional edits to 
      SCA calculation from G. Ramstein  
      **05B_calcSCA.R**  
      function to calculate SCA estimates for multiple environments and pooling  
      **05C_calcRecip.R**  
      function to calculate reciprocal effect estimates for multiple environments and pooling  
      **05D_pool_Griffing_Method3.R**  
      Griffing's ANOVA using method III (without parents) to estimate Baker's ratio
**06_GGE_biplots.R**  
GGE biplots for diallel phenotypes & parents following Frutos et al. (2014) 
**07_BayesDiallel_fulldiallelanalyze.R**  
applies and stores AFD objects from DiallelAnalyzer function in BayesDiallel  
**08_read_AFD_objects.R**  
reads in AFD objects for subsequent analyses in BayesDiallel  
**09_BayesDiallel_plots.R**  
plots observed vs. expected phenotypes, highest posterior density (HPD) intervals, and strawplots 
(see Lenarcic et al. 2012 and BayesDiallel documentation)   
**10_create_psq_df.R**  
exports table of posterior PSq values (diallel variance projection [VarP] in Crowley et al. 2014) 
for all traits  
**11_VarP_plot.R**  
plots relative contribution of inhertiance classes to VarP  
**12_degree_of_dominance.R**
uses BayesDiallel AFD results to estimate the degree of dominance and the dominance index for crosses in a diallel (see Maurizio et al. 2018)  
**13_BayesDiallel_fulldiallelanalyze_by_env.R**  
applies and stores AFD objects for each environment  
**14_read_AFD_by_env.R**  
reads in AFD objects by environment
**15_ranks_by_environment.R**  
estimates and hybrid rankings in each environment
**16_plot_ranks_by_environment.R**  
plots hybrid rankings by environment


## References
Crowley JJ, Kim Y, Lenarcic AB, Quackenbush CR, Barrick C, Adkins DE, Shaw GS, Miller DR, Pardo Manuel de Villena F, Sullivan PF, 
Valdar W (2014) Genetics of adverse reactions to haloperidol in a mouse diallel: A drug-placebo experiment and Bayesian causal analysis. 
*Genetics* 196(1):321-47.  

Frutos E, Purificación Galindo M (2014) An interactive biplot implementation in R for modeling genotype-by-environment interaction. *Stoch Environ Res Risk Assess* 28:1629-1641.

Griffing B (1956) Concept of general and specific combining ability in relation to diallel crossing systems. *Aust. J. Biol. Sci.* 9:463-493.   

Lenarcic AB, Svenson KL, Churchill GA, Valdar W (2012) A general Bayesian approach to analyzing diallel crosses of inbred strains. 
*Genetics* 190:413-435. doi: 10.1534/genetics.111.132563   

Rosyara U (2014) plantbreeding: Analysis and visualization of data from plant breeding and genetics experiments. 
http://R-Forge.R-project.org/projects/plantbreeding/   

van Buuren S and Groothuis-Outshoorn K (2011) mice: multivariate imputation by chained equations in R. *J. Stat. Softw.* 45:1-67  

BayesDiallel website: http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html

