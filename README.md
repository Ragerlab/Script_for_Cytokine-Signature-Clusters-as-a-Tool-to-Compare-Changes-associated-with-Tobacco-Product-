# Cytokine Compartment Analysis
> Exploratory analyses to determine association between smoking status/ cotinine concentrations and cytokine concentrations.

> All cytokine compartment analyses in this folder are designated by their figure number or table number in the manuscript in parantheses.


## Descriptive Statistics (Table S1)
- Calculated mean, median, min, max, and standard deviation for the overall cohort, non-smokers, cigarette smokers, and e-cigarette smokers by compartment.

## Cytokine Distribution Comparison (Table S2, Table S3)
- Using Wilcoxon Rank Sum tests to compare baseline & cigarette smoker distributions or baseline & e-cigarette smokers by compartment (Table S2).
- Using Wilcoxon Rank Sum tests to compare baseline cytokine concentrations across compartments (not in manuscript). 
- Used ANOVA test as a crude model (serving as similar function as original Wilcoxon Rank Sum tests) to compare ANCOVA results. ANCOVA was run to control for age (Table S3).

## Cluster Distribution Comparison (Table S5)
- Used ANOVA test as a crude model (serving as similar function as original Wilcoxon Rank Sum tests) to compare ANCOVA results. ANCOVA was run to control for age (Table S3). This is the same analysis in the cytokine distribution comparison, but now we are comparing clusters of cytokines rather than individual cytokines.
   
## Baseline Analysis (Figure 1, Table 2, and Figure S1)
- Calculating and visualizing mean and standard deviation of baseline (non-smokers) concentrations for each cytokine (Figure 1). 
- Running Shapiro-Wilk's test for normality and Spearman rank correlation tests to determine correlation between each compartment (NLF, ELF, Sputum, or Serum) at baseline (Table 2). 
- Calculating and visualizing median and standard deviation of baseline (non-smokers) concentrations for each cytokine (Figure S1). 

## Simulation Analysis (Figure 2)
- Simulated 500 concentration data points to see if randomly generated data would yield similar statistically significant and highly correlated results between each across compartmental analysis. 

## Demographics Analysis (Table 1)
Calculated demographics statistics including fisher's exact and ANOVA

## Cytokine Demographics Distribution Comparison (not in manuscript)
Running Wilcoxon Rank Sum test and stratifying subjects based on the demographic variable of interest (race, ethnicity, sex, age, or bmi) further stratifying compartment to determine how cytokine distributions varied for each variable. 

## Individual Cytokine Rankings (not in manuscript)
- Gave each cigarette or e-cigarette smoker a score based on the amount of deviation from baseline for each cytokine .
- Each smoker was then ranked from 1-42 with the top 5 deviators used for the top deviators visualization.

## Top Deviators Analysis (not in manuscript)
- Calculating and visualizing ('Individual_Smoker_Rankings') mean and standard deviation of smokers (either cigarette or e-cigarette smokers) that were the top 5 deviators amongst the smokers. 
- Baseline mean and standard deviation were also plotted for comparison and the plot was stratified by compartment. 

## Cotinine Analysis (not in manuscript)
- Running Spearman correlation tests to determine association between cytokine and cotinine concentrations across compartments.
