# SetDesign_reproduce
Code to reproduce simulation results, figures, and real data analysis results from the paper "Large Impact of Genetic Data Processing Steps on Stability and Reproducibility of Set-Based Analyses in Genome-Wide Association Studies" by Naishu Kui, Yao Yu, Jaihee Choi, Xihao Li, Zachary R. McCaw, Chad Huff and Ryan Sun.

Please first install the SetDesign package, for example through the command devtools::install_github("naishukui/SetDesign"), then follow the next steps to generate the simulation results and figures.
### Step 1

Submit batch jobs(bsub <simu.lsf) for the following scripts to generate simulation results for both simulated and anlytical power:

- `snpC.R`: Linear outcomes without correlation.
- `snpCcor.R`: Linear outcomes with correlation.
- `snpD.R`: Binary outcomes without correlation.
- `snpDcor.R`: Binary outcomes with correlation.

Key parameters for these scripts:
- `n`: Total number of subjects.
- `k`: Total number of SNPs.
- `rho`: Correlation coefficient between SNPs.
- `list1`: Vector of genetic effect sizes for the first alternative allele in tri-allelic position.
- `list2`: Vector of genetic effect sizes for the second alternative allele in tri-allelic position.
 
### Step 2

Run `results.R` to:
- Calculate analytical power using the new method developed in this package.
- Generate plots comparing the four different power calculation methods across various effect size combinations.


