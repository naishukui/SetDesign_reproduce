# SetDesign_reproduce
The simulation figures can be generated following the steps here.


### Step 1

Submit batch jobs for the following scripts to generate simulation results for both simulated and theoretical power:

- `snpC.R`: Linear outcomes without correlation.\
- `snpCcor.R`: Linear outcomes with correlation.\
- `snpD.R`: Binary outcomes without correlation.\
- `snpDcor.R`: Binary outcomes with correlation.\

Key parameters for these scripts:
- `n`: Total number of subjects.\
- `k`: Total number of SNPs.\
- `rho`: Correlation coefficient between SNPs.\
- `list1`: Vector of genetic effect sizes for the first alternative allele.\
- `list2`: Vector of genetic effect sizes for the second alternative allele.\
 
### Step 2

Run `results.R` to:
- Calculate theoretical power using the new method developed in this package.
- Generate plots comparing the four different power calculation methods across various effect size combinations.

