<p align="center">
<img src="https://www.youwish.nl/wp-content/uploads/2023/03/Tutti-Frutti-Recovered.jpg" width="300" height="300">
</p>

# TuttiFrutti is a collection of miscellaneous R functions and R scripts.  
You are free to use these [R functions](https://github.com/bonifazi/R_utils/edit/main/README.md#list-of-r-functions) and [R scripts](https://github.com/bonifazi/R_utils/edit/main/README.md#list-of-r-functions), for more on their usage, see the MIT [License](https://github.com/bonifazi/R_utils/blob/main/LICENSE) documentation.  
Please acknowledge its source, i.e., this repository (https://github.com/bonifazi/R_utils), so that it may be useful to others as well.
Your suggestions on improving these functions are very welcomed. Send me an email at (renzo.bonifazi@outlook.it) or open an issue here on GitHub.

## List of R functions:
### Description on usage:
The R functions are all documented with the R [docstring](https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html) package. Read the documentation to know more how to use them, and what options are available. Where possible, I provide some examples to show their functionalities.  
To view the functions' documentation, first load the `docstring` package in R, then view the documentation of the function by running `docstring(fun = "<functionname.R>")` in the console.
### List:
* [plot convergence of MiXBLUP](https://github.com/bonifazi/R_utils/blob/main/PlotConvergeneMiXBLUP.R). `MiXBLUP` has a gnuplot code that automatically generates a convergence graph. When gnuplot is not available on your system, you can generate the same graph using this function. The output is a `ggplot2` R object.
* [Convert a MiX99 parameter file into (co)variance matrices](https://github.com/bonifazi/R_utils/blob/main/meltParfile.R). `MiX99` parameter file for (co)variance components has a lower triangular 'long' format as "`effect_number,i,j,covar_value`". This function converts it into full symmetric (co)variance matrices.
* [rebase EBV](https://github.com/bonifazi/R_utils/blob/main/rebase_ebv.R). A function to rebase vector(s) of EBVs given a list of animals to be considered as the base population.
* [compute LR method statistics](https://github.com/bonifazi/R_utils/blob/main/compute_LR_stats.R). A function that, taking EBVs from a partial and a whole evaluation, computes LR method statistics following [Legarra and Reverter, 2018, GSE, 50:53](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0426-6). Several options are available. Main ones are that: 1) the statistics can be computed with and without providing a 'focal group' of individuals, i.e., a validation group; 2) Plotting can be used to investigate better if the validation group is homogenous. 3) Bootstrapping with replacement can be used to compute standard errors of the estimated statistics.
* ... [new functions will be added here]

## List of R scripts:
### Description on usage:
These Rscripts are to be run in a command-line style, e.g. `Rscript --vanilla script.R --option1 data --option2 output.csv`. They have an help option which you can call using: `Rscript --vanilla script.R -h` or `Rscript --vanilla script.R --help`. 
### List:
* [Analyse Plink ROH](https://github.com/bonifazi/R_utils/blob/main/Analyse_Plink_ROH.R). A command-line Rscript to analyse genomic inbreeding from Runs of Homozygosity (ROH) obtained from Plink using the `detectRUNS` R package. This script produces .csv and .pdf files on ROH inbreeding. To run it:  
`Rscript --vanilla Analyse_Plink_ROH.R --plink_files plinkcleaned --plink_roh ROH.hom --group geno_BRD --pedigree ped.ped --output results_dir 2>&1 | tee logfile.log`  
To see a description of each argument use: `Rscript Analyse_Plink_ROH.R --help.`  
Note that the `detectRuns` package groups the results based on the number of groups in the first column of the ROH files, which I guess can be used if you want to define sub-populations. For now, the Rscript internally overrides the ‘group’ column of the `--plink_roh` file with that given in the `--group` label.
* ... [new Rscripts will be added here]

## Why TuttiFrutti? 
TuttiFrutti is sometimes used for products that have a mix of different colours and flavours, such as candies. My Swedish collegue often said to me 'TuttiFrutti!' to show off his italian skills. In Italian, 'Tutti Frutti' translates to 'all flavours', which is similar to what this repository is intended to be: a mix of different R things (and maybe more later on).
