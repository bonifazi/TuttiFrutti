<p align="center">
<img src="https://github.com/bonifazi/TuttiFrutti/assets/74569672/448b4952-b723-4f35-9ac0-fafb162cbade" width="300" height="300">
</p>

# TuttiFrutti is a collection of miscellaneous R functions and R scripts.  

## Licence
You are free to use these [R functions](https://github.com/bonifazi/R_utils/edit/main/README.md#list-of-r-functions) and [R scripts](https://github.com/bonifazi/R_utils/edit/main/README.md#list-of-r-functions).
This project is licensed under the MIT license - see the [License](https://github.com/bonifazi/R_utils/blob/main/LICENSE) file for details.

## Citation
Please acknowledge the original referenced papers for their methodology as indicated within the specific R functions or R scripts under `@references`. Moreover, if you use this code in your research or find it helpful, please consider acknowledging this repository (https://github.com/bonifazi/TuttiFrutti) by citing it as:  
_TuttiFrutti: A Collection of Miscellaneous R Functions and R Scripts. Bonifazi R. 2023. GitHub Repository. Available online: https://github.com/bonifazi/TuttiFrutti_

```bibtext
@misc{TuttiFrutti,
title={TuttiFrutti: A Collection of Miscellaneous R Functions and R Scripts},
author={Renzo Bonifazi},
year={2023},
howpublished={\url{https://github.com/bonifazi/TuttiFrutti}},
}
```

## Contact
Your suggestions on improving these functions and scripts are very welcome. If you have any suggestions, questions, or need support, feel free to contact renzo.bonifazi@outlook.it, or open an issue here on GitHub.

## List of R functions:
### Description on usage:
The R functions are all documented with the R [docstring](https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html) package. Read the documentation to know more about how to use them, and what options are available. Where possible, I provide some examples to show their functionalities.  
To view the functions' documentation, first load the [`docstring`](https://github.com/Dasonk/docstring) package in R, then view the documentation of the function by running `docstring(fun = "<functionname.R>")` in the console.
### List:
* [Plot convergence of MiXBLUP](https://github.com/bonifazi/R_utils/blob/main/PlotConvergeneMiXBLUP.R). `MiXBLUP` has a gnuplot code that automatically generates a convergence graph. When gnuplot is not available on your system, you can generate the same graph using this function. The output is a `ggplot2` R object.
* [Convert a MiX99 parameter file into (co)variance matrices](https://github.com/bonifazi/R_utils/blob/main/meltParfile.R). `MiX99` parameter file for (co)variance components has a lower triangular 'long' format as "`effect_number,i,j,covar_value`". This function converts it into full symmetric (co)variance matrices.
* [Rebase (G)EBVs](https://github.com/bonifazi/R_utils/blob/main/rebase_ebv.R). A function to rebase vector(s) of EBVs given a list of animals to be considered as the base population.
* [Compute LR method statistics](https://github.com/bonifazi/R_utils/blob/main/compute_LR_stats.R). A function that, taking EBVs from a partial and a whole evaluation, computes LR method statistics following [Legarra and Reverter, 2018, GSE, 50:53](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0426-6). Several options are available. The main ones are that: 1) the statistics can be computed with and without providing a 'focal group' of individuals, i.e., a validation group; 2) Plotting can be used to investigate better if the validation group is homogenous. 3) Bootstrapping with replacement can be used to compute standard errors of the estimated statistics.  
* [Compute general validation statistics](https://github.com/bonifazi/TuttiFrutti/blob/main/Validation_stats.R). Function to compute general validation statistics (to be documented).
* [Functions for the integration of international (G)EBVs into national evaluations](https://github.com/bonifazi/Integration_EBV_and_GEBV). R functions for integrating (genomic) estimated breeding values ((G)EBVs) into national evaluations following [Bonifazi et al., 2023, GSE](https://doi.org/10.1186/s12711-023-00813-2). See the [README](https://github.com/bonifazi/Integration_EBV_and_GEBV/blob/main/README.md) for more information.
* ... [new functions will be added here]

## List of R scripts:
### Description on usage:
These Rscripts are to be run in a command-line style, e.g. `Rscript --vanilla script.R --option1 data --option2 output.csv`. They have a help option which you can call using: `Rscript --vanilla script.R -h` or `Rscript --vanilla script.R --help`. 
### List:
* [Analyse Plink ROH](https://github.com/bonifazi/R_utils/blob/main/Analyse_Plink_ROH.R). A command-line Rscript to analyse genomic inbreeding from Runs of Homozygosity (ROH) obtained from Plink using the `detectRUNS` R package. This script produces .csv and .pdf files on ROH inbreeding. To run it:  
`Rscript --vanilla Analyse_Plink_ROH.R --plink_files plinkcleaned --plink_roh ROH.hom --group geno_BRD --pedigree ped.ped --output results_dir 2>&1 | tee logfile.log`  
To see a description of each argument use: `Rscript Analyse_Plink_ROH.R --help.`  
Note that the `detectRuns` package groups the results based on the number of groups in the first column of the ROH files, which I guess can be used if you want to define sub-populations. For now, the Rscript internally overrides the ‘group’ column of the `--plink_roh` file with that given in the `--group` label.  
* [Extract EBV and REL from asreml .sln file](https://github.com/bonifazi/TuttiFrutti/blob/main/ExtractAsremlSolutions.R). A command-line Rscript to extract EBV and REL from asreml .sln file. This script produces a .csv file with ID, EBV, REL, and (user-provided) VAR(A) for each trait in the asreml solution file (.sln). To run it:  
`Rscript --vanilla ExtractAsremlSolutions.R --file my_path/asreml.sln --effect_name effect5 --trait_names "Trait1, Trait2, Trait3" --varA "varA_trait1, varA_trait2, varA_trait3" --output myoutput.csv`  
To see a description of each argument use: `Rscript ExtractAsremlSolutions.R --help.`  
* [Plot postgibbsf90 convergence](https://github.com/bonifazi/TuttiFrutti/blob/main/plot_postgibbsf90.R). Rscript to plot postgibbsf90 covariances. The script produces a .pdf file. See --help for usage and read the details section of the Rscript for more information on input and settings. To run it:  
`Rscript --vanilla plot_postgibbsf90.R --file my_path/postgibbs_samples --output postgibbs_plots.pdf`
or with more control on output, e.g., adding trait names and making two main groups:  
`Rscript --vanilla plot_postgibbsf90.R -f my_path/postgibbs_samples -o postgibbs_plots.pdf -g 2 -t "BRD1_AWW, BRD2_AWW, BRD1_CE, BRD2_CE" `
* ... [new Rscripts will be added here]

## Why TuttiFrutti? 
TuttiFrutti is used for products that mix different colours and flavours, such as candies. My Swedish colleague often told me 'TuttiFrutti!' to show off his Italian skills. In Italian, 'Tutti Frutti' translates to 'all flavours', which is similar to what this repository is intended to be: a mix of different R things (and maybe more later on).
