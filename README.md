# Collection of miscellaneous R functions.  
You are free to use these R functions, for more on their usage, see the MIT [License](https://github.com/bonifazi/R_utils/blob/main/LICENSE) documentation.  
Please acknowledge its source, i.e., this repository (https://github.com/bonifazi/R_utils), so that it may be useful to others as well.
Your suggestions on improving these functions are very welcomed. Send me an email at (renzo.bonifazi@outlook.it) or open an issue here on GitHub.

## Description:
These functions are all documented with the R [docstring](https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html) package. Read the documentation to know more how to use them, and what options are available. Where possible, I provide some examples to show their functionality.  
To view the R functions' documentation, first load the `docstring` package, then view the documentation of the function by running in the console `docstring(fun = "<functionname.R>")`.

## List of functions:
* [plot convergence of MiXBLUP](https://github.com/bonifazi/R_utils/blob/main/PlotConvergeneMiXBLUP.R). `MiXBLUP` has a gnuplot code that automatically generates a convergence graph. When gnuplot is not available on your system, you can generate the same graph using this function. The output is a `ggplot2` R object.
* [Convert a MiX99 parameter file into (co)variance matrices](https://github.com/bonifazi/R_utils/blob/main/meltParfile.R). `MiX99` parameter file for (co)variance components has a lower triangular 'long' format as "`effect_number,i,j,covar_value`". This function converts it into full symmetric (co)variance matrices.
* [rebase EBV](https://github.com/bonifazi/R_utils/blob/main/rebase_ebv.R). A function to rebase vector(s) of EBVs given a list of animals to be considered as the base population.
* [compute LR method statistics](https://github.com/bonifazi/R_utils/blob/main/compute_LR_stats.R). A function that, taking EBVs from a partial and a whole evaluation, computes LR method statistics following [Legarra and Reverter, 2018, GSE, 50:53](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0426-6). Several options are available. Main ones are that: 1) the statistics can be computed with and without providing a 'focal group' of individuals, i.e., a validation group; 2) Plotting can be used to investigate better if the validation group is homogenous. 3) Bootstrapping with replacement can be used to compute standard errors of the estimated statistics.
* ... [new functions will be added here]
