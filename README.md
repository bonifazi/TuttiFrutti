# Collection of miscellaneous R functions.  
You are free to use these R functions, for more on their usage, see the MIT [License](https://github.com/bonifazi/R_utils/blob/main/LICENSE) documentation.  
Please acknowledge its source, i.e., this repository (https://github.com/bonifazi/R_utils), so that it may be useful to others as well.
Your suggestions on improving these functions are very welcomed. Send me an email at (renzo.bonifazi@outlook.it) or open an issue here on GitHub.

## Description:
These functions are documented with the R [docstring](https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html) package.  
To view the documentation, first load the `docstring` package, then view the documentation of the function by running in the console `docstring(fun = "<functionname.R>")`.

## List of functions:
* A function to [plot convergence of MiXBLUP](https://github.com/bonifazi/R_utils/blob/main/PlotConvergeneMiXBLUP.R). `MiXBLUP` has a gnuplot code that automatically generates a convergence graph. When gnuplot is not available on your system, you can generate the same graph using this function. The output is a `ggplot2` R object.
* [Convert a MiX99 parameter file into (co)variance matrices](https://github.com/bonifazi/R_utils/blob/main/meltParfile.R). `MiX99` parameter file for (co)variance components has a lower triangular 'long' format as "`effect_number,i,j,covar_value`". This function converts it into full symmetric (co)variance matrices.
* A function to [rebase EBV](https://github.com/bonifazi/R_utils/blob/main/rebase_ebv.R) provided a list of animals in the base population.
* A function to [compute LR method statistics](https://github.com/bonifazi/R_utils/blob/main/compute_LR_stats.R). This function takes EBV from a partial and a whole evaluaton and computes LR method statistics following [Legarra and Reverter, GSE, 2018,50:53](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0426-6). Several options are available. The statistics can be computed with and without providing a 'focal group' of individuals, i.e., a validation group. Plotting can be used to investigate better if the validation group is homogenous.
* ... [new functions will be added here]
