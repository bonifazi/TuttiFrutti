# R_utils
Collection of miscellaneous R functions.  
These functions are documented with the R [docstring](https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html) package.  
To view the documentation, first load the `docstring` package, then view the documentation of the function by running in the console `docstring(fun = "<functionname.R>")`.

List of functions:
* A function to [plot convergence of MiXBLUP](https://github.com/bonifazi/R_utils/blob/main/PlotConvergeneMiXBLUP.R). `MiXBLUP` has a gnuplot code that automatically generates a convergence graph. When gnuplot is not available on your system, you can generate the same graph using this function. The output is a `ggplot2` R object.
* [Melt a MiX99 parameter file into (co)variance matrices](https://github.com/bonifazi/R_utils/blob/main/meltParfile.R). `MiX99` parameter file for (co)variance components has a lower triangular 'long' format as `effect,i,j,covar`. This function converts it into full symmetric (co)variance matrices.
* ...
