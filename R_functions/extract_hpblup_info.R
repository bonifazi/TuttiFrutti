extract_hpblup_info <- function(hpblup_logfile) {
#' Extract Information from HPBLUP Log File
#'
#' This function extracts key information from an HPBLUP log file, such as the
#' number of iterations, min and max eigenvalues, spectral condition number,
#' and relative error in the computed solution.
#' 
#' @param hpblup_logfile A character string specifying the path to the HPBLUP log file.
#'
#' @return A data frame containing the following columns:
#'   \item{iterations_n}{The number of iterations (integer).}
#'   \item{min_eigen}{The minimum eigenvalue (numeric).}
#'   \item{max_eigen}{The maximum eigenvalue (numeric).}
#'   \item{spect_nr}{The spectral condition number (numeric).}
#'   \item{rel_err_val1}{The first value of the relative error in the computed solution (numeric).}
#'   \item{rel_err_val2}{The second value of the relative error in the computed solution (numeric).}
#'
#' @examples
#' # Assuming you have a valid HPBLUP log file path
#' hpblup_logfile <- "path/to/hpblup_log.txt"
#' result <- extract_hpblup_info(hpblup_logfile)
#' print(result)
#'
#' @importFrom stringr str_extract str_subset str_split
#' @importFrom utils readLines
#'
#' @export
############# END DOCUMENTATION ############################################
  library(stringr)
  if(file.exists(hpblup_logfile)) {
    lines <- readLines(hpblup_logfile)
  } else {
    cat("\n File ", hpblup_logfile, "does not exists.")
  }

  orig_wd <- getwd()
  setwd(dirname(hpblup_logfile))
  all_iter <- as.integer(str_extract(str_subset(lines, "All iterations \\("), "(?<=\\().+?(?=\\))"))
  eigen <- as.numeric(unlist(str_extract_all(str_subset(lines, "Min\\|max eigenvalue"), "(?:\\d+\\.|\\.\\d)\\d*(?:E[+-]?\\d+)?")))
  spectn <- as.numeric(str_extract(str_subset(lines, "^Spectral condition number"), "(?<=:\\s).*"))
  relat_err <- str_extract(str_subset(lines, "^Relative error in the computed solution"), "(?<=\\().*?(?=\\))")
  relat_err <- as.numeric(str_split(relat_err, ",")[[1]])

  values <- data.frame(
    iterations_n = all_iter, # num. of all iterations
    min_eigen =  eigen[1],   # min eigenvalue
    max_eigen =  eigen[2],   # max eigenvalue
    spect_nr = spectn,       # spectral number
    rel_err_val1 = relat_err[1], # relative error in computed solution (first value)
    rel_err_val2 = relat_err[2] # relative error in computed solution (second value)

  )
  setwd(orig_wd)
  return(values)
}