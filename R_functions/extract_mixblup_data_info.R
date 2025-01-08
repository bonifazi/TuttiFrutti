extract_mixblup_data_info <- function(inpfile, keywords) {
#' Extract file information from a MiXBLUP instruction file
#'
#' This function extracts file paths and counts the number of lines in the files associated with given keywords
#' from a MiXBLUP instruction file. The result is returned as a data frame containing the provided keywords,
#' file paths, and the corresponding number of lines in each file.
#'
#' @param inpfile A character string specifying the path to the MiXBLUP instruction file.
#' @param keywords A character vector containing the keywords to search for in the MiXBLUP instruction file. 
#'
#' @return A data frame with three columns:
#'   \describe{
#'     \item{keyword}{A character vector containing the keywords provided in the `keywords` argument.}
#'     \item{file}{A character vector with the file paths extracted for each keyword}
#'     \item{lines}{A numeric vector containing the number of lines in each file}
#'   }
#' 
#' @details
#' The function searches for each keyword in the file, extracts the corresponding file path, and counts the number of lines in the files. 
#' The function assumes that provided keywords are defined in the instruction file and they are followed by a file path 
#' in the format `keyword <file_path>`.
#' For example, for the keyword "DATAFILE", the instruction file should contain a line such as `DATAFILE <path_to_data_file>`.
#'
#' @examples
#' \dontrun{
#' # Assuming 'mixblup_instructions.txt' is the MiXBLUP instruction file
#' # and it contains keywords like 'DATAFILE', 'PEDFILE', etc.
#' keywords <- c("DATAFILE", "PEDFILE", "ERMFILE")
#' result <- extract_mixblup_data_info("mixblup.inp", keywords)
#' print(result)
#' }
#'
#' @importFrom stringr str_extract str_subset
#' @export
############# END DOCUMENTATION ############################################
  library(stringr)
  # Check if the file exists
  if (file.exists(inpfile)) {
    lines <- readLines(inpfile)  # Read lines from the input file
  } else {
    cat("\n File ", inpfile, "does not exist.")  # Print error message if file doesn't exist
  }
  
  orig_wd <- getwd() # save and move to inst file wd (to account for relative paths)
  setwd(dirname(inpfile))

  # Initialize a list to store the results
  file_info <- list()

  # Iterate through each keyword in the provided vector
  for (keyword in keywords) {
    # Create the regular expression for the given keyword
    pattern <- paste0("(?<=^", keyword, "\\s)([^\\s]+)")
    
    # Extract the file path for the current keyword
    file_path <- str_extract(str_subset(lines, paste0("^", keyword)), pattern)
    
    # If file exists, count the number of lines
    if (length(file_path) > 0) {
      file_lines <- as.integer(system(paste0("wc -l < ", file_path), intern = TRUE))
    } else {
      stop("ERROR: file for keyword ", keyword, " not found! Check your keyword and MiXBLUP instruction file.")
    }
    
    # store file path and line count
    file_info[[keyword]] <- list(file = file_path, lines = file_lines) 
  }
  
  # Convert the result list into a data frame
  results_df <- data.frame(
    keyword = rep(keywords, each = 1),
    file = sapply(file_info, function(x) x$file),
    lines = sapply(file_info, function(x) x$lines),
    row.names = NULL
  )
  
  setwd(orig_wd) # go back to original wd
  return(results_df)
}