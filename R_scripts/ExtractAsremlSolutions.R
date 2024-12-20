#!/usr/bin/env Rscript
# Purpose: Extract asreml solutions and s.e. and convert them to EBV nad REL.
# Author: Renzo Bonifazi - 12-06-2023
# Usage:
# # Rscript --vanilla ExtractAsremlSolutions.R --file my_path/asreml.sln --effect_name effect5 --trait_names "Trait1, Trait2, Trait3" --varA "varA_trait1, varA_trait2, varA_trait3" --output myoutput.csv
#
# Example:
# # Rscript --vanilla ExtractAsremlSolutions.R --file Renzo/Documents/asreml.sln --effect_name Tr.pedID --trait_names "Trait1, Trait2, Trait3" --varA "1.0, 12.5, 1.12" --output myoutput.csv
#
# Details:
# Extracts the solutions (as EBV and REL) of the provided effect name (effect_name above) from the .sln file of asreml, and format it into a tabular form so that
# each animal is one row and each column is one trait EBV with corresponding REL (and genetic variance provided - varA).
# The output file is generated into the same directory as input --file and is a .csv file with an header.
#
# Versions:
# * 0.0 - initial script.
#
# Notes:
# In priciple this script could be quickly adapted to extract all info for any effect estimates in .sln file.
# I only need to change the behaviour to keep the s.e. instead of converting them to REL values and add an extra flag to switch on and off this option on specific effects.
###############################################################################################################################################################################################################
cat("Program Started"); Sys.time() # print starting time

# 0. Load packages --------
if (!require("pacman")) {install.packages("pacman", quiet = T)}
pacman::p_load(optparse, data.table, dplyr, tidyr, stringr, tidylog)

# 1. Read cmd line arguments ---------------------------------------------------------------
option_list <- list(
  make_option(c("--file"),
    type = "character", default = NULL,
    help = "Input .sln Asreml file name", metavar = "character"
  ),
   make_option(c("--effect_name"),
    type = "character", default = NULL,
    help = "Name of the effect to be tabulated, i.e. the string name as in the asreml.sln file first column", metavar = "character"
  ),
   make_option(c("--trait_names"),
    type = "character", default = NULL,
    help = "String list of trait names separated by a comma, e.g. 'Trait1, Trait2, Trait3'", metavar = "character"
  ),
    make_option(c("--varA"),
    type = "character", default = NULL,
    help = "String of genetic variance(s) for defined trait(s) separated by a comma. It MUST follow the order of --trait_names provided, e.g. '1, 20, 0.10'", metavar = "character"
  ),
  make_option(c("--output"),
    type = "character", default = NULL,
    help = "Output file name", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list) # process list of args above
opt <- parse_args(opt_parser) # make args accessible in an object

# [test usage only] opt$file <- "asreml.sln"
# [test usage only] opt$effect_name <- "Tr.pedID"
# [test usage only] opt$trait_names <- "trt1, trt2, trt3"
# [test usage only] opt$varA <- "1.0, 12.5, 1.12"
# [test usage only] opt$output <- "outfile.csv"
# [test usage only] setwd("C:/Users/Documents/mycodes/")

# check all opts are provided
if (
  any(sapply(
    c("file", "effect_name", "trait_names", "varA", "output"),
    function(opt_arg) is.null(opt[[opt_arg]]))) 
  ) {
  print_help(opt_parser)
  stop(cat("All arguments must be supplied! \n Current path is:", getwd(), "\n See the help above on required options."), call. = FALSE)
}

# 2. Vars and set working directory --------
# ___ 2.1 Set variables and wd

infile <- opt$file
effect_name <- opt$effect_name
trait_names <- strsplit(opt$trait_names, split = ", ", fixed= T)
variances <- strsplit(opt$varA, split = ", ", fixed= T)
outfile <- opt$output
setwd(dirname(infile)) # set wd to input filename

# ___ 2.2 Fixed vars
col_names <- c("Trait", "ID", "EBV", "se")

# 3. Read files --------
sln <- fread(infile,
  header = F, stringsAsFactors = F,
  strip.white = T, col.names = col_names,
  colClasses = list(character = c(1,2), double = c(3:4)) # first 2 cols as character, last 2 as double/numeric
)

# 4. Edit files or any other operations --------
tabled <- sln %>%
filter(Trait == opt$effect_name) %>% # retain only the effect requested
separate(ID, into = c("trait", "IDcode"), sep = "\\.", remove = FALSE) %>% 
mutate(trait = as.integer(trait)) %>% # convert to integer
mutate(IDcode= str_remove(IDcode, "^0+")) %>% # remove leading 0's from the IDcode
mutate(varA = as.numeric(variances[[1]][trait])) %>% # associate trait variance to each trait
mutate(trait = trait_names[[1]][trait]) %>% # associate trait names to each trait
mutate(REL = (1-((se^2) / varA)) ) %>% # compute REL as REL = 1-(se^2/varA)
pivot_wider(id_cols = IDcode, names_from = trait, values_from = c(varA, EBV, REL)) %>% # tabulate the varA, EBV, and REL for each animal ID
rename(ID = IDcode) # rename IDcode to ID to avoid possible confusion

# print some summary stats
summary(tabled)

# 6. Save output  ------
fwrite(x = tabled, file = outfile, quote = F, row.names = F, append = F, sep = ",")

cat("Program ended"); Sys.time() # print ending time
