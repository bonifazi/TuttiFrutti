#!/usr/bin/env Rscript
#
# Purpose: Sort a pedigre file from oldest to youngest animals.
# Author: Renzo Bonifazi - 20-12-2024
# Usage:
# # Rscript --vanilla sort_ped.R --ped ped.txt [--sep SEP] --output sorted_ped.txt (default: --sep " " (empty space))
#
# Details:
# The Rscript sort the pedigree so that partents appear before the offspring.
# For now, no date of birth or other checks for now.
# Versions:
# * 0.0 - initial script.
#
# Details:
# Input file: a pedigree file without header having
# animal ID, sire ID, and dam ID, in col1, col2, col3, respectevly.
# missing code is assumed to be 0.
# Parents that do not appear as animals will be added with missing code = 0.
# Other cols after col3 may be present and will be deleted
#
# Notes:
# TO DO:
# can add the option to return renumbered 1:n pedigree form purgeR pkg (now keeping original IDs)
###################################################################################
Sys.time() # print starting time

# 0. Load packages --------
if (!require("pacman")) {install.packages("pacman", quiet = T)}
pacman::p_load(dplyr, data.table, optparse, purgeR)
# [test usage only] pacman::p_load(tidylog)

# 1. Read cmd line arguments ---------------------------------------------------------------
# args parser # https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
option_list <- list(
  make_option(c("-p", "--ped"),
    type = "character", default = NULL,
    help = "Input file name", metavar = "character"
  ),
  make_option(c("-s", "--sep"),
    type = "character", default = " ",
    help = "Output field separator", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Output file name", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list) # process list of args above
opt <- parse_args(opt_parser) # make args accessible in an object
if (is.null(opt$ped) | is.null(opt$output)) {
  print_help(opt_parser)
  stop(cat("Input and output file name must be supplied! \n Current path is:", getwd(), "\n"), call. = FALSE)
}

# 2. Vars and set working directory --------
# ___ 2.1 Set variables and wd
infile   <- opt$ped
file_sep <- opt$sep
outfile  <- opt$output

# [test usage only]
# infile <- "PBD_fam.ped"
# outfile <- "test_sort.ped"
# file_sep <- " "

# ___ 2.2 Fixed vars
# col_names <- c("id", "sire", "dam")
'%not_in%' <- Negate('%in%')

# 3. Read files --------
file_proc <- fread(infile, stringsAsFactors = F, header = F) %>%
  as.data.frame() %>% .[, c(1:3)]

# find missing sires and dams
sires_to_add <- file_proc[, 2] %>% 
  as.data.frame() %>% 
  distinct() %>% 
  filter(. %not_in% file_proc[, 1]) %>%  # Check that Sire is not already in Animal ID
  filter(. != 0)

dams_to_add <- file_proc[, 3] %>% 
  as.data.frame() %>% 
  distinct() %>% 
  filter(. %not_in% file_proc[, 1]) %>%  # Check that Dam is not already in Animal ID
  filter(. != 0)

# Combine both missing sires and dams into a single data frame
parents_to_add <- rbind(
  cbind(sires_to_add[,1], 0, 0),
  cbind(dams_to_add[,1],  0, 0)
)
# append parents to add to the ped file
ped_to_sort <- rbind(
  parents_to_add,
  file_proc
) %>% distinct()

# 4. Sort ped --------
pacman::p_load(purgeR)

order <- ped_sort(
 as.data.frame(ped_to_sort),
 id   = colnames(ped_to_sort)[1],
 sire = colnames(ped_to_sort)[2],
 dam  = colnames(ped_to_sort)[3],
 keep_names = T)

sorted_ped <- ped_to_sort[match(order$names, ped_to_sort[,1]), ]

# 5. Save output  ------
fwrite(x = sorted_ped, file = outfile,
  quote = F, append = F,
  col.names = F, row.names = F,
  sep = file_sep)

Sys.time() # print ending time