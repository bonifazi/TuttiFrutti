#!/usr/bin/env Rscript
# Purpose: This scripts produces .csv and .pdf report for Plink ROH.
# Author: Renzo Bonifazi - 26-04-2023
# Usage:
# Rscript --vanilla Analyse_Plink_ROH.R --plink_files plinkcleaned --plink_roh ROH.hom --group geno_BRD --pedigree ped.ped --output results_dir
#
# Example:
# # Rscript --vanilla Make_pedigree.R --file Renzo/Documents/file.csv --output Output.csv
#
# Details:
#
# Versions:
# * 0.0 - initial script.
#
# Notes:
#
###################################################################################
Sys.time() # print starting time

# 0. Load and install any missing package --------
check_pkgs <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = T)) {
      cat("Installing the ", pkg, ' package, using "install.packages(', pkg, ")")
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
  }
}
check_pkgs(pkg_list = c("optparse", "dplyr", "detectRUNS", "ggplot2", "tools", "tidyr"))

# [test usage only] pacman::p_load(optparse, dplyr, detectRUNS, ggplot2, tools, tidyr)
# [test usage only] pacman::p_load(tidylog)

# 1. Read cmd line arguments ---------------------------------------------------------------

# args parser # https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
option_list <- list(
  make_option(c("--plink_files"),
              type = "character", default = NULL,
              help = "Plink filename for .fam and .bim files (without providing the .fam and .bim suffix)", metavar = "character"
  ),
  make_option(c("--plink_roh"),
              type = "character", default = NULL,
              help = "Plink ROH file", metavar = "character"
  ),
  make_option(c("--group"),
              type = "character", default = NULL,
              help = "Plink ROH file", metavar = "character"
  ),
  make_option(c("--pedigree"),
              type = "character", default = NULL,
              help = "Pedigree file", metavar = "character"
  ),
  make_option(c("--output_dir"),
              type = "character", default = NULL,
              help = "Output directory name", metavar = "character"
  ),
  # optionally user-defined interval for YOB plotting (default range is taken from data)
  make_option(c("--min_YOB_plot"),
              type = "integer", default = NULL,
              help = "OPTIONAL: Minimum Year of Birth used in plots. Default range taken from data.", metavar = "character"
  ),
  make_option(c("--max_YOB_plot"),
              type = "integer", default = NULL,
              help = "OPTIONAL: Maximum Year of Birth used in plots. Default range taken from data.", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = T, 
                           usage = "\nRscript --vanilla Analyse_Plink_ROH.R --plink_files plinkcleaned --plink_roh ROH.hom --group geno_BRD --pedigree ped.ped --output results_dir") # process list of args above
opt <- parse_args(opt_parser) # make args accessible in an object
check_opt <- c("plink_files", "plink_roh", "group", "pedigree", "output_dir")
if (any(unlist(lapply(opt[check_opt], is.null)))) {
  print_help(opt_parser)
  stop("\n\n Check that all required inputs are provided using the correct flags. Required flags are:", paste0(" --", check_opt, " "))
}

# 2. Vars and set working directory --------
# ___ 2.1 Set variables and wd
plink_files <- opt$plink_files
roh_file <- opt$plink_roh
group_label <- opt$group
pedigree <- opt$pedigree
output_dir <- opt$output_dir

# [test usage only] plink_files <- "plink_detectRUNS"
# [test usage only] roh_file <-    "ROH.hom"
# [test usage only] pedigree <-    "ped_Retriever.txt"
# [test usage only] output_dir <-  "C:/Users/bonif002/Downloads/ROH"
# [test usage only] setwd(output_dir)

# extend plink files to .fam and .bim
fam_file_name <- paste0(plink_files, ".fam")
bim_file_name <- paste0(plink_files, ".bim")

# set wd to input filename
setwd(dirname(output_dir))

# ___ 2.2 Fixed vars

# Define output filenames
ROH_pdf_file <- file.path(output_dir, "ROH_report.pdf")
F_ROH_chromo_wide_csv_file <- file.path(output_dir, "FROH_chromo_wide.csv")
F_ROH_genome_wide_csv_file <- file.path(output_dir, "FROH_genome_wide.csv")

# 3. Read input files ----------------------------------------------------------
# NOTE: {detectRUNS} reads by itself the file, you only provide the name string
# read .bim file. For .bim file {detectRUNS} expects only the first 4 out of the 6 cols
if( ncol(read.table(bim_file_name))>4 ) {
  stop("Number of columns in ", bim_file_name, "are > 4, this will create an error because For .bim file {detectRUNS} expects only the first 4 out of the 6 cols. \n Manually retain only col1 to 4.")
}
# bim_file <- read.table(bim_file)[,c(1:4)]

# NOTE: {detectRUNS} reads by itself the file, you only provide the name string
# read .fam file. For .fam file {detectRUNS} expects the grouping col in the 1st column.
# fam_file <- read.table(fam_file_name,
#                        col.names = c("group", "id", "parent1", "parent2", "col5", "col6"),
#                        colClasses = rep("character", 6),
#                        stringsAsFactors = F)

# n. ch from .bim file
n_ch <- max(read.table(bim_file_name)[, 1])


# read pedigree
pedigree <- read.table(pedigree,
                       sep = " ", header = F,
                       col.names = c("ID", "sire", "dam", "YOB"),
                       colClasses = rep("character", 4), stringsAsFactors = F
)

# Read ROHs from plink
runs <- readExternalRuns(roh_file, program = "plink")
# edit the group label
cat("\n There are :", length(unique(runs$group)), "groups in the ROH file provided.")
warning("\nBy default the group is set to the '--group' labelled provided --> ", group_label, ". This assumes there is only 1 group (i.e., only one population)")
# the group col is used for grouping ROH results. Assume there is only one group
runs <- runs %>% mutate(group = group_label)
cat(length(unique(runs$group)))
cat("\n There are :", length(unique(runs$group)), "groups in the ROH file provided.")
cat("---- Finished reading input files ----")

# 4. Run calculations ----------------------------------------------------------
# ___ 4.1 Run calculations -----------------------------------------------------
cat("---- Start ROH stats calculations in detectRUNS ----")
summaryList <- summaryRuns(
  runs = runs, mapFile = bim_file_name, genotypeFile = fam_file_name,
  Class = 6, snpInRuns = TRUE
)

summaryList$summary_ROH_count # ROH counts
summaryList$summary_ROH_mean_chr # ROH size stats
# head(summaryList$SNPinRun) # proportion of times an SNP falls inside a run in any given pop/group

# ___ 4.2 Make plots -----------------------------------------------------------
pdf(file = ROH_pdf_file)
# Plot ROHs
plot_InbreedingChr(runs = runs, mapFile = bim_file_name, groupSplit = T, style = "All")
plot_DistributionRuns(runs = runs, mapFile = bim_file_name, groupSplit = T, style = "All")
# Manhattan plot ROHs
plot_manhattanRuns(runs, genotypeFile = fam_file_name, mapFile = bim_file_name)
plot_Runs(runs)
plot_StackedRuns(runs)
for (chr in 1:n_ch) {
  print(paste("Chromosome", chr))
  plot_SnpsInRuns(
    runs = runs[runs$chrom == chr, ],
    genotypeFile = fam_file_name, mapFile = bim_file_name
  )
}
dev.off()

# ___ 4.3) Calculate F_ROHs ----------------------------------------------------
# ___ 4.3.1) Calculate F_ROHs chromosome-wide ----------------------------------
F_ROH_ChomosomeWide <- Froh_inbreeding(runs = runs, mapFile = bim_file_name, genome_wide = F)
# F_ROH_ChomosomeWide <- Froh_inbreedingClass(runs = runs, mapFile = mapfile)
write.csv(x = F_ROH_ChomosomeWide, F_ROH_chromo_wide_csv_file, quote = F, row.names = F)
# dim(F_ROH_ChomosomeWide)
# F_ROH_ChomosomeWide[1:10, 1:10]

# ___ 4.3.2) Calculate F_ROHs genome-wide --------------------------------------
# Calculate F_ROHs genome-wide
F_ROH_GenomeWide <- Froh_inbreeding(runs = runs, mapFile = bim_file_name, genome_wide = T)
write.csv(x = F_ROH_GenomeWide, file = F_ROH_genome_wide_csv_file, quote = F, row.names = F)

# ___ 4.3.3) Calculate F_ROHs genome-wide and add averages per YOB and delta_F  ----------------------------------
# extract min and max YOB from pediree
min_YOB_plot <- min_YOB <- as.integer(min(pedigree$YOB))
max_YOB_plot <- max_YOB <- as.integer(max(pedigree$YOB))

F_ROH_GenomeWide_YOB_grouped <- F_ROH_GenomeWide %>%
  left_join(pedigree %>% select(ID, YOB), by = c("id" = "ID")) %>% # add info on YOB
  group_by(group, YOB) %>% # take average per YOB
  summarise(
    mean_fROH = mean(Froh_genome), # average F_roh per year-group
    n_animals = n() # count animals per year-group
  ) %>%
  mutate(LN_1_minus_F = log(1 - mean_fROH)) %>% # LN(1-F)
  mutate(YOB = as.integer(YOB)) %>%
  complete(YOB = full_seq(min_YOB:max_YOB, period = 1)) %>% # add all years in range 1920-2020, misisng years simply get NAs
  relocate(group, YOB, n_animals, mean_fROH, LN_1_minus_F) # order columns

write.csv(
  x = F_ROH_GenomeWide_YOB_grouped,
  file = paste0(file_path_sans_ext(F_ROH_genome_wide_csv_file), "_YOB_grouped.csv"),
  quote = F, row.names = F
)
# head(F_ROH_GenomeWide)
# head(F_ROH_ChomosomeWide)

# ___ 4.3.4) make R_ROH plots by YOB
pdf(file.path(output_dir, "F_ROH_plots.pdf"),width = 13)
for(column in colnames(F_ROH_GenomeWide_YOB_grouped)[c(3:4)]){ # n of animals and F_ROH
  p <- F_ROH_GenomeWide_YOB_grouped %>%
    select(group, YOB, any_of(column)) %>%
    pivot_longer(cols = -c(group, YOB),
                 names_to = "F_param", values_to = "value") %>%
    filter(YOB > min_YOB_plot & YOB < max_YOB_plot) %>%
    ggplot(aes(x = YOB, y = value)) +
    geom_line() + geom_point()+
    labs(title = column) +
    scale_x_continuous(breaks = seq(min_YOB_plot, max_YOB_plot, 2),minor_breaks = waiver()) +
    theme(axis.text.x = element_text(angle = +90))
  print(p)
}
dev.off()

cat("\n\n ++++++++++++++++++++++++++++++\n")
Sys.time() # print ending time
cat("+++++ Program finished  +++++ \n")