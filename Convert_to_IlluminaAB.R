#!/usr/bin/env Rscript
#
# Purpose: 
# This script processes genotype data from a specified input file, converts genotype information into a standardized Illumina AB format, 
# and outputs the processed data to a new file. It also handles the conversion of genotyping from BOT to TOP strand if necessary, 
# and allows for custom processing based on user input options such as skipping lines in the manifest file, using multiple CPUs, 
# and formatting the output in a long or wide format.
#
# Author: Renzo Bonifazi - 10-11-2024
#
# Usage:
# # Rscript --vanilla Convert_to_IlluminaAB.R --geno_file "path/to/genotype_file.lgen" --manifest "path/to/manifest.csv" --output "output_file.txt" --longformat TRUE --cpus 4 --skip 5 --TOPonly TRUE
# 
# See 'Rscript --vanilla Convert_to_IlluminaAB.R --manual' for more information on usage.
#
# Example:
# # Rscript --vanilla Convert_to_IlluminaAB.R --geno_file "genofile.lgen" --manifest "BovineSNP50_v3_A2.csv" --output "genofile_ABcoded.lgen" --longformat TRUE --cpus 0 --skip 7 --TOPonly TRUE
#
# Details:
# This script takes as input a genotype file (.lgen format) and a manifest file (CSV format). The script merges the data from 
# these two sources, and processes the genotypes by converting them to an Illumina AB format or a 0125 format.
# Additionally, it handles the conversion of BOT genotypes to TOP genotypes if only TOP genotypes are available. The script 
# can run with different options such as using multiple CPUs for parallel processing, skipping lines in the manifest file, 
# or formatting the output in a specific way.
#
# Versions:
# * 0.0 - Initial script, functionality for merging genotype data with manifest and converting genotypes.
#
# Notes:
# - Column names are expected in input files. See --manual for more info.
# - All SNP_names must be present in the manifest file.
# - Only SNPs with a TOP or BOT IlmnStrand are kept.
# - The script uses `toupper()` to standardize SNP names and manifest names for consistency.
# - There is a placeholder for handling errors when encountering missing data in TOP/BOT genotype coding.
# - Missing genotypes in .lgen are marked as '00' in the input file and '00' in the AB coding.
#
# todo: see chatgpt suggestions on how to improve the code
###################################################################################
Sys.time() # print starting time

# 0. Load packages ----------------------------------------------------------------
check_pkgs <- function(pkg_list){
      for(pkg in pkg_list) {
        if(!require(pkg, character.only = T, quietly = T)){
            cat('Installing the ', pkg, ' package, using "install.packages(', pkg, ')')
            install.packages(pkg, repos = "http://cran.us.r-project.org")
        }
      }
    }
suppressMessages(check_pkgs(pkg_list = c("optparse", "dplyr", "data.table", "dtplyr")))
# [test usage only] if (!require("pacman")) {install.packages("pacman", quiet = T)}
# [test usage only] pacman::p_load(optparse, dplyr, data.table, dtplyr)
# [test usage only] pacman::p_load(tidylog)

# 1. Read cmd line arguments ------------------------------------------------------
option_list <- list(
  make_option(c("-g", "--geno_file"),
    type = "character", default = NULL,
    help = "Genotype file in .lgen format", metavar = "character"
  ),
  make_option(c("-m","--manifest"),
    type = "character", default = NULL,
    help = "Manifest file", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Output file name", metavar = "character"
  ),
  # optional options
  make_option(c("--TOPonly"), # optional, default if FALSE (both TOP and BOT coding are present in genotype file)
    type = "logical", default = FALSE,
    help = "Logical [TRUE or FALSE] to indicate if genotypes have only TOP genotype coding. [default \"%default\"]", metavar = "logical"
  ),
  make_option(c("--longformat"), # optional, default is TRUE, should all info be returned or only AB and 0125 cols?
    type = "logical", default = T,
    help = "Logical [TRUE or FALSE]. Should all computed columns be retuned  or only the AB and 0125 coding [default \"%default\"]", metavar = "logical"
  ),
  make_option(c("--skip"), # optional, default is none
    type = "integer", default = 0,
    help = "Number of lines to skip in manifest file. [default \"%default\"]", metavar = "integer"
  ),
  make_option(c("-c", "--cpus"), # optional, default is all
    type = "integer", default = 0, # Use 0 to automatically detect all cores [default]
    help = "Number of CPUs to use. [default \"%default\"]", metavar = "integer"
  ),
  make_option(c("--version"),
    action="store_true", default=FALSE, 
    help="Version: 0.0'."
  ),
  make_option(c("--manual"),
    action="store_true", default=FALSE, 
    help="Show manual of this Rscript by running 'Rscript <Rscript_name.R> --manual'."
  )
)

opt_parser <- OptionParser(option_list = option_list) # process list of args above
opt <- parse_args(opt_parser) # make args accessible in an object
# 1.1. Manaul -------
if (opt$manual) {
  cat("
      Usage:
      # Rscript --vanilla Convert_to_IlluminaAB.R --geno_file \"path/to/genotype_file.lgen\" --manifest \"path/to/manifest.csv\" --output \"output_file.txt\" [--longformat TRUE] [--cpus 4] [--skip 5] [--TOPonly TRUE]

      Example:
      # Rscript --vanilla Convert_to_IlluminaAB.R --geno_file \"genofile.lgen\" --manifest \"BovineSNP50_v3_A2.csv\" --output \"genofile_ABcoded.lgen\" --longformat TRUE --cpus 0 --skip 7 --TOPonly TRUE

      In the above command:
       --geno_file: [mandatory] Path to the genotype file in .lgen format. See below for expected column names.
       --manifest: [mandatory] Path to the manifest file containing SNP details. See below for expected column names.
       --output: [mandatory] Path to the output file where results will be saved.
       --longformat: [optional, default TRUE] Whether to return all computed columns or only AB and 0125 coding.
       --cpus: [optional, default 0] Number of CPUs to use for computation. 0 will detect automatically.
       --skip: [optional, default 0] Number of lines to skip in the manifest file.
       --TOPonly: [optional, default FALSE] Whether the genotype data contains only TOP coding.

      Expected column names:
      NOTE: column names are case-sensitive!
      1. Genotype file (.lgen) must contain the following columns:
         - 'SNP_name': Name identifier for each SNP (e.g., HAPMAP43437-BTA-101873 or rsID).
           NOTE: 'SNP_name' column is used to match with the 'Name' column in the manifest file.
           All SNP_names must be present in the manifest file.
         - 'genotype': The genotype for each individual at the SNP position (e.g., AA, AB, BB).

      2. Manifest file (CSV) must contain the following columns:
         - 'Name': Name identifier for each SNP (e.g., HAPMAP43437-BTA-101873 or rsID - see previous 'NOTE').
         - 'IlmnStrand': Strand information (e.g., 'TOP' or 'BOT').
         - 'Chr': Chromosome number.
         - 'MapInfo': SNP map position on the chromosome.
         - 'SNP': SNP identifier (tipically in the format of e.g., [A/G], [T/C]) which denotes the alleles at the SNP site.

      Details:
      This script converts genotype data to Illumina AB format. The genotyping data is merged with a manifest file containing SNP information. 
      The genotypes are then recoded to the AB format and optionally converted into the 0125 format.
      
      Example files:
      1. Genotype file
      BREED ID SNP_name genotype A1 A2
      Breed1 ID000001 HAPMAP43437-BTA-101873 AG A G
      Breed1 ID000002 ARS-BFGL-NGS-98142 GG G G
      Breed1 ID000003 HAPMAP34944-BES1_CONTIG627_1906 CC C C
      Breed1 ID000004 ARS-BFGL-NGS-66449 AA A A
      Breed1 ID000005 ARS-BFGL-BAC-32770 GG G G
      Breed1 ID000006 ARS-BFGL-NGS-105096 GG G G
      Breed1 ID000007 ARS-BFGL-NGS-16466 AG A G
      Breed1 ID000008 ARS-BFGL-NGS-51647 CG C G
      Breed1 ID000009 HAPMAP53946-RS29015852 AA A A

      2. Manifest
      Illumina, Inc.,,,,,,,,,,,,,,,,,
      [Heading],,,,,,,,,,,,,,,,,,
      Descriptor File Name,BovineSNP50_v3_A1.bpm,,,,,,,,,,,,,,,,,
      Assay Format,Infinium HTS,,,,,,,,,,,,,,,,,
      Date Manufactured,1/14/2016,,,,,,,,,,,,,,,,,
      Loci Count ,53218,,,,,,,,,,,,,,,,,
      [Assay],,,,,,,,,,,,,,,,,,
      IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand
      ABCA12_r2-1_T_F_2277749139,ABCA12,TOP,[A/G],0059616496,CTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC,,,UMD_3.1,2,103548215,diploid,Bos taurus,UMD3.1,1,TOP,ACTCTGGTGGATGGTTCATAATCTGCTAAGATGAATAAGTTACTGGGGAAACTGGTGCATTTATTTTAAATATAAATTATATAGTCTGTAAGATATAAAGACTGCCTAATTTATTTGAACACCATACTGATCTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC[A/G]CTCCCTTCCAGCTTACCTCAACAGCCTGAATAATTTCCTCCTGCGAGTTAACATGTCAAAATATGATGCTGCCCGACATGGTAAAGTTATTTACATAGGAGCTCCTTGTATTGAAACTCTTGCTACTCTCCATGTGAAAATATACATTAGACCCCATTTTCCTCCCTGTGGCAGCTAT,ACTCTGGTGGATGGTTCATAATCTGCTAAGATGAATAAGTTACTGGGGAAACTGGTGCATTTATTTTAAATATAAATTATATAGTCTGTAAGATATAAAGACTGCCTAATTTATTTGAACACCATACTGATCTTGTCTTCTTTTGGAATGTTACAGGTATGGTATGATCCAGAAGGCTATC[A/G]CTCCCTTCCAGCTTACCTCAACAGCCTGAATAATTTCCTCCTGCGAGTTAACATGTCAAAATATGATGCTGCCCGACATGGTAAAGTTATTTACATAGGAGCTCCTTGTATTGAAACTCTTGCTACTCTCCATGTGAAAATATACATTAGACCCCATTTTCCTCCCTGTGGCAGCTAT,1241,3,-
      APAF1_dup-1_B_F_2327661418,APAF1,BOT,[T/C],0041654401,ATATTGTGCAACTGGGCCTCTGTGAACTGGAAACTTCAGAGGTTTATCGG,,,UMD_3.1,5,63150400,diploid,Bos taurus,UMD3.1,1,BOT,CCATTTCCTAATATTGTGCAACTGGGCCTCTGTGAACTGGAAACTTCAGAGGTTTATCGG[T/C]AAGCTAAGCTGCAGGCCAAGCAGGAGGTCGATAACGGAATGCTTTACCTGGAGTGGGTGT,ACACCCACTCCAGGTAAAGCATTCCGTTATCGACCTCCTGCTTGGCCTGCAGCTTAGCTT[A/G]CCGATAAACCTCTGAAGTTTCCAGTTCACAGAGGCCCAGTTGCACAATATTAGGAAATGG,1241,3,+
      ")
  quit(save = "no")  # Exits the script after printing help
}

# 1.2. Check input files -------
if (is.null(opt$geno_file) | is.null(opt$manifest)) {
  print_help(opt_parser)
  stop(cat("Input genotype file and manifest file must be supplied! \n Current path is:", getwd(), "\n"), call. = FALSE)
}

# 3. Functions -------------------------------------------------------------------
# Function to convert from BOT to TOP manifest SNPs
convert_BOT_to_TOP <- function(line) {
  tryCatch({
    if(line$IlmnStrand=="BOT") {
      line[4:11] <- illumina_df[
          as.integer(illumina_df[illumina_df$SNP==line$SNP & illumina_df$TOP_BOT=="BOT",'BOTtoTOPindex']) # index row conversion
          ,] # all other cols (SNP TOP_BOT A B G0 G1 G2 BOTtoTOPindex)
    return(line)
    } else {line}
 }, 
  error = function(e) {
    message("An error occurred: ", e$message)
    print(line)  # Print the line causing the error
    return(NULL)  # Optionally handle the error
  },
  warning = function(w) {
    message("A warning occurred: ", w$message)
    print(line) # Print the line causing the warning
    return(NULL) 
  })
}

# 4. Set vars -------------------------------------------------------------------
cat("\n >>> Setting vars <<< \n")
geno <- opt$geno_file
manifest_file <- opt$manifest
outfile <- opt$output
longformat <- opt$longformat
setDTthreads(opt$cpus)

# [test usage only] geno <- "genofile.lgen"
# [test usage only] manifest_file <- "BovineSNP50_v3_A2.csv"
# [test usage only] outfile <- "test.txt"
# [test usage only] longformat <- TRUE
# [test usage only] opt$cpus <- 0
# [test usage only] opt$skip <- 7
# [test usage only] opt$TOPonly <- TRUE
# [test usage only] setwd("C:/Documents/mycodes/")

# Create the data frame using a matrix to fill it row-wise
illumina_df <- data.frame(
  matrix(
    c(
        #   SNP  TOP/BOT      A      B      G0     G1     G2   BOTtoTOPindex
        "[A/C]",   "TOP",   "A",   "C",   "AA",  "AC",  "CC",   NA, 
        "[C/A]",   "TOP",   "A",   "C",   "AA",  "AC",  "CC",   NA,
        "[A/G]",   "TOP",   "A",   "G",   "AA",  "AG",  "GG",   NA,
        "[G/A]",   "TOP",   "A",   "G",   "AA",  "AG",  "GG",   NA,
        "[T/C]",   "BOT",   "T",   "C",   "TT",  "TC",  "CC",    3,
        "[C/T]",   "BOT",   "T",   "C",   "TT",  "TC",  "CC",    3,
        "[T/G]",   "BOT",   "T",   "G",   "TT",  "TG",  "GG",    1,
        "[G/T]",   "BOT",   "T",   "G",   "TT",  "TG",  "GG",    1,
        "[A/T]",   "TOP",   "A",   "T",   "AA",  "AT",  "TT",   NA,
        "[T/A]",   "TOP",   "A",   "T",   "AA",  "AT",  "TT",   NA,
        "[G/C]",   "BOT",   "G",   "C",   "GG",  "CG",  "CC",   15,
        "[C/G]",   "BOT",   "G",   "C",   "GG",  "CG",  "CC",   15,
        "[A/T]",   "BOT",   "T",   "A",   "TT",  "AT",  "AA",    9,
        "[T/A]",   "BOT",   "T",   "A",   "TT",  "AT",  "AA",    9,
        "[G/C]",   "TOP",   "C",   "G",   "CC",  "CG",  "GG",   NA,
        "[C/G]",   "TOP",   "C",   "G",   "CC",  "CG",  "GG",   NA
    ),
    ncol = 8,
    byrow = TRUE
  )
)
colnames(illumina_df) <- c("SNP", "TOP_BOT", "A", "B", "G0", "G1", "G2", "BOTtoTOPindex")
cat("Illumina conversion table: \n"); print(illumina_df)

# 5. Read data -----------
cat("\n >>> reading data <<< \n")
cat("\n >>> reading manifest <<< \n")
manifest  <- read.csv(manifest_file,
                        skip = if(!is.null(opt$skip)) {opt$skip} else {0},
                        header = T)
head(manifest[1:10,1:10])
cat("\n Column names of manifest file:", names(manifest))
cat("\n Dimension of manifest file:", dim(manifest))

cat("\n >>> reading .lgen <<< \n")
lgeno <- fread(geno, header = T, showProgress = T)
head(lgeno)
cat("\n Dimension of genotype file:", dim(lgeno))

# 6. Merge geno and manifest -------------------------------------------------------------------
# 6.1 Expand the manifest with Illumina Table
cat("\n >>> Adding Illumina codes to manifest <<< \n")
manifest_exp <- manifest %>%
  mutate(SNP_name = toupper(Name)) %>% # for consistency in matching
  select(SNP_name, Chr, MapInfo, SNP, IlmnStrand) %>%
  left_join(
    illumina_df,
    by = c("SNP" = "SNP", "IlmnStrand" = "TOP_BOT")
    )
# head(manifest_exp)

cat("\n TOP/BOT combinations in manifest:\n")
manifest_exp %>% 
  filter(IlmnStrand %in% c("TOP", "BOT")) %>% 
  count(SNP, IlmnStrand)

# 6.2 if only TOP genotypes are present, convert BOT coding to TOP coding in table ---
if(opt$TOPonly == TRUE) {
  cat("\n >>> Converting BOT to TOP in manifest <<< \n")
  manifest_exp <- do.call(rbind, lapply(1:nrow(manifest_exp), function(i) {
    convert_BOT_to_TOP(manifest_exp[i, ])}
  ))
}
if(opt$TOPonly == TRUE & any(manifest_exp$IlmnStrand=="BOT")) {
  stop("'--TOPonly' option is set to 'TRUE' but there are still BOT variants in manifest:\n",
  "BOT variants count: ", print(nrow(filter(manifest_exp, IlmnStrand=="BOT"))
  ))
}

# 6.3 Merge geno with expanded manifest ---
cat("\n >>> Merging .lgen and manifest info <<< \n")
merged <- lgeno %>% 
    mutate(SNP_name = toupper(SNP_name)) %>% 
    left_join(
        manifest_exp,
        by = "SNP_name")

# 6.4. Print some summary stats on provided data ---
cat("\n Number of SNPs per Chromosome: \n")
print(count(merged, Chr))
cat("\n Number of SNPs per SNP-coding and TOP/BOT Strand")
print(count(merged, IlmnStrand, SNP))

# 7. Recode genotypes as AB or 0125 -------------------------------------------------------------
cat("\n >>> Recoding genotypes to AB and 0125 <<< \n")
# 7.1 Check that all SNPs information are linked. If any is missing, return error with some info ------
check_info <- merged %>%
  filter(is.na(A) | is.na(B) | is.na(G0) | is.na(G1) | is.na(G2)) 
if(nrow(check_info) > 0) {
  print("\n SNP information absent/present:\n")
  print(
    merged %>%
      count(is.na(A), is.na(B), is.na(G0), is.na(G1), is.na(G2), name = "line_count")
  )
  error_message <- paste0(
    "There are missing information to convert the genotype for ", nrow(check_info), " rows.\n",
    "Number of unique SNPs affected: ", length(unique(check_info$SNP_name)), "\n",
    "First 10 SNPs with missing info:\n", paste0(head(unique(check_info$SNP_name), 10), collapse = "\n"), "\n",
    "First 10 lines with missing info:\n", 
    paste0(capture.output(head(check_info)), collapse = "\n"), "\n",
    "Execution halted."
  )
  stop(error_message)
}

# 7.2 Convert genotypes into AB format --------------------------------------------------------------
cat("\n >>> Recoding genotypes to AB <<< \n")
setDT(merged)
merged[, AB_genotype := {
  fcase( 
    genotype == paste0(A, A), "AA", # If genotype is AA (both alleles are A)
    genotype == paste0(B, B), "BB", # If genotype is BB (both alleles are B)
    genotype == paste0(A, B), "AB", # If genotype is AB (A and B in any order)
    genotype == paste0(B, A), "AB", # If genotype is BA (same as AB, order doesn't matter)
    default = "00"                  # If genotype doesn't match the above patterns
  )
}]

# 7.3 Convert AB coding into 0125 -----------------------------------------------------------------------
cat("\n >>> Recoding genotypes to 0125 <<< \n")
merged[, geno_0125 := fcase(
  AB_genotype == "AA", 0,  # If AB_genotype is AA, set 0125_genotype to 0
  AB_genotype == "BB", 2,  # If AB_genotype is BB, set 0125_genotype to 2
  AB_genotype == "AB", 1,  # If AB_genotype is AB, set 0125_genotype to 1
  default = 5    # If none of the above, set to 5
)]

# 7.4 Print some count stats  -----------------------------------------------------------------------
cat("\n Count of genotypes and their AB and 0125 coding: \n\n")

merged %>% count(genotype, IlmnStrand, AB_genotype, geno_0125)

# 8. Write to disk --------------------------------------------------------------------------------------
cat("\n >>> Writing to disk <<< \n")
if(longformat == T) {
  fwrite(merged, outfile, quote = F, col.names = T, row.names = F, sep = " ")
} else if (longformat == F) {
  fwrite(
    merged %>% select(all_of(c(names(lgeno), "AB_genotype", "geno_0125"))),
    outfile, quote = F, col.names = T, row.names = F, sep = " ")
}

Sys.time() # print end time