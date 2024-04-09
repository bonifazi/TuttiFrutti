#!/usr/bin/env Rscript
#
# Purpose: Plot convergence of postgibbsf90.
# Author: Renzo Bonifazi - 08-04-2024
# Usage:
# # Rscript --vanilla plot_postgibbsf90.R --file my_path/postgibbs_samples --output postgibbs_plots.pdf --traits BRD1_AWW, BRD2_AWW, BRD1_CE, BRD2_CE
#
# Details:
# * you first have to run postgibbsf90 to produce the input file 'postgibbs_samples'
# * '--traits' are provided in a comma seperated list
# * by default output is in the same working directory as the input file
# * output is a .pdf file with one page for variances and one for covariances 
#   unless 'trait_groups' is provided. In this case, (co)variances are further divided in blocks
#   and blocks of (co)variances are plotted in separted pages. This is useful with large (co)variances.
#   The (co)variance matrix is divided in NxN blocks, with N = --trait_groups.
# 
# Versions:
# * 0.0 - initial script.
#
###################################################################################
Sys.time() # print starting time

# 0. Load packages --------
if (!require("pacman")) {install.packages("pacman", quiet = T)}
pacman::p_load(optparse, dplyr, tidyr, data.table, stringr, ggplot2)
# [test usage only] pacman::p_load(tidylog)

# 1. Read cmd line arguments ---------------------------------------------------------------
# args parser # https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
option_list <- list(
  make_option(c("-f", "--file"),
    type = "character", default = 'postgibbs_samples',
    help = "Input file name [default: 'postgibbs_samples']", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = 'postgibbs_samples_plots.pdf',
    help = "Output file name [default: 'postgibbs_samples_plots.pdf']", metavar = "character"
  ),
  make_option(c("-t", "--traits"),
    type = "character", default = NULL,
    help = "Comma separated name for traits [default 't1, t2 .. tn']", metavar = "character"
  ),
  make_option(c("-g", "--trait_groups"),
    type = "integer", default = NULL,
    help = "Number of traits groups to form [ default is none]", metavar = "integer"
  )
)

opt_parser <- OptionParser(option_list = option_list) # process list of args above
opt <- parse_args(opt_parser) # make args accessible in an object

# [test usage only] Set variables --------
# [test usage only] opt <- list()
# [test usage only] opt$file <- "postgibbs_samples"
# [test usage only] opt$output <- "postgibbs_samples.pdf"
# [test usage only] opt$traits <-  "T1_B1, T1_B2, T1_B3, T2_B1, T2_B2, T2_B3"
# [test usage only] opt$trait_groups <- 2
# [test usage only] setwd("C:/Users/Documents/mycodes/")

# 2. Check input files, set wd, dead files --------
setwd(dirname(opt$file)) # set wd to same as input file
if (!file.exists(opt$file)) {
  print_help(opt_parser)
  stop(cat("Input file does NOT exists! \n Current path is:", getwd(), "\n"), call. = FALSE)
}

file_proc <- read.table(opt$file)
# [test usage only] file_proc[, 25:33] <- NULL

# extract num of cols with (co)vars from input file
n_cols <- ncol(file_proc)-3
# extract num of traits (covar is a squared upper-tri with diag)
n_traits <- (-1 + sqrt(1 + 8 * n_cols)) / 2

# trait labels
if(is.null(opt$traits)) {
    tr <- paste0("trt_", 1:n_traits) # default labels
} else {
    tr <- strsplit(opt$traits, ",")[[1]]                 # user-provided labels
}

# 3. Define variables ------------------------------------
# [deprecated]
# # find the column position for variances
# vars_positions <- vector()
# for (i in 1:n_traits) {
#     if(i==1){
#         vars_positions[i] <- 1
#         } else {
#         vars_positions[i] <- vars_positions[i-1]+(n_traits- (i-2))
#     }
# }

# define labels for plots
labels <- data.frame()
for(i in 1:n_traits){
    for(j in 1:n_traits){
        if(j>=i) {
        # cat(i , "-",  j, "\n")
        if(j==i) {
            labels_n <- data.frame(
                label = paste0("VAR_", tr[i]),
                index_i = i,
                index_j = j
            )
        }
        if(j!=i) {
            labels_n <- data.frame(
                label = paste0("COVAR_", tr[i], "-", tr[j]),
                index_i = i,
                index_j = j
            )
        }
        labels <- rbind(labels, labels_n)
        }
    }
}

if(!is.null(opt$trait_groups)) {
    if(n_traits %% opt$trait_groups != 0) {
        stop ("--trait_groups is not a multiple of number of traits, n_traits:", n_traits, "--trait_groups = ", opt$trait_groups)
    }
    block_size <- n_traits/ opt$trait_groups
    index_blocks <-  as.vector(which(1:n_traits %% block_size== 0))

    # define block size
    range_block <- data.frame()
    for(i in 1:length(index_blocks)) {
        if(i == 1) {
        range_block_n = 1:index_blocks[i]
        } else {
        range_block_n = (index_blocks[i-1]+1) : index_blocks[i]
        }
        range_block <- rbind(range_block, range_block_n)
    }
    
    # assign blocks's value
    block_indeces = data.frame()
    counter=0
    for(i in 1:nrow(range_block)){
        for(j in 1:nrow(range_block)){
            if(j>=i){
            counter = counter+1
            block_n = cbind(range_block[i,], range_block[j,], counter)
            block_indeces = rbind(block_indeces, block_n)
            
            }
        }
    }

    # add block count to the labels
    labels$block <- NA
    for(i in 1:nrow(labels)) {
        row <- intersect(
        as.data.frame(which(labels[i,]$index_i == block_indeces[, 1:((ncol(block_indeces)-1)/2)], arr.ind = TRUE))[,1],
        as.data.frame(which(labels[i,]$index_j == block_indeces[, (((ncol(block_indeces)-1)/2)+1):(ncol(block_indeces)-1)], arr.ind = TRUE))[,1]
    )
        labels[i,]$block <- block_indeces[row,]$counter
    }
    
} else {
    labels$block <- 1
}

# 4. Plot and save pdf --------
plot_data <- file_proc[, -c(1,3)] %>% # remove not-needed cols
    setnames(c("iteration", labels$label)) %>% 
    pivot_longer(cols  = -iteration, values_to = "covar_value", names_to = "label") %>% 
    left_join(labels %>% select(label, block), by="label") %>% # add the block
    mutate(main_page = paste0(str_extract(label, "^[^_]+"), "_", block))

pdf(opt$output, width = 14)
for(page in unique(plot_data$main_page)) {
    
    p <- plot_data %>% 
    filter(main_page == page) %>% 
    ggplot(aes(x = iteration, y = covar_value, colour = label)) +
        geom_line() +
        theme_bw() +
        labs(y="") +
        labs(colour = "") +  # Remove legend title
        theme(legend.position = "bottom")  # Position legend at the bottom
    print(p)

}
dev.off()

Sys.time() # print ending time
