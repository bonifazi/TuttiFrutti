# Purpose: general function to compute validation statistics, plot results (and convergence for NextGP). 
#
# Author: Renzo Bonifazi
# Date: 07-08-2023
# TO DO: add documentation
############################################################
# 1) Load packages ---------------
if (!require("pacman", character.only = T)) {
          stop('Please install the pacman package, using install.packages("pacman"),
                or change the way packages are loaded in this Rscript', call. = F)
        }
pacman::p_load(boot, parallel, ggplot2, patchwork, dplyr, tidylog)

# 2) Define functions ------------

# function to return prediction results (including bootstrapping for SE)
compute_prediction_results <- function(
    data_x_axis,
    data_y_axis,
    VG_label, # VG_label for result (e.g. "validation" and "training")
    bootstrap = T,
    boot_samples = 10000L,
    parallel = "no", # provide "multicore" (on Windows) or "snow" (on Linux)
    ncpus = 1L) {

# cbind dataset 
data <- as.data.frame(cbind(data_x_axis, data_y_axis))

  # define core function for bootstrapping
  core_calcs <- function(data, indices, ...) {
          stats <- list() # list to store results
          sampled_data <- data[indices, ] # subset data based on sampled_ID
          # correlation
          stats$correlation <- cor(x = sampled_data[,1], y = sampled_data[,2])
          # slope and intercept
          reg <- lm(sampled_data[,2] ~ sampled_data[,1], data = sampled_data)
          stats$slope     <- as.numeric(coef(reg)[2])
          stats$intercept <- as.numeric(coef(reg)[1])
          # mean difference as mean_x - mean_y
          stats$mean_diff <- mean(sampled_data[,1]) - mean(sampled_data[,2])
          # stats$size_samples <- nrow(sampled_data) # if you want to check size of sample used = nrow(data analysed)
          return(unlist(stats))
  }

  # without bootstrapping
  if (bootstrap == F) {
        
        results <- as.data.frame(core_calcs(data = data))
        colnames(results) <- VG_label  # add a VG_label for Validation Group
        
  # with bootstrapping
      } else if (bootstrap == T) {
          # print clustering info
          cat(
            "\nRunning bootsrapping, using:",
            boot_samples, " samples.\n",
            "Other options are:\n parallel = ", parallel, "\n",
            "n. of CPUs:", ncpus, "\n"
          )
          # make a cluster and export in each node needed objects
          cl <- makeCluster(ncpus)
          clusterExport(
            cl = cl,
            varlist = c("data"),
            envir = environment(core_calcs)
          )
          # run bootstrap
          b <- boot(
            data = data,
            # VG_label = VG_label,
            statistic = core_calcs,
            R = boot_samples,
            parallel = parallel,
            ncpus = ncpus,
            cl = cl # use the created cluster
          )
          # close cluster
          stopCluster(cl)
           # convert to data.frame and compute SE
          results <- data.frame(
            value = b$t0,
            SE = apply(b$t, MARGIN = 2, sd)
          )
          # append n. samples to results
          results <- rbind(results, data.frame(value = boot_samples, SE = NA, row.names = "boot_samples"))
          colnames(results) <- paste0(VG_label, "_", colnames(results))
      }
  return(results)
}

# function to make scatterplots predictions
plot_GP <- function(data, xcol, ycol, title, subtitle, slope, intercept ) {
                    # [NOT WORKING ]  SE_slope, SE_intercept) {
  p <- ggplot(data, aes(x = !!sym(xcol), y = !!sym(ycol), color = Bron2)) + 
        geom_point(aes(color = Bron2, shape = Bron2))+
        geom_abline(intercept = intercept, slope = slope, color = "blue")+
        geom_abline(slope = 1, color = "gray") +
        # [NOT WORKING ] # Add shaded regions for SE of slope and intercept
        # geom_ribbon(aes(ymin = intercept - SE_intercept, ymax = intercept + SE_intercept), 
        #             fill = "lightblue", alpha = 0.5) +
        # geom_ribbon(aes(ymin = slope - SE_slope, ymax = slope + SE_slope), 
        #             fill = "lightgreen", alpha = 0.5) +
        # [geom_smooth can be used I think that not the same SE as in boostrapping are visualized here]
        # geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
        theme_bw() +
                labs(
                  x = element_text(paste0("GEBV_", xcol)), y = ycol,
                  title = title,
                  subtitle = subtitle
                ) +
                theme(aspect.ratio = 1)
  return(p)
}

# function to plot convergence of parameters
plot_conv <- function(file, x_axis) {

  data_p <- read.table(paste0("output/", file), header = T) # read file
  data_p <- cbind(x_axis_info, data_p) # attach iterations for x-axis plotting

  p <- data_p %>% 
          ggplot(aes(x= x_axis_info, y=  .data[[colnames(data_p)[2]]]) )+
          geom_line()+
          xlab("iterations")+
          theme_bw()+
          theme(axis.text.x=element_text(angle=+90))+
          labs(
            title = file
          )
  return(p)
}

# 2) read data -------------------------------------
# ___2.1) read predictions
# predictions on training
training_pred <- read.table("training_predictions.csv", sep = ",", header= T)
head(training_pred)
# predictions on validation
validation_pred <- read.table("validation_predictions.csv", sep = ",", header= T)
head(validation_pred)

#  ___2.2) read and add info on DRP
training_set   <- read.table("training_set.txt", header = T)   %>% rename(mean_DRP_orig = mean_DRP)
validation_set <- read.table("validation_set.txt", header = T) %>% rename(mean_DRP_orig = mean_DRP)

# merge DRP info with predictions
validation_pred_ext <- cbind(validation_pred, validation_set)
training_pred_ext   <- cbind(training_pred,     training_set)

# Implement an error check if mean_DRP in input files do not match against the Julia mean_DRP output files
if( !(all(validation_pred_ext$mean_DRP == validation_pred_ext$mean_DRP_orig)==T && 
      all(training_pred_ext$mean_DRP   ==   training_pred_ext$mean_DRP_orig)==T )
) {stop("ERROR: order of mean_DRP in input files and output files do not match: there is likely a mistmatch or re-ordering!")}

# 3) compute validation metrics -----------------
results_list <- list() # list to store results

for(col_to_validate in colnames(validation_pred)) {
  cat("Now doing:", col_to_validate, "\n")

  if(col_to_validate != "mean_DRP") {

    results <- cbind(
      # validation group
      compute_prediction_results(
        data_x_axis = validation_pred_ext[col_to_validate],
        data_y_axis = validation_pred_ext$mean_DRP,
        VG_label = paste0("validation_", col_to_validate),
        bootstrap = T, boot_samples = 10000, parallel = "snow", ncpus = 10),
      
      # training group
      compute_prediction_results(
        data_x_axis = training_pred_ext[col_to_validate],
        data_y_axis = training_pred_ext$mean_DRP,
        VG_label = paste0("training_", col_to_validate),
        bootstrap = T, boot_samples = 10000, parallel = "snow", ncpus = 10)
    )
  results_list[[col_to_validate]] <- results
  rm(results)
  }
}

cat("\n\n Prediction results: \n\n"); names(results_list)
lapply(results_list, head) %>% as.data.frame

# write results to disk
write.csv(x = lapply(results_list, head) %>% as.data.frame,
              file = "prediction_results.txt", quote = F)

# 4) Make plots ---------------------------------
# __4.1) make training and validation plots -----
pdf("prediction_plots.pdf", height = 9, width = 14)
  for(col_to_validate in colnames(validation_pred)) {
  cat("Now plotting:", col_to_validate, "\n")
  if(col_to_validate != "mean_DRP") {
      # validation
      p_val <- plot_GP(data = validation_pred_ext, title="validation set",
              subtitle=col_to_validate,
              xcol = col_to_validate,
              ycol = "mean_DRP",
              slope =     results_list[[col_to_validate]]["slope",     paste0("validation_", col_to_validate, "_value")],
              intercept = results_list[[col_to_validate]]["intercept", paste0("validation_", col_to_validate, "_value")]
              ) 
              # SE_slope = results["slope", "validation_SE"],    SE_intercept = results["intercept", "validation_SE"])
      
      # training
      p_tra <- plot_GP(data = training_pred_ext, title="training set",
              subtitle=col_to_validate,
              xcol = col_to_validate,
              ycol = "mean_DRP",
              slope =     results_list[[col_to_validate]]["slope",     paste0("training_", col_to_validate, "_value")],
              intercept = results_list[[col_to_validate]]["intercept", paste0("training_", col_to_validate, "_value")]
              )
    p_comb <- p_val + p_tra
    print(p_comb)
    }
  }
dev.off()

# __4.2) make convergence plots

# __4.2.1) get info on runLMEM
# Define the grep command to match lines starting with "runLMEM" and not starting with "#"
grep_command <- "grep '^runLMEM' *.jl | grep -v '^#'"
# Use system() to execute the grep command and capture the output
julia_settings <- system(grep_command, intern = TRUE)

julia_numbers <- as.numeric(unlist(
                  strsplit(
                    # Extract runLMEM numbers from the output using regular expressions
                    gsub(".*runLMEM\\(f,pheno_training,([0-9]+),([0-9]+),([0-9]+).*",
                          "\\1 \\2 \\3",
                          julia_settings),
                      " ") # Split the extracted numbers into a numeric vector
                 ))

x_axis_info <- seq(from=julia_numbers[2]+julia_numbers[3], to=julia_numbers[1], by=julia_numbers[3])

# __4.2.2) read all files to plot (all that start with b and var, but exclude betaM - too many effects)
files_to_plot <- list.files("output/") %>%
                  as.data.frame %>%
                  filter(grepl("b|var|pi", .)) %>% 
                  filter(!grepl("^betaM", .)) # exclude betaM (too many, i.e. = # SNPs)

# __4.2.3) write output in a pdf, each page is a parameter
pdf("convergence.pdf")
  for(f in 1:nrow(files_to_plot)) {
    file = files_to_plot[f, ]
    p <- plot_conv(file = file, x_axis = x_axis_info)
    print(p)
  }
dev.off()
