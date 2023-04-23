compute_LR_stats <-
  function(EBV_whole,
           EBV_partial,
           val_groupIDs = NA,
           val_subgroup = NA,
           VAR_A = NA,
           average_F = NA,
           bootstrap = F,
           plot = F,
           plot_verbose = F,
           plot_subgroups = F) {
    # ------------------------- START DOCUMENTATION ----------------------------
    #' Compute LR method validation statistics
    #'
    #' @description
    #' Compute LR method (Legarra and Reverter, GSE, 2018) statistics and plot
    #' the comparison. LR statistics are computed by comparing EBV from a 'whole' evaluation with EBV from a 'partial'
    #' evaluation. Whole and partial EBV are provided by the user (see below for more information).
    #' LR statistics are computed for a validation group of individuals, called 'focal group'.
    #' The IDs of the focal group are provided in 'val_groupIDs'.
    #' If no ID are provided, the function computes the statistics using all individuals.
    #' The function can return a scatter plot for visual inspection.
    #'
    #' @param EBV_whole a data.frame with animals' IDs in first column and EBV from
    #' whole evaluations in the second column.
    #' @param EBV_partial a data.frame with animals' IDs in first column and EBV from
    #' whole evaluations in the second column.
    #' @param val_groupIDs a data.frame with animals' IDs in first column. These IDs
    #' are used to compute LR method statistics.
    #' @param val_subgroup #to be implemented. split statistics per each sub-group provided
    #' @param VAR_A additive genetic variance (a scalar). When provided, level bias is also returned
    #' scaled in genetic standard deviations (see @details for formula).
    #' @param average_F average inbreeding level (a scalar) of the validation group. When provided,
    #' it is used to return the accuracy of partial evaluations (see @details for formula).
    #' @param bootstrap [in development] logical (default = FALSE). When '= T', bootstrapping will be used to
    #' estimate Standard Errors of LR statistics (see @details for more information).
    #' @param plot logical (default = FALSE). Return a scatter plot of EBV partial and EBV whole.
    #' @param plot_verbose logical (default = FALSE). Add more details to the returned plot
    #' (for now, adds the number of animals used in a subtitle).
    #' @param plot_subgroups logical (default = FALSE). When TRUE, it reads a second column in val_groupIDs
    #' that is considered a "group column". This column has characters that label IDs in val_groupIDs in 'n-groups'.
    #' IDs will be coloured by 'group' in the returned plot. Color sub-groups within the val_groupIDs
    #' can be useful for detecting heterogeneous groups.
    #'
    #' @return a list with: level bias, level bias in genetic standard deviations (level_bias_GSD)
    #' (only if VAR_A is provided), dispersion bias, ratio of accuracies (rho), increase in accuracies (inc_acc),
    #' accuracy of partial EBV (accuracy_partial) (only if VAR_A and average_F are provided), and
    #' a scatter plot (if plot = T).
    #'
    #' @usage compute_LR_stats(EBV_whole, EBV_partial, val_groupIDs, ...)
    #'
    #' @importFrom assertthat assert_that is.scalar
    #' @importFrom ggplot2 ggplot
    #' @importFrom methods hasArg
    #'
    #' @author Renzo Bonifazi
    #' 
    #' @references Legarra and Reverter, Genet Sel Evol, 2018. Semi-parametric estimates of population accuracy and bias of
    #' predictions of breeding values and future phenotypes using the LR method.
    #' https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0426-6
    #'
    #' @details
    #' Validation metrics are implemented following Legarra and Reverter (GSE, 2018) as:
    #' level_bias = mean EBV partial - mean EBV full;
    #' level bias in genetic standard deviations (level_bias_GSD) = level_bias / sqrt(VAR_A);
    #' level bias in GSD (if VAR_A is provided);
    #' dispersion = cov(partial EBV, whole EBV)/var(partial EBV);
    #' accuracy of partial EBV = sqrt(cov(partial EBV, whole EBV) / ((1 - average_F) * VAR_A));
    #' ratio of accuracies (rho) = cor(partial EBV, whole EBV);
    #' increases in population accuracies (inc_acc) = 1/rho
    #'
    #' @examples
    #' # generate some data
    #' data <- as.data.frame(matrix(nrow = 10000, ncol = 3))
    #' data[, 1] <- 1:10000
    #' data[, 2] <- rnorm(10000, mean =   1, 12)
    #' data[, 3] <- rnorm(10000, mean = 0.8, 12)
    #' val_groupIDs <- data.frame(ID = c(1:10), group = c(rep("A", 5), rep("B", 5)))
    #' colnames(data) <- c("AID", "whole", "partial")
    #'
    #' # not providing ID to subset will use all ID in 'data'
    #' results <- compute_LR_stats(EBV_whole = data[,c(1,2)], EBV_partial = data[,c(1,3)])
    #' # providing IDs to subset will use only these ID from 'data' in val_groupIDs
    #' results <- compute_LR_stats(EBV_whole = data[,c(1,2)], EBV_partial = data[,c(1,3)], val_groupIDs = data.frame(c(1:10)))
    #' # providing 1 ID that is not present in 'data' will return an error
    #' results <- compute_LR_stats(EBV_whole = data[,c(1,2)], EBV_partial = data[,c(1,3)], val_groupIDs = data.frame(c("-100")))
    #' # providing a mix of IDs present in 'data' and one ID that is not present in 'data' will return an error
    #' results <- compute_LR_stats(EBV_whole = data[,c(1,2)], EBV_partial = data[,c(1,3)], val_groupIDs = data.frame(c("-100", 1:10)))
    #' # using all features (plot etc.)
    #' compute_LR_stats(EBV_whole = data[,c(1,2)], EBV_partial = data[,c(1,3)], val_groupIDs = val_groupIDs, average_F = 0.01, VAR_A = 20.5, plot = T, plot_verbose = T, plot_subgroups = T)
    #' 
    #' @section
    #' TODO:
    #' [done]
    #' - return average_F and var_A (can be useful for checks)
    #' - check on having loaded pkgs loaded (implement warning on loaded ggplot2 and assertthat pkgs)
    #' 
    #' [Higher priority]
    #' - make a verbose args with T as default to print out some useful info (# animals read or used and samples used in boot)
    #' e.g. "no IDS provided, using all animals in EBV_whole' and 'EBV_partial'.
    #' - bootstrap for SE (R boot pkg with 10K default samples, using drawing with replacment)
    #' - average_F take average from df with [ID, F] cols (with checks implemented as well)
    #' - acc_p = check on warnings on NaN it may be the numerator as well. Implmenet a check?
    #' 
    #' [Lower priority]
    #' - better example dataset: make a multi-variate sample with correlation (Chol can be used, but maybe there is a base R func): would be usefel to simulate real level bias, dispersion, and at least correlation
    #' - clarify that the acc. partial EBV is different from that in the original paper and more like Macedo paper / or my paper.
    #' - implement val_subgroups
    #' 
    # --------------------- END DOCUMENTATION ----------------------------------
    # 0) Check on pkgs  --------------------------------------------------------
    check_pkgs <- function(...){
      for (pkg in c(...)) {
        if(!require(pkg, character.only = T)){
          stop('Please install the ', pkg, ' package, using "install.packages(', pkg, ')"', call. = F)
        }
      }
    }
    check_pkgs("assertthat", "ggplot2")
    
    # 1) Checks on given input args --------------------------------------------
    # checks on EBV_partial and EBV_whole
    assert_that(is.data.frame(EBV_partial),
                is.data.frame(EBV_whole),
                msg = "Given args EBV_whole and/or EBV_partial are not data.frame. Check your input args structure."
    )
    # checks on val_groupIDs
    assert_that(is.data.frame(val_groupIDs) || is.na(val_groupIDs),
                msg = "Given args 'val_groups' is not data.frame. Check your input arg structure."
    )
    # checks on VAR_A, and average_F
    assert_that(is.scalar(VAR_A) || is.na(VAR_A),
                is.scalar(average_F) || is.na(average_F),
                msg = "Given args 'VAR_A' and/or 'average_F' are not scalars. Check your input args structure."
    )
    
    # check on validation IDs (if user provided them)
    if (hasArg(val_groupIDs)) {
      # check that >0 IDs can be found in EBV_partial and EBV_whole
      if (length(EBV_partial[EBV_partial[, 1] %in% val_groupIDs[, 1], ][, 1]) == 0 ||
          length(EBV_whole[EBV_whole[, 1] %in% val_groupIDs[, 1], ][, 1]) == 0) {
        stop("No IDs provided in 'val_groupIDs' was found in 'EBV_partial' or 'EBV_whole'. Check that the provided IDs are into these data.frames.")
      }
      # check that the same IDs can be found in both EBV_partial and EBV_whole
      if (!all(
        (EBV_partial[EBV_partial[, 1] %in% val_groupIDs[, 1], ][, 1]) %in% (EBV_whole[EBV_whole[, 1] %in% val_groupIDs[, 1], ][, 1])
      )) {
        stop("'Not all validation group IDs were found in both 'EBV_partial' and 'EBV_whole' args. Check that both provided data.frame contains the IDs provided in 'val_groupIDs'")
      }
    }
    
    # 2) Data preparation ------------------------------------------------------
    
    # subset EBV_partial and EBV_whole for the given IDs in 'val_groupIDs'
    if (hasArg(val_groupIDs)) {
      EBV_p <- EBV_partial[EBV_partial[, 1] %in% val_groupIDs[, 1], ][, c(1, 2)] # carry-on only ID and EBV cols
      EBV_w <- EBV_whole[EBV_whole[, 1] %in% val_groupIDs[, 1], ][, c(1, 2)] # carry-on only ID and EBV cols
    } else { # otherwise keep all
      EBV_p <- EBV_partial
      EBV_w <- EBV_whole
    }
    
    # merge EBV_p and EBV_w to make sure that EBV_p and EBV_w are assigned to same ID
    # it may be that IDs in EBV_p and EBV_w are shuffled
    # The first col of EBV_p and EBV_w is used as index for merging
    merged <- merge(EBV_p[, c(1, 2)], EBV_w[, c(1, 2)], by.x = colnames(EBV_p)[1], by.y = colnames(EBV_w)[1])
    
    # check that all IDs were found in either EBV_partial and EBV_whole, otherwise stop (I have this as default behavior for now)
    if (hasArg(val_groupIDs)) {
      if (nrow(merged) < nrow(val_groupIDs)) {
        stop("Not all IDs provided in 'val_groupIDs' were sub-set from the 'EBV_partial' and 'EBV_whole' data.frames. Make sure to provide in 'val_groupIDs' args only IDs that exist in data.frames 'EBV_partial' and 'EBV_whole'")
      }
    }
    
    # 3) Compute LR statistics -------------------------------------------------
    
    res <- list() # initialize a list to store results
    
    # level bias (as mean EBV partial - mean EBV full)
    res$level_bias <- mean(merged[, 2]) - mean(merged[, 3])
    if (!is.na(VAR_A)) {
      # level bias expressed in genetic standard deviation
      res$level_bias_in_GSD <- res$level_bias / sqrt(VAR_A)
    } else {
      res$level_bias_in_GSD <- NULL # to keep returned list dimension consistent 
    }
    # dispersion bias
    res$dispersion_bias <- cov(merged[, 2], merged[, 3]) / var(merged[, 2])
    if (!is.na(average_F) & !is.na(VAR_A)) {
      # accuracy of partial evaluation
      res$accuracy_partial <- sqrt(
        cov(merged[, 2], merged[, 3]) / ((1 - average_F) * VAR_A)
      )
      res$average_F <- average_F
      res$VAR_A <- VAR_A
      if (is.na(res$accuracy_partial)) {
        warning("NaN produced for accuracy partial, likely due to computing the sqrt of a negative value in '(1 - average_F * VAR_A)'.")
      }
    } else {
      res$accuracy_partial <- NULL # to keep returned list dimension consistent 
      res$average_F <- NULL
      res$VAR_A <- NULL
    }
    # ratio of accuracies (rho)
    res$rho <- cor(merged[, 2], merged[, 3], method = "pearson")
    # increase in accuracies (1/rho)
    res$inc_acc <- 1 / res$rho
    
    # 4) Plot ------------------------------------------------------------------
    # plot whole on partial
    if (plot == T) {
      if (plot_subgroups == F) { # plot without subgroups
        p <- ggplot(data = merged, aes_string(
          x = colnames(merged)[2], # partial
          y = colnames(merged)[3] )) + # whole
          geom_abline(slope = 1, color = "gray") +
          geom_point()
      } else if (plot == T & plot_subgroups == T) { # add subgroups to plot
        # add group IDs to merged object
        merged_subgrp <-  merge(merged, val_groupIDs[,c(1,2)], by.x = colnames(merged)[1], by.y = colnames(val_groupIDs)[1])
        
        p <- ggplot(data = merged_subgrp, aes_string(
          x = colnames(merged_subgrp)[2],     # partial
          y = colnames(merged_subgrp)[3],     # whole
          color = colnames(merged_subgrp)[ncol(merged_subgrp)] )) + # group
          geom_abline(slope = 1, color = "gray") +
          geom_point()
      }
      # add common elements to plot
      p <- p +
        coord_fixed(ratio = 1) +
        theme_bw() +
        labs(
          x = element_text("EBV Partial"), y = element_text("EBV Whole")
        )
      # if verbose -> add more details into plot
      if (plot_verbose == T) { 
        # add subtitle
        p <- p + labs(subtitle = element_text(paste0("N. animals in plot = ", nrow(merged))))
        # will add here new details to plot as function develops
      }
      res$plot <- p
    } else {
      res$plot <- NULL # to keep returned list dimension consistent 
    }
    # return results
    return(res)
  }

# function(data, sample_ID) {subdata <- data[sample_ID,]; res_1 <- cor(subdata[,2], subdata[,3]); res_2 <- mean(subdata[,2]) - mean(subdata[,3]); unlist(list(cor = res_1, lb = res_2))}
# b <- boot(data = data, statistic = f_cor, 1000)
# str(b)
