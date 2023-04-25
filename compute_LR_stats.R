compute_LR_stats <-
  function(EBV_partial,
           EBV_whole,
           val_groupIDs = NA,
           val_subgroup = NA,
           VAR_A = NA,
           average_F = NA,
           inbreeding = NULL,
           plot = F,
           plot_verbose = F,
           plot_subgroups = F,
           bootstrap = F,
           boot_samples = NA,
           parallel = "no",
           ncpus = 1L,
           verbose = F) {
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
    #' @param inbreeding data.frame with animal's ID in the first column and
    #' inbreeding coefficients in the second column. The data.frame must include entries for animals in the validation group.
    #' Other animals may be included as well. This argument is required when 'bootstrap = T', and it is used to compute the average
    #' inbreeding value of validation animals sampled within each bootstrap sample.
    #' @param plot logical (default = FALSE). Return a scatter plot of EBV partial and EBV whole.
    #' If 'VAR_A' is provided then EBVs are scaled by the genetic standard deviation.
    #' @param plot_verbose logical (default = FALSE). Add more details to the returned plot:
    #' adds n. animals used in subtitle and computed stats in a caption. If 'VAR_A' is provided,
    #' level bias in caption is in GSD. 
    #' @param plot_subgroups logical (default = FALSE). When TRUE, it reads a second column in val_groupIDs
    #' that is considered a "group column". This column has characters that label IDs in val_groupIDs in 'n-groups'.
    #' IDs will be coloured by 'group' in the returned plot. Color sub-groups within the val_groupIDs
    #' can be useful for detecting heterogeneous groups.
    #' @param bootstrap logical (default = FALSE). When '= T', bootstrapping will be used to
    #' estimate standard errors of each LR statistics (see @details for more information).
    #' @param boot_samples integer (default = NA). Number of samples to use for bootsrapping (see @details for more information).
    #' @param parallel character (default = "no"). Defines the parallel option for {boot} package.
    #' Usually, on Windows use 'parallel = "snow"', on Linux use 'parallel = "multicore"'.
    #' (see @details for more information).
    #' @param ncpus integer (default = 1). Number of CPUs to use for parallelization when bootstrapping.
    #' (see @details for more information).
    #' @param verbose logical (default = F). When set to TRUE, prints more details of internal input/output.
    #'
    #'
    #' @return the function returns a list with these elements:
    #' 1) 'stats', a data.frame with two columns named "value" and "SE". SE provided only if 'bootstrap = TRUE'.
    #' 'value' are the estimates of the LR statistics for: level bias, level bias in genetic standard deviations (level_bias_GSD)
    #' (only if VAR_A is provided), dispersion bias, ratio of accuracies (rho), increase in accuracies (inc_acc),
    #' accuracy of partial EBV (accuracy_partial) (only if 'VAR_A' and 'average_F' OR 'inbreeding' are provided)
    #' 2) 'plot' is a scatter plot of whole EBV on partial EBV. Plot is returned only when 'plot = TRUE'.
    #'
    #' @usage compute_LR_stats(EBV_whole, EBV_partial, val_groupIDs, ...)
    #'
    #' @importFrom assertthat assert_that is.scalar
    #' @importFrom ggplot2 ggplot
    #' @importFrom parallel makeCluster clusterExport
    #' @importFrom boot boot
    #' @importFrom methods hasArg
    #'
    #' @author Renzo Bonifazi
    #'
    #' @references
    #' Github repo: https://github.com/bonifazi/R_utils
    #'
    #' Referenced papers:
    #' Legarra and Reverter, Genet Sel Evol, 2018. Semi-parametric estimates of population accuracy and bias of
    #' predictions of breeding values and future phenotypes using the LR method.
    #' https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0426-6
    #'
    #'
    #' @details
    #' \subsection{Validation metrics estimated}{
    #'   Validation metrics are implemented following Legarra and Reverter (GSE, 2018) as:
    #'   level_bias = mean EBV partial - mean EBV full;
    #'   level bias in genetic standard deviations (level_bias_GSD) = level_bias / sqrt(VAR_A);
    #'   level bias in GSD (if VAR_A is provided);
    #'   dispersion = cov(partial EBV, whole EBV)/var(partial EBV);
    #'   accuracy of partial EBV = sqrt(cov(partial EBV, whole EBV) / ((1 - average_F) * VAR_A));
    #'   ratio of accuracies (rho) = cor(partial EBV, whole EBV);
    #'   increases in population accuracies (inc_acc) = 1/rho.  
    #'   Validation metrics are rounded to 6 decimals.
    #' }
    #' 
    #' \subsection{Obtaining standard errors of estimated metrics}{
    #'   Standard Errors (SE) are obtained using bootstrapping with re-sampling via the {boot} package.
    #'   For more details on bootrapping see the {boot} package manual (https://stat.ethz.ch/R-manual/R-devel/library/boot/html/boot.html)
    #' }
    #' \subsection{Plot}{
    #'   Plot shows 2 regression lines: a grey one with slope = 1 used as an indicator,
    #'   and a blue one which is made by regressing whole EBV on partial EBV for 'val_groupIDs'
    #' }
    #' 
    #'
    #' @examples
    #' # generate some data
    #' data <- as.data.frame(matrix(nrow = 10000, ncol = 3))
    #' data[, 1] <- 1:10000 # animals
    #' data[, 3] <- rnorm(10000, mean =   0, 12) # generate whole EBV
    #' data[, 2] <- 5+(0.80*data[,3])+rnorm(n = nrow(data), mean = 0, sd = 2) # generate partial EBV
    #' val_groupIDs <- data.frame(ID = c(1:100), group = c(rep("A", 50), rep("B", 50)))
    #' colnames(data) <- c("AID", "partial", "whole")
    #' inbreeding <- data.frame(ID = data$AID, F = rnorm(n = nrow(data), mean = 0.05, sd = 0.01))
    #'
    #' # not providing ID to subset will use all ID in 'data'
    #' results <- compute_LR_stats(EBV_partial = data[,c(1,2)], EBV_whole = data[,c(1,3)])
    #' # providing IDs to subset will use only these ID from 'data' in val_groupIDs
    #' results <- compute_LR_stats(EBV_partial = data[,c(1,2)], EBV_whole = data[,c(1,3)], val_groupIDs = data.frame(c(1:10)))
    #' # providing 1 ID that is not present in 'data' will return an error
    #' results <- compute_LR_stats(EBV_partial = data[,c(1,2)], EBV_whole = data[,c(1,3)], val_groupIDs = data.frame(c("-100")))
    #' # providing a mix of IDs present in 'data' and one ID that is not present in 'data' will return an error
    #' results <- compute_LR_stats(EBV_partial = data[,c(1,2)], EBV_whole = data[,c(1,3)], val_groupIDs = data.frame(c("-100", 1:10)))
    #' # using plot
    #' results <- compute_LR_stats(EBV_partial = data[,c(1,2)], EBV_whole = data[,c(1,3)], val_groupIDs = val_groupIDs, average_F = 0.01, VAR_A = 20.5, plot = T, plot_verbose = T, plot_subgroups = T)
    #' # using plot and bootstrapping
    #' results <- compute_LR_stats(EBV_partial = data[,c(1,2)], EBV_whole = data[,c(1,3)], val_groupIDs = val_groupIDs, VAR_A = 20.5, inbreeding = inbreeding, plot = T, plot_verbose = T, plot_subgroups = T, bootstrap = T, boot_samples = 100, parallel = "snow", ncpus = 4)
    #' # view the LR stats results
    #' results$stats
    #' # view the plot
    #' results$plot
    #'
    #' @section
    #' TODO:
    #' [done]
    #' - fixed error with parallel on windows by definin the env of the cluster varlist,
    #' following https://stackoverflow.com/questions/25728178/r-cluster-export-error-object-not-found
    #'
    #' [Higher priority]
    #' - allow for providing inbreeding data.frame instead of average_F
    #'   use something like this in the core_calc function or outside.
    #'   implement a check as well upon input?
    #'   acc_p_no_boot <- !is.na(average_F) & !is.na(VAR_A) | hasArg(inbreeding) & !is.na(VAR_A) # without bootstrapping
    #' - change = NA to  = NULL or just remove the NA.
    #'    > https://stackoverflow.com/questions/7964830/test-if-an-argument-of-a-function-is-set-or-not-in-r
    #'    > https://stackoverflow.com/questions/9877271/how-to-check-existence-of-an-input-argument-for-r-functions
    #' - expand the verbose args to print out some useful info (# animals read or used and samples used in boot)
    #' e.g. "no IDS provided, using all animals in EBV_whole' and 'EBV_partial'.
    #' - acc_p = check on warnings on NaN it may be the numerator as well. Implement a check?
    #'
    #' [Lower priority]
    #' - get rid of this ggplot warning:
    #'   5: `aes_string()` was deprecated in ggplot2 3.0.0.
    #'   ℹ Please use tidy evaluation idioms with `aes()`.
    #'   ℹ See also `vignette("ggplot2-in-packages")` for more information.
    #'   This warning is displayed once every 8 hours.
    #'   Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    #'   generated.
    #' - clarify that the acc. partial EBV is different from that in the original paper and more like Macedo paper / or my paper.
    #' - implement val_subgroups
    #' - add progress bar when bootsrapping
    #' - look into boot.ci: https://www.youtube.com/watch?v=ZHspLhlw97k
    #' - Chol can be used for better simulated data and show that stats returna values as simulated, but maybe there is a base R func): would be usefel to simulate real level bias, dispersion, and at least correlation
    #'
    # --------------------- END DOCUMENTATION ----------------------------------
    # 0) Check on pkgs  --------------------------------------------------------
    check_pkgs <- function(...) {
      for (pkg in c(...)) {
        if (!require(pkg, character.only = T)) {
          stop("Please install the ", pkg, ' package, using "install.packages(', pkg, ')"', call. = F)
        }
      }
    }
    check_pkgs("assertthat", "ggplot2")
    if (bootstrap == T) {
      check_pkgs("parallel", "boot")
    }
    
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
    # without bootstrapping
    if (bootstrap == F) {
      assert_that(is.scalar(VAR_A) || is.na(VAR_A),
                  is.scalar(average_F) || is.na(average_F),
                  msg = "Given args 'VAR_A' and/or 'average_F' are not scalars. Check your input args structure."
      )
    } else if (
      # with bootstrapping
      bootstrap == T) {
      assert_that(is.scalar(VAR_A) || is.na(VAR_A),
                  msg = "Given args 'VAR_A' and/or 'average_F' are not scalars. Check your input args structure."
      )
      assert_that(is.data.frame(inbreeding) || is.null(inbreeding),
                  msg = "Given args 'inbreeding' is not a data.frame. Check your input args structure."
      )
      # check if average_F was mistakenly put instead of inbreeding when requiring bootstrap
      if (bootstrap == T & !is.na(average_F)) {
        stop("When 'bootstrap = T'  provide the 'inbreeding' arg instead of 'average_F'. Check also that you are not providing both 'inbreeding' and 'average_F'.")
      }
    }
    # checks on required args when bootsrapping = T
    if (bootstrap == T & !hasArg(boot_samples)) {
      stop("Argument 'boot_samples' is missing, while requesting bootstrap = T. Please provide the number of samples to use.")
    }
    
    # checks on validation IDs (if user provided them)
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
        stop("Not all validation group IDs were found in both 'EBV_partial' and 'EBV_whole' args. Check that both provided data.frame contains the IDs provided in 'val_groupIDs'")
      }
      # if bootstrap required & inbreeding file provided, check that:
      #  validation IDs can be found and that there are no NAs
      if (bootstrap == T & hasArg(inbreeding)) {
        if (!all(val_groupIDs[, 1] %in% inbreeding[, 1])) {
          stop("Not all validation group IDs were found in the 'inbreeding' data.frame. Check that col. 2 of 'inbreeding' data.frame has the IDs provided in 'val_groupIDs'")
        }
        if (any(
          is.na(inbreeding[inbreeding[, 1] %in% val_groupIDs[, 1], 2])
        )
        ) {
          stop("NAs were found in 'inbreeding' data.frame. Check it.")
        }
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
    # The first col of EBV_p and EBV_w are used as index for merging
    merged <- merge(EBV_p[, c(1, 2)], EBV_w[, c(1, 2)], by.x = colnames(EBV_p)[1], by.y = colnames(EBV_w)[1])
    
    # check that all IDs were found in either EBV_partial and EBV_whole, otherwise stop (I have this as default behavior for now)
    if (hasArg(val_groupIDs)) {
      if (nrow(merged) < nrow(val_groupIDs)) {
        stop("Not all IDs provided in 'val_groupIDs' were sub-set from the 'EBV_partial' and 'EBV_whole' data.frames. Make sure to provide in 'val_groupIDs' args only IDs that exist in data.frames 'EBV_partial' and 'EBV_whole'")
      }
    }
    
    # 3) Compute LR statistics -------------------------------------------------
    # 3.1) Function for computing statistics -----------------------------------
    # core_calcs <- function(data, sampled_IDs, ..., VAR_A, average_F = average_F, inbreeding = inbreeding, bootstrap = bootstrap) {
    core_calcs <- function(data, sampled_IDs, ...) {
      stats <- list() # list to store sample results
      sampled_data <- data[sampled_IDs, ] # subset data based on sampled_ID
      # 1. level bias (as mean EBV partial - mean EBV full)
      stats$level_bias <- mean(sampled_data[, 2]) - mean(sampled_data[, 3])
      # 2. level bias expressed in genetic standard deviation
      if (!is.na(VAR_A)) {
        stats$level_bias_in_GSD <- stats$level_bias / sqrt(VAR_A)
      } else {
        stats$level_bias_in_GSD <- NULL # to keep returned list dimension consistent
      }
      # 3. dispersion bias # cov(partial EBV, whole EBV)/var(partial EBV)
      stats$dispersion_bias <- cov(sampled_data[, 2], sampled_data[, 3]) / var(sampled_data[, 2])
      # 4. accuracy of partial evaluation
      # two possible conditions to go ahead with acc_p calculation
      acc_p_no_boot <- !is.na(average_F) & !is.na(VAR_A) # without bootstrapping
      acc_p_yes_boot <- bootstrap == T & hasArg(inbreeding) & !is.na(VAR_A) # with bootstrapping
      acc_p_yes_boot <- bootstrap == T & !is.na(VAR_A) # with bootstrapping
      # if any of the two = T, go ahead
      if (acc_p_no_boot == T | acc_p_yes_boot == T) {
        # if bootstrapping active, compute average_F of the bootstrap sampled IDs
        if (acc_p_yes_boot == T) {
          average_F <- mean(inbreeding[sampled_IDs, 2])
        }
        # acc_p = sqrt(cov(partial EBV, whole EBV) / ((1 - average_F) * VAR_A)
        stats$accuracy_partial <- sqrt(
          cov(sampled_data[, 2], sampled_data[, 3]) / ((1 - average_F) * VAR_A)
        )
        stats$average_F <- average_F
        stats$VAR_A <- VAR_A
        if (is.na(stats$accuracy_partial)) {
          warning("NaN produced for accuracy partial, likely due to computing the sqrt of a negative value in '(1 - average_F * VAR_A)'.")
        }
      }
      # 5. ratio of accuracies (rho)
      stats$rho <- cor(sampled_data[, 2], sampled_data[, 3], method = "pearson")
      # 6. increase in accuracies (1/rho)
      stats$inc_acc <- 1 / stats$rho
      return(unlist(stats))
    }
    
    # 3.2) Run core function ---------------------------------------------------
    results <- list() # initialize list to store results
    
    # ___ 3.2.A) No bootstrapping ----------------------------------------------
    if (bootstrap == F) {
      stats <- core_calcs(data = merged, sampled_IDs = seq_along(merged[, 1]))
      stats <- data.frame(value = stats)
    } else if (bootstrap == T) {
      # ___ 3.2.B) With bootstrapping ------------------------------------------
      if (verbose == T) {
        cat(
          "\nRunning 'compute_LR_stats' with bootsrapping, using:",
          boot_samples, " samples.\n",
          "Other options are:\nparallel = ", parallel, "\n",
          "n. of CPUs:", ncpus, "\n"
        )
      }
      # make a cluster and export in each node needed objects
      cl <- makeCluster(ncpus)
      clusterExport(
        cl = cl,
        varlist = c(
          "inbreeding", "val_groupIDs",
          "VAR_A", "average_F", "bootstrap"
        ),
        envir = environment(core_calcs)
      )
      if (verbose == T) {
        print(cl)
      }
      # run bootstrap
      b <- boot(
        data = merged,
        statistic = core_calcs,
        R = boot_samples,
        parallel = parallel,
        ncpus = ncpus,
        cl = cl # use the created cluster
      )
      stopCluster(cl) # close cluster
      # convert to data.frame and compute SE
      stats <- data.frame(
        value = b$t0,
        SE = apply(b$t, MARGIN = 2, sd)
      )
      # append n. samples to results
      stats <- rbind(stats, data.frame(value = boot_samples, SE = NA, row.names = "boot_samples"))
    }
    results$stats <- round(stats, 6)
    
    # 4) Plot ------------------------------------------------------------------
    # plot whole on partial
    if (plot == T) {
      # geoms layers that apply on all conditions
      geoms_layers <- list(
        geom_abline(slope = 1, color = "gray"),
        geom_point()
      )
      # make a reg slope that fits plotted data (with plotted data as defined below)
      # note that lm is fitted because level bias is not the same as the intercept of a regression slope (see Legarra and Reverter, 2018)
      make_reg_line <- function(plot_data) {
        coef_plot <- coef(lm(plot_data[, 3] ~ plot_data[, 2]))
        reg_line <- geom_abline(intercept = coef_plot[1], slope = coef_plot[2], color = "blue")
        return(reg_line)
      }
      # scale x- and y- axis by trait GSD if VAR_A is provided
      if (!is.na(VAR_A)) {
        plot_data <- merged
        plot_data[, 2] <- plot_data[, 2] / sqrt(VAR_A) # scale partial
        plot_data[, 3] <- plot_data[, 3] / sqrt(VAR_A) # scale whole
        reg_line <- make_reg_line(plot_data) # reg line will use level bias in GSD
      } else {
        plot_data <- merged # if no VAR_A then no scaling applied
        reg_line <- make_reg_line(plot_data) # reg line will use level bias
      }
      # plot without subgroups
      if (plot_subgroups == F) {
        p <- ggplot(data = plot_data, aes_string(
          x = colnames(plot_data)[2], # partial
          y = colnames(plot_data)[3] # whole
        )) +
          geoms_layers + reg_line
      } else if (plot == T & plot_subgroups == T) { # add subgroups to plot
        # add group IDs to merged object
        merged_subgrp <- merge(plot_data, val_groupIDs[, c(1, 2)], by.x = colnames(plot_data)[1], by.y = colnames(val_groupIDs)[1])
        
        p <- ggplot(data = merged_subgrp, aes_string(
          x = colnames(merged_subgrp)[2], # partial
          y = colnames(merged_subgrp)[3], # whole
          color = colnames(merged_subgrp)[ncol(merged_subgrp)] # group
        )) +
          geoms_layers + reg_line
      }
      # add common elements to plot
      p <- p +
        coord_fixed(ratio = 1) +
        theme_bw() +
        labs(
          x = element_text("EBV Partial"), y = element_text("EBV Whole")
        ) +
        theme(aspect.ratio = 1)
      # if verbose -> add more details into plot
      if (plot_verbose == T) {
        if (!is.na(VAR_A)) {
          details <- paste0("LB (GSD): ", round(results$stats["level_bias_in_GSD", 1], 3))
        } else {
          details <- paste0("LB: ", round(results$stats["level_bias", 1], 3))
        }
        details <- paste0(
          details,
          "; disp.: ", round(results$stats["dispersion_bias", 1], 3),
          "; rho: ", round(results$stats["rho", 1], 3)
        )
        # add subtitle
        p <- p + labs(
          subtitle = element_text(paste0("N. animals in plot = ", nrow(plot_data))),
          caption = details
        )
      }
      results$plot <- p
    } else {
      results$plot <- NULL # to keep returned list dimension consistent
    }
    # return results
    return(results)
  }
