rebase_ebv <- function(data,
                       base_pop = NULL,
                       constant_value = NULL,
                       decimals = 12,
                       plot = F,
                       verbose = T) {
  # ------------------------- START DOCUMENTATION ------------------------------
  #' Rebase estimated breeding values
  #'
  #' @description
  #' This function re-bases estimated breeding values (EBV) in 'data' for a set of IDs provided in 'base_pop'.
  #' The base_pop IDs are considered the base population and will have a mean EBV of 0.
  #' All other IDs in 'data' will have their EBV re-based by removing the mean EBV of 'base_pop'.
  #'
  #' This re-basing operation is done for each column from column 2 onward.
  #' For instance, you can provide 'data' with structure "ID, EBV1, EBV2, EBV3, ..., EBVn",
  #' where each EBV is for a specific trait.
  #'
  #' Another alternative usage is to provide a fixed EBV value against which all EBVs are going to be expressed (arg `constant_value`).
  #' Be careful with this usage and the interpretation of the results. All EBVs across all columns are going to be
  #' "shifted" and expressed against this constant value.
  #'
  #' @param data a data.frame with individuals' ID in first column, and EBV in col2, col3 .. coln.
  #' @param base_pop a data.frame with individuals' ID in first column to be considered as base population.
  #' @param constant_value an integer (default = NULL). EBVs will be expressed as deviations from this value. It cannot be used if base_pop is provided.
  #' @param decimals numeric. Number of decimals used internally for checking that new
  #' mean EBV for base population IDs is 0 (default = 12).
  #' @param verbose logic (default = T). Print log details to user.
  #' @param plot logic (default = F). Scatterplot of EBV before vs after rebasing. When 'plot = T' a list is
  #' returned (see @retrun for more information). Black line is a regression line with intercept = 0 and slope = 1,
  #' Blue line is regression line of EBV after rebasing on EBV before rebasing.
  #' Intercept of the blue line should be the negation of the mean EBV of the base pop and the intercept should be 1.
  #'
  #' @return
  #' if 'plot = F' (default), returns a data.base as provided in 'data' with IDs in col1, re-based EBV in col2, col3, .. coln.
  #' if 'plot = T', returns a list with:
  #' First element named 'rebased_ebv' and having the results and having same structure as when 'plot = F'),
  #' Second element named 'plot' with as many sub-elements named after the column names of 'data' (excluding col1);
  #' and with each sub-element being a ggplot2 object. Thus, to access the plots use e.g. 'output$plots$col1_name'.
  #'
  #' @usage rebase_ebv(data = data, base_pop = base_pop)
  #'
  #' @importFrom assertthat assert_that
  #' @importFrom ggplot2 ggplot (only when 'plot = T')
  #'
  #' @author Renzo Bonifazi
  #'
  #' @examples
  #' # generate some data
  #' data <- as.data.frame(matrix(nrow = 10000, ncol = 3))
  #' data[, 1] <- 1:10000
  #' data[, 2] <- rnorm(10000, mean =   1, 15)
  #' data[, 3] <- rnorm(10000, mean =  17, 13)
  #' base_pop <- data.frame(ID = c(1:10))
  #' colnames(data) <- c("ID", "ebv1", "ebv2")
  #'
  #' # run the function without plotting
  #' output <- rebase_ebv(data = data, base_pop = base_pop)
  #' head(output) # check results
  #'
  #' # run the function with plotting
  #' output <- rebase_ebv(data = data, base_pop = base_pop, plot = T)
  #' # check output 'str' (now has 2 levels)
  #' str(output, max.level =  2, give.attr = F)
  #'
  #' head(output$rebased_ebv) # check results for
  #' output$plot$ebv1 # view scatterplot for column 1 named 'ebv1'
  #' output$plot$ebv2 # view scatterplot for column 1 named 'ebv2'
  #'
  #' # re-base EBV against a constant value
  #' output <- rebase_ebv(data = data, constant_value = 400, plot = T)
  #' output$plot$ebv1 # view scatterplot for column 1 named 'ebv1'
  #' output$plot$ebv2 # view scatterplot for column 1 named 'ebv2'
  #'
  #' @section
  #' TODO:
  #' [Mid priority]
  #' - (possible): add option to provide as many constant_values as number of cols in data,
  #' e.g. constant_value = c(val_1, val_2, val_3)
  #' [Lower priority]
  #' - (possible): remove dependency on assertthat
  #'
  # ------------------------- END DOCUMENTATION --------------------------------
  # 0) Check on pkgs
  if (!require("assertthat")) {
    stop('Please install the "assertthat" package, using \n:
                                    "install.packages("assertthat")"')
  }

  # 1) Checks on given input args ----------------------------------------------
  # check that at least one between constant_value and base_pop are provided
  if(is.null(constant_value) & is.null(base_pop)) {
    stop('Please provide `base_pop` or `constant_value` as arguments (only one of the two)')
  }

  # check that if constant_value is provided base_pop is not provided
  if(!is.null(constant_value) & !is.null(base_pop)) {
    stop('Please provide only one between `base_pop` and `constant_value`')
  }

  # Check the str of input args when using base_pop
  if(!is.null(base_pop)){
    assert_that(is.data.frame(data), is.data.frame(base_pop), is.numeric(decimals))
    # Check all IDs in 'base_pop' must be found in col1 of 'data' data.frame
    if (!all(base_pop[, 1] %in% data[, 1])) {
      stop("Not all IDs in 'base_pop' were found in 'data'. Check you input args.")
    }
  }

  # Check the str of input args when using constant_value
  if(!is.null(constant_value)) {
    assert_that(is.data.frame(data), is.numeric(constant_value), is.numeric(decimals))
  }

  # 2) rebase all EBV in 'ebv' -------------------------------------------------

  # __ 2.1.A) if base_pop is provided
  if(is.null(constant_value) & !is.null(base_pop)) {
    # subset data to retain only base_pop IDs
    data_bp <- merge(data, base_pop, by.x = colnames(data)[1], by.y = colnames(base_pop)[1])
    # make an empty copy of data & fill it with NA's
    new_data <- data
    new_data[, -1] <- NA
  } else if(
    # __ 2.1.B) if constant_value is provided
    !is.null(constant_value) & is.null(base_pop)) {
    # make an empty copy of data & fill it with NA's
    new_data <- data
    new_data[, -1] <- NA
  }

  if (verbose == T & is.null(constant_value)) {
    cat(
      " N. individuals in data:", nrow(data), "\n",
      "N. individuals in base population:", nrow(data_bp)
    )
  }

  # __2.2) loop for all columns besides ID col (== col1)
  for (ebv_col in colnames(data)[-1]) {

    # __2.2.A) when base_pop is provided
    if(is.null(constant_value) & !is.null(base_pop)){
      # take mean base pop
      mean_basepop <- mean(data_bp[, ebv_col])
      # subtract from each EBV col of 'data' the mean EBV of the base_pop group of animals for the corresponding column

      new_data[, ebv_col] <- data[, ebv_col] - mean_basepop
      new_mean <- mean(merge(new_data, base_pop, by.x = colnames(new_data)[1], by.y = colnames(base_pop)[1])[, ebv_col])

      if (verbose == T) {
        cat(
          "\n\n Mean EBV of base pop. individulas for col. ", ebv_col, " before rebasing = ", round(mean_basepop, 4), "\n",
          "Mean EBV of base pop. individulas for col. ", ebv_col, " after rebasing = ", round(new_mean, 4), "\n"
        )
      }

      # check that the mean re-based EBV for base_pop is now 0.
      #   some rounding is involved, therefore, decimals are provided next to error.
      if (round(new_mean, decimals) != 0) {
        stop(
          "The new mean ebv of base_pop animals is not 0 (value was rounded to ", decimals, " decimals. Current new mean EBV for this group is: ", new_mean,
          "Check your data for possible errors or decrease the rounding precision in 'decimals' arg."
        )
      }
    } else if(
      # __ 2.2.B) when constant_value is provided
      !is.null(constant_value) & is.null(base_pop)) {
      mean_basepop <- constant_value
      new_data[, ebv_col] <- data[, ebv_col] - mean_basepop
    }

    # 2.3) check that correlation and slope between re-based EBV and provided (before re-basing) EBV is 1
    test_cor <- cor(data[, ebv_col], new_data[, ebv_col])
    reg <- lm(new_data[, ebv_col] ~ data[, ebv_col])
    reg_coef <- as.numeric(coef(reg))
    if (round(test_cor, decimals) != 1) {
      stop("Error: correlation of EBV before and after rebasing is not 1 for column: ", ebv_col, "\n\n Correlation of EBV is: ", test_cor)
    }
    if (round(reg_coef[2], decimals) != 1) {
      stop("Error: slope of EBV before and after rebasing is not 1 for column: ", ebv_col, "\n\n Slope of EBV rebased on EBV before rebasing is: ", reg_coef[2])
    }
    # check that intercept has changed by value of - mean of base POP
    if (round(reg_coef[1], 4) != -round(mean_basepop, 4)) { # round to 4 because of how precision of lm coeff.
      stop("Error: rebased EBV intercept from regressing EBV after on EBV before reabsing is not equal to - the mean EBV of the base POP.\  Intercept of EBV rebased on EBV before rebasing is:", round(reg_coef[1], 4))
    }
  }

  # if plot is not requested (default), results is a data.frame
  if (plot == F) {
    return(new_data)
  } else if (plot == T) { # if plot is requested, results are a list.
    # make a plot of EBVs before and after re-basing
    if (!suppressWarnings(require("ggplot2"))) {
      stop('Please install the "ggplot2" package, using \n:
    "install.packages("ggplot2")"')
    }
    results <- list(rebased_ebv = new_data)

    for (ebv_col in colnames(data)[-1]) { # loop for all columns besides ID col (== col1)
      reg <- lm(new_data[, ebv_col] ~ data[, ebv_col])
      pdat <- as.data.frame(cbind(data[, ebv_col], new_data[, ebv_col]))
      results[["plot"]][[ebv_col]] <- ggplot(pdat, aes(x = V1, y = V2)) +
        geom_point(colour = "black") +
        geom_abline(intercept = 0, slope = 1, col = "black") + # line passing through 0
        geom_abline(intercept = coef(reg)[1], slope = coef(reg)[2], col = "blue") +
        theme_bw()+
        xlab("original EBV") +
        ylab("rebased EBV") +
        labs(title = ebv_col)

      if(!is.null(constant_value)) { # add a line on y-axis showing the shift by constant_value
        results[["plot"]][[ebv_col]] <- results[["plot"]][[ebv_col]] +
          geom_hline(yintercept = mean(round(reg_coef[1], 4)), col = "red", linetype = "dotted")
      }
    }
    return(results)
  }
}
