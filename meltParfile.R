meltParfile <- function(parfile, storeVcov = NA, storeCorr = NA) {
  #################### DOCUMENTATION ###########################################
  #' Melt a parameter file (e.g. MiX99) in covariances and correlation matrices.
  #'
  #' @description This function re-arrange a parameter file in the format of MiX99
  #' in matrices of covariances and correlations. See below for MiX99 parameter file
  #' format.
  #'
  #' @param parfile data.frame or matrix in MiX99 parameter file format:
  #'             - col1: effect number (integer)
  #'             - col2: row number    (integer)
  #'             - col3: column number (integer)
  #'             - col4: (co)variance  (numeric)
  #' @param storeVcov If you want to store the covariances matrices, provide a string (character). Default is NA.
  #' @param storeCorr If you want to store the correlation matrices, provide a string (character). Default is NA.
  #'
  #' @usage
  #' meltParfile(parfile=parfile, storeVcov="covar_effectN_")
  #' 
  #' meltParfile(parfile=parfile, storeVcov="covar_effectN_", storeCorr="corr_effectN_")
  #'
  #' @returns
  #' An R list with covariances and correlation symmetric matrices (one per effect)
  #'
  #' @references
  #' This is an adapted function from Mohammad Ali Nilforooshan's github script
  #' (https://gist.github.com/nilforooshan/874ec62d2ccafc26b66a89d0bed1d70c)
  #'
  ##############################################################################
  list_res <- list() # store results in a list
  for (k in 1:nlevels(as.factor(parfile$V1))) { # melt parfile by effect number
    tmp <- split(parfile, parfile$V1)[[k]]
    size <- max(tmp[, 2:3])
    vcov <- matrix(0, nrow = size, ncol = size)
    for (j in 1:size) {
      for (i in 1:size) {
        if (length(tmp$V4[tmp$V2 == i & tmp$V3 == j]) == 1) {
          vcov[i, j] <- vcov[j, i] <- tmp$V4[tmp$V2 == i & tmp$V3 == j]
        }
      }
    }
    corr <- cov2cor(vcov)
    if (!is.na(storeVcov)) {
      list_res[[paste(storeVcov, k, sep = "")]] <- vcov
    }
    if (!is.na(storeCorr)) {
      list_res[[paste(storeCorr, k, sep = "")]] <- corr
    }
  }
  return(list_res)
}
