PlotConvergeneMiXBLUP <- function(ConvFilePath, solver = "solver", title = NA) {
  #################### DOCUMENTATION ###########################################
  #' Plot MiXBLUP model convergence
  #'
  #' @description This function returns a ggplot with different convergences of MiXBLUP.
  #'
  #' @param ConvFilePath (character) pathway where convergence file is located.
  #' @param solver (character) Either "solver" (default) or "hpblup. "The default solver is MiXBLUP default solver,i.e. solver.exe with associated convergence file
  #' named "Conlog.txt". Solver can be set to "hpblup" for hpblup.exe which has associated convergence file
  #' named "convergence.dat"
  #' @param title (character) title for convergence plot. Default is NA.
  #' 
  #' @usage
  #' PlotConvergeneMiXBLUP(ConvFilePath = "RUNS/PBLUP")
  #' PlotConvergeneMiXBLUP(ConvFilePath = "RUNS/ssSNPBLUP", solver="hpblup", title="single-step convergence")
  #'
  #' @returns
  #' An R ggplot2 plot object
  #' 
  #' @importFrom ggplot2 dplyr tidyr
  #' 
  #' @details 
  #' Column names are assumed to be fixed.
  #' Color vectors are manually defined.
  #' 
  #' @note  
  #' TODO:
  #' - colors can be better manipulated, now is quite descriptive,
  #' - legend lines can be sorted by custom order (e.g. CD criteria first etc.)
  #' - better way to read column names for CD, CK, ecc. from files
  ##############################################################################
  if (solver == "hpblup") {   # mixblup solver is: hpblup.exe
    ConvFileExt <- "/convergence.dat"
    data <- read.delim(paste0(ConvFilePath, ConvFileExt), sep = "")
    data <- data[, c(2, 5, 8:10, 12)]
    colnames(data)[1:4]<- c("Iteration", "tau", "CR^2", "CD^2")
    col_vec <- c("log(CR^2)" = "purple", "log(CD^2)" = "forestgreen", "log(CK)" = "dodgerblue", "log(CM)" = "orange",  "log(tau)" = "gold") # line colors
  
    } else if (solver == "solver") { # mixblup solver is: solver.exe
    ConvFileExt <- "/Conlog.txt"
    data <- read.table(paste0(ConvFilePath, ConvFileExt), skip = 10, fill = NA) # skip all header part
    data <- data[, c(1:5)]
    colnames(data) <- c("Iteration", "CA", "CR", "CM", "CD")
    col_vec <- c("log(CR)" = "purple", "log(CD)" = "forestgreen", "log(CA)" = "dodgerblue", "log(CM)" = "orange") # line colors
  }
  
  data <- data %>%
    mutate(across(.cols = -1, ~ log(.x)/log(10), .names = "log({.col})")) %>% # convert to log10 values
    select("Iteration", contains("log")) %>%
    pivot_longer(cols = -1, names_to = "Criteria", values_to = "value") # data in longer format for ggplot

  # plot
  p <- data %>%
    ggplot() +
    geom_line(aes(x = Iteration, y = value, color=Criteria)) +
    scale_color_manual(values = col_vec) +
    theme_bw()+
    labs(title = title, x = "Iterations", y = "Convergence") +
    theme(
      axis.title = element_text(color = "black", size = 12),
      axis.text = element_text(color = "black", size = 10),
      legend.position = "bottom"
    )
  
  return(p)
}
