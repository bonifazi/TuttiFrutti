PlotConvergeneMiXBLUP <- function(ConvFilePath, solver = "solver", title = NA) {
  
  if (solver == "hpblup") {   # mixblup solver is: hpblup.exe
    ConvFileExt <- "/convergence.dat"
    data <- read.delim(paste0(ConvFilePath, ConvFileExt), sep = "")
    data <- data[, c(2, 5, 8:10, 12)]
    colnames(data)[1:4]<- c("Iteration", "tau", "CR^2", "CD^2")
    col_vec <- c("log(CR^2)" = "purple", "log(CD^2)" = "forestgreen", "log(CK)" = "dodgerblue", "log(CM)" = "orange",  "log(tau)" = "gold") # line colors
  
    } else { # mixblup solver is: solver.exe
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
