#Single dataset simulations install

install_benchmark_packages <- function() {
  cran_pkgs <- c(
    "ggplot2", "dplyr", "tidyr", "purrr", "tibble",
    "Matrix", "pls", "pROC", "mclust", "irlba", "glmpca", "reticulate"
  )
  
  missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
  if (length(missing_cran) > 0) {
    install.packages(missing_cran)
  }
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  bioc_pkgs <- c("MOFA2", "sva", "RUVSeq", "mixOmics", "limma", "DESeq2", "edgeR")
  
  missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc)
  }
}

install_benchmark_packages()
