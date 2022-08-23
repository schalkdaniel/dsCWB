#' Read the heart disease data.
#' @param file (`character(1L)`) Path to the csv file.
#' @param add_source (`logical(1L)`) Indicator whether the site should be added as variable.
#' @param add_id (`logical(1L)`) Add a column with ids.
#' @param rm_pcols (`logical(1L)`) Remove problematic columns num, thal, ca, slope, fbs, and chol.
#' @param add_sim_col (`logical(1L)`) Add a column with random noise.
readData = function(file, add_source = FALSE, add_id = FALSE, rm_pcols = TRUE, add_sim_col = FALSE) {
  if (grepl("reprocessed", file))
    tmp = read.csv(file, sep = " ", header = FALSE, na.strings = c("?", "-9", "-9.0"))
  else
    tmp = read.csv(file, header = FALSE, na.strings = c("?", "-9", "-9.0"))

  cnames = c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg", "thalach",
    "exang", "oldpeak", "slope", "ca", "thal", "num")

  fvals = c("sex", "cp", "fbs", "restecg", "exang", "slope", "ca", "thal", "num")
  colnames(tmp) = cnames

  #for (fv in cnames) {
    #var = tmp[[fv]]
    #idx_msg = (var == "?") | (var == "-9") | (var == "-9.0")
    #var[idx_msg] = NA
    #tmp[[fv]] = var
  #}
  for (fv in fvals) {
    tmp[[fv]] = as.factor(as.integer(as.character(tmp[[fv]])))
  }
  if (add_source) tmp$source = strsplit(file, "[.]")[[1]][2]
  if (add_id) tmp$id = seq_len(nrow(tmp))

  tmp$heart_disease = ifelse(tmp$num == "1", "yes", "no")

  if (rm_pcols) {
    tmp$num = NULL
    tmp$thal = NULL # 486 missings in total
    tmp$ca = NULL # 611 missings in total
    tmp$slope = NULL # 309 missings in total
    tmp$fbs = NULL # Too many missings for switzerland
    tmp$chol = NULL # Only zeros for switzerland ... :-(
  }

  if (add_sim_col) {
    x = sort(runif(nrow(tmp), 0, 100))
    xn = rep(0, nrow(tmp))
    idx_yes = tmp$heart_disease == "yes"
    idx_yes = ifelse(rbinom(nrow(tmp), 1, 0.1), ! idx_yes, idx_yes)
    xn[idx_yes] = tail(x, sum(idx_yes))
    xn[! idx_yes] = head(x, sum(! idx_yes))
    tmp$simcol = xn
  }
  return(tmp)
}


