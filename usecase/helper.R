#' Read the heart disease data.
#' @param file (`character(1L)`) Path to the csv file.
#' @param add_source (`logical(1L)`) Indicator whether the site should be added as variable.
#' @param add_id (`logical(1L)`) Add a column with ids.
#' @param rm_pcols (`logical(1L)`) Remove problematic columns num, thal, ca, slope, fbs, and chol.
#' @param add_sim_col (`logical(1L)`) Add a column with random noise.
#' @param process (`logical(1L)`) Remove 'suspicious' data.
readData = function(file, add_source = FALSE, add_id = FALSE, rm_pcols = TRUE, add_sim_col = FALSE,
                    process = TRUE) {
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

  if (process) {
    idx_remove_trestpbs = tmp$trestbps != 0
    tmp = tmp[idx_remove_trestpbs, ]
  }

  return(tmp)
}

#' Create a data skeleton with varying `feature` values `x` and constant values for all other
#' variables in `dat` (numeric = mean, ordingal = mode).
#' @param x (`vector()`) Varying values.
#' @param dat (`data.frame()`) Data frame with all variables (used to calculate constants).
#' @param feature (`character(1L)`) Name of the varying variable.
createDFSkeleton = function(x, dat, feature) {
  checkmate::assertVector(x)
  checkmate::assertDataFrame(dat)
  checkmate::assertChoice(feature, colnames(dat))

  n = length(x)
  idx = seq_len(n)
  d0 = do.call(data.frame, lapply(dat, function(ff) {
    if (is.numeric(ff)) {
      rep(mean(ff), n)
    } else {
      tx = table(ff)
      max_idx = which.max(tx)
      rep(names(tx)[max_idx], n)
    }
  }))
  d0[[feature]] = x
  return(d0)
}

#' Create a data.frame with partial effects of a `gam` model produced by mgcv. Column
#' names are `feature`, pred, method (= "mgcv"), and if `! is.null(sitevar)` server.
#' @param mod_mgcv (`gam`) GAM model from `mgcv`.
#' @param feature (`character(1L)`) Name of the varying variable.
#' @param sitevar (`character(1L)`) Name of the variable containing the site.
#' @param bpattern (`character(1L)`) Regexp to filter for a specific base learner.
#' @param x (`vector()`) Varying values.
getMGCVPE = function(mod_mgcv, feature, sitevar = NULL, bpattern = NULL, x = NULL) {
  checkmate::assertClass(mod_mgcv, "gam")
  checkmate::assertChoice(feature, colnames(mod_mgcv$model))
  checkmate::assertChoice(sitevar, colnames(mod_mgcv$model), null.ok = TRUE)
  checkmate::assertCharacter(bpattern, len = 1L, null.ok = TRUE)
  checkmate::assertVector(x, null.ok = TRUE)

  if (is.null(x)) {
    d0 = mod_mgcv$model
  } else {
    d0 = createDFSkeleton(x, mod_mgcv$model, feature)
  }
  if (! is.null(sitevar)) {
    snames = unique(mod_mgcv$model[[sitevar]])
    d0 = do.call(rbind, lapply(snames, function(sn) {
      d1 = d0
      d1[[sitevar]] = sn
      return(d1)
    }))
  }
  p = predict(mod_mgcv, type = "terms", newdata = d0)

  cnames = colnames(p)
  if ((! is.numeric(d0[[feature]])) && is.null(sitevar)) {
    fidx = which(cnames == feature)
  } else {
    if (is.null(sitevar)) {
#      bpattern = "s[(]"
      fidx = which(grepl(feature, cnames) & grepl("s[()]", cnames))
    } else {
      cnames0 = gsub(" ", "", cnames)
      site_pattern = paste0(",", sitevar)
      fidx = which(grepl(feature, cnames) & grepl(site_pattern, cnames0))
    }
  }
  pe = unname(p[, fidx])

  pedf = data.frame(x = d0[[feature]], pred = pe, method = "mgcv")
  if (! is.null(sitevar)) {
    pedf$server = d0[[sitevar]]
  }
  names(pedf)[1] = feature
  return(pedf)
}

#' Create a data.frame with partial effects of a CWB model produced by compboost. Column
#' names are `feature`, pred, method (= "compboost"), and if `! is.null(sitevar)` server.
#' @param mod_compboost (`Compboost`) CWB model from `compboost`.
#' @param feature (`character(1L)`) Name of the varying variable.
#' @param sitevar (`character(1L)`) Name of the variable containing the site.
#' @param bpattern (`character(1L)`) Regexp to filter for a specific base learner.
#' @param x (`vector()`) Varying values.
getCompboostPE = function(mod_compboost, feature, sitevar = NULL, bpattern = NULL, x = NULL) {
  checkmate::assertR6(mod_compboost, "Compboost")
  checkmate::assertChoice(feature, colnames(mod_compboost$data))
  checkmate::assertChoice(sitevar, colnames(mod_compboost$data), null.ok = TRUE)
  checkmate::assertCharacter(bpattern, len = 1L, null.ok = TRUE)
  checkmate::assertVector(x, null.ok = TRUE)

  xnew = x
  if (is.null(x)) {
    xnew  = mod_compboost$data[[feature]]
  }
  d0 = data.frame(x = x)
  names(d0) = feature

  if (! is.null(sitevar)) {
    snames = unique(mod_compboost$data[[sitevar]])
    d0 = do.call(rbind, lapply(snames, function(sn) {
      d1 = d0
      d1[[sitevar]] = sn
      return(d1)
    }))
  }
  xnew = suppressWarnings(mod_compboost$prepareData(d0))
  blnames = mod_compboost$getBaselearnerNames()

  bln = blnames[grep(feature, blnames)]
  if (is.null(bpattern)) {
    if (is.null(sitevar)) {
      bln = bln[! grepl("_tensor", bln)]
    } else {
      bln = bln[grepl("_tensor", bln)]
    }
  } else {
    bln = bln[grepl(bpattern, bln)]
  }

  pe = try(as.vector(mod_compboost$model$predictFactoryNewData(bln, xnew)), silent = TRUE)
  if (inherits(pe, "try-error")) pe = rep(0, length(x))


  pedf = data.frame(x = d0[[feature]], pred = pe, method = "compboost")
  if (! is.null(sitevar)) {
    s = xnew[[sitevar]]$getRawData()
    pedf$server = d0[[sitevar]]
  }
  names(pedf)[1] = feature
  return(pedf)
}

#### TODO
getDistCWBPE = function(mod_dsCWB, feature, sitevar = NULL, bpattern = NULL, x = NULL) {
  getCompboostPE(mod_dsCWB, feature, sitevar, bpattern, x)
}

#' Create a data.frame with partial effects of a CWB model produced by compboost. Column
#' names are `feature`, pred, method (= "compboost"), and if `! is.null(sitevar)` server.
#' @param mod_compboost (`Compboost`) CWB model from `compboost`.
#' @param feature (`character(1L)`) Name of the varying variable.
#' @param sitevar (`character(1L)`) Name of the variable containing the site.
#' @param bpattern (`character(1L)`) Regexp to filter for a specific base learner.
#' @param x (`vector()`) Varying values.
addEffects = function(dshared, dsites, key, keys, name_add_val = "pred") {
  checkmate::assertDataFrame(dshared)
  checkmate::assertChoice(key, colnames(dsites))
  checkmate::assertCharacter(keys, min.len = 1L)
  checkmate::assertDataFrame(dsites, nrows = nrow(dshared) * length(keys))
  if (! is.factor(dsites[[key]])) dsites[[key]] = as.factor(dsites[[key]])
  nuisance = lapply(keys, checkmate::assertChoice, choices = levels(dsites[[key]]))

  dagg = do.call(rbind, lapply(keys, function(k) {
    dsub = dsites[dsites[[key]] == k, ]
    dsub[[name_add_val]] = dsub[[name_add_val]] + dshared[[name_add_val]]
    return(dsub)
  }))
  return(dagg)
}

#' Create a data.frame with partial effects for the three methods "distributed CWB" (`dsCWB`),
#' "CWB" (`compboost`), and "GAMM" (`mgcv`) . Names are `feature`, pred, method, and
#' if `! is.null(sitevar)` server.
#' @param feature (`character(1L)`) Name of the varying variable.
#' @param mod_dsCWB (`dsCWB`) Distributed CWB model from `dsCWB`.
#' @param mod_compboost (`Compboost`) CWB model from `compboost`.
#' @param mod_mgcv (`gam`) GAM model from `mgcv`.
#' @param site (`logical(1L)`) Indicator whether the PEs are site specific or not.
#' @param add_effects (`logical(1L)`) Should shared and site-specific effects added up to the overall effect?
#' @param num_points (`integer(1L)`) Number of points used to create the grid for numerical features.
fVizData = function(feature, mod_dsCWB, mod_compboost, mod_mgcv, site = FALSE, add_effects = FALSE, num_points = 100L) {
  checkmate::assertChoice(feature, colnames(mod_mgcv$model))
  if (! missing(mod_dsCWB)) checkmate::assertR6(mod_dsCWB, "dsCWB")
  checkmate::assertR6(mod_compboost, "Compboost")
  checkmate::assertClass(mod_mgcv, "gam")
  checkmate::assertLogical(site, len = 1L)
  checkmate::assertLogical(add_effects, len = 1L)
  checkmate::assertCount(num_points)

  dat = mod_mgcv$model
  x = dat[[feature]]
  if (is.numeric(x)) {
    xnew = seq(min(x), max(x), length.out = num_points)
  } else {
    xnew = unique(x)
  }
  d0 = createDFSkeleton(x, dat, feature)

  svar_mgcv = "src"
  svar_compboost = "source"
  if (add_effects) {
    pe_mgcv_shared = getMGCVPE(mod_mgcv, feature, sitevar = NULL, x = xnew)
    pe_mgcv_site = getMGCVPE(mod_mgcv, feature, sitevar = svar_mgcv, x = xnew)

    skeys = levels(pe_mgcv_site$server)

    pe_mgcv = addEffects(pe_mgcv_shared, pe_mgcv_site, key = "server", keys = skeys)

    pe_compboost_shared = getCompboostPE(mod_compboost, feature, sitevar = NULL, x = xnew)
    pe_compboost_site = getCompboostPE(mod_compboost, feature, sitevar = svar_compboost, x = xnew)
    pe_compboost = addEffects(pe_compboost_shared, pe_compboost_site, key = "server", keys = skeys)

  } else {
    if (! site) {
      svar_mgcv = NULL
      svar_compboost = NULL
    }
    pe_mgcv = getMGCVPE(mod_mgcv, feature, sitevar = svar_mgcv, x = xnew)
    pe_compboost = getCompboostPE(mod_compboost, feature, sitevar = svar_compboost, x = xnew)
  }
  return(rbind(pe_mgcv, pe_compboost))
}

#' Create a plot for the partial effects for the three methods "distributed CWB" (`dsCWB`),
#' "CWB" (`compboost`), and "GAMM" (`mgcv`). A ggplot is returned, if `plot = FALSE` the
#' partial effects data (see `fVizData`) are returned.
#' @param feature (`character(1L)`) Name of the varying variable.
#' @param mod_dsCWB (`dsCWB`) Distributed CWB model from `dsCWB`.
#' @param mod_compboost (`Compboost`) CWB model from `compboost`.
#' @param mod_mgcv (`gam`) GAM model from `mgcv`.
#' @param site (`logical(1L)`) Indicator whether the PEs are site specific or not.
#' @param add_effects (`logical(1L)`) Should shared and site-specific effects added up to the overall effect?
#' @param num_points (`integer(1L)`) Number of points used to create the grid for numerical features.
fViz = function(feature, mod_dsCWB, mod_compboost, mod_mgcv, site = FALSE, add_effects = FALSE,
  num_points = 100L, plot = TRUE)  {

  checkmate::assertChoice(feature, colnames(mod_mgcv$model))
  if (! missing(mod_dsCWB)) checkmate::assertR6(mod_dsCWB, "dsCWB")
  checkmate::assertR6(mod_compboost, "Compboost")
  checkmate::assertClass(mod_mgcv, "gam")
  checkmate::assertLogical(site, len = 1L)
  checkmate::assertLogical(add_effects, len = 1L)
  checkmate::assertCount(num_points)
  checkmate::assertLogical(plot, len = 1L)

  pe_data = fVizData(feature, mod_dsCWB, mod_compboost, mod_mgcv, site, add_effects, num_points)

  if (! plot) return(invisible(pe_data))
  if (is.numeric(pe_data[[feature]])) {
    if (site || add_effects) {
      gg = ggplot(pe_data, aes_string(x = feature, y = "pred", color = "method", linetype = "method")) +
        geom_line() +
        facet_wrap(~ server, ncol = 2)
    } else {
      gg = ggplot(pe_data, aes_string(x = feature, y = "pred", color = "method", linetype = "method")) +
        geom_line()
    }
  } else {
    if (site || add_effects) {
      gg = ggplot(pe_data, aes_string(x = feature, y = "pred", color = "method", shape = "method")) +
        geom_point(position = position_dodge(0.2)) +
        facet_wrap(~ server, ncol = 2)
    } else {
      gg = ggplot(pe_data, aes_string(x = feature, y = "pred", color = "method", shape = "method")) +
        geom_point(position = position_dodge(0.2))
    }
  }
  attr(gg, "pe_data") = pe_data
  return(gg)
}

mytheme = function(base_size = 14) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(size = rel(1), face = "bold", margin = margin(0,0,5,0), hjust = 0),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size = rel(0.85), face = "bold"),
      axis.text = element_text(size = rel(0.70), face = "bold"),
      axis.line = element_line(color = "black", arrow = arrow(length = unit(0.3, "lines"), type = "closed")),
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.70), face = "bold"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(1.5, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "#17252D", color = "#17252D"),
      strip.text = element_text(size = rel(0.85), face = "bold", color = "white", margin = margin(5,0,5,0))
    )
}
