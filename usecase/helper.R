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

peSkeleton = function(dat, feature, servers = NULL, npoints = 100L) {
  x = dat[[feature]]
  d0 = data.frame(val = seq(min(x), max(x), length.out = npoints), pred = NA)
  names(d0) = c(feature, "pred")
  if (! is.null(servers)) {
    d0 = do.call(rbind, lapply(servers, function(s) cbind(d0, server = s)))
  }
  return(d0)
}

createDFSkeleton = function(x, dat, feature) {
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

getMGCVPE = function(mod_mgcv, feature, sitevar = NULL, bpattern = NULL, x = NULL) {
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
  if (is.null(bpattern)) {
    if (is.null(sitevar)) {
      bpattern = "s[(]"
    } else {
      bpattern = "ti[(]"
    }
  }
  fidx = which(grepl(feature, cnames) & grepl(bpattern, cnames))
  if ((! is.numeric(d0[[feature]])) && is.null(sitevar)) fidx = which(cnames == feature)
  pe = unname(p[, fidx])

  pedf = data.frame(x = d0[[feature]], pred = pe, method = "mgcv")
  if (! is.null(sitevar)) {
    pedf$server = d0[[sitevar]]
  }
  names(pedf)[1] = feature
  return(pedf)
}

getCompboostPE = function(mod_compboost, feature, sitevar = NULL, bpattern = NULL, x = NULL) {
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


getDistCWBPE = function(mod_dsCWB, feature, sitevar = NULL, bpattern = NULL, x = NULL) {
  getCompboostPE(mod_dsCWB, feature, sitevar, bpattern, x)
}

addEffects = function(dshared, dsites, key, keys, name_add_val = "pred") {
  dagg = do.call(rbind, lapply(keys, function(k) {
    dsub = dsites[dsites[[key]] == k, ]
    dsub[[name_add_val]] = dsub[[name_add_val]] + dshared[[name_add_val]]
    return(dsub)
  }))
  return(dagg)
}

fVizData = function(feature, mod_dsCWB, mod_compboost, mod_mgcv, site = FALSE, add_effects = FALSE, num_points = 100L) {
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

fViz = function(feature, mod_dsCWB, mod_compboost, mod_mgcv, site = FALSE, add_effects = FALSE,
  num_points = 100L, plot = TRUE)  {

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

fVizOld = function(feature, mod_dsCWB, mod_compboost, mod_mgcv, site = FALSE, dat, add_effects = FALSE, svar = "source")  {


  cboost_names = mod_compboost$getBaselearnerNames()

  fnum = is.numeric(dat[[feature]])
  if (site || add_effects) {
    if (missing(mod_dsCWB)) {
      d0 = mod_compboost$data
      pred = peSkeleton(d0, feature, unique(d0[[svar]]))
    } else {
      if (fnum) {
        pred = siteFEDataNum(mod_dsCWB, feature)
      } else {
        pred = siteFEDataCat(mod_dsCWB, feature)
      }
    }

    pred$method = "dsCWB"
    dnew = pred[, c(feature, "server")]
    names(dnew)[2] = "source"
    xnew = suppressWarnings(mod_compboost$prepareData(dnew))

    bln = cboost_names[grep(feature, cboost_names)]
    bln = bln[grepl("_tensor", bln)]
    pred_compboost = try({
      y_compboost = mod_compboost$model$predictFactoryNewData(bln, xnew)
      #if (inherits(y_compboost, "try-error")) y_compboost = 0
      if (fnum) {
        x = xnew[[feature]]$getData()
      } else {
        x = xnew[[1]]$getRawData()
      }
      s = xnew[["source"]]$getRawData()
      pred_compboost = data.frame(x = x, y = y_compboost, server = s, method = "compboost")

      names(pred_compboost) = names(pred)

      pred_compboost
    }, silent = TRUE)
    if (inherits(pred_compboost, "try-error")) {
      pred_compboost = pred
      pred_compboost$pred = 0
      pred_compboost$method = "compboost"
    }

    # FIXME FOR CATEGORICAL FEATURES
    pred_mgcv = getMGCVPE(mod_mgcv, feature = feature, site = TRUE)
    names(pred_mgcv) = names(pred_compboost)

    pred = na.omit(rbind(pred, pred_compboost, pred_mgcv))
    if (add_effects) psite = pred

    if (! add_effects) {
      gg = ggplot(pred, aes_string(x = feature, y = "pred", color = "method", linetype = "method")) +
        geom_line() +
        facet_wrap(~ server, ncol = 2)
    }

  }

  if ((! site) || add_effects) {
    if (missing(mod_dsCWB)) {
      d0 = mod_compboost$data
      pred = peSkeleton(d0, feature)
    } else {
      if (fnum) {
        pred = sharedFEDataNum(mod_dsCWB, feature)
      } else {
        pred = sharedFEDataCat(mod_dsCWB, feature)
      }
    }
    pred$method = "dsCWB"
    xnew = suppressWarnings(mod_compboost$prepareData(pred[, feature, drop = FALSE]))

    pred_compboost = try({
      bln = cboost_names[grep(feature, cboost_names)]
      bln = bln[! grepl("_tensor", bln)]
      y_compboost = mod_compboost$model$predictFactoryNewData(bln, xnew)
      if (fnum) {
        x = xnew[[1]]$getData()
      } else {
        x = xnew[[1]]$getRawData()
      }
      pred_compboost = data.frame(x = x, y = y_compboost, method = "compboost")
      names(pred_compboost) = names(pred)

      pred_compboost
    }, silent = TRUE)

    if (inherits(pred_compboost, "try-error")) {
      pred_compboost = pred
      pred_compboost$pred = 0
      pred_compboost$method = "compboost"
    }

    # FIXME FOR CATEGORICAL FEATURES
    pred_mgcv = getMGCVPE(mod_mgcv, feature = feature, FALSE)
    names(pred_mgcv) = names(pred_compboost)

    pred = na.omit(rbind(pred, pred_compboost, pred_mgcv))

    if (! add_effects) {
      gg = ggplot(pred, aes_string(x = feature, y = "pred", color = "method", linetype = "method")) +
      geom_line()
    }

    if (add_effects) pshared = pred
  }
  if (add_effects) {
    psite_cwb = psite[psite$method != "mgcv", ]
    pshared_cwb = pshared[pshared$method != "mgcv", ]

    psi_mgcv = psite[psite$method == "mgcv", ]
    psh_mgcv = pshared[pshared$method == "mgcv", ]

    for (s in unique(psite$server)) {
      psite_cwb$pred[psite_cwb$server == s] = psite_cwb$pred[psite_cwb$server == s] + pshared_cwb$pred
    }
    psi_mgcv$pred = psi_mgcv$pred + psh_mgcv$pred

    psite = rbind(psite_cwb, psi_mgcv)

    gg = ggplot(psite, aes_string(x = feature, y = "pred", color = "method", linetype = "method")) +
      geom_line() +
      facet_wrap(~ server, ncol = 2)
  }
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


fVizCategorical = function(feature, mod_dsCWB, mod_compboost, mod_mgcv) {
  cds = mod_dsCWB$coef

}



