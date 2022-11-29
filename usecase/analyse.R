################################################################################
##                                                                            ##
##                             DISTRIBUTED CWB                                ##
##                                                                            ##
################################################################################

USE_REAL_TEST_ENV = FALSE
GGSAVE = TRUE

## SETUP DataSHIELD ENVIRONMENT
## ========================================================================= ##

if (USE_REAL_TEST_ENV) {

  ## Real DataSHIELD server:
  ## =================================================================
  # This code uploads the data to the DataSHIELD test server. The
  # advantage is that the code is tested on a real DataSHIELD instance.
  # The drawback is the time it takes to fit a model.

  source(here::here("usecase/upload-data.R"))

  library(DSI)
  library(DSOpal)
  library(dsBaseClient)

  builder = newDSLoginBuilder()

  surl     = "https://opal-demo.obiba.org/"
  username = "administrator"
  password = "password"

  sources = c("cleveland", "hungarian", "switzerland", "va")
  for (i in seq_along(sources)) {
    builder$append(
      server   = sources[i],
      url      = surl,
      user     = username,
      password = password,
      table    = paste0("SLDS-TEST.", sources[i])
    )
  }

  ## Get data of the servers:
  connections = datashield.login(logins = builder$build(), assign = TRUE)
} else {

  ## DSLite test server:
  ## =================================================================
  # This code allocates a local DataSHIELD test instance with DSLite.
  # It is much faster than the real server but behaves slightly different
  # than a real DataSHIELD instance because of the JAVA parser in the
  # middle.

  source(here::here("usecase/setup-dslite.R"))

  library(dsBaseClient)
}

datashield.symbols(connections)

## Data dimensions per server:
(ddim = ds.dim("D"))

## TRAIN THE MODEL
## ========================================================================= ##


#library(dsCWB)
devtools::load_all()

df = 2.2
df_ri = 3
n_knots = 10L
mstop = 100000L
pat = 5L
lr = 0.1

symbol = "D"
target = "heart_disease"
feature_names = setdiff(dsBaseClient::ds.names("D", connections)[[1]], target)

cwb = dsCWB(connections, symbol, target, feature_names, mstop = mstop,
  val_fraction = 0.2, patience = pat, seed = 31415L, positive = "yes",
  df = df, learning_rate = lr, random_intercept_df = df_ri, nknots = n_knots)


## FIGURES
## ========================================================================== ##

if (FALSE) {
  ## Values for the paper:
  pCommand = function(cmd, value, digits = 4)
    paste0("\\newcommand{", cmd, "}{", round(value, digits), "\\xspacekkkj}\n")

  # Get log data with
  # - selected base learner
  # - train and test risk
  # - effect type (site or shared)
  l = cwb$getLog()
  table(l$bl)
  table(l$effect_type)
  table(paste0(l$bl, "-", l$effect_type))

  tab = table(l$effect_type[-1])

  cat(pCommand("\\mstop", max(l$iteration)))
  cat(pCommand("\\nselshared", tab["shared"]))
  cat(pCommand("\\nselsite", tab["site"]))
  cat(pCommand("\\patience", pat))
  cat(pCommand("\\dfMain", df))
  cat(pCommand("\\dfRI", df_ri))
  cat(pCommand("\\learningrate", lr))
  cat(pCommand("\\nknots", n_knots))

  library(ggplot2)
  library(patchwork)

  ## Base learner traces:
  ## =================================================================

  gg_traces_effects = plotBaselearnerTraces(cwb, n_legend = 4L, add_effect_type = TRUE, pretty_names = TRUE) +
    ggsci::scale_color_jama() +
    ggsci::scale_fill_jama() +
    mytheme() +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    ylab("Proportion of\nadded base learners") +
    guides(fill = guide_legend(ncol = 1))

  gg_traces_effects

  ## Risk traces:
  ## =================================================================

  gg_risk = ggplot(l, aes(x = iteration)) +
    geom_line(aes(y = risk_train, color = "Train risk")) +
    geom_line(aes(y = risk_val, color = "Validation risk")) +
    ggsci::scale_color_jama() +
    mytheme() +
    theme(legend.title = element_blank()) +
    ylab("Risk") +
    xlab("Iteration")

  gg_risk

  if (GGSAVE) {
    ggsave(gg_risk, filename = here::here("usecase/figures/gg-risk.pdf"), width = 13,
      height = 4.2, units = "cm")
  }

  ## Feature importance:
  ## =================================================================

  vip = function(risk, features) {
    agg = aggregate(risk, by = list(features), FUN = sum)
    ord = order(agg$x, decreasing = TRUE)
    agg = agg[ord, ]
    return(data.frame(feature = agg[[1]], vip = agg[[2]]))
  }
  getFName = function(fnames) {
    sapply(fnames, function(fname) {
      strsplit(fname, "-")[[1]][1]
    }, USE.NAMES = FALSE)
  }

  #vi = vip(l$risk_train[-1], paste0(l$bl[-1], "-", l$effect_type))
  vi = vip(l$risk_train[-1], getFName(l$bl[-1]))#, "-", l$effect_type))
  gg_vip = ggplot(vi, aes(x = reorder(feature, vip), y = vip)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    mytheme() +
    xlab("") +
    ylab("VIP") +
    theme(legend.position = "bottom")

  gg_traces_vip = gg_traces_effects + guide_area() + gg_vip +
    plot_layout(widths = c(2, 1.2, 2), guides = "collect")
  gg_traces_vip

  if (GGSAVE) {
    ggsave(gg_traces_vip, filename = here::here("usecase/figures/gg-traces-vip4.pdf"),
      width = 20, height = 5.3, units = "cm")
  }

  ## Partial feature effects:
  ## =================================================================

  peDCWBNum = function(feat) {
    fe_op_sh = sharedFEDataNum(cwb, feat)
    ggshared = ggplot(fe_op_sh, aes_string(x = feat, y = "pred")) +
      geom_line() +
      mytheme() +
      ylab("Partial feature\neffect") +
      ggtitle("Shared effect")

    fe_op_si = siteFEDataNum(cwb, feat)
    ggsite = ggplot(fe_op_si, aes_string(x = feat, y = "pred", color = "server")) +
      geom_line() +
      mytheme() +
      ggsci::scale_color_jama() +
      ylab("Partial feature\neffect") +
      ggtitle("Site effects") +
      theme(legend.title = element_blank())

    fe_op = fe_op_si
    fe_op$pred = fe_op$pred + rep(fe_op_sh$pred, 4)
    ggagg = ggplot(fe_op, aes_string(x = feat, y = "pred", color = "server")) +
      geom_line() +
      mytheme() +
      ggsci::scale_color_jama() +
      ylab("Partial feature\neffect") +
      ggtitle("Shared + site effects") +
      theme(legend.title = element_blank())

    gg = ggsite + ggagg #& theme(legend.position = "bottom")
    gg = gg + plot_layout(guides = "collect")

    gg = ggshared + gg + plot_layout(widths = c(1, 2))
    gg
  }

  peDCWBCat = function(feat) {
    fe_sh = sharedFEDataCat(cwb, feat)
    fe_sh$xnum = as.integer(as.factor(fe_sh[[feat]]))
    fe_op = siteFEDataCat(cwb, feat)

    snames = unique(fe_op$server)
    fe_op = do.call(rbind, lapply(snames, function(s) {
      d0 = fe_op[fe_op$server == s, ]
      d0$predadd = d0$pred + fe_sh$pred
      d0$pred = fe_sh$pred
      d0$xnum = as.integer(as.factor(d0[[feat]]))
      return(d0)
    }))

    xo = 0.3
    fe_op$xshift = seq(-xo, xo, length.out = 5)[-3][as.integer(as.factor(fe_op$server))]

    fe_sh$server = "shared"
    ggplot() +
      geom_hline(yintercept = 0, color = "dark gray") +
      geom_segment(data = fe_sh, aes(x = xnum, xend = xnum, y = 0, yend = pred, color = server), size = 0.3) +
      geom_segment(data = fe_sh, aes(x = xnum - xo, xend = xnum + xo, y = pred, yend = pred, color = server), size = 0.3) +
      geom_segment(data = fe_op, aes(x = xnum + xshift, xend = xnum + xshift,
        y = pred, yend = predadd, color = server), size = 0.3) +
      geom_point(data = fe_op, aes(x = xnum + xshift, y = predadd, color = server, shape = server), show.legend = FALSE) +
      coord_flip() +
      xlab(feat) +
      ylab("Coefficient") +
      mytheme() +
      ggsci::scale_color_jama() +
      labs(color = "") +
      scale_x_continuous(breaks = fe_sh$xnum, labels = fe_sh[[feat]])
  }



  gg_oldpeak2 = peDCWBNum("oldpeak")
  gg_oldpeak2

  if (GGSAVE) {
    ggsave(gg_oldpeak2, filename = here::here("usecase/figures/gg-oldpeak2.pdf"),
      width = 25, height = 4.5, units = "cm")
  }

  fnums = c("age", "trestbps", "thalach", "oldpeak")
  fcats = c("sex", "cp", "restecg", "exang")

  ggs = lapply(fnums, peDCWBNum)
  gg_nums = ggs[[1]] / ggs[[2]] / ggs[[3]] / ggs[[4]]

  if (GGSAVE) {
    ggsave(gg_nums, filename = here::here("usecase/figures/gg-pes-num.pdf"),
      width = 25, height = 29, units = "cm")
  }

  ggs = lapply(fcats, peDCWBCat)
  gg_cats = (ggs[[1]] + ggs[[2]]) / (ggs[[3]] + ggs[[4]]) & theme(legend.position = "bottom")
  gg_cats = gg_cats + plot_layout(guides = "collect")

  if (GGSAVE) {
    ggsave(gg_cats, filename = here::here("usecase/figures/gg-pes-cat.pdf"),
      width = 15, height = 17, units = "cm")
  }


  ## Individual contribution to prediction score:
  ## =================================================================

  newdata = data.frame(age = 60, sex = "1", cp = "4", trestbps = 90,
    restecg = "2", oldpeak = 3, thalach = 160, exang = "1")
  pinds = cwb$predictIndividual(newdata)


  preds = unlist(pinds$shared)
  dfinds = rbind(
    data.frame(pred = unname(preds), po = 0, bl = names(preds), server = "shared"),
    do.call(rbind, lapply(names(pinds$site), function(s) {
      preds = unlist(pinds$site[[s]])
      po = unlist(pinds$shared)
      data.frame(pred = unname(preds) + po, po = po, bl = names(preds), server = s)
    })))

  dfinds$bli = as.integer(as.factor(dfinds$bl))
  xo = 0.5
  dfinds$xshift = seq(-xo, xo, length.out = 5)[as.integer(as.factor(dfinds$server))]

  dfs = do.call(rbind, lapply(names(preds), function(pn) {
    df0 = dfinds[dfinds$server == "shared", ]
    i = df0$bl == pn
    fn = strsplit(df0[i, "bl"], "-")[[1]][1]
    flab = paste0(fn, " (", newdata[[fn]], ")")
    data.frame(pred = df0[i, "pred"], bli = df0[i, "bli"], xshift = xo,
      server = df0[i, "server"], feature = flab)
  }))
  gg_pind = ggplot(dfinds) +
    geom_segment(data = dfs, aes(x = bli - xshift, xend = bli + xshift,
      y = pred, yend = pred, color = server), size = 0.3) +
    geom_segment(data = dfinds,
      aes(x = bli + xshift, xend = bli + xshift,
        y = po, yend = pred, color = server), size = 0.3) +
    geom_point(data = dfinds[dfinds$server != "shared", ],
      aes(x = bli + xshift, y = pred, color = server, shape = server), show.legend = FALSE) +
    coord_flip() +
    mytheme() +
    ggsci::scale_color_jama() +
    xlab("") +
    ylab("Individual Contribution") +
    scale_x_continuous(breaks = dfs$bli, labels = dfs$feature) +
    labs(color = "")

  if (GGSAVE) {
    ggsave(gg_pind, filename = here::here("usecase/figures/gg-ind.pdf"), width = 13,
      height = 6, units = "cm")
  }
  sapply(names(cwb$coef$site), function(s) cwb$predict(newdata, s))


}


datashield.logout(connections)
