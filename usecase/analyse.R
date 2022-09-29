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

  gg_traces_effects = plotBaselearnerTraces(cwb, n_legend = 7L, add_effect_type = TRUE) +
    ggsci::scale_color_jama() +
    ggsci::scale_fill_jama() +
    mytheme() +
    theme(legend.title = element_blank(), legend.position = "bottom") +
    ylab("Relative frequency\nof included base learner")

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
      height = 6, units = "cm")
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

  gg_traces_vip = (gg_traces_effects + gg_vip) / guide_area() +
    plot_layout(widths = c(2, 1), heights = c(4, 1), guides = "collect")
  gg_traces_vip

  if (GGSAVE) {
    ggsave(gg_traces_vip, filename = here::here("usecase/figures/gg-traces-vip.pdf"),
      width = 25, height = 10, units = "cm")
  }

  ## Partial feature effects:
  ## =================================================================

  fe_op_sh = sharedFEDataNum(cwb, "oldpeak")
  gg_oldpeak_shared = ggplot(fe_op_sh, aes(x = oldpeak, y = pred)) +
    geom_line() +
    mytheme() +
    ylab("Partial feature effect") +
    ggtitle("Shared effect")

  fe_op_si = siteFEDataNum(cwb, "oldpeak")
  gg_oldpeak_site = ggplot(fe_op_si, aes(x = oldpeak, y = pred, color = server)) +
    geom_line() +
    mytheme() +
    ggsci::scale_color_jama() +
    ylab("Partial feature effect") +
    ggtitle("Site effects") +
    theme(legend.title = element_blank())

  fe_op = fe_op_si
  fe_op$pred = fe_op$pred + rep(fe_op_sh$pred, 4)
  gg_oldpeak_agg = ggplot(fe_op, aes(x = oldpeak, y = pred, color = server)) +
    geom_line() +
    mytheme() +
    ggsci::scale_color_jama() +
    ylab("Partial feature effect") +
    ggtitle("Shared + site effects") +
    theme(legend.title = element_blank())

  gg_oldpeak_site2 = gg_oldpeak_site + gg_oldpeak_agg & theme(legend.position = "bottom")
  gg_oldpeak_site2 = gg_oldpeak_site2 + plot_layout(guides = "collect")

  gg_oldpeak2 = gg_oldpeak_shared + gg_oldpeak_site2 + plot_layout(widths = c(1, 2))
  gg_oldpeak2

  if (GGSAVE) {
    ggsave(gg_oldpeak2, filename = here::here("usecase/figures/gg-oldpeak2.pdf"),
      width = 25, height = 8, units = "cm")
  }
}


datashield.logout(connections)
