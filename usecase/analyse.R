################################################################################
##                                                                            ##
##                             DISTRIBUTED CWB                                ##
##                                                                            ##
################################################################################

library(ggplot2)
library(patchwork)

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

symbol = "D"
target = "heart_disease"
feature_names = setdiff(dsBaseClient::ds.names("D", connections)[[1]], target)

cwb = dsCWB(connections, symbol, target, feature_names, mstop = 20000L,
  val_fraction = 0.2, patience = 3L, seed = 31415L, positive = "yes",
  df = 5, learning_rate = 0.05, random_intercept_df = 2)


## FIGURES
## ========================================================================== ##

if (FALSE) {
  # Get log data with
  # - selected base learner
  # - train and test risk
  # - effect type (site or shared)
  l = cwb$getLog()

  table(l$bl)
  table(l$effect_type)
  table(paste0(l$bl, "-", l$effect_type))


  ## Base learner traces:
  ## =================================================================

  gg_traces_effects = plotBaselearnerTraces(cwb, n_legend = 7L, add_effect_type = TRUE) +
    ggsci::scale_color_jama() +
    ggsci::scale_fill_jama() +
    mytheme() +
    theme(legend.title = element_blank()) +
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
    ggsave(gg_risk, filename = here::here("usecase/gg-risk.pdf"), width = 13,
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
    ylab("VIP")

  gg_traces_vip = gridExtra::grid.arrange(gg_traces_effects, gg_vip,
    layout_matrix = cbind(1, 1, 2))

  gg_traces_vip

  if (GGSAVE) {
    ggsave(gg_traces_vip, filename = here::here("usecase/gg-traces-vip.pdf"),
      width = 25, height = 8, units = "cm")
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

#  gg_oldpeak = gridExtra::grid.arrange(gg_oldpeak_shared, gg_oldpeak_site,
#    layout_matrix = cbind(1, 2))
#
#  if (GGSAVE) {
#    ggsave(gg_oldpeak, filename = here::here("usecase/gg-oldpeak.pdf"),
#      width = 25, height = 8, units = "cm")
#  }

  gg_oldpeak_site2 = gg_oldpeak_site + gg_oldpeak_agg & theme(legend.position = "bottom")
  gg_oldpeak_site2 = gg_oldpeak_site2 + plot_layout(guides = "collect")

  gg_oldpeak2 = gg_oldpeak_shared + gg_oldpeak_site2 + plot_layout(widths = c(1, 2))
  gg_oldpeak2

  if (GGSAVE) {
    ggsave(gg_oldpeak2, filename = here::here("usecase/gg-oldpeak2.pdf"),
      width = 25, height = 8, units = "cm")
  }

  ## TODO: Categorical features
  #ggplot(sharedFEDataNum(cwb, "trestbps"), aes(x = trestbps, y = pred)) + geom_line()

  #ggplot(siteFEDataNum(cwb, "age"), aes(x = age, y = pred, color = server)) + geom_line()

  #ggplot(siteFEDataCat(cwb, "cp"), aes(x = "", y = pred, color = server, shape = server)) +
    #geom_point(size = 4, position = position_dodge2(width = 0.75)) +
    #geom_boxplot(show.legend = FALSE) +
    #xlab("") +
    #facet_wrap(cp ~ ., ncol = 4)

  #ggplot(sharedFEDataCat(cwb, "cp"), aes(x = cp, y = pred, color = cp, shape = cp)) +
    #geom_point(size = 4, position = position_dodge2(width = 0.75)) +
    #geom_boxplot(show.legend = FALSE) +
    #xlab("")

}

################################################################################
##                                                                            ##
##                          COMPARISON: COMPBOOST                             ##
##                                                                            ##
################################################################################

## PREPARE DATA
## ========================================================================== ##

## Load again, preparing the server is removing all custom functions.
source(here::here("usecase/helper.R"))

#  ## Get train indices
#  ## ===================================================
#  # To get comparable results, we train compboost
#  # on exact the same samples as the distributed
#  # version.
#
#  # cms = getDSLiteData(connections, "cm")
#  # val_idx = lapply(cms, function(cm) cm$getTrainValIndex()$test)
#
#  datasets = list.files(here::here("usecase/data"), full.names = TRUE, pattern = ".data")
#
#  ll_dfs = lapply(datasets, readData, add_source = TRUE)
#  df_full = do.call(rbind, lapply(ll_dfs, function(df) {
#    s = df$source[1]
#    idx = val_idx[[s]]
#    df = na.omit(df)
#    df$val = FALSE
#    df$val[idx] = TRUE
#    return(df)
#  }))
#  val_idx = which(df_full$val)
#  df_full$val = NULL

# save(df_full, file = here::here("usecase/data/df-full.Rda"))
# save(val_idx, file = here::here("usecase/data/val-idx.Rda"))

## TRAIN MODEL
## ========================================================================== ##

# remotes::install_github("schalkdaniel/compboost", ref = "dev")
library(compboost)

## Load again, preparing the server is removing all custom functions.
source(here::here("usecase/helper.R"))

load(file = here::here("usecase/data/df-full.Rda"))
load(file = here::here("usecase/data/val-idx.Rda"))

df = 5
df_random_intercept = 2
anisotrop = TRUE

site_var = "source"
target = "heart_disease"

fnums = c("age", "trestbps", "thalach", "oldpeak")
fcats = c("sex", "cp", "restecg", "exang")

cboost = Compboost$new(df_full, target = target, loss = LossBinomial$new(),
  optimize = OptimizerCoordinateDescent$new(4L), learning_rate = 0.05,
  positive = "yes", test_idx = val_idx, use_early_stopping = TRUE,
  stop_args = list(patience = 3L))#, oob_fraction = 0.2)

cboost$addBaselearner(site_var, "one-hot", BaselearnerCategoricalRidge, df = df_random_intercept)

for (f in fnums) {
  cboost$addBaselearner(f, "spline", BaselearnerPSpline, df = df)
  cboost$addTensor(f, site_var, df1 = df, df2 = df_random_intercept, anisotrop)
}

for (f in fcats) {
  uvals = unique(df_full[[f]])
  df0 = ifelse(length(uvals) < df, length(uvals), df)
  cboost$addBaselearner(f, "one-hot", BaselearnerCategoricalRidge, df = df0)
  cboost$addTensor(f, site_var, df1 = df0, df2 = df_random_intercept, anisotrop)
}

cboost$train(20000)

## COMPARISON OF LOGS
## ========================================================================== ##
if (FALSE) {
  l_dcwb = cwb$getLog()
  l_cwb = cboost$getLoggerData()

  l = data.frame(
    bl_dsCWB = paste0(l_dcwb$bl, "-", l_dcwb$effect_type),
    bl_compboost = l_cwb$baselearner,
    risk_dsCWB = l_dcwb$risk_train,
    risk_compboost = l_cwb$train_risk,
    val_dsCWB = l_dcwb$risk_val,
    val_compboost = l_cwb$oob_risk,
    iteration = seq_len(nrow(l_dcwb)) - 1)

  l
}

################################################################################
##                                                                            ##
##                            COMPARISON: MGCV                                ##
##                                                                            ##
################################################################################

library(mgcv)

df_mgcv = df_full
df_mgcv$src = as.factor(df_mgcv$source)
df_mgcv$source = NULL
df_mgcv$heart_disease = ifelse(df_mgcv$heart_disease == "yes", 1, 0)

fmgcv = heart_disease ~
  #sex + cp + restecg + exang +
  #s(age, bs = "ps", m = c(3, 2)) +
  #s(trestbps, bs = "ps", m = c(3, 2)) +
  #s(thalach, bs = "ps", m = c(3, 2)) +
  #s(oldpeak, bs = "ps", m = c(3, 2)) +

  ti(age, src, bs = c('ps', 're')) +
  ti(sex, src, bs = c('re', 're')) +
  ti(cp, src, bs = c('re', 're')) +
  ti(trestbps, src, bs = c('ps', 're')) +
  ti(restecg, src, bs = c('re', 're')) +
  ti(thalach, src, bs = c('ps', 're')) +
  ti(exang, src, bs = c('re', 're')) +
  ti(oldpeak, src, bs = c('ps', 're')) +

  #s(age, by = src, bs = "ps") +
  #s(sex, by = src, bs = "re") +
  #s(cp, by = src, bs = "re") +
  #s(trestbps, by = src, bs = "ps") +
  #s(restecg, by = src, bs = "re") +
  #s(thalach, by = src, bs = "re") +
  #s(exang, by = src, bs = "re") +
  #s(oldpeak, by = src, bs = "ps") +

  s(src, bs = "re")

mod_mgcv = gam(fmgcv, family = binomial(), data = df_mgcv)
#pe_mgcv = plot(mod_mgcv)


################################################################################
##                                                                            ##
##                          COMPARISON FIGURES                                ##
##                                                                            ##
################################################################################

source(here::here("usecase/helper.R"))

library(ggplot2)
library(patchwork)

# Without dsCWB:
fViz("oldpeak", mod_compboost = cboost, mod_mgcv = mod_mgcv, add_effects = TRUE, dat = df_full)
fViz("oldpeak", mod_compboost = cboost, mod_mgcv = mod_mgcv, site = FALSE, dat = df_full)
fViz("oldpeak", mod_compboost = cboost, mod_mgcv = mod_mgcv, site = TRUE, dat = df_full)

## Example for age:
## ================================================

ggleft = fViz("age", cwb, cboost, mod_mgcv, site = FALSE, dat = df_full) +
  ggtitle("age", "Shared feature effect")

ggright = fViz("age", cwb, cboost, mod_mgcv, site = TRUE, dat = df_full) +
  ggtitle("", "Site-effects")

ggleft + ggright

fViz("age", cwb, cboost, mod_mgcv, site = FALSE, dat = df_full, add_effects = TRUE) +
  ggtitle("age", "Shared feature effect")


## All numerical figures:
## ================================================

fnums = c("age", "trestbps", "thalach", "oldpeak")
fcats = c("sex", "cp", "restecg", "exang")

ggs = lapply(fnums, function(fname) {
  ggleft = fViz(fname, cwb, cboost, mod_mgcv, site = FALSE, dat = df_full) +
    ggtitle(fname, "Shared feature effect")
  ggright = fViz(fname, cwb, cboost, mod_mgcv, site = TRUE, dat = df_full) +
    ggtitle("", "Site-effects")

  gg = ggleft + ggright &
    theme(legend.position = "bottom") &
    ggsci::scale_color_jama() &
    mytheme()

  return(gg + plot_layout(guides = "collect"))
})
gg_nums = eval(parse(text = paste(paste0("ggs[[", seq_along(ggs), "]]"), collapse = " / ")))

if (GGSAVE) {
  ggsave(gg_nums, filename = here::here("usecase/fe-nums.pdf"), width = 20, height = 40, unit = "cm")
}

ggs = lapply(fnums, function(fname) {
  gg = fViz(fname, cwb, cboost, mod_mgcv, dat = df_full, add_effects = TRUE)
  if (fname == "age") {
    gg = gg +
      ggtitle("Aggregated site-specific feature effects", fname)
  } else {
    gg = gg +
      ggtitle("",fname)
  }
  gg = gg +
    theme(legend.position = "bottom") +
    ggsci::scale_color_jama() +
    mytheme() +
    ylab("Partial feature effect") +
    theme(legend.position = "bottom") +
    labs(color = "", linetype = "")

  return(gg)
})
gg_all = (ggs[[1]] + ggs[[2]]) / (ggs[[3]] + ggs[[4]]) & theme(legend.position = "bottom")
gg_all = gg_all + plot_layout(guides = "collect")

gg_all

if (GGSAVE) {
  ggsave(gg_all, filename = here::here("usecase/fe-nums-add-effects.pdf"), width = 22, height = 20, unit = "cm")

  ggsave(ggs[[4]] + ggtitle("Aggregated site-specific feature effects", "oldpeak"),
    filename = here::here("usecase/fe-oldpeak-mgcv.pdf"), width = 11, height = 10, unit = "cm")
}

datashield.logout(connections)
