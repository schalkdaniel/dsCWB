################################################################################
##                                                                            ##
##                             DISTRIBUTED CWB                                ##
##                                                                            ##
################################################################################

USE_REAL_TEST_ENV = FALSE
GGSAVE = TRUE
MAX_MSTOP = 100000L
RERUN_DIST_CWB = FALSE

## SETUP DataSHIELD ENVIRONMENT
## ========================================================================= ##
if (RERUN_DIST_CWB) {

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

mod_dsCWB = dsCWB(connections, symbol, target, feature_names, mstop = MAX_MSTOP,
  positive = "yes", df = 2.2, learning_rate = 0.1, random_intercept_df = 3,
  nknots = 10L)

save(mod_dsCWB, file = here::here("usecase/data/mod-dsCWB.Rda"))

} else {

devtools::load_all()
load(file = here::here("usecase/data/mod-dsCWB.Rda"))

}

################################################################################
##                                                                            ##
##                          COMPARISON: COMPBOOST                             ##
##                                                                            ##
################################################################################

## PREPARE DATA
## ========================================================================== ##

## Load again, preparing the server is removing all custom functions.
# source(here::here("usecase/helper.R"))

#  ## Get train indices
#  ## ===================================================
#  # To get comparable results, we train compboost
#  # on exact the same samples as the distributed
#  # version.
#
#  cms = getDSLiteData(connections, "cm")
#  val_idx = lapply(cms, function(cm) cm$getTrainValIndex()$test)
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
#library(compboost)
devtools::load_all("~/repos/compboost")

## Load again, preparing the server is removing all custom functions.
source(here::here("usecase/helper.R"))

load(file = here::here("usecase/data/df-full.Rda"))

df = 2.2
df_random_intercept = 3
isotrop = FALSE
n_knots = 10L

site_var = "source"
target = "heart_disease"

fnums = c("age", "trestbps", "thalach", "oldpeak")
fcats = c("sex", "cp", "restecg", "exang")

mod_compboost = Compboost$new(df_full, target = target, loss = LossBinomial$new(),
  optimize = OptimizerCoordinateDescent$new(4L), learning_rate = 0.1,
  positive = "yes")

mod_compboost$addBaselearner(site_var, "one-hot", BaselearnerCategoricalRidge, df = df_random_intercept)

for (f in fnums) {
  mod_compboost$addBaselearner(f, "spline", BaselearnerPSpline, df = df, n_knots = n_knots)
  mod_compboost$addTensor(f, site_var, df1 = df, df2 = df_random_intercept, isotrop, n_knots = n_knots)
}

for (f in fcats) {
  uvals = unique(df_full[[f]])
  df0 = ifelse(length(uvals) < df, length(uvals), df)
  mod_compboost$addBaselearner(f, "one-hot", BaselearnerCategoricalRidge, df = df0)
  mod_compboost$addTensor(f, site_var, df1 = df0, df2 = df_random_intercept, isotrop)
}

mod_compboost$train(MAX_MSTOP)


## COMPARISON OF LOGS
## ========================================================================== ##

if (FALSE) {
  l_dcwb = mod_dsCWB$getLog()
  l_cwb = mod_compboost$getLoggerData()

  nmax = nrow(l_cwb)
  nmin = nrow(l_dcwb)

  l = data.frame(
    bl_dsCWB = paste0(l_dcwb$bl, "-", l_dcwb$effect_type),
    bl_compboost = l_cwb$baselearner,
    risk_dsCWB = l_dcwb$risk_train,
    risk_compboost = l_cwb$train_risk,
    iteration = seq_len(nrow(l_cwb)) - 1)

  table(l$bl_dsCWB, l$bl_compboost)
  #write.csv(l, file = here::here("usecase/data/log.csv"))
  #l
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
  sex + cp + restecg + exang +

  s(sex, src, bs = c('re')) +
  s(cp, src, bs = c('re')) +
  s(restecg, src, bs = c('re')) +
  s(exang, src, bs = c('re')) +

  s(age, bs = "ps", m = c(3, 2)) +
  s(trestbps, bs = "ps", m = c(3, 2)) +
  s(thalach, bs = "ps", m = c(3, 2)) +
  s(oldpeak, bs = "ps", m = c(3, 2)) +

  ti(age, src, bs = c('ps', 're')) +
  ti(trestbps, src, bs = c('ps', 're')) +
  ti(thalach, src, bs = c('ps', 're')) +
  ti(oldpeak, src, bs = c('ps', 're')) +

  s(src, bs = "re")

mod_mgcv = gam(fmgcv, family = binomial(), data = df_mgcv)

# Calculate risk
d0 = mod_compboost$data
d0$src = d0$source
p = predict(mod_mgcv, newdata = d0)
y = as.vector(mod_compboost$response$getResponse())

mean(log(1 + exp(-y * p)))
#> [1] 0.4441
min(mod_compboost$model$getRiskVector())
#> [1] 0.4245
min(mod_dsCWB$getLog()$risk_train)
#> [1] 0.4245

################################################################################
##                                                                            ##
##                          COMPARISON FIGURES                                ##
##                                                                            ##
################################################################################

source(here::here("usecase/helper.R"))

library(ggplot2)
library(patchwork)

## All numerical figures:
## ================================================

fnums = c("age", "trestbps", "thalach", "oldpeak")
fcats = c("sex", "cp", "restecg", "exang")

ggs_num = lapply(fnums, function(fname) {
  gg = fViz(fname, mod_dsCWB = mod_dsCWB, mod_compboost = mod_compboost,
    mod_mgcv = mod_mgcv, add_effects = TRUE)

  #gg = fViz(fname, mod_compboost = mod_compboost, mod_mgcv = mod_mgcv, add_effects = TRUE, dummy_dist_cwb = TRUE)

  if (fname == "age") {
    gg = gg +
      ggtitle("Aggregated site-specific feature effects", fname)
  } else {
    gg = gg +
      ggtitle("", fname)
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
gg_all_num = (ggs_num[[1]] + ggs_num[[2]]) / (ggs_num[[3]] + ggs_num[[4]]) & theme(legend.position = "bottom")
gg_all_num = gg_all_num + plot_layout(guides = "collect")

gg_all_num

if (GGSAVE) {
  ggsave(gg_all_num, filename = here::here("usecase/figures/fe-nums-add-effects.pdf"), width = 22, height = 20, unit = "cm")

  ggsave(ggs_num[[4]] + ggtitle("Aggregated site-specific feature effects", "oldpeak"),
    filename = here::here("usecase/figures/fe-oldpeak-mgcv.pdf"), width = 11, height = 10, unit = "cm")
}


ggs_cat = lapply(fcats, function(fname) {
  gg = fViz(fname, mod_dsCWB = mod_dsCWB, mod_compboost = mod_compboost,
    mod_mgcv = mod_mgcv, add_effects = TRUE)

  if (fname == "age") {
    gg = gg +
      ggtitle("Aggregated site-specific feature effects", fname)
  } else {
    gg = gg +
      ggtitle("", fname)
  }
  gg = gg +
    theme(legend.position = "bottom") +
    ggsci::scale_color_jama() +
    mytheme() +
    ylab("Partial feature effect") +
    theme(legend.position = "bottom") +
    labs(color = "", shape = "")

  return(gg)
})
gg_all_cat = (ggs_cat[[1]] + ggs_cat[[2]]) / (ggs_cat[[3]] + ggs_cat[[4]]) & theme(legend.position = "bottom")
gg_all_cat = gg_all_cat + plot_layout(guides = "collect")

gg_all_cat

if (GGSAVE) {
  ggsave(gg_all_cat, filename = here::here("usecase/figures/fe-cats-add-effects.pdf"), width = 22, height = 20, unit = "cm")

  ggsave(ggs_cat[[2]] + ggtitle("Aggregated site-specific feature effects", "cp"),
    filename = here::here("usecase/figures/fe-cp-mgcv.pdf"), width = 11, height = 10, unit = "cm")
}


datashield.logout(connections)
