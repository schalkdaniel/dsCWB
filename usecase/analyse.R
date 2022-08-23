source(here::here("usecase/upload-data.R"))

library(DSI)
library(DSOpal)
library(dsBaseClient)

#library(dsCWB)
devtools::load_all()

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
datashield.symbols(connections)

## Data dimensions per server:
(ddim = ds.dim("D"))


symbol = "D"
target = "heart_disease"
feature_names = setdiff(dsBaseClient::ds.names("D", connections)[[1]], target)

cwb = dsCWB(connections, symbol, target, feature_names, mstop = 100L,
  val_fraction = 0.2, patience = 3L, seed = 31415L, positive = "yes",
  df = 5, learning_rate = 0.01, derivs = 3L)

# Get hyperparameters:
# cwb_site_hps = datashield.aggregate(connections, "getBlHyperpars('cm')")

datashield.logout(connections)

#save(cwb, cwb_site_hps, file = here::here("usecase/model-lr01.Rda"))
#load(here::here("usecase/model-lr01.Rda"))

l = cwb$getLog()

library(ggplot2)

plotBaselearnerTraces(cwb, n_legend = 7L, add_effect_type = FALSE)
plotBaselearnerTraces(cwb, n_legend = 7L, add_effect_type = TRUE)

table(l$bl)
table(l$effect_type)
table(paste0(l$bl, "-", l$effect_type))

ggplot(l, aes(x = iteration)) +
  geom_line(aes(y = risk_train, color = "Train risk")) +
  geom_line(aes(y = risk_val, color = "Validation risk"))

vip = function(risk, features) {
  agg = aggregate(risk, by = list(features), FUN = sum)
  ord = order(agg$x, decreasing = TRUE)
  agg = agg[ord, ]
  return(data.frame(feature = agg[[1]], vip = agg[[2]]))
}
vi = vip(l$risk_train, paste0(l$bl, "-", l$effect_type))
ggplot(vi, aes(x = reorder(feature, vip), y = vip)) +
  geom_bar(stat = "identity") +
  coord_flip()

ggplot(sharedFEDataNum(cwb, "oldpeak"), aes(x = oldpeak, y = pred)) + geom_line()
ggplot(sharedFEDataNum(cwb, "trestbps"), aes(x = trestbps, y = pred)) + geom_line()

ggplot(siteFEDataNum(cwb, "oldpeak"), aes(x = oldpeak, y = pred, color = server)) + geom_line()
ggplot(siteFEDataNum(cwb, "age"), aes(x = age, y = pred, color = server)) + geom_line()

ggplot(siteFEDataCat(cwb, "cp"), aes(x = "", y = pred, color = server, shape = server)) +
  geom_point(size = 4, position = position_dodge2(width = 0.75)) +
  geom_boxplot(show.legend = FALSE) +
  xlab("") +
  facet_wrap(cp ~ ., ncol = 4)

ggplot(sharedFEDataCat(cwb, "cp"), aes(x = cp, y = pred, color = cp, shape = cp)) +
  geom_point(size = 4, position = position_dodge2(width = 0.75)) +
  geom_boxplot(show.legend = FALSE) +
  xlab("")

# ---------------------------------------------------------------------------------------- #
#                                Comparison to compboost                                   #
# ---------------------------------------------------------------------------------------- #


source(here::here("usecase/helper.R"))

datasets = list.files(here::here("usecase/data"), full.names = TRUE)

ll_dfs = lapply(datasets, readData, add_source = TRUE)
df_full = do.call(rbind, lapply(ll_dfs, na.omit))

#remotes::install_github("schalkdaniel/compboost", ref = "dev")
library(compboost)

df = 3
anistrop = TRUE

site_var = "source"
target = "heart_disease"

fnum = c("age", "trestbps", "thalach", "oldpeak")
fcat = c("sex", "cp", "restecg", "exang")

set.seed(10L)
cboost = Compboost$new(df_full, target = target, loss = LossBinomial$new(),
  learning_rate = 0.01, oob_fraction = 0.2, use_early_stopping = TRUE,
  stop_args = list(eps_for_break = 0, patience = 3))

for (f in fnum) {
  cboost$addBaselearner(f, "spline", BaselearnerPSpline, df = df)
  cboost$addTensor(f, site_var, df1 = df, df2 = 1, anistrop)
}

for (f in fcat) {
  uvals = unique(df_full[[f]])
  df0 = ifelse(length(uvals) < df, length(uvals), df)
  cboost$addBaselearner(f, "one-hot", BaselearnerCategoricalRidge, df = df0)
  cboost$addTensor(f, site_var, df1 = df0, df2 = 1, anistrop)
}

# Get list of base learners:
cboost$getBaselearnerNames()

cboost$train(10000L)

table(cboost$getSelectedBaselearner())

plotRisk(cboost)
plotBaselearnerTraces(cboost)

plotPEUni(cboost, "oldpeak")
plotTensor(cboost, "oldpeak_source_tensor")
plotTensor(cboost, "cp_source_tensor")


## Manually calculate penalties =============================================== ##

(ptensor = cboost$baselearner_list$oldpeak_source_tensor$factory$getPenalty())
pspline = cboost$baselearner_list$oldpeak_spline$factory$getPenalty()


dmat = cboost$baselearner_list$oldpeak_source_tensor$factory$getData()
xtx = dmat %*% t(dmat)
pmat = cboost$baselearner_list$oldpeak_source_tensor$factory$getPenaltyMat()

cpsp::demmlerReinsch(xtx, pmat, df)
ptensor


xtx_source = diag(table(df_full$source))
pmat_source = diag(ncol(xtx_source))
psource = cpsp::demmlerReinsch(xtx_source, pmat_source, 1)

penaltyKron = function(A, B) {
  Pa = diag(ncol(A))
  Pb = diag(ncol(B))

  return(kronecker(A, Pb) + kronecker(Pa, B))
}

a = rnorm(1)
b = rnorm(1)

pp = penaltyKron(a * pmat, b * pmat_source)

psource * pspline
ptensor
pmat[1:10, 1:10]

## Manually calculate penalties =============================================== ##





# Compare to the distributed algorithm:
devtools::load_all()
load(here::here("usecase/model-lr01.Rda"))

l = cwb$getLog()
max(l$iteration)

cwb_site_hps$cleveland[["oldpeak-spline"]]
bl = cwb$bls[["oldpeak-spline"]]
bl$getHyperpars()$penalty

