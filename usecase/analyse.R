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

cwb = dsCWB(connections, symbol, target, feature_names, mstop = 2000L,
  val_fraction = 0.2, patience = 3L, seed = 31415L, positive = "yes",
  df = 5, learning_rate = 0.01, derivs = 3L)

#save(cwb, file = "usecase/model-lr01.Rda")
#load("usecase/model-lr01.Rda")

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

readData = function(file, add_source = FALSE, add_id = FALSE, rm_pcols = TRUE) {
  if (grepl("reprocessed", file))
    tmp = read.csv(file, sep = " ", header = FALSE, na.strings = c("?", "-9", "-9.0"))
  else
    tmp = read.csv(file, header = FALSE, na.strings = c("?", "-9", "-9.0"))

  cnames = c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg", "thalach",
    "exang", "oldpeak", "slope", "ca", "thal", "num")

  fvals = c("sex", "cp", "fbs", "restecg", "exang", "slope", "ca", "thal", "num")
  colnames(tmp) = cnames

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
  }
  return(tmp)
}

datasets = list.files(here::here("usecase/data"), full.names = TRUE)

ll_dfs = lapply(datasets[-3], readData, rm_pcols = TRUE, add_source = TRUE)
df_full = do.call(rbind, lapply(ll_dfs, na.omit))

datashield.logout(connections)
