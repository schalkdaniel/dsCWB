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
    server   = paste0("ds", i),
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

#mstop = 100L
#learning_rate = 0.1
#val_fraction = 0.1
#patience = 5L
#seed = 31415L
#df = 5
#nknots = 20L
#ord = 3L
#derivs = 2L

#eps_for_break = NULL
#positive = "yes"


devtools::load_all()

cwb = dsCWB(connections, symbol, target, feature_names, mstop = 1000L,
  val_fraction = 0.05, patience = 3L, seed = 31415L, positive = "yes",
  df = 5, force_shared_iters = 100L)

l = cwb$getLog()
table(l$bl)
table(l$effect_type)

library(ggplot2)

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

sharedFEDataNum = function(cwb, feature) {
  blns = names(cwb$coef$shared)
  bln = blns[grepl(feature, blns)]
  coefs = cwb$coef$shared[[bln]]

  bl = cwb$bls[[bln]]

  lms = bl$getKnotRange()
  ndat = data.frame(x = seq(lms[1], lms[2], len = 100L))
  colnames(ndat) = feature

  df_fe = data.frame(val = ndat[[feature]], pred = bl$predictNewdata(ndat, par = coefs))
  colnames(df_fe)[1] = feature

  return(df_fe)

}
siteFEDataNum = function(cwb, feature) {
  ss = names(cwb$coef$site)
  blns = names(cwb$coef$site[[1]])
  bln = blns[grepl(feature, blns)]
  coefs = lapply(ss, function(s) cwb$coef$site[[s]][[bln]])
  names(coefs) = ss

  bl = cwb$bls[[bln]]

  lms = bl$getKnotRange()
  ndat = data.frame(x = seq(lms[1], lms[2], len = 100L))
  colnames(ndat) = feature

  df_fe = do.call(rbind, lapply(ss, function(s) {
    pred = bl$predictNewdata(ndat, par = coefs[[s]])
    out = data.frame(val = ndat[[feature]], pred = pred, server = s)
    return(out)
  }))
  colnames(df_fe)[1] = feature

  return(df_fe)
}
siteFEDataCat = function(cwb, feature) {
  ss = names(cwb$coef$site)
  blns = names(cwb$coef$site[[1]])
  bln = blns[grepl(feature, blns)]
  bl = cwb$bls[[bln]]
  coefs = lapply(ss, function(s) {
    cout = cwb$coef$site[[s]][[bln]]
    names(cout) = names(bl$getDictionary())
    cout
  })
  names(coefs) = ss

  df_fe = do.call(rbind, lapply(ss, function(s) {
    pred = coefs[[s]]
    out = data.frame(val = names(pred), pred = pred, server = s)
    return(out)
  }))
  colnames(df_fe)[1] = feature

  return(df_fe)
}

ggplot(sharedFEDataNum(cwb, "oldpeak"), aes(x = oldpeak, y = pred)) + geom_line()

ggplot(siteFEDataNum(cwb, "oldpeak"), aes(x = oldpeak, y = pred, color = server)) + geom_line()
ggplot(siteFEDataNum(cwb, "age"), aes(x = age, y = pred, color = server)) + geom_line()
ggplot(siteFEDataNum(cwb, "chol"), aes(x = chol, y = pred, color = server)) + geom_line()
ggplot(siteFEDataNum(cwb, "trestbps"), aes(x = trestbps, y = pred, color = server)) + geom_line()

ggplot(siteFEDataCat(cwb, "cp"), aes(x = "", y = pred, color = server, shape = server)) +
  geom_point(size = 4, position = position_dodge2(width = 0.75)) +
  geom_boxplot(show.legend = FALSE) +
  xlab("") +
  facet_wrap(cp ~ ., ncol = 4)

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

library(compboost)

datashield.logout(connections)
