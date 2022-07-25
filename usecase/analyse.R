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

mstop = 100L
learning_rate = 0.1
val_fraction = 0.1
patience = 5L
seed = 31415L
df = 5
nknots = 20L
ord = 3L
derivs = 2L

eps_for_break = NULL
positive = "yes"


cwb = dsCWB(connections, symbol, target, feature_names, mstop = 100L,
  val_fraction = 0.1, patience = 5L, seed = 31415L, positive = "yes")

l = cwb$getLog()


datashield.logout(connections)
