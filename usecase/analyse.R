source(here::here("usecase/upload-data.R"))


library(DSI)
library(DSOpal)
library(dsBaseClient)

library(dsCWB)

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
conn = datashield.login(logins = builder$build(), assign = TRUE)
datashield.symbols(conn)

## Data dimensions per server:
(ddim = ds.dim("D"))





datashield.logout(conn)
