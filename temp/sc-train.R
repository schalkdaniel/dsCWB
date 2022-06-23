q()
R
library(testthat)
devtools::load_all()

# Install package on DataSHIELD test server:
surl     = "https://opal-demo.obiba.org/"
username = "administrator"
password = "password"

opal = opalr::opal.login(username = username, password = password, url = surl)
opalr::dsadmin.install_github_package(opal = opal, pkg = "dsCWB", username = "schalkdaniel", ref = "main")
opalr::dsadmin.publish_package(opal = opal, pkg = "dsCWB")
opalr::opal.logout(opal, save = FALSE)

# Establish connection to the DataSHIELD servers:
library(DSI)
library(DSOpal)

builder = newDSLoginBuilder()

builder$append(
  server   = "ds-test-server-dummy1",
  url      = surl,
  user     = username,
  password = password
)
builder$append(
  server   = "ds-test-server-dummy2",
  url      = surl,
  user     = username,
  password = password
)
connections = datashield.login(logins = builder$build(), assign = TRUE)
datashield.assign(connections, "dat", quote(iris))

symbol = "dat"
target = "Sepal.Length"
feature_names = setdiff(names(iris), target)
mstop = 100L
learning_rate = 0.1
df = 5
nknots = 20L
ord = 3L
derivs = 2L
val_fraction = NULL
patience = NULL
eps_for_break = NULL
positive = NULL
seed = NULL

datashield.logout(connections)

