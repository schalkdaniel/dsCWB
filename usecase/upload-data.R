## ================================================== ##
##                  PROCESS RAW DATA                  ##
## ================================================== ##

## Download files from https://archive.ics.uci.edu/ml/datasets/heart+disease

sources = c("cleveland", "hungarian", "switzerland", "va")
dnames = paste0("data/processed.", sources, ".data")

#  1. #3 (age) years
#  2. #4 (sex) 1 male, 0 female
#  3. #9 (cp) chest pain type
#             -- Value 1: typical angina
#             -- Value 2: atypical angina
#             -- Value 3: non-anginal pain
#             -- Value 4: asymptomatic
#  4. #10 (trestbps) resting blood pressure (in mm Hg on admission to the hospital)
#  5. #12 (chol) serum cholestoral in mg/dl
#  6. #16 (fbs) (fasting blood sugar > 120 mg/dl) (1 = true; 0 = false)
#  7. #19 (restecg) resting electrocardiographic results
#                   -- Value 0: normal
#                   -- Value 1: having ST-T wave abnormality (T wave inversions and/or ST elevation or depression of > 0.05 mV)
#                   -- Value 2: showing probable or definite left ventricular hypertrophy by Estes' criteria
#  8. #32 (thalach) maximum heart rate achieved
#  9. #38 (exang) exercise induced angina (1 = yes; 0 = no)
# 10. #40 (oldpeak) oldpeak = ST depression induced by exercise relative to rest
# 11. #41 (slope)  the slope of the peak exercise ST segment
#                  -- Value 1: upsloping
#                  -- Value 2: flat
#                  -- Value 3: downsloping
# 12. #44 (ca) number of major vessels (0-3) colored by flourosopy
# 13. #51 (thal) 3 = normal; 6 = fixed defect; 7 = reversable defect
# 14. #58 (num) diagnosis of heart disease (angiographic disease status, target variable)
#               -- Value 0: < 50% diameter narrowing
#               -- Value 1: > 50% diameter narrowing

readData = function(file, add_source = FALSE) {
  tmp = read.csv(file, header = FALSE)

  cnames = c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg", "thalach",
    "exang", "oldpeak", "slope", "ca", "thal", "num")

  fvals = c("sex", "cp", "fbs", "restecg", "exang", "slope", "ca", "thal", "num")
  colnames(tmp) = cnames

  for (fv in fvals) {
    tmp[[fv]] = as.factor(tmp[[fv]])
  }
  if (add_source) a$source = strsplit(file, "[.]")[[1]][2]
  tmp$id = seq_len(nrow(tmp))
  return(tmp)
}

## ================================================== ##
##                UPLOAD TO DATASHIELD                ##
## ================================================== ##

surl     = "https://opal-demo.obiba.org/"
username = "administrator"
password = "password"

opal = opalr::opal.login(username = username, password = password, url = surl)

opalr::opal.project_delete(opal, project = "SLDS-TEST")
opalr::opal.project_create(opal, project = "SLDS-TEST", database = "mongodb")

datasets = list.files(here::here("usecase/data"), full.names = TRUE)

for (d in datasets) {

  opalr::opal.table_save(opal,
    tibble = tibble::as_tibble(readData(d)),
    project = "SLDS-TEST",
    table = strsplit(d, "[.]")[[1]][2],
    overwrite = TRUE,
    force = TRUE)
}

## ================================================== ##
##                   INSTALL PACKAGE                  ##
## ================================================== ##

ch1 = opalr::dsadmin.install_github_package(opal = opal, pkg = "dsCWB", username = "schalkdaniel", ref = "main")

if (ch1) {
  opalr::dsadmin.publish_package(opal = opal, pkg = "dsCWB")
} else {
  cat("Could not install package! :-(")
}

opalr::opal.logout(opal)

rm(list = ls())
