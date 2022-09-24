tcheck = function(cwb) {
  checkmate::assertR6(cwb, classes = "HostModel")
  if (is.null(cwb$getLog()))
    stop("Model was not trained yet or at least no log is available!")
}

#' @title Get data for a numerical shared effect
#' @param cwb ()\cr
#' @param feature (`character(1L)`)\cr
#'   Name of the feature.
#' @param npoints (`integer(1L)`)\cr
#'   Number of feature points (equidistand grid) used for predicting.
#' @return `data.frame()` with feature values and the prediction.
#' @export
sharedFEDataNum = function(cwb, feature, npoints = 100L) {
  tcheck(cwb)
  checkmate::assertChoice(feature, choices = cwb$getFeatureNames())
  checkmate::assertCount(npoints)

  blns = names(cwb$coef$shared)
  bln = blns[grepl(feature, blns)]
  coefs = cwb$coef$shared[[bln]]

  bl = cwb$bls[[bln]]

  lms = bl$getKnotRange()
  ndat = data.frame(x = seq(lms[1], lms[2], len = npoints))
  colnames(ndat) = feature

  df_fe = data.frame(val = ndat[[feature]], pred = bl$predictNewdata(ndat, par = coefs))
  colnames(df_fe)[1] = feature

  return(df_fe)

}

#' @title Get data for a categorical shared effect
#' @param cwb ()\cr
#' @param feature (`character(1L)`)\cr
#'   Name of the feature.
#' @return `data.frame()` with feature values and the prediction.
#' @export
sharedFEDataCat = function(cwb, feature) {
  tcheck(cwb)
  checkmate::assertChoice(feature, choices = cwb$getFeatureNames())

  blns = names(cwb$coef$shared)
  bln = blns[grepl(feature, blns)]
  bl = cwb$bls[[bln]]
  coefs = cwb$coef$shared[[bln]]

  df_fe = data.frame(
    val = names(bl$getDictionary()),
    pred = coefs)

  colnames(df_fe)[1] = feature

  return(df_fe)
}

#' @title Get data for a numerical site effect
#' @param cwb ()\cr
#' @param feature (`character(1L)`)\cr
#'   Name of the feature.
#' @param npoints (`integer(1L)`)\cr
#'   Number of feature points (equidistand grid) used for predicting.
#' @return `data.frame()` with feature values, the prediction, and the server.
#' @export
siteFEDataNum = function(cwb, feature, npoints = 100L) {
  tcheck(cwb)
  checkmate::assertChoice(feature, choices = cwb$getFeatureNames())
  checkmate::assertCount(npoints)

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
    if (is.null(coefs[[s]])) {
      pred = 0
    } else {
      pred = bl$predictNewdata(ndat, par = coefs[[s]])
    }
    out = data.frame(val = ndat[[feature]], pred = pred, server = s)
    return(out)
  }))
  colnames(df_fe)[1] = feature

  return(df_fe)
}

#' @title Get data for a categorical site effect
#' @param cwb ()\cr
#' @param feature (`character(1L)`)\cr
#'   Name of the feature.
#' @return `data.frame()` with feature values, the prediction, and the server.
#' @export
siteFEDataCat = function(cwb, feature) {
  tcheck(cwb)
  checkmate::assertChoice(feature, choices = cwb$getFeatureNames())

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
