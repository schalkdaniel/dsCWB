% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client-funs.R
\name{updateClientPenalty}
\alias{updateClientPenalty}
\title{Update the penalties of the client models.}
\usage{
updateClientPenalty(
  model_symbol,
  ll_pen_binary,
  anisotrop = TRUE,
  simple = FALSE,
  update_random_intercept = FALSE
)
}
\arguments{
\item{model_symbol}{(`character(1L)`)\cr
The name of the model at the server.}

\item{ll_pen_binary}{(`character(1L)`)\cr
Named list with new penalty terms encoded by `encodeObject`.}

\item{anisotrop}{(`logical(1L)`)\cr
Flag indicating whether the penalty should be done anisotrop or isotrop.}

\item{simple}{(`logical(1L)`)\cr
Flag indicating that just the penalty is updated but no tensor operation is applied.}

\item{update_random_intercept}{(`logical(1L)`)\cr
Flag indicating whether the random intercept penalty should also get updated or not.}
}
\description{
Update the penalties of the client models.
}
