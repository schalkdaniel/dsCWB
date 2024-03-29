% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client-funs.R
\name{getClientSSE}
\alias{getClientSSE}
\title{Get the SSE from the client models.}
\usage{
getClientSSE(model_symbol = "cwb", par_list_binary = NULL)
}
\arguments{
\item{model_symbol}{(`character(1L)`)\cr
The name of the model at the server.}

\item{par_list_binary}{(`list()`)\cr
Encoded list of parameter values. Each entry must have a name
that corresponds to a base learner.}
}
\value{
SSE values of the of the registered base learners.
}
\description{
Get the SSE from the client models.
}
