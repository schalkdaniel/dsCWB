% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/client-model.R
\name{ClientModel}
\alias{ClientModel}
\title{Client model}
\description{
This class creates a CWB model container for the sites.
It collects all base learner and also the data (data frame,
pseudo residuals, loss function, etc.) per site.
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{loss}}{(`Loss`)\cr
Loss object used for modelling.}

\item{\code{bls}}{(`list(Baselearner)`)\cr
List of base learner.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{ClientModel$new()}}
\item \href{#method-getLossOptimalConstant}{\code{ClientModel$getLossOptimalConstant()}}
\item \href{#method-initConstantModel}{\code{ClientModel$initConstantModel()}}
\item \href{#method-getSSE}{\code{ClientModel$getSSE()}}
\item \href{#method-getRisk}{\code{ClientModel$getRisk()}}
\item \href{#method-getXtX}{\code{ClientModel$getXtX()}}
\item \href{#method-getXty}{\code{ClientModel$getXty()}}
\item \href{#method-getPseudoResids}{\code{ClientModel$getPseudoResids()}}
\item \href{#method-getBlNames}{\code{ClientModel$getBlNames()}}
\item \href{#method-getBlParams}{\code{ClientModel$getBlParams()}}
\item \href{#method-updateCWBParts}{\code{ClientModel$updateCWBParts()}}
\item \href{#method-update}{\code{ClientModel$update()}}
\item \href{#method-updatePenalty}{\code{ClientModel$updatePenalty()}}
\item \href{#method-addBaselearners}{\code{ClientModel$addBaselearners()}}
\item \href{#method-getFeatureNames}{\code{ClientModel$getFeatureNames()}}
\item \href{#method-getPrediction}{\code{ClientModel$getPrediction()}}
\item \href{#method-getOffset}{\code{ClientModel$getOffset()}}
\item \href{#method-getTarget}{\code{ClientModel$getTarget()}}
\item \href{#method-getTrainValObs}{\code{ClientModel$getTrainValObs()}}
\item \href{#method-getTrainValIndex}{\code{ClientModel$getTrainValIndex()}}
\item \href{#method-clone}{\code{ClientModel$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Creates a new instance of this [R6][R6::R6Class] class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$new(
  symbol = "D",
  target,
  feature_names = NULL,
  learning_rate = 0.1,
  df = 5L,
  nknots = 20L,
  ord = 3L,
  derivs = 2L,
  val_fraction = NULL,
  positive = NULL,
  seed = NULL,
  envir = parent.frame(),
  random_intercept = TRUE,
  random_intercept_df = 2
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{symbol}}{(`character(1L)`)\cr
Character containing the name of the data.}

\item{\code{target}}{(`character(1L)`)\cr
Character containing the name of the target variable.}

\item{\code{feature_names}}{(`character()`)\cr
Character vector of all target variables.}

\item{\code{learning_rate}}{(`numeric(1L)`)\cr
Learning rate.}

\item{\code{df}}{(`numeric(1L)`)\cr
Degrees of freedom.}

\item{\code{nknots}}{(`integer(1L)`)\cr
Number of inner knots.}

\item{\code{ord}}{(`integer(1L)`)\cr
Polynomial degree of basis functions.}

\item{\code{derivs}}{(`integer(1L)`)\cr
Number of penalized differences.}

\item{\code{val_fraction}}{(`numeric(1L)`)\cr
Fraction of observations used for validation.}

\item{\code{positive}}{(`integer(1L)`)\cr
Character indicating the positive class in the binary classification setting.}

\item{\code{seed}}{(`integer(1L)`)\cr
Seed used for drawing the validation data.}

\item{\code{envir}}{(`environment()`)\cr
Environment of the data.}

\item{\code{random_intercept}}{(`logical(1L)`)\cr
Indicator whether random intercepts should be added or not.}

\item{\code{random_intercept_df}}{(`integer(1L)`)\cr
Penalty of the random intercept.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getLossOptimalConstant"></a>}}
\if{latex}{\out{\hypertarget{method-getLossOptimalConstant}{}}}
\subsection{Method \code{getLossOptimalConstant()}}{
Get loss optimal constant.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getLossOptimalConstant()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-initConstantModel"></a>}}
\if{latex}{\out{\hypertarget{method-initConstantModel}{}}}
\subsection{Method \code{initConstantModel()}}{
Initialize the loss optimal constant.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$initConstantModel(const)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{const}}{(`numeric(1L)`)\cr
Loss optimal constant.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getSSE"></a>}}
\if{latex}{\out{\hypertarget{method-getSSE}{}}}
\subsection{Method \code{getSSE()}}{
Get the SSE for each base learner.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getSSE(par_list = NULL, par_list_binary = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{par_list}}{(`list()`)\cr
List of parameter values. Each entry must have a name
that corresponds to a base learner.}

\item{\code{par_list_binary}}{(`character(1L)`)\cr
Parameter list as encoded binary object.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getRisk"></a>}}
\if{latex}{\out{\hypertarget{method-getRisk}{}}}
\subsection{Method \code{getRisk()}}{
Get the risk.


Get target from the data from the symbol.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getRisk(type = "train")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{(`character(1L)`)\cr
Type which target shold be returned. Choices are
`full` (for the whole vector), `train` (for the vector used for training),
and `val` (for the validation vector).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getXtX"></a>}}
\if{latex}{\out{\hypertarget{method-getXtX}{}}}
\subsection{Method \code{getXtX()}}{
Return XtX per base learner.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getXtX()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getXty"></a>}}
\if{latex}{\out{\hypertarget{method-getXty}{}}}
\subsection{Method \code{getXty()}}{
Get Xty per base learner.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getXty()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getPseudoResids"></a>}}
\if{latex}{\out{\hypertarget{method-getPseudoResids}{}}}
\subsection{Method \code{getPseudoResids()}}{
Get pseudo residuals.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getPseudoResids(type = "full")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{(`character(1L)`)\cr
Type which target shold be returned. Choices are
`full` (for the whole vector), `train` (for the vector used for training),
and `val` (for the validation vector).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getBlNames"></a>}}
\if{latex}{\out{\hypertarget{method-getBlNames}{}}}
\subsection{Method \code{getBlNames()}}{
Get base learner names.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getBlNames()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getBlParams"></a>}}
\if{latex}{\out{\hypertarget{method-getBlParams}{}}}
\subsection{Method \code{getBlParams()}}{
Get base learner parameters.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getBlParams()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-updateCWBParts"></a>}}
\if{latex}{\out{\hypertarget{method-updateCWBParts}{}}}
\subsection{Method \code{updateCWBParts()}}{
Update the model by adding a base learner.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$updateCWBParts(blname, par = NULL, par_binary = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{blname}}{(`character(1L)`)\cr
Character containing the base learner name.}

\item{\code{par}}{(`numeric()`)\cr
Parameter used to update pseudo residuals and predictions.}

\item{\code{par_binary}}{(`character(1L)`)\cr
Binary encoded parameter used to update pseudo residuals and predictions.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-update"></a>}}
\if{latex}{\out{\hypertarget{method-update}{}}}
\subsection{Method \code{update()}}{
Update the model by adding a base learner.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$update(blname)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{blname}}{(`character(1L)`)\cr
Character containing the base learner name.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-updatePenalty"></a>}}
\if{latex}{\out{\hypertarget{method-updatePenalty}{}}}
\subsection{Method \code{updatePenalty()}}{
Re-set the penalty terms.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$updatePenalty(
  ll_pen,
  anisotrop = TRUE,
  simple = FALSE,
  update_random_intercept = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ll_pen}}{(`list()`)\cr
Named list with new penalty terms.}

\item{\code{anisotrop}}{(`logical(1L)`)\cr
Flag indicating whether the penalty should be done anisotrop or isotrop.}

\item{\code{simple}}{(`logical(1L)`)\cr
Flag indicating that just the penalty is updated but no tensor operation is applied.}

\item{\code{update_random_intercept}}{(`logical(1L)`)\cr
Flag indicating whether the random intercept penalty should also get updated or not.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addBaselearners"></a>}}
\if{latex}{\out{\hypertarget{method-addBaselearners}{}}}
\subsection{Method \code{addBaselearners()}}{
Add base learners to the model.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$addBaselearners(ll_init)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ll_init}}{(`list()`)\cr
List containing the initialization parameters.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getFeatureNames"></a>}}
\if{latex}{\out{\hypertarget{method-getFeatureNames}{}}}
\subsection{Method \code{getFeatureNames()}}{
Get feature names
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getFeatureNames()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getPrediction"></a>}}
\if{latex}{\out{\hypertarget{method-getPrediction}{}}}
\subsection{Method \code{getPrediction()}}{
Get prediction.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getPrediction(type = "full")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{(`character(1L)`)\cr
Type which target shold be returned. Choices are
`full` (for the whole vector), `train` (for the vector used for training),
and `val` (for the validation vector).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getOffset"></a>}}
\if{latex}{\out{\hypertarget{method-getOffset}{}}}
\subsection{Method \code{getOffset()}}{
Return offset of the model.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getOffset()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getTarget"></a>}}
\if{latex}{\out{\hypertarget{method-getTarget}{}}}
\subsection{Method \code{getTarget()}}{
Get target from the data from the symbol.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getTarget(type = "full")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{(`character(1L)`)\cr
Type which target shold be returned. Choices are
`full` (for the whole vector), `train` (for the vector used for training),
and `val` (for the validation vector).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getTrainValObs"></a>}}
\if{latex}{\out{\hypertarget{method-getTrainValObs}{}}}
\subsection{Method \code{getTrainValObs()}}{
Get the number of rows per train and validation set.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getTrainValObs()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getTrainValIndex"></a>}}
\if{latex}{\out{\hypertarget{method-getTrainValIndex}{}}}
\subsection{Method \code{getTrainValIndex()}}{
Get List of train and test indices.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$getTrainValIndex()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ClientModel$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
