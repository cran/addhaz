\name{MultAddHaz}
\alias{MultAddHaz}

\title{Fit Multinomial Additive Hazards Models}

\description{
This function fits multinomial additive hazards models subject to linear inequality constraints using the function \code{constrOptim} in the \code{stats} package for multinomial (multi-category) outcomes. It also calculates the cause-specific contributions to the disability prevalence for each category of the response variable based on the extension of the attribution method, as proposed by Yokota, Van Oyen, and Molenberghs (2015).
}

\usage{
MultAddHaz(formula, data, subset, weights, na.action, model = TRUE,
           contrasts = NULL, start.val, attrib = FALSE, attrib.var,
           type.attrib = "abs", set.seed = FALSE, seed, bootstrap = FALSE,
           nbootstrap, parallel = FALSE, type.parallel = "snow", ncpus = 4,...)
}

\arguments{
  \item{formula}{A formula expression of the form response ~ predictors, similar to other regression models. In case of \code{attrib = TRUE}, the first predictor in the formula should be the attrib.var. See example.}
  \item{data}{an optional data frame or matrix containing the variables in the model. If not found in \code{data}, the variables are taken from environment(formula), typically the environment from which \code{MultAddHaz} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process. All observations are included by default.}
  \item{weights}{an optional vector of survey weights to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. The default is set by the \code{na.action} setting of options, and is \code{na.fail} if that is unset. The 'factory-fresh' default is \code{na.omit}.}
  \item{model}{logical. If \code{TRUE}, the model frame is included as a component of the returned object.}
  \item{contrasts}{an optional list to be used for some or all of the factors appearing as variables in the model formula.}
  \item{start.val}{vector of the starting values. Starting values are required and they should lie inside the parameter space.}
  \item{attrib}{logical. Should the attribution of chronic diseases/conditions for each disability level be calculated? Default is \code{FALSE}.}
  \item{attrib.var}{character indicating the name of the attribution variable to be used if \code{attrib = TRUE}. If missing, the attribution variable is the intercept (no stratification of the attribution results). The attribution variable must be the first variable included in the linear predictor in \code{formula}. See example.}
  \item{type.attrib}{type of attribution to be calculated. The options are \code{"abs"} for  absolute contribution, \code{"rel"} for relative contribution, or \code{"both"} for both absolute and relative contributions. Default is \code{"abs"}.}
  \item{set.seed}{logical. Should a random seed be set? Default is \code{FALSE}.}
  \item{seed}{integer indicating the random seed. Only valid if \code{set.seed = TRUE}.}
  \item{bootstrap}{logical. Should bootstrap percentile confidence intervals be calculated for the parameters and attributions? Default is \code{FALSE}. See details.}
  \item{nbootstrap}{integer. Number of bootstrap replicates.}
  \item{parallel}{logical. Should parallel calculations be used to obtain the  bootstrap percentile confidence intervals? Default is \code{FALSE}.}
  \item{type.parallel}{type of parallel operation to be used (if \code{parallel = TRUE}), with options: \code{"multicore"} and \code{"snow"}. Default is \code{"snow"}. See details.}
  \item{ncpus}{integer. Number of processes to be used in parallel operation: typically one would choose this to the number of available CPUs. Default is 4.}
  \item{...}{other arguments passed to or from the other functions.}
}

\details{
  The model is a generalized linear model with a non-canonical link function, which requires a restriction on the linear predictor (\eqn{\eta \ge 0}) to produce valid probabilities. This restriction is implemented in the optimization procedure, with an adaptive barrier algorithm, using the function \code{constrOptim} in the \code{stats} package.

 The variance-covariance matrix is based on the observed information matrix.

 This version of the the package only allows the calculation of non-parametric bootstrap percentile confidence intervals (CI), which are based on the 2.5 and 97.5 percentiles. Stratified bootstrap is applied to each category of the outcome. Also,  the function gives the user the option to do parallel calculation of the bootstrap CI. The \code{snow} parallel option is available for all operating systems (Windows, Linux, and Mac OS) while the \code{multicore} option is  only available for Linux and Mac OS systems. These two calculations are done by calling the \code{boot} function in the \code{boot} package. For more details see the documentation of the \code{boot} package.

 More information about the multinomial additive hazards model and the calculation of the contribution of chronic conditions to the disability prevalence can be found in the references.
}

\value{
  A list with arguments:
  \item{coefficients}{column matrix with the regression coefficients.}
  \item{ci}{matrix with confidence intervals calculated via bootstraping (if \code{bootstrap = TRUE}) or as the inverse of the observed information matrix.}
  \item{resDeviance}{residual deviance.}
  \item{df}{degrees of freedom.}
  \item{pvalue}{column matrix of p-values based on the Wald test. Only provided if \code{bootstrap = FALSE}.}
  \item{stdError}{column matrix with the standard errors for the parameter              estimates based on the inverse of the observed information matrix. Only provided if \code{bootstrap = FALSE}.}
  \item{vcov}{variance-covariance matrix (inverse of the observed information matrix). Only provided if \code{bootstrap = FALSE}.}
  \item{contribution}{for \code{type.attrib = "abs"} or \code{"rel"}, a matrix is provided. For \code{type.attrib = "both"}, a list with two matrices (\code{"abs"} and \code{"rel"}) is provided. This represents the contribution of each predictor in the model (usually diseases or conditions) to the disability prevalence. Percentile bootstrap confidence intervals are provided if \code{bootstrap = TRUE}.}
  \item{bootsRep}{matrix with the bootstrap replicates.}
  \item{bootstrap}{logical. Was bootstrap CI requested?}
  \item{call}{the matched call.}
}

\references{
Yokota, R.T.C., Van Oyen, H., Molenberghs, G. (2015). Use of a multinomial additive hazards model to assess the disability burden using cross-sectional data. Available at <\url{https://uhdspace.uhasselt.be/dspace/bitstream/1942/19407/1/11321252014093.pdf}>.
}

\author{
Renata T C Yokota. This function is based on the R code developed by Caspar W N Looman and Wilma J Nusselder for the binomial additive hazards model with modifications and adaptations for the multinomial case.
}

\seealso{
\code{\link{BinAddHaz}}
}

\examples{
  data(disabData)

  ## Model without bootstrap percentile confidence intervals and no attribution
\dontrun{
  fit1 <- MultAddHaz(dis.mult ~ diab + arth + stro , data = disabData, weights = wgt,
                     start = TRUE, start.val = c(rep(0.5, 4), rep(0.3, 4)))
  summary(fit1)

  ## Model with bootstrap percentile confidence intervals and attribution without
  ## stratification, no parallel calculation
  # Warning message due to the low number of bootstrap replicates

  fit2 <- MultAddHaz(dis.mult ~ diab + arth + stro , data = disabData, weights = wgt,
                     start = TRUE, start.val = c(rep(0.5, 4), rep(0.3, 4)), attrib = TRUE,
                     type.attrib = "both", set.seed = TRUE, seed = 111, bootstrap = TRUE,
                     nbootstrap = 5)
  summary(fit2)

  ## Model with bootstrap percentile confidence intervals with attribution based
  ## on the age variable and parallel calculation of bootstrap CI
  # Warning message due to the low number of bootstrap replicates

  diseases <- as.matrix(disabData[,c("diab", "arth", "stro")])

  fit3 <- MultAddHaz(dis.mult ~ factor(age) -1 + diseases: factor(age),
                     data = disabData, weights = wgt, start = TRUE,
                     start.val = c(rep(0.5, 8), rep(0.3, 8)), attrib.var = age,
                     attrib = TRUE, type.attrib = "both", set.seed = TRUE,
                     seed = 111, bootstrap = TRUE, nbootstrap = 5, parallel = TRUE,
                     type.parallel = "snow", ncpus = 4)
  summary(fit3)}
}

\keyword{models}
\keyword{regression}
\keyword{survey}
