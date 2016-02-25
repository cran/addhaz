\name{BinAddHaz}
\alias{BinAddHaz}

\title{Fit Binomial Additive Hazards Models}

\description{
This function fits binomial additive hazards models subject to linear inequality constraints using the function \code{constrOptim} in the \code{stats} package for binary outcomes. Additionally, it calculates the cause-specific contributions to the disability prevalence based on the attribution method, as proposed by Nusselder and Looman (2004).
}

\usage{
BinAddHaz(formula, data, subset, weights, na.action, model = TRUE,
          contrasts = NULL, start = FALSE, start.val, attrib = FALSE,
          attrib.var, type.attrib = "abs", set.seed = FALSE, seed,
          bootstrap = FALSE, nbootstrap, parallel = FALSE,
          type.parallel = "snow", ncpus = 4,...)
}

\arguments{
  \item{formula}{A formula expression of the form response ~ predictors, similar to other regression models. In case of \code{attrib = TRUE}, the first predictor in the formula should be the \code{attrib.var}. See example.}
  \item{data}{an optional data frame or matrix containing the variables in the model. If not found in \code{data}, the variables are taken from environment(formula), typically the environment from which \code{BinAddHaz} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process. All observations are included by default.}
  \item{weights}{an optional vector of survey weights to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. The default is set by the \code{na.action} setting of options, and is \code{na.fail} if that is unset. The 'factory-fresh' default is \code{na.omit}.}
  \item{model}{logical. If \code{TRUE}, the model frame is included as a component of the returned object.}
  \item{contrasts}{an optional list to be used for some or all of the factors appearing as variables in the model formula.}
  \item{start}{logical. Will the user provide starting values? Default is \code{FALSE}, in which case the starting values will be the absolute values of the coefficients for a Poisson glm.}
  \item{start.val}{vector of the starting values. They should be provided if \code{start = TRUE} and lie inside the parameter space.}
  \item{attrib}{logical. Should the attribution of disability to chronic diseases/conditions be calculated? Default is \code{FALSE}.}
  \item{attrib.var}{character indicating the name of the attribution variable to be used if \code{attrib = TRUE}. If missing, the attribution variable is the intercept (no stratification of the attribution results). The attribution variable must be the first variable included in the linear predictor in \code{formula}. See example.}
  \item{type.attrib}{type of attribution to be calculated. The options are \code{"abs"} for  absolute contribution, \code{"rel"} for relative contribution, or \code{"both"} for both absolute and relative contributions. Default is \code{"abs"}.}
  \item{set.seed}{logical. Should a random seed be set? Default is \code{FALSE}.}
  \item{seed}{integer indicating the random seed. Only valid if \code{set.seed = TRUE}.}
  \item{bootstrap}{logical. Should bootstrap percentile confidence intervals be calculated for the parameters and attributions? Default is \code{FALSE}. See details.}
  \item{nbootstrap}{integer. Number of bootstrap replicates.}
  \item{parallel}{logical. Should parallel calculations be used to obtain the bootstrap percentile confidence intervals? Only valid if \code{bootstrap = TRUE}. Default is \code{FALSE}.}
  \item{type.parallel}{type of parallel operation to be used (if \code{parallel = TRUE}), with options: \code{"multicore"} and \code{"snow"}. Default is \code{"snow"}. See details.}
  \item{ncpus}{integer. Number of processes to be used in parallel operation: typically one would choose this to be the number of available CPUs. Default is 4.}
  \item{...}{other arguments passed to or from the other functions.}
}

\details{
 The model is a generalized linear model with a non-canonical link function, which requires a restriction on the linear predictor (\eqn{\eta \ge 0}) to produce valid probabilities. This restriction is implemented in the optimization procedure, with an adaptive barrier algorithm, using the function \code{constrOptim} in the \code{stats} package.

 The variance-covariance matrix is based on the observed information matrix.

 This version of the the package only allows the calculation of non-parametric bootstrap percentile confidence intervals (CI), which are based on the 2.5 and 97.5 percentiles. Also,  the function gives the user the option to do parallel calculation of the bootstrap CI. The \code{snow} parallel option is available for all operating systems (Windows, Linux, and Mac OS) while the \code{multicore} option is  only available for Linux and Mac OS systems. These two calculations are done by calling the \code{boot} function in the \code{boot} package. For more details, see the documentation of the \code{boot} package.

 More information about the binomial additive hazards model and the calculation of the contribution of chronic conditions to the disability prevalence can be found in the references.
}

\value{
  A list with arguments:
  \item{coefficients}{numeric vector with the regression coefficients.}
  \item{ci}{confidence intervals calculated via bootstraping (if \code{bootstrap = TRUE}) or as the inverse of the observed information matrix.}
  \item{resDeviance}{residual deviance.}
  \item{df}{degrees of freedom.}
  \item{pvalue}{numerical vector of p-values based on the Wald test. Only provided if  \code{bootstrap = FALSE}.}
  \item{stdError}{Numerical vector with the standard errors for the parameter estimates based on the inverse of the observed information matrix. Only provided if  \code{bootstrap = FALSE}.}
  \item{vcov}{variance-covariance (inverse of the observed information matrix). Only provided if bootstrap = FALSE.}
  \item{contribution}{for  \code{type.attrib = "abs"} or  \code{"rel"}, a matrix is provided. For  \code{type.attrib = "both"}, a list with two matrices ( \code{"abs"} and  \code{"rel"}) is provided. This represents the contribution of each predictor in the model (usually diseases or conditions) to the disability prevalence. Percentile bootstrap confidence intervals are provided if  \code{bootstrap = TRUE}.}
  \item{bootsRep}{matrix with the bootstrap replicates.}
  \item{bootstrap}{logical. Was bootstrap CI requested?}
  \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictor by the inverse of the link function.}
  \item{residuals}{difference between the response variable and the fitted values.}
  \item{call}{the matched call.}
}

\references{
Nusselder, W.J., Looman, C.W.N. (2004). Decomposition of differences in health expectancy by cause. Demography, 41(2), 315-334.

Nusselder, W.J., Looman, C.W.N. (2010). WP7: Decomposition tools: technical report on attribution tool. European Health Expectancy Monitoring Unit (EHEMU). Available at
<\url{http://www.eurohex.eu/pdf/Reports_2010/2010TR7.2_TR\%20on\%20attribution\%20tool.pdf}>.
}

\author{
Renata T C Yokota. This function is based on the R code developed by Caspar W N Looman and Wilma J Nusselder for non R-users, with modifications. Original R code is available upon request to Wilma J Nusselder (w.nusselder@erasmusmc.nl).
}

\seealso{
\code{\link{MultAddHaz}}
}

\examples{
  data(disabData)

  ## Model without bootstrap percentile confidence intervals and no attribution

  fit1 <- BinAddHaz(dis.bin ~ diab + arth + stro , data = disabData, weights = wgt)
  summary(fit1)

  ## Model with bootstrap percentile confidence intervals and attribution without
  ## stratification, no parallel calculation
  # Warning message due to the low number of bootstrap replicates
\dontrun{
  fit2 <- BinAddHaz(dis.bin ~ diab + arth + stro , data = disabData, weights = wgt,
                    attrib = TRUE, type.attrib = "both", set.seed = TRUE, seed = 111,
                    bootstrap = TRUE, nbootstrap = 5)
  summary(fit2)

  ## Model with bootstrap percentile confidence intervals with attribution based on
  ## the age variable and parallel calculation of bootstrap CI
  # Warning message due to the low number of bootstrap replicates

  diseases <- as.matrix(disabData[,c("diab", "arth", "stro")])

  fit3 <- BinAddHaz(dis.bin ~ factor(age) -1 + diseases:factor(age), data = disabData,
                    weights = wgt, attrib.var = age, attrib = TRUE, type.attrib = "both",
                    set.seed = TRUE, seed = 111, bootstrap = TRUE, nbootstrap = 10,
                    parallel = TRUE, type.parallel = "snow", ncpus = 4)
  summary(fit3)}
}

\keyword{models}
\keyword{regression}
\keyword{survey}
