\name{survival}
\alias{survival}
\docType{data}
\title{Example data set (time-to-event endpoint)}
\description{
Example data set based on a trial with a time-to-event endpoint.
}
\usage{data(survival)}
\format{
  A data set with 359 observations and 12 variables:
  \describe{
    \item{\code{outcome}}{Time-to-event outcome variable. A larger value of the outcome variable (longer survival) indicates a beneficial effect.}
    \item{\code{outcome_censor}}{Binary censoring variable (0 or 1). The value of 1 indicates censoring.}    
    \item{\code{treatment}}{Binary treatment variable (0 or 1). The values of 0 and 1 correspond to the control and experimental treatments, respectively.}
    \item{\code{cont1}}{Continuous covariate used in the primary analysis model.}
    \item{\code{cont2}}{Continuous covariate used in the primary analysis model.}
    \item{\code{class1}}{Class/categorical covariate used in the primary analysis model.}
    \item{\code{class2}}{Class/categorical covariate used in the primary analysis model.}  
    \item{\code{biomarker1}}{Numeric biomarker.}  
    \item{\code{biomarker2}}{Numeric biomarker.}  
    \item{\code{biomarker3}}{Numeric biomarker.}  
    \item{\code{biomarker4}}{Numeric biomarker.}  
    \item{\code{biomarker5}}{Nominal biomarker.}  
    }
}
\keyword{datasets}
