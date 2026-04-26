\name{rsides-package}
\alias{rsides-package}
\alias{rsides}
\docType{package}
\title{Subgroup search}
\description{
The package implements a family of subgroup search algorithms based on the SIDES (Subgroup Identification based on Differential Effect Search) method for clinical trials with normally distributed, binary and time-to-event endpoints. The package supports complex analysis models with an adjustment for continuous and categorical covariates (analysis of covariance models, logistic regression models, Cox proportional hazards models).
}
\details{
\tabular{ll}{
Package: \tab rsides\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2024-05-27\cr
License: \tab GPL-2\cr
}
Key functions included in the package:
\itemize{
\item \link{SubgroupSearch}: Perform a SIDES-based subgroup search.
\item \link{GenerateReport}: Generate a detailed summary of subgroup search results in a Microsoft Word format.
}
The package comes with three example data sets:
\itemize{
\item \link{continuous}: Data set based on a trial with a continuous endpoint.
\item \link{binary}: Data set based on a trial with a binary endpoint.
\item \link{survival}: Data set based on a trial with a time-to-event endpoint. 
}
Three case studies are included in this manual to illustrate subgroup identification in clinical trials:
\itemize{
\item \link{Example1}: Subgroup search in a clinical trial with a continuous endpoint.
\item \link{Example2}: Subgroup search in a clinical trial with a binary endpoint.
\item \link{Example3}: Subgroup search in a clinical trial with a time-to-event endpoint. 
}
}

\references{
Lipkovich, I., Dmitrienko, A., Denne, J., Enas, G. (2011). Subgroup Identification based on Differential Effect Search (SIDES): A recursive partitioning method for establishing response to treatment in patient subpopulations. \emph{Statistics in Medicine}. 30, 2601-2621. 

Lipkovich, I., Dmitrienko A.  (2014). Strategies for identifying predictive biomarkers and subgroups with enhanced treatment effect in clinical trials using SIDES. \emph{Journal of Biopharmaceutical Statistics}. 24, 130-153.

Lipkovich, I., Dmitrienko, A. (2014). Biomarker identification in clinical trials. \emph{Clinical and Statistical Considerations in Personalized Medicine}. Carini, C., Menon, S., Chang, M. (editors). Chapman and Hall/CRC Press, New York.

Lipkovich, I., Dmitrienko, A., D'Agostino, R.B. (2017). Tutorial in Biostatistics: Data-driven subgroup identification and analysis in clinical trials. \emph{Statistics in Medicine}. 36, 136-196.
}
\keyword{datasets}
