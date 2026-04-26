\name{SubgroupSearch}
\alias{SubgroupSearch}
\alias{print.SubgroupSearchResults}
\title{
Perform a SIDES-based subgroup search  
}
\description{
This function performs a SIDES-based subgroup search for clinical trials with normally distributed, binary and time-to-event endpoints. The function implements the following subgroup search procedures:
\itemize{
\item SIDES procedure: Basic subgroup search procedure (Lipkovich et al., 2011).
\item Fixed and Adaptive SIDEScreen procedures: Two-stage subgroup search procedure with biomarker selection (Lipkovich and Dmitrienko, 2014).
}

}
\usage{
SubgroupSearch(parameters)
}
\arguments{
  \item{parameters}{List defining the subgroup search's parameters. The list includes three sublists: \code{endpoint_parameters}, \code{data_set_parameters}, \code{algorithm_parameters}. The parameters that need to be defined with each of these lists are defined below.
    \itemize{

      \item \code{endpoint_parameters}: List defining the parameters of the primary outcome variable and analysis method. The following parameters need to be specified:
      
        \itemize{
          \item \code{outcome_variable}: Character value defining the name of the outcome variable in the data set specified in \code{data_set_parameters}. 
          \item \code{label}: Character value defining the outcome variable's label.
          \item \code{type}: Character value defining the outcome variable's type: \cr \code{"continuous"} if the outcome variable is a continuous endpoint, \cr \code{"binary"} if the outcome variable is a binary endpoint, \cr \code{"survival"} if the outcome variable is a time-to-event endpoint.
          \item \code{outcome_censor_variable}: Character value defining the name of the censoring variable in the data set. This argument is required only if the outcome variable is a time-to-event endpoint.
          \item \code{outcome_censor_value}: Character value defining the value of the censoring variable that corresponds to censored outcomes. This argument is required only if the outcome variable is a time-to-event endpoint.
          \item \code{direction}: Numeric value defining the direction of a beneficial effect: \cr \code{1} if a higher value of the outcome variable indicates a beneficial effect, \cr \code{-1} if a lower value of the outcome variable indicates a beneficial effect.
          \item \code{analysis_method}: Character value defining the analysis method for the outcome variable: \cr \code{"T-test"}, \code{"Z-test for proportions"} or \code{"Log-rank test"} for continuous, binary and  time-to-event endpoint without covariate adjustment, respectively, \cr (\code{"ANCOVA"}, \code{"Logistic regression"} or \code{"Cox regression"} for continuous, binary and  time-to-event endpoint with covariate adjustment, respectively. The covariates to be included in the model are specified using \cr \code{cont_covariates} and \code{class_covariates}. 
          \item \code{cont_covariates}: Vector of character values defining the names of the continuous covariates to be included in the model if \code{analysis_method} is set to \code{"ANCOVA"}, \code{"Logistic regression"} or \code{"Cox regression"}. This argument is not required if \code{analysis_method} is set to \code{"T-test"}, \cr \code{"Z-test for proportions"} or \code{"Log-rank test"}.
          \item \code{class_covariates}: Vector of character values defining the names of the class/categorical covariates to be included in the model if \code{analysis_method} is set to \code{"ANCOVA"}, \code{"Logistic regression"} or \code{"Cox regression"}. This argument is not required if \code{analysis_method} is set to \code{"T-test"}, \cr  \code{"Z-test for proportions"} or \code{"Log-rank test"}.
        }

      \item \code{data_set_parameters}: List defining the data set and its characteristics. The following parameters need to be specified:
  
        \itemize{
          \item \code{data_set}: Character value defining the name of the clinical trial data set. The package comes with three data sets that are used in the examples:\cr
                                  \link{continuous}: Data set based on a trial with a continuous endpoint.\cr
                                  \link{binary}: Data set based on a trial with a binary endpoint.\cr
                                  \link{survival}: Data set based on a trial with a time-to-event endpoint. 
          \item \code{treatment_variable_name}: Character value defining the name of the treatment variable in the data set. Only two-arm trials are supported.
          \item \code{treatment_variable_control_value}: Character value defining the value of the treatment variable that corresponds to the control arm. 
          \item \code{biomarker_names}: Vector of character values defining the names of the candidate biomarkers.
          \item \code{biomarker_types}: Vector of character values defining the types of the candidate biomarkers:\cr
                                  \code{"numeric"} if the biomarker is a continuous variable,  \cr 
                                  \code{"nominal"} if the biomarker is a nominal variable.
          \item \code{biomarker_levels}: Vector of numeric values defining the first subgroup search level at which each biomarker is introduced. For example, if the level is 1, the biomarker will be used in all subgroups and, if the level is 2, the biomarker will be used only in the second-level and deeper subgroups. By default, the level is set to 1 for each biomarker.
        }

        \item \code{algorithm_parameters}: List of the subgroup search algorithm's parameters. The following parameters need to be specified:
        
          \itemize{
            \item \code{subgroup_search_algorithm}: Character value defining the name of the subgroup search algorithm: \cr \code{"SIDES procedure"}: Basic subgroup search procedure, \cr \code{"Fixed SIDEScreen procedure"}: SIDEScreen procedure with a fixed number of biomarkers selected for the second stage, \cr \code{"Adaptive SIDEScreen procedure"}: SIDEScreen procedure with a data-driven number of biomarkers selected for the second stage. 
            \item \code{depth}: Integer value defining the subgroup search depth. The default value is 2.
            \item \code{width}: Integer value defining the subgroup search width. The default value is 2.
            \item \code{gamma}: Vector of numeric values defining the complexity parameters (also known as the child-to-parent ratios). The complexity parameters must be between 0 and 1 (unless no complexity control is applied at a certain search level in which case the complexity parameter at this level is set to NA) and the vector's length must be equal to the search depth. The default value is 1 at each search level, i.e., by default \code{gamma} is equal to \code{c(1, 1)} if the depth parameter is set to 2.
            \item \code{n_perms_mult_adjust}: Integer value defining the number of permutations for computing multiplicity-adjusted treatment effect p-values within the promising subgroups. The default value is 1000.
            \item \code{ncores}: Integer value defining the number of processor cores that will be used for computing multiplicity-adjusted treatment effect p-values. The default value is 1.
            \item \code{nperc}: Integer value defining the minimum number of unique values for continuous biomarkers for applying a percentile transformation. The default value is 20.
            \item \code{min_subgroup_size}: Integer value defining the minimum total number of patients in a promising subgroup. The default value is 30.
            \item \code{n_top_biomarkers}: Integer value defining the number of best biomarkers selected for the second stage of the SIDEScreen procedure. This argument is only required for the Fixed SIDEScreen procedure. The default value is 3.
            \item \code{multiplier}: Numeric value defining the multiplier in the data-driven biomarker selection rule for the second stage of the SIDEScreen procedure. This argument is only required for the Adaptive SIDEScreen procedure. The default value is 1.
            \item \code{n_perms_vi_score}: Numeric value defining the number of permutations used in the data-driven biomarker selection rule for the second stage of the SIDEScreen procedure. This argument is only required for the Adaptive SIDEScreen procedure. The default value is 100.
            \item \code{random_seed}: Integer value defining the random seed that will be used for computing permutation-based multiplicity-adjusted treatment effect p-values. The default value is 49291.
          }

    }
  }

}
\value{
The function returns an object of class \samp{SubgroupSearchResults}. This object is a list with the following components:
  \item{parameters}{a list containing the user-specified parameters, i.e., endpoint, data set and algorithm parameters.}
  \item{patient_subgroups}{a list containing the subgroup search results, in particular, a summary of the subgroup effects, a variable importance summary and a brief summary of the algorithm's parameters. The summary of subgroup effects provides information on the treatment effect in the overall population and promising subgroups identified by the selected algorithm. The summary includes the number of patients in each subgroup by trial arm, treatment effect estimate as well as raw and multiplicity-adjusted p-values. For a continuous primary endpoint, the treatment effect estimate is defined as the sample mean difference or the mean difference computed from the ANCOVA model. For a binary primary endpoint, the treatment effect estimate is defined as the sample difference in proportions if the Z-test for proportions is carried out or the odds ratio computed from the logistic regression model. Finally, if the primary endpoint is a time-to-event endpoint, the treatment effect estimate is defined as the hazard ratio based on an exponential distribution assumption if the analysis is based on the log-rank test or the hazard ratio computed from the Cox proportional hazards model if a model-based analysis is employed.}
  A detailed summary of the subgroup search results can be generated using the \link{GenerateReport} function.  
}

\seealso{
  \code{\link{GenerateReport}}
}
\references{
Lipkovich, I., Dmitrienko, A., Denne, J., Enas, G. (2011). Subgroup Identification based on Differential Effect Search (SIDES): A recursive partitioning method for establishing response to treatment in patient subpopulations. \emph{Statistics in Medicine}. 30, 2601-2621. 

Lipkovich, I., Dmitrienko A.  (2014). Strategies for identifying predictive biomarkers and subgroups with enhanced treatment effect in clinical trials using SIDES. \emph{Journal of Biopharmaceutical Statistics}. 24, 130-153.

Lipkovich, I., Dmitrienko, A., D'Agostino, R.B. (2017). Tutorial in biostatistics: data-driven subgroup identification and analysis in clinical trials. \emph{Statistics in Medicine}.36, 136-196. 
}

