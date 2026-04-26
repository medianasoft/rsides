\name{GenerateReport}
\alias{GenerateReport}
\title{
Generate a Word-based summary of subgroup search results  
}
\description{
This function creates a detailed summary of subgroup search results in a Microsoft Word format.
}
\usage{
GenerateReport(results, report_title, report_filename)
}
\arguments{
  \item{results}{Object of class \code{SubgroupSearchResults} created by the \code{SubgroupSearch} function.  
  }
  \item{report_title}{Character value defining the report's title.
  }
  \item{report_filename}{Character value defining the report's filename. The report is saved in the current working directory.
  }
}
\value{No return value, called for side effects}
\seealso{
  \code{\link{SubgroupSearch}}
}
\examples{
  \donttest{
    # Example of a subgroup search in a trial with a binary endpoint
    # This example used an example data set (binary) that comes with the package
    # Primary endpoint parameters

    endpoint_parameters = list(outcome_variable = "outcome", 
                              type = "binary",    
                              label = "Outcome",
                              analysis_method = "Z-test for proportions", 
                              direction = 1)

    ##############################################################################

    # Data set parameters

    # Set of candidate biomarkers
    biomarker_names = c("biomarker1", "biomarker2", 
                        "biomarker3", "biomarker4", 
                        "biomarker5")

    # Biomarker type 
    biomarker_types = c(rep("numeric", 4), "nominal")

    # Data set parameters
    data_set_parameters = list(data_set = binary,
                              treatment_variable_name = "treatment",
                              treatment_variable_control_value = "0",
                              biomarker_names = biomarker_names,
                              biomarker_types = biomarker_types)

    ##############################################################################

    # Algorithm parameters for the basic SIDES procedure

    # Algorithm
    subgroup_search_algorithm = "SIDES procedure"

    # Number of permutations to compute multiplicity-adjusted treatment 
    # effect p-values within promising subgroups
    n_perms_mult_adjust = 10

    # Number of processor cores (use less or equal number of CPU cores on the current host)
    ncores = 1

    # Default values for the search depth (2), search width (2), 
    # maximum number of unique values for continuous biomarkers (20)

    # Algorithm parameters
    algorithm_parameters = list(
      n_perms_mult_adjust = n_perms_mult_adjust,
      min_subgroup_size = 60,
      subgroup_search_algorithm = subgroup_search_algorithm,
      ncores = ncores,
      random_seed = 3011)

    # Perform subgroup search

    # List of all parameters
    parameters = list(endpoint_parameters = endpoint_parameters,
                      data_set_parameters = data_set_parameters,
                      algorithm_parameters = algorithm_parameters)

    # Perform subgroup search
    results = SubgroupSearch(parameters)

    # Simple summary of subgroup search results
    results

    # Generate a detailed Word-based report with a summary of subgroup search results
    GenerateReport(
      results,
      report_title = "Subgroup search report", 
      report_filename = tempfile(
        "Subgroup search report.docx", 
        fileext=".docx"
      )
    ) 
  }
}
