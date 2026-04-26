# Primary endpoint parameters

# Analysis of the binary endpoint without accounting for any covariates
endpoint_parameters = list(outcome_variable = "outcome", 
  type = "binary",
  label = "Outcome", 
  analysis_method = "Z-test for proportions", 
  direction = 1)

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

# Algorithm parameters for the basic SIDES procedure

# Algorithm
subgroup_search_algorithm = "SIDES procedure"

# Number of permutations to compute multiplicity-adjusted treatment effect p-values within promising subgroups
n_perms_mult_adjust = 1

# Number of processor cores 
ncores = 1

# Default values for the search depth (2), search width (2), maximum number of unique values for continuous biomarkers (20)

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

test_that("It works!!!", {

  results = SubgroupSearch(parameters)

  patient_subgroups = results$patient_subgroups

  expect_equal(results$patient_subgroups$Subgroups$size, c(359, 119, 268, 123, 78))
  expect_equal(round(results$patient_subgroups$`Variable Importance`$vi, 3), c(0.427, 0.407, 0.102, 0.076, 0))

})
