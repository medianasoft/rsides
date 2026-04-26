test_that("Example 1-b works", {
  ##############################################################################
  
  # Algorithm parameters for the Fixed SIDEScreen procedure
  
  # Algorithm
  subgroup_search_algorithm = "Fixed SIDEScreen procedure"
  
  # Number of permutations to compute multiplicity-adjusted treatment 
  # effect p-values within promising subgroups
  n_perms_mult_adjust = 10
  
  # Number of processor cores (use less or equal number of CPU cores on the current host)
  ncores = 1
  
  # Number of biomarkers selected for the second stage in the Fixed SIDEScreen algorithm
  n_top_biomarkers = 3
  
  # Default values for the search depth (2), search width (2), 
  # maximum number of unique values for continuous biomarkers (20)
  
  # Algorithm parameters
  algorithm_parameters = list(
    n_perms_mult_adjust = n_perms_mult_adjust,
    min_subgroup_size = 60,
    subgroup_search_algorithm = subgroup_search_algorithm,
    ncores = ncores,
    n_top_biomarkers = n_top_biomarkers,
    random_seed = 3011)
  
  # Perform subgroup search
  
  # List of all parameters
  parameters = list(endpoint_parameters = endpoint_parameters_e1_s2,
                    data_set_parameters = data_set_parameters_e1,
                    algorithm_parameters = algorithm_parameters)
  
  results = SubgroupSearch(parameters)
  
  GenerateReport(results,
    report_title = "Subgroup search report", 
    report_filename = tempfile("Continuous endpoint (Fixed SIDEScreen).docx", fileext = ".docx")
  )

  # Test patient subgroups
  patient_subgroups = results$patient_subgroups
  expect_equal(results$patient_subgroups$Subgroups$size, c(359, 292, 185, 209, 290, 205, 198))
  expect_equal(round(results$patient_subgroups$`Variable Importance`$vi, 3), c(0, 0, 0, 0, 0))

  # # Reset value for expect_snapshot
  # results$patient_subgroups$`Elapsed time (seconds)` = 0.1
  # # Check all values
  # expect_snapshot(results)
})
