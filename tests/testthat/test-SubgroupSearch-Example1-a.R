test_that("Example 1-a works", {
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
  parameters = list(endpoint_parameters = endpoint_parameters_e1_s2,
                    data_set_parameters = data_set_parameters_e1,
                    algorithm_parameters = algorithm_parameters)
  
  results = SubgroupSearch(parameters)

  GenerateReport(results,
    report_title = "Subgroup search report",
    report_filename = tempfile("Continuous endpoint (SIDES).docx", fileext = ".docx")
  )

  # ***************************************
  # Subgroup search results
  # ***************************************
  #            Subgroup Total size (Control, Treatment)           Estimate (SE)
  #            Overall population                  359 (182, 177)         0 (0)
  #                 biomarker1>5:                  292 (143, 149)         0 (0)
  #    biomarker1>5:biomarker1>7:                    185 (89, 96)         0 (0)
  #   biomarker1>5:biomarker2>-2:                  209 (100, 109)         0 (0)
  #                biomarker2>-3:                  290 (149, 141)         0 (0)
  #   biomarker2>-3:biomarker1>6:                   205 (97, 108)         0 (0)
  #  biomarker2>-3:biomarker2>-1:                   198 (103, 95)         0 (0)
  #  P-value (Adjusted p-value)
  #                      1 (NA)
  #                       1 (1)
  #                       1 (0)
  #                       1 (0)
  #                       1 (0)
  #                       1 (0)
  #                       1 (0)

  # Test patient subgroups
  patient_subgroups = results$patient_subgroups
  expect_equal(results$patient_subgroups$Subgroups$size, c(359, 292, 185, 209, 290, 205, 198))
  expect_equal(round(results$patient_subgroups$`Variable Importance`$vi, 3), c(0.000, 0.000, 0.000, 0.000, 0.000))

  # # Reset value for expect_snapshot
  # results$patient_subgroups$`Elapsed time (seconds)` = 0.5
  # # Check all values
  # expect_snapshot(results)
})
