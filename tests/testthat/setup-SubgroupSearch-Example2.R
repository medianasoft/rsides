##############################################################################

# Primary endpoint parameters

# Analysis strategy 1: Analysis of the binary endpoint without 
# accounting for any covariates
endpoint_parameters_e2_s1 = list(outcome_variable = "outcome", 
  type = "binary",
  label = "Outcome", 
  analysis_method = "Z-test for proportions", 
  direction = 1)

# Analysis strategy 2: Analysis of the continuous endpoint using an ANCOVA 
# model that accounts for two continuous covariates (cont1, cont2) and 
# two class/categorical covariates (class1, class2)
endpoint_parameters_e2_s2 = list(outcome_variable = "outcome", 
  type = "binary",
  label = "Outcome", 
  analysis_method = "Logistic regression", 
  cont_covariates = "cont1, cont2", 
  class_covariates = "class1, class2", 
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
data_set_parameters_e2 = list(data_set = binary,
  treatment_variable_name = "treatment",
  treatment_variable_control_value = "0",
  biomarker_names = biomarker_names,
  biomarker_types = biomarker_types)

cat("setup-SubgroupSearch-Example2.R executed.\n")