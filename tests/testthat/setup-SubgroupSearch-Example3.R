##############################################################################

# Primary endpoint parameters

# Analysis strategy 1: Analysis of the continuous endpoint without 
# accounting for any covariates
endpoint_parameters_e3_s1 = list(
  outcome_variable = "outcome", 
  outcome_censor_variable = "outcome_censor",
  outcome_censor_value = "1",
  type = "survival",
  label = "Outcome", 
  analysis_method = "Log-rank test", 
  direction = 1
)

# Analysis strategy 2: Analysis of the continuous endpoint using a Cox model 
# that accounts for two continuous covariates (cont1, cont2) and 
# two class/categorical covariates (class1, class2)
endpoint_parameters_e3_s2 = list(
  outcome_variable = "outcome", 
  outcome_censor_variable = "outcome_censor",
  outcome_censor_value = "1",
  type = "survival",
  label = "Outcome", 
  analysis_method = "Cox regression", 
  cont_covariates = "cont1, cont2", 
  class_covariates = "class1, class2", 
  direction = 1
)

##############################################################################

# Data set parameters

# Set of candidate biomarkers
biomarker_names = c("biomarker1", "biomarker2", 
                    "biomarker3", "biomarker4", 
                    "biomarker5")

# Biomarker type 
biomarker_types = c(rep("numeric", 4), "nominal")

# Data set parameters
data_set_parameters_e3 = list(data_set = survival,
  treatment_variable_name = "treatment",
  treatment_variable_control_value = "0",
  biomarker_names = biomarker_names,
  biomarker_types = biomarker_types)

cat("setup-SubgroupSearch-Example3.R executed.\n")