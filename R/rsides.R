require(Rcpp)
require(RcppEigen)
require(xml2)
require(foreach) 
require(doParallel)
require(doRNG)
require(officer)
require(flextable)
require(parallel)

trim0<-function(s)
{
  return(gsub("(^ +)|( +$)", "", s))
}

# Remove all whitespaces from a string
RemoveWhitespaces = function (string) gsub(" ", "", string, fixed = TRUE)

# Remove the leading or trailing whitespace from a string
trim = function (string) gsub("^\\s+|\\s+$", "", string)

# guesses type of variable in dt data frame
guesstype<-function(dt)
{ctype<-rep('',ncol(dt))
 for (i in 1:ncol(dt))
 {if (is.numeric(dt[1,i])) ctype[i]<-"numeric" else ctype[i]<-"char"}
 return (ctype)
}

convert2num<-function(dt,index.nom)
  # recodes values in columns of dt indexed by index.nom into integers 1,2,3, and saves dictionary
  
{
  cnames<-colnames(dt)[index.nom]
  #  dt<-as.data.frame(dt[,index.nom])
  dt.coded<-as.data.frame(dt)
  #  colnames(dt.coded)<-colnames(dt)
  
  for (j.col in 1:length(index.nom)) # walks only over variables that need to be recoded
  {
    mycol.factor<-as.factor(dt[,index.nom[j.col]])
    
    lev<-levels(mycol.factor)
    dictionary<-as.data.frame(matrix(0,length(lev),2))
    colnames(dictionary)<-c("original","coded")
    for (i in 1:nrow(dictionary))
    {dictionary[i,1]<-lev[i]
     dictionary[i,2]<-i
    }
    
    # applies dictionary to code data
    repl.tmp<-rep(NA,nrow(dt))
    for (i in 1:nrow(dictionary))
    {repl.tmp[mycol.factor==dictionary[i,1]]<-dictionary[i,2]
    }
    dt.coded[,index.nom[j.col]]<-as.numeric(repl.tmp)
    if (j.col==1) 
    {dic.list<-list(dictionary)} else {dic.list<-append(dic.list,list(dictionary))}
  }
  names(dic.list)<-cnames
  all<-list(dt.coded,dic.list)
  names(all)<-c("data","dictionary")
  return(all)
}

# generates xml code for SIDES from sides parameters
AQ<-function(s) {
  return(paste0("\"",s,"\""))
}

# Perform data standardization and create a data frame with standardized data as well as a data dictionary
StandardizedSIDESParameters = function(parameters) {
  
  df = parameters$data_set_parameters$data_set
  if (is.null(df)) stop("Missing data set", call. = FALSE)
  if (!is.data.frame(df)) stop("The data set must be passed as a data frame", call. = FALSE)

  vnames=colnames(df)

  biomarker_names = parameters$data_set_parameters$biomarker_names
  if (is.null(biomarker_names)) {    
    stop("Biomarker names must be specified", call. = FALSE)
  }

  n_biomarkers=length(biomarker_names)

  biomarker_types = parameters$data_set_parameters$biomarker_types
  if (is.null(biomarker_types)) {    
    stop("Biomarker types must be specified", call. = FALSE)
  }
  
  if (n_biomarkers != length(biomarker_types)) stop("The biomarker type must be specified for all biomarkers", call. = FALSE)
  
  biomarker_types_standardized = rep(0, length(biomarker_types))
  
  for (i in 1:length(biomarker_types)) {
    if (!is.na(biomarker_types[i])) {
      if (!(biomarker_types[i] %in% c("numeric", "nominal"))) stop("The biomarker type must be either numeric or nominal", call. = FALSE)
      # Parameter standardization
      if (biomarker_types[i] == "numeric") biomarker_types_standardized[i] = 1
      if (biomarker_types[i] == "nominal") biomarker_types_standardized[i] = 2
      
    }
  }

  biomarker_levels = parameters$data_set_parameters$biomarker_levels
  if (is.null(biomarker_levels)) {    
    biomarker_levels = rep(1, n_biomarkers)
  }
  parameters$data_set_parameters$biomarker_levels = biomarker_levels

  if (sum(is.na(match(biomarker_names, vnames)))>0) stop("Some specified biomarkers are not found in the data set", call. = FALSE)

  if (n_biomarkers != length(biomarker_levels)) stop("Biomarker levels are not specified correctly", call. = FALSE)
  
  treatment_variable_name = parameters$data_set_parameters$treatment_variable_name
  if (is.null(treatment_variable_name)) {    
    stop("Treatment variable name names must be specified", call. = FALSE)
  }

  if (!is.null(parameters$data_set_parameters$treatment_variable_control_value)) {    
    treatment_variable_control_value = parameters$data_set_parameters$treatment_variable_control_value
  } else {
    stop("Value identifying the control arm must be specified", call. = FALSE)
  }

  if (!is.null(parameters$endpoint_parameters$outcome_variable)) {    
    outcome_variable_name = parameters$endpoint_parameters$outcome_variable
    parameters$data_set_parameters$outcome_variable_name = outcome_variable_name
  } else {
    stop("Outcome variable name must be specified", call. = FALSE)
  }

  if (!is.null(parameters$endpoint_parameters$direction)) {    
    outcome_variable_direction = parameters$endpoint_parameters$direction
    parameters$data_set_parameters$outcome_variable_direction = outcome_variable_direction
  } else {
    stop("Outcome variable direction must be specified", call. = FALSE)
  }

  if (!is.null(parameters$endpoint_parameters$type)) {    
    outcome_variable_type = parameters$endpoint_parameters$type
    parameters$data_set_parameters$outcome_variable_type = outcome_variable_type
  } else {
    stop("Outcome type must be specified", call. = FALSE)
  }

  if (!(outcome_variable_direction %in% c(- 1, 1))) stop("The outcome variable direction must be -1 or 1", call. = FALSE)
  
  if (outcome_variable_type=="binary" & outcome_variable_direction==-1) df$outcome=1-df$outcome

  outcome_censor_name = parameters$endpoint_parameters$outcome_censor_variable
  parameters$data_set_parameters$outcome_censor_name = outcome_censor_name
  if (!is.null(outcome_censor_name)) {
    censor.index=match(outcome_censor_name,vnames)
    if (is.na(censor.index)) censor.index = NULL
  } else {
    censor.index = NULL
  }
  
  outcome_censor_value = parameters$endpoint_parameters$outcome_censor_value
  parameters$data_set_parameters$outcome_censor_value = outcome_censor_value
  if (outcome_variable_type=="survival" & (is.null(outcome_censor_name) | is.null(censor.index))) stop("Outcome censoring name must be specified", call. = FALSE)
  if (outcome_variable_type=="survival" & is.null(outcome_censor_value)) stop("Censoring value must be specified", call. = FALSE)
  
  outcome.index=match(outcome_variable_name, vnames)
  if (is.na(outcome.index) | is.null(outcome_variable_name)) stop("Outcome variable is not found in the data set", call. = FALSE)

  trt.index=match(treatment_variable_name,vnames)
  if (is.na(trt.index)| is.null(treatment_variable_name)) stop("Treatment variable is not found in the data set", call. = FALSE)
  
  treatment = (df[,trt.index]!=treatment_variable_control_value)*1

  outcome = df[, outcome.index]

  if (is.null(censor.index)) 
  {outcome_censor=rep(0,nrow(df))} else
  {outcome_censor = df[,censor.index]}
  
  if (outcome_variable_type=="survival")
  {outcome_censor=(outcome_censor==outcome_censor_value)*1}

  if (any(is.na(treatment))) stop("The treatment variable may not have missing values", call. = FALSE)
  if (any(is.na(outcome))) stop("The outcome variable may not have missing values", call. = FALSE)
  if (outcome_variable_type=="survival" & any(is.na(outcome_censor))) stop("The censoring variable may not have missing values", call. = FALSE)
  
  biomarker_df = as.data.frame(df[,match(biomarker_names,colnames(df))])
  
  index.char=NULL
  coded.all=NULL
  if (length(biomarker_names[biomarker_types=="nominal"]) > 0) {
    index.nom<-match(biomarker_names[biomarker_types=="nominal"], biomarker_names) # index nominal vars among biomarkers
    # of all nominal variables, if any find character types that need to be recoded
    gtypes<-guesstype(as.data.frame(biomarker_df[,index.nom]))
    charcov<-biomarker_names[index.nom][gtypes=="char"]
    if (length(charcov)>0)
    {index.char=match(charcov,biomarker_names)}  #finds index of character nominal variables in the biomarker data set 
    
  }
  
  if (!is.null(index.char)) {
    coded.all<-convert2num(biomarker_df,index.char)
    biomarker_df=as.data.frame(coded.all[[1]]) # replaces the original biomarker data with decoded 
  }
  
  df=cbind(biomarker_df,treatment,outcome,outcome_censor)
  colnames(df)=c(biomarker_names,"treatment","outcome","outcome_censor")

  if (!is.na(outcome_variable_type)) {
    if (!(outcome_variable_type %in% c("continuous", "binary", "survival"))) stop("The outcome type must be continuous, binary or survival", call. = FALSE)
    # Parameter standardization   
    if (outcome_variable_type == "continuous") outcome_type_standardized = 1
    if (outcome_variable_type == "binary") outcome_type_standardized = 2
    if (outcome_variable_type == "survival") outcome_type_standardized = 3
  }

  analysis_method = parameters$endpoint_parameters$analysis_method 

  # Default values  
  analysis_method_standardized = NA 
  cont_covariates = as.matrix(1:length(treatment))
  class_covariates = as.matrix(1:length(treatment))
  n_cont_covariates = 0
  n_class_covariates = 0 

  # Create a matrix of covariates
  if (outcome_type_standardized == 1) {

    if (!(analysis_method %in% c("T-test", "ANCOVA"))) stop("The analysis method must be T-test or ANCOVA", call. = FALSE)    

    # Parameter standardization   
    if (analysis_method == "T-test") analysis_method_standardized = 1
    if (analysis_method == "ANCOVA") {

      cont_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$cont_covariates), ","))

      class_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$class_covariates), ","))

      n_cont_covariates = length(cont_names)
      n_class_covariates = length(class_names)

      # Analysis method
      # analysis_method = 1: Basic test without covariates 
      # analysis_method = 2: ANCOVA with 1 continuous covariate
      # analysis_method = 3: ANCOVA with 2 continuous covariates
      # analysis_method = 4: ANCOVA with 2 continuous covariates and 1 class covariate
      # analysis_method = 5: ANCOVA with 3 continuous covariates and 1 class covariate
      # analysis_method = 6: ANCOVA with 4 continuous covariates and 1 class covariate

      if (n_cont_covariates == 1 & n_class_covariates == 0) analysis_method_standardized = 2
      if (n_cont_covariates == 2 & n_class_covariates == 0) analysis_method_standardized = 3
      if (n_cont_covariates == 2 & n_class_covariates == 1) analysis_method_standardized = 4
      if (n_cont_covariates == 3 & n_class_covariates == 1) analysis_method_standardized = 5
      if (n_cont_covariates == 4 & n_class_covariates == 1) analysis_method_standardized = 6

      if (is.na(analysis_method_standardized)) stop("The number of continuous and class covariates in the model is not correctly defined.", call. = FALSE)
   
      # Extract continuous covariates
      if (n_cont_covariates > 0) {
        # cont_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$cont_covariates), ","))
        cont_covariates = as.matrix(parameters$data_set_parameters$data_set[,match(cont_names, colnames(parameters$data_set_parameters$data_set))])
      } 

      # Extract class covariates
      if (n_class_covariates > 0) {
        # class_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$class_covariates), ","))
        class_covariates = as.matrix(parameters$data_set_parameters$data_set[,match(class_names, colnames(parameters$data_set_parameters$data_set))])
      } 

    }

  }

  # Create a matrix of covariates
  if (outcome_type_standardized == 2) {

    if (!(analysis_method %in% c("Z-test for proportions", "Logistic regression"))) stop("The analysis method must be Z-test for proportions or Logistic regression", call. = FALSE)    

    # Parameter standardization   
    if (analysis_method == "Z-test for proportions") analysis_method_standardized = 1
    if (analysis_method == "Logistic regression") {

      cont_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$cont_covariates), ","))

      class_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$class_covariates), ","))

      n_cont_covariates = length(cont_names)
      n_class_covariates = length(class_names)

      if (n_cont_covariates > 0 | n_class_covariates > 0) analysis_method_standardized = 2

      if (is.na(analysis_method_standardized)) stop("The number of continuous and class covariates in the model is not correctly defined", call. = FALSE)

      # Extract continuous covariates
      if (n_cont_covariates > 0) {
        cont_covariates = as.matrix(parameters$data_set_parameters$data_set[,match(cont_names, colnames(parameters$data_set_parameters$data_set))])
      } 

      # Extract class covariates
      if (n_class_covariates > 0) {
        class_covariates = as.matrix(parameters$data_set_parameters$data_set[,match(class_names, colnames(parameters$data_set_parameters$data_set))])
      }        

    }

  }

  # Create a matrix of covariates
  if (outcome_type_standardized == 3) {

    if (!(analysis_method %in% c("Log-rank test", "Cox regression"))) stop("The analysis method must be Log-rank test or Cox regression", call. = FALSE)    

    # Parameter standardization   
    if (analysis_method == "Log-rank test") analysis_method_standardized = 1
    if (analysis_method == "Cox regression") {

    if (any(duplicated(outcome))) stop("There must be no ties for the time-to-event outcome variable", call. = FALSE)

      cont_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$cont_covariates), ","))

      class_names = unlist(strsplit(RemoveWhitespaces(parameters$endpoint_parameters$class_covariates), ","))

      n_cont_covariates = length(cont_names)
      n_class_covariates = length(class_names)

      if (n_cont_covariates > 0 | n_class_covariates > 0) analysis_method_standardized = 2

      if (is.na(analysis_method_standardized)) stop("The number of continuous and class covariates in the model is not correctly defined", call. = FALSE)

      # Extract continuous covariates
      if (n_cont_covariates > 0) {
        cont_covariates = as.matrix(parameters$data_set_parameters$data_set[,match(cont_names, colnames(parameters$data_set_parameters$data_set))])
      } 

      # Extract class covariates
      if (n_class_covariates > 0) {
        class_covariates = as.matrix(parameters$data_set_parameters$data_set[,match(class_names, colnames(parameters$data_set_parameters$data_set))])
      }        

    }    
  }

  # Create the list of covariate to be passed to the SIDES library
  covariate_object = list(outcome = outcome,
                          censor = outcome_censor,
                          treatment = treatment,
                          cont_covariates = cont_covariates,
                          class_covariates = class_covariates,
                          n_cont_covariates = n_cont_covariates,
                          n_class_covariates = n_class_covariates)

  if (!is.null(parameters$algorithm_parameters$depth)) {    
    depth = parameters$algorithm_parameters$depth
    if (depth <=0) stop("Search depth must be a positive integer", call. = FALSE)
  } else {
    # Default value
    depth = 2
  }

  parameters$algorithm_parameters$depth = depth

  if (!is.null(parameters$algorithm_parameters$width)) { 
    width = parameters$algorithm_parameters$width
    if (width <=0) stop("Search depth must be a positive integer", call. = FALSE)
  } else {
    # Default value
    width = 2
  }

  parameters$algorithm_parameters$width = width
  
  if (!is.null(parameters$algorithm_parameters$min_subgroup_size)) { 
    min_subgroup_size = parameters$algorithm_parameters$min_subgroup_size
    if (min_subgroup_size <=0) stop("Minimum number of patients must be a positive integer", call. = FALSE)
  } else {
    # Default value
    min_subgroup_size = 30
  }

  # Error check for criterion_type
  criterion_type = parameters$algorithm_parameters$criterion_type
  if (is.null(criterion_type) || is.na(criterion_type)) {
    criterion_type = 1
  } else {
    # Supported splitting criteria:
    # criterion_type = 1: Differential criterion 
    # criterion_type = 2: Maximum criterion 
    # criterion_type = 3: Directional criterion 
    if (!(criterion_type %in% c(1, 2, 3))) stop("The splitting criterion type must be 1, 2 or 3", call. = FALSE)
  }

  if (!is.null(parameters$algorithm_parameters$ncores)) { 
    ncores = parameters$algorithm_parameters$ncores
  } else {
    ncores = 1
  }

  # Account for the maximum number of cores
  ncores = min(ncores, detectCores())
  
  # Error check for gamma
  gamma_standardized = rep(1, depth)

  if (!is.null(parameters$algorithm_parameters$gamma)) { 
    gamma = parameters$algorithm_parameters$gamma
    if (length(gamma) != depth) stop("The length of the vector of complexity parameters must be equal to the depth parameter", call. = FALSE)
    
    for (i in 1:depth) {
      if (!is.na(gamma[i])) {
        if (gamma[i] <= 0) stop("Complexity parameters must be positive", call. = FALSE)
      }
      # Parameter standardization
      if (is.na(gamma[i])) {
        gamma_standardized[i] = -1
      } else {
        gamma_standardized[i] = gamma[i]    
      }
    } 
  } 

  parameters$algorithm_parameters$gamma_standardized = gamma_standardized

  # Error check for pvalue_max
  if (!is.null(parameters$algorithm_parameters$pvalue_max)) {
    
    pvalue_max = parameters$algorithm_parameters$pvalue_max
    # ensure that pvalue_max is a number between 0 and 1
    pvalue_max =min(max(pvalue_max,0),1)
    
  } else {
    # Do not apply restrictions on subgroup p-values
    pvalue_max = 1
  }
  
  # Error check for local_mult_adj
  if (!is.null(parameters$algorithm_parameters$local_mult_adj)) {
    
    local_mult_adj = parameters$algorithm_parameters$local_mult_adj
    # ensure that local_mult_adj is a number between 0 and 1
    local_mult_adj =min(max(local_mult_adj,0),1)
    
  } else {
    # Apply local multiplicity adjustment by default
    local_mult_adj = 1
  }
    
  # Error check for n_perms_mult_adjust
  if (!is.null(parameters$algorithm_parameters$n_perms_mult_adjust)) {
    
    n_perms_mult_adjust = round(parameters$algorithm_parameters$n_perms_mult_adjust)
    if (n_perms_mult_adjust <= 0) stop("Number of permutations must be a positive integer", call. = FALSE)
    
    
  } else {
    # Default value
    n_perms_mult_adjust = 1000
  }

  parameters$algorithm_parameters$n_perms_mult_adjust = n_perms_mult_adjust
  
  # Error check for perm_type
  if (!is.null(parameters$algorithm_parameters$perm_type)) {

    if (!(perm_type %in% 1)) stop("The permutation type must be 1", call. = FALSE)

  } else {
    # Default value
    perm_type = 1
  }


  subgroup_search_algorithm = parameters$algorithm_parameters$subgroup_search_algorithm
  
  # Error check for subgroup_search_algorithm
  if (is.null(subgroup_search_algorithm)) {
    subgroup_search_algorithm = "SIDES procedure"
    subgroup_search_algorithm_standardized = 1
  } else {
    
    if (!(subgroup_search_algorithm %in% c("SIDES procedure", "Fixed SIDEScreen procedure", "Adaptive SIDEScreen procedure"))) stop("The subgroup search procedure type must be SIDES procedure, Fixed SIDEScreen procedure or Adaptive SIDEScreen procedure", call. = FALSE)
    # Parameter standardization   
    if (subgroup_search_algorithm == "SIDES procedure") subgroup_search_algorithm_standardized = 1
    if (subgroup_search_algorithm == "Fixed SIDEScreen procedure") subgroup_search_algorithm_standardized = 2
    if (subgroup_search_algorithm == "Adaptive SIDEScreen procedure") subgroup_search_algorithm_standardized = 3
  }
  
  # Error check for n_top_biomarkers
  if (!is.null(parameters$algorithm_parameters$n_top_biomarkers)) {
    
    # ensure that n_top_biomarkers is a positive integer
    # ensure that n_top_biomarkers is no greater than the total number of biomarkers
    n_top_biomarkers = round(parameters$algorithm_parameters$n_top_biomarkers)
    n_top_biomarkers = min(max(1,n_top_biomarkers),n_biomarkers)  
  } else {
    
    n_top_biomarkers = min(3,n_biomarkers)
  }

  parameters$algorithm_parameters$n_top_biomarkers = n_top_biomarkers
  
  # Error check for multiplier
  if (!is.null(parameters$algorithm_parameters$multiplier)) {
    
    if (!is.numeric(parameters$algorithm_parameters$multiplier)){ 
      multiplier =1
    } else {multiplier = parameters$algorithm_parameters$multiplier}
    
    
  } else {
    
    multiplier = 1
  }
  
  parameters$algorithm_parameters$multiplier = multiplier


  # Error check for n_perms_vi_score
  if (!is.null(parameters$algorithm_parameters$n_perms_vi_score)) {
    
    # ensure that n_perms_vi_score is a positive integer
    n_perms_vi_score = round(parameters$algorithm_parameters$n_perms_vi_score)
    if (n_perms_vi_score <= 0) n_perms_vi_score = 100
    
  } else {
    # Default value
    n_perms_vi_score = 100
  }

  parameters$algorithm_parameters$n_perms_vi_score = n_perms_vi_score

  # Error check for nperc
  if (!is.null(parameters$algorithm_parameters$nperc)) {
    
    # ensure that nperc is a positive integer
    nperc = round(parameters$algorithm_parameters$nperc)
    if (nperc < 0) nperc = 20
    
  } else {
    # Default value
    nperc = 20
  }

  parameters$algorithm_parameters$nperc = nperc

  # Error check for depth_second_stage
  if (!is.null(parameters$algorithm_parameters$depth_second_stage)) {
    
    # ensure that depth_second_stage is a positive integer
    depth_second_stage = round(parameters$algorithm_parameters$depth_second_stage)
    if (depth_second_stage <=0) depth_second_stage = round(depth)
    
  } else {
    # Default value if the parameter is not specified
    depth_second_stage = round(depth)
  }
  
  # Error check for width_second_stage
  if (!is.null(parameters$algorithm_parameters$width_second_stage)) {
    
    # ensure that width_second_stage is a positive integer
    width_second_stage = round(parameters$algorithm_parameters$width_second_stage)
    if (width_second_stage <= 0) width_second_stage = round(parameters$algorithm_parameters$width)
    
  } else {
    # Default value if the parameter is not specified
    width_second_stage = round(width)
  }
  
  # Error check for gamma_second_stage
  if (!is.null(parameters$algorithm_parameters$gamma_second_stage)) {
    
    # Error check for gamma
    gamma_second_stage = parameters$algorithm_parameters$gamma_second_stage
    
    if (length(gamma_second_stage) != depth_second_stage) stop("The length of the vector of complexity parameters in the second stage must be equal to the depth parameter in the second stage", call. = FALSE)
    
    gamma_second_stage_standardized = rep(0, depth_second_stage)
    
    for (i in 1:depth_second_stage) {
      if (!is.na(gamma_second_stage[i])) {
        if (gamma_second_stage[i] <= 0) stop("Complexity parameters in the second stage must be positive", call. = FALSE)
      }
      # Parameter standardization
      if (is.na(gamma_second_stage[i])) {
        gamma_second_stage_standardized[i] = -1
      } else {
        gamma_second_stage_standardized[i] = gamma_second_stage[i]    
      }
    } 
    
  } else {
    # Default value if the parameter is not specified
    gamma_second_stage_standardized = gamma_standardized
  }
  
  # Error check for pvalue_max_second_stage
  if (!is.null(parameters$algorithm_parameters$pvalue_max_second_stage)) {
    
    # ensure that pvalue_max_second_stage is between 0 and 1
    pvalue_max_second_stage =min(max(0,parameters$algorithm_parameters$pvalue_max_second_stage),1)
    
  } else {
    # Default value if the parameter is not specified
    pvalue_max_second_stage = pvalue_max
  }
  
  # Error check for random_seed
  if (!is.null(parameters$algorithm_parameters$random_seed)) {
    
    # ensure that random_seed is an integer
    random_seed = round(parameters$algorithm_parameters$random_seed)
    
  } else {

    # Default value if the parameter is not specified
    random_seed = 49291
  }

  parameters$algorithm_parameters$random_seed = random_seed

  # Precision parameter used in model fitting
  if (!is.null(parameters$algorithm_parameters$precision)) {
    
    precision = parameters$algorithm_parameters$precision
    
  } else {

    # Default value if the parameter is not specified
    precision = 0.001
  }

  parameters$algorithm_parameters$precision = precision

  # Number of interations used in model fitting
  if (!is.null(parameters$algorithm_parameters$max_iter)) {
    
    max_iter = parameters$algorithm_parameters$max_iter
    
  } else {

    # Default value if the parameter is not specified
    max_iter = 10
  }

  parameters$algorithm_parameters$max_iter = max_iter

  parlist<-list(  
    min_subgroup_size, 
    outcome_variable_direction, 
    outcome_type_standardized,
    analysis_method_standardized, 
    biomarker_types_standardized,
    biomarker_levels, 
    criterion_type, 
    width, 
    depth, 
    gamma_standardized, 
    pvalue_max,
    local_mult_adj,
    n_perms_mult_adjust,
    perm_type,
    subgroup_search_algorithm_standardized,
    n_top_biomarkers,
    multiplier,
    n_perms_vi_score,
    nperc,
    width_second_stage,
    depth_second_stage,
    gamma_second_stage_standardized,
    pvalue_max_second_stage,
    random_seed,
    precision,
    max_iter
  )
  
  names(parlist)<-c(  
    "min_subgroup_size", 
    "outcome_variable_direction", 
    "outcome_type_standardized", 
    "analysis_method",
    "biomarker_types_standardized", 
    "biomarker_levels", 
    "criterion_type", 
    "width", 
    "depth", 
    "gamma_standardized", 
    "pvalue_max",
    "local_mult_adj",
    "n_perms_mult_adjust",
    "perm_type",
    "subgroup_search_algorithm_standardized",
    "n_top_biomarkers",
    "multiplier",
    "n_perms_vi_score",
    "nperc",
    "width_second_stage",
    "depth_second_stage",
    "gamma_second_stage_standardized",
    "pvalue_max_second_stage",
    "random_seed",
    "precision",
    "max_iter"
  )
  
  return(list(parlist,df,coded.all,covariate_object,parameters))


}


CreateProjectFile<-function(datafname, parameters, project_filename="sides_project.xml", exdata_filename="exdata_std.csv") {  
  # creates xml file sides_project.xml
  
  rescode<-0

  # Compute standardized SIDES parameters
  tryCatch(sidesParamStd<-StandardizedSIDESParameters(parameters), 
           error=function(e) {message(e); rescode<<-1})
  if (rescode!=0) stop("Exiting due to errors", call. = FALSE)
  
  param_std<-sidesParamStd[[1]]
  df<-sidesParamStd[[2]]
  dic<-sidesParamStd[[3]]

  # Create the list of covariate to be passed to the SIDES library
  covariate_object = sidesParamStd[[4]]

  # List of parameters with default values
  parameters = sidesParamStd[[5]]
  
  # stdfilename="exdata_std.csv"
  write.csv(df,exdata_filename,row.names=FALSE, na=".", quote=FALSE)
  
  biomarker_types=parameters$data_set_parameters$biomarker_types
  biomarker_names=parameters$data_set_parameters$biomarker_names
  biomarker_levels=parameters$data_set_parameters$biomarker_levels
  dataset=parameters$data_set_parameters$data_set
  treatment_variable_name = parameters$data_set_parameters$treatment_variable_name
  treatment_variable_control_value = parameters$data_set_parameters$treatment_variable_control_value
  outcome_censor_name = parameters$data_set_parameters$outcome_censor_name
  outcome_censor_value = parameters$data_set_parameters$outcome_censor_value
  outcome_variable_name = parameters$data_set_parameters$outcome_variable_name
  outcome_variable_type = parameters$data_set_parameters$outcome_variable_type
  analysis_method = parameters$data_set_parameters$analysis_method
  nperc = parameters$algorithm_parameters$nperc
  outcome_variable_direction= parameters$data_set_parameters$outcome_variable_direction
  
  criterion_type = param_std$criterion_type
  
  vnames=colnames(dataset)
  ncov=length(biomarker_names)
  
  myxml<-rep("",1000)
  myxml[1]=  "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
  myxml[2]= "<head>"
  myxml[3]=  paste0("<data skipfirstrow=",AQ(1)," ncol=",AQ(ncol(df))," nrow=", AQ(nrow(df)),
                    " file=",AQ(datafname)," stddata=",AQ(exdata_filename)," />")
  myxml[4]= "<structure>"
  myxml[5]= "<biomarkers>"
  i.xml=5
  
  for (i in 1:ncov) {
    myxml[i.xml+i]=paste0("<biomarker name=",AQ(biomarker_names[i])," column=",AQ(match(biomarker_names[i],vnames))," numeric=",AQ(1*(biomarker_types[i]=="numeric"))," level=",AQ(biomarker_levels[i]),"/>")
  }
  i.xml=i.xml+ncov
  
  myxml[i.xml+1]= "</biomarkers>"
  
  if (outcome_variable_type=="binary") {
    myxml[i.xml+2]=paste0("<outcome name=",AQ(outcome_variable_name)," column=",AQ(match(outcome_variable_name,vnames))," type=",AQ(outcome_variable_type)," direction=",AQ("larger")," desirable=",AQ(1)," />")
    i.xml=i.xml+2         
  } 
  if (outcome_variable_type=="survival") {
    if (outcome_variable_direction==1) {
      myxml[i.xml+2]=paste0("<outcome name=",AQ(outcome_variable_name)," column=",AQ(match(outcome_variable_name,vnames))," type=",AQ("time")," direction=",AQ("larger")," desirable=",AQ("")," />")
    } else {
      myxml[i.xml+2]=paste0("<outcome name=",AQ(outcome_variable_name)," column=",AQ(match(outcome_variable_name,vnames))," type=",AQ("time")," direction=",AQ("smaller")," desirable=",AQ("")," />")          
    }
    myxml[i.xml+3]=paste0("<outcome_censor name=",AQ(outcome_censor_name)," column=",AQ(match(outcome_censor_name,vnames)),"  value=",AQ(outcome_censor_value),"/>")
    i.xml=i.xml+3
  } 
  if (outcome_variable_type=="continuous") {
    if (outcome_variable_direction==1) {
      myxml[i.xml+2]=paste0("<outcome name=",AQ(outcome_variable_name)," column=",AQ(match(outcome_variable_name,vnames))," type=",AQ(outcome_variable_type)," direction=",AQ("larger")," desirable=",AQ("")," />")  
    } else {
      myxml[i.xml+2]=paste0("<outcome name=",AQ(outcome_variable_name)," column=",AQ(match(outcome_variable_name,vnames))," type=",AQ(outcome_variable_type)," direction=",AQ("smaller")," desirable=",AQ("")," />")
    }
    i.xml=i.xml+2    
  }
  if (outcome_variable_type=="ancova") {
    if (outcome_variable_direction==1) {
      myxml[i.xml+2]=paste0("<outcome name=",AQ(outcome_variable_name)," column=",AQ(match(outcome_variable_name,vnames))," type=",AQ(outcome_variable_type)," direction=",AQ("larger")," desirable=",AQ("")," />")  
    } else {
      myxml[i.xml+2]=paste0("<outcome name=",AQ(outcome_variable_name)," column=",AQ(match(outcome_variable_name,vnames))," type=",AQ(outcome_variable_type)," direction=",AQ("smaller")," desirable=",AQ("")," />")
    }
    i.xml=i.xml+2    
  }
  myxml[i.xml+1]=paste0("<treatment name=",AQ(treatment_variable_name)," column=",AQ(match(treatment_variable_name,vnames))," value=",AQ(treatment_variable_control_value)," />")
  myxml[i.xml+2]="</structure>"
  myxml[i.xml+3]="<parameters>"
  myxml[i.xml+4]=paste0("<criterion_type value=",AQ(criterion_type)," />")
  myxml[i.xml+5]=paste0("<width value=",AQ(param_std$width)," />")
  myxml[i.xml+6]=paste0("<depth value=",AQ(param_std$depth)," />")
  myxml[i.xml+7]="<complexity_control>"
  i.xml=i.xml+7
  for (i in 1:param_std$depth) {
    myxml[i.xml+i]=paste0("<gamma id=",AQ(i)," value=",AQ(param_std$gamma_standardized[i])," />")
  }
  i.xml=i.xml+param_std$depth
  myxml[i.xml+1]="</complexity_control>"
  myxml[i.xml+2]=paste0("<min_subgroup_size value=",AQ(param_std$min_subgroup_size)," />")
  myxml[i.xml+3]=paste0("<subgroup_search_algorithm value=",AQ(param_std$subgroup_search_algorithm_standardized)," />")
  myxml[i.xml+4]=paste0("<local_mult_adj value=",AQ(param_std$local_mult_adj)," />")
  myxml[i.xml+5]=paste0("<n_perms_mult_adjust value=",AQ(param_std$n_perms_mult_adjust)," />")
  myxml[i.xml+6]=paste0("<perm_type value=",AQ(param_std$perm_type)," />")
  myxml[i.xml+7]=paste0("<n_top_biomarkers value=",AQ(param_std$n_top_biomarkers)," />")
  myxml[i.xml+8]=paste0("<n_perms_vi_score value=",AQ(param_std$n_perms_vi_score)," />")
  myxml[i.xml+9]=paste0("<multiplier value=",AQ(param_std$multiplier)," />")
  myxml[i.xml+10]=paste0("<analysis_method value=",AQ(param_std$analysis_method)," />")
  myxml[i.xml+11]=paste0("<nperc value=",AQ(param_std$nperc)," />")
  myxml[i.xml+12]=paste0("<precision value=",AQ(param_std$precision)," />")
  myxml[i.xml+13]=paste0("<max_iter value=",AQ(param_std$max_iter)," />")
  myxml[i.xml+14]="</parameters>"
  myxml[i.xml+15]="</head>"
  length(myxml)=i.xml+15
  # write to the file
  fileConn<-file(project_filename)  
  writeLines(myxml, fileConn)
  close(fileConn)
  
  return(list(colnames(df),dic, covariate_object, parameters))
}

read_SIDESxml<-function(fname="sides_output.xml") {
  xml2ver<-packageVersion("xml2") 
  if (substr(xml2ver,1,1)>=1) {ver.offset<-0} else {ver.offset<-1}
  myxml<-as_list(read_xml(fname))
  
  ## extract subgroups
  
  #sub.fields<-c("subgroup","bio1.name","bio1.sign","bio1.value","bio2.name","bio2.sign","bio2.value","bio3.name","bio3.sign","bio3.value", "subgroup_size","splitting_criterion","splitting_criterion_log_p_value","p_value","adjusted_p_value")
  sub.fields<-c("subgroup","bio1.name","bio1.sign","bio1.value","bio2.name","bio2.sign","bio2.value","bio3.name","bio3.sign","bio3.value", "size","size_control","size_treatment", "prom_estimate","prom_sderr","prom_sd","p_value","adjusted_p_value")

  nfields<-length(sub.fields)
  subgroup.frame<-as.data.frame(matrix(NA,1,nfields))
  colnames(subgroup.frame)<-sub.fields
  empty.rec<-subgroup.frame
  
  myxml = myxml$head

  nsub<-length(myxml$subgroups) 
  i.sub<-0

  for (i in 1:nsub)
  {
    if (length(myxml$subgroups[[i]])>1) {
      i.sub<-i.sub+1
      subgroup.frame$subgroup[i.sub]<-attr(myxml$subgroups[[i]]$definition$component,"description")
      if (substr(xml2ver,1,1)>=1) {nbio<-length(myxml$subgroups[[i]]$definition)} else {nbio<-(length(myxml$subgroups[[i]]$definition)-1)/2} 
      
      if (nbio>0 & attr(myxml$subgroups[[i]]$definition[[1+ver.offset]],"biomarker") !="") {    
        for (j in 1:nbio) {
          subgroup.frame[i.sub,1+(j-1)*3+1]<-paste0("biomarker",attr(myxml$subgroups[[i]]$definition[[(1+ver.offset)*j]],"biomarker"))

          if (!is.null(attr(myxml$subgroups[[i]]$definition[[(1+ver.offset)*j]],"sign"))) {
             subgroup.frame[i.sub,1+(j-1)*3+2]<-attr(myxml$subgroups[[i]]$definition[[(1+ver.offset)*j]],"sign")
          }
          if (!is.null(attr(myxml$subgroups[[i]]$definition[[(1+ver.offset)*j]],"value"))) {
               subgroup.frame[i.sub,1+(j-1)*3+3]<-attr(myxml$subgroups[[i]]$definition[[(1+ver.offset)*j]],"value")
           }
        } 
      } 
      # read in parameters
      subgroup.frame$size_control[i.sub]<-attr(myxml$subgroups[[i]]$parameters,"size_control")
      subgroup.frame$size_treatment[i.sub]<-attr(myxml$subgroups[[i]]$parameters,"size_treatment")
      subgroup.frame$size[i.sub]<-as.numeric(subgroup.frame$size_control[i.sub])+as.numeric(subgroup.frame$size_treatment[i.sub])
      subgroup.frame$prom_estimate[i.sub]<-attr(myxml$subgroups[[i]]$parameters,"prom_estimate")
      subgroup.frame$prom_sderr[i.sub]<-attr(myxml$subgroups[[i]]$parameters,"prom_sderr")
      subgroup.frame$prom_sd[i.sub]<-attr(myxml$subgroups[[i]]$parameters,"prom_sd")
      subgroup.frame$p_value[i.sub]<-attr(myxml$subgroups[[i]]$parameters,"p_value")
      subgroup.frame$adjusted_p_value[i.sub]<-attr(myxml$subgroups[[i]]$parameters,"adjusted_p_value")
      
      subgroup.frame<-rbind(subgroup.frame,empty.rec)
      
    }
  }
  subgroup.frame<-subgroup.frame[-nrow(subgroup.frame),]

  
  ## extract variable importance
  
  vi.fields<-c("biomarker","vi")
  nfields<-length(vi.fields)
  vi.frame<-as.data.frame(matrix(NA,1,nfields))
  colnames(vi.frame)<-vi.fields
  empty.rec<-vi.frame
  
  if (length(myxml$vi_scores) >= 1) {
    i.vi<-0
    for (i in 1:length(myxml$vi_scores)) {
      if (is.list(myxml$vi_scores[[i]])) {
        i.vi<-i.vi+1
        vi.frame$biomarker[i.vi]<- paste0("biomarker",attr(myxml$vi_scores[[i]],"biomarker"))
        vi.frame$vi[i.vi]<- attr(myxml$vi_scores[[i]],"value")
        vi.frame<-rbind(vi.frame,empty.rec)
      }
    }
  }
  
  vi.frame<-vi.frame[-nrow(vi.frame),]

# extract vi threshold if present
  if (!is.null(myxml$vi_threshold)) {
    vi_threshold=list()
    vi_threshold$mean<-attr(myxml$vi_threshold,"mean")
    vi_threshold$sd<-attr(myxml$vi_threshold,"sd")
    vi_threshold$threshold<-attr(myxml$vi_threshold,"threshold")
  } else {
    vi_threshold<-NULL
  } 

# extract the Type I error rate 
 if (!is.null(myxml$error_rate)) {
   error_rate<-attr(myxml$error_rate,"value")
 } else error_rate<-NULL

  return(list(subgroup.frame,vi.frame,vi_threshold,error_rate))
 
}  


recode_subgroups <-function(sides.res,codes) {
  # recode subgroup table
  vnames.std<-codes[[1]]
  dic<-codes[[2]]
  
  subgr.table <-sides.res[[1]]
  
  for (i in 1:nrow(subgr.table)) {
    signature="" 
    for (j in 1:3) {
      b.col=1+(j-1)*3+1
      b.sign =1+(j-1)*3+2
      b.val=1+(j-1)*3+3
      if (!is.na(subgr.table[i,b.col])) {
        b.ind=as.numeric(substr(subgr.table[i,b.col],10,nchar(subgr.table[i,b.col]))) 
        subgr.table[i,b.col]<-vnames.std[b.ind]
      }
      if (!is.na(subgr.table[i,b.sign]) & !is.null(dic)) {
        which=match(vnames.std[b.ind],names(dic[[2]]))
        # replace value(s) from the dictionary
        if (!is.na(which)) {
          dictionary=as.data.frame(dic[[2]][which])
          spt=strsplit(subgr.table[i,b.val]," ")[[1]]
          value=NULL
          for (iw in 1:length(spt)) {
            w<-dictionary[as.numeric(spt[iw]),1]
            if (is.null(value)) value=paste(value,w,sep="") else value=paste(value,w,sep=",")
          }
          subgr.table[i,b.val]<-value
        } 
      }
      if (!is.na(subgr.table[i,b.col])) {
        signature <- paste0(signature,subgr.table[i,b.col],subgr.table[i,b.sign],subgr.table[i,b.val], ":")
      }
    }
    if (signature !="") {subgr.table$subgroup[i]<-signature}
  }
  
  #replace names in variable importance table 
  
  vi.table <-sides.res[[2]]
  for (i in 1: nrow(vi.table))
  {b.ind=as.numeric(substr(vi.table[i,1],10,nchar(vi.table[i,1])))
   vi.table[i,1]=vnames.std[b.ind]
  }
  vi.table[,2]<-as.numeric(vi.table[,2])
  res<-list(subgr.table,vi.table,sides.res[[3]],sides.res[[4]])
  return(res)
}

SubgroupSearch <- function(parameters) { 

    # Start time 
    start_time = Sys.time()

    # Temporary files
    project_filename = tempfile(pattern = "sides_project", fileext = ".xml")
    exdata_filename  = tempfile(pattern = "exdata_std", fileext = ".csv")
    output_filename  = tempfile(pattern = "sides_output",  fileext = ".xml")

    # Save the number of permutations  
    n_perms_mult_adjust = parameters$algorithm_parameters$n_perms_mult_adjust

    # Run SIDES without a multiplicity adjustment
    parameters$algorithm_parameters$n_perms_mult_adjust = 1  

    # Random seed
    random_seed = parameters$algorithm_parameters$random_seed

    # Create the project file
    datafname = "placeholder.csv"
    rescode = 0
    tryCatch(codes<-CreateProjectFile(datafname, parameters, project_filename, exdata_filename), 
            error=function(e) {message(e); rescode<<-1})
    if (rescode!=0) stop("Exiting due to errors", call. = FALSE)

    # Create the list of covariate to be passed to the SIDES library
    covariate_object = codes[[3]]

    # List of parameters with default values
    parameters = codes[[4]]

    # SIDES(ancova)
    SIDES(covariate_object$outcome, 
          covariate_object$censor,
          covariate_object$treatment, 
          covariate_object$cont_covariates, 
          covariate_object$class_covariates, 
          covariate_object$n_cont_covariates, 
          covariate_object$n_class_covariates,
          project_filename,
          output_filename)    

    sides.res<- read_SIDESxml(output_filename) 

    sides.res<- recode_subgroups(sides.res,codes) 
    subgr<-sides.res[[1]]
    subgr.short<-subgr[,-seq(2,10)]
    n_subgroups = dim(subgr.short)[1] - 1
    vi<-sides.res[[2]]
    ord<-order(vi[,2],decreasing=TRUE)
    vi<-vi[ord,]
    if (!is.null(sides.res[[3]])) vi.thresh<-as.numeric(sides.res[[3]]) else vi.thresh<-NA

    error_rate = NA

    if (!is.null(parameters$algorithm_parameters$ncores)) { 
      ncores = parameters$algorithm_parameters$ncores
    } else {
      ncores = 1
    }

    # Account for the maximum number of cores
    ncores = min(ncores, detectCores())

    parameters$algorithm_parameters$n_perms_mult_adjust = ceiling(n_perms_mult_adjust / ncores)  

    # Run SIDES with a multiplicity adjustment if at least one subgroup is found
    if (n_subgroups >= 1) {
        counter = NULL

        # Create the project file
        rescode<-0
        tryCatch(codes<-CreateProjectFile(datafname, parameters, project_filename, exdata_filename), 
                error=function(e) {message(e); rescode<<-1})
        if (rescode!=0) stop("Exiting due to errors", call. = FALSE)

        if (parameters$algorithm_parameters$subgroup_search_algorithm == "SIDES procedure") {  
          
          if (ncores==1) {
            null_pvalue = SIDESAdjP(covariate_object$outcome, 
                        covariate_object$censor,
                        covariate_object$treatment, 
                        covariate_object$cont_covariates, 
                        covariate_object$class_covariates, 
                        covariate_object$n_cont_covariates, 
                        covariate_object$n_class_covariates,
                        random_seed + 100,
                        project_filename)
          } else {
            cl = makeCluster(ncores)
            registerDoParallel(cl)
            null_pvalue = foreach(counter = 1:ncores, .combine=c, .packages = c("rsides")) %dorng% { 
              SIDESAdjP(covariate_object$outcome, 
                        covariate_object$censor,
                        covariate_object$treatment, 
                        covariate_object$cont_covariates, 
                        covariate_object$class_covariates, 
                        covariate_object$n_cont_covariates, 
                        covariate_object$n_class_covariates,
                        random_seed + 100 * counter,
                        project_filename)
            }
            stopCluster(cl)
          }
        }

        if (parameters$algorithm_parameters$subgroup_search_algorithm == "Fixed SIDEScreen procedure") {  
          
          if (ncores==1) {
            null_pvalue = FixedSIDEScreenAdjP(covariate_object$outcome, 
                        covariate_object$censor,
                        covariate_object$treatment, 
                        covariate_object$cont_covariates, 
                        covariate_object$class_covariates, 
                        covariate_object$n_cont_covariates, 
                        covariate_object$n_class_covariates,
                        random_seed + 100,
                        project_filename)
          } else {
            cl = makeCluster(ncores)
            registerDoParallel(cl)
            null_pvalue = foreach(counter = 1:ncores, .combine=c, .packages = c("rsides")) %dorng% { 
              FixedSIDEScreenAdjP(covariate_object$outcome, 
                        covariate_object$censor,
                        covariate_object$treatment, 
                        covariate_object$cont_covariates, 
                        covariate_object$class_covariates, 
                        covariate_object$n_cont_covariates, 
                        covariate_object$n_class_covariates,
                        random_seed + 100 * counter,
                        project_filename)
            }
            stopCluster(cl)
          }
        }

        if (parameters$algorithm_parameters$subgroup_search_algorithm == "Adaptive SIDEScreen procedure") {  
          
          if (ncores==1) {
            null_pvalue = AdaptiveSIDEScreenAdjP(covariate_object$outcome, 
                        covariate_object$censor,
                        covariate_object$treatment, 
                        covariate_object$cont_covariates, 
                        covariate_object$class_covariates, 
                        covariate_object$n_cont_covariates, 
                        covariate_object$n_class_covariates,
                        random_seed + 100,
                        project_filename)
          } else {
            cl = makeCluster(ncores)
            registerDoParallel(cl)
            null_pvalue = foreach(counter = 1:ncores, .combine=c, .packages = c("rsides")) %dorng% { 
              AdaptiveSIDEScreenAdjP(covariate_object$outcome, 
                        covariate_object$censor,
                        covariate_object$treatment, 
                        covariate_object$cont_covariates, 
                        covariate_object$class_covariates, 
                        covariate_object$n_cont_covariates, 
                        covariate_object$n_class_covariates,
                        random_seed + 100 * counter,
                        project_filename)
            }
            stopCluster(cl)
          }
       }

        null_pvalue = as.numeric(null_pvalue)
        p_value = as.numeric(subgr.short[["p_value"]][2:(1+n_subgroups)])
        adj_pvalue = rep(1, 1 + n_subgroups)
        adj_pvalue[1] = NA
        for (i in 1:n_subgroups) adj_pvalue[i + 1] = mean(null_pvalue <= p_value[i])
        subgr.short[["adjusted_p_value"]] = adj_pvalue
        error_rate = mean(null_pvalue <= 0.025) 

print(p_value)
print(null_pvalue)

    }

    end_time = Sys.time()

    elapsed_time = difftime(end_time, start_time, units = "secs")
    elapsed_time = round(as.numeric(elapsed_time), 1)

    # Delete temporary files
    temp_files = c("exdata_std.csv", "sides_output.xml", "sides_project.xml")
    for (i in 1:length(temp_files)) {
      if (file.exists(temp_files[i])) file.remove(temp_files[i])
    }
    
    # Save the results and assumptions
    results = list()

    patient_subgroups = list(subgr.short, vi, vi.thresh, error_rate, elapsed_time)
    names(patient_subgroups)<-c("Subgroups","Variable Importance","VI Threshold","Type I error rate","Elapsed time (seconds)")

    results$patient_subgroups = patient_subgroups
    results$parameters = parameters

    class(results) = "SubgroupSearchResults"

    return(results)
}

print.SubgroupSearchResults = function (x, ...) {

    results = x

    endpoint_parameters = results$parameters$endpoint_parameters
    subgroup_search_parameters = results$parameters
    patient_subgroups = results$patient_subgroups 

    cat("***************************************\n\n")
    cat("Subgroup search results\n\n")
    cat("***************************************\n\n")

    x = results$patient_subgroups$Subgroups
    if (!is.na(x[1, 9])) {
      if (x[1, 9] == -1) x[1, 9] = NA
    }
    n_subgroups = dim(x)[1]
    y = as.data.frame(matrix(0, n_subgroups, 4))

    if (endpoint_parameters$analysis_method == "T-test") {

      for (i in 1:n_subgroups) {
        y[i, 1] = x[i, 1] 
        y[i, 2] = paste0(x[i, 2], " (", x[i, 3], ", ", x[i, 4], ")")
        y[i, 3] = paste0(round(as.numeric(x[i, 5]), 2), " (", round(as.numeric(x[i, 7]), 2), ")")
        y[i, 4] = paste0(round(as.numeric(x[i, 8]), 4), " (", round(as.numeric(x[i, 9]), 4), ")")
      }
      colnames(y) = c("Subgroup", "Total size (Control, Treatment)", "Estimate (SD)", "P-value (Adjusted p-value)")   
      print(y, row.names = FALSE)

    } 

    if (endpoint_parameters$analysis_method %in% c("Z-test for proportions", "Log-rank test")) {

      for (i in 1:n_subgroups) {
        y[i, 1] = x[i, 1] 
        y[i, 2] = paste0(x[i, 2], " (", x[i, 3], ", ", x[i, 4], ")")
        y[i, 3] = round(as.numeric(x[i, 5]), 3)
        y[i, 4] = paste0(round(as.numeric(x[i, 8]), 4), " (", round(as.numeric(x[i, 9]), 4), ")")
      }
      colnames(y) = c("Subgroup", "Total size (Control, Treatment)", "Estimate", "P-value (Adjusted p-value)")   
      print(y, row.names = FALSE)

    } 

    if (endpoint_parameters$analysis_method %in% c("ANCOVA", "Logistic regression", "Cox regression")) {

      for (i in 1:n_subgroups) {
        y[i, 1] = x[i, 1] 
        y[i, 2] = paste0(x[i, 2], " (", x[i, 3], ", ", x[i, 4], ")")
        y[i, 3] = paste0(round(as.numeric(x[i, 5]), 2), " (", round(as.numeric(x[i, 6]), 2), ")")
        y[i, 4] = paste0(round(as.numeric(x[i, 8]), 4), " (", round(as.numeric(x[i, 9]), 4), ")")
      }
      colnames(y) = c("Subgroup", "Total size (Control, Treatment)", "Estimate (SE)", "P-value (Adjusted p-value)")   
      print(y, row.names = FALSE)

    }

    vi = results$patient_subgroups$`Variable Importance`
    vi[, 2] = round(vi[, 2], 3)
    vi = vi[vi[, 2] > 0, ]

    if(dim(vi)[1] > 0) {
  
      cat("\n***************************************\n\n")
      cat("Variable importance\n\n")
      cat("***************************************\n\n")

      colnames(vi) = c("Biomarker", "Variable importance score")
      print(vi, row.names = FALSE)

    }

    cat("\n***************************************\n\n")
    cat("Algorithm's parameters\n\n")
    cat("***************************************\n\n")

    y = as.data.frame(matrix(0, 2, 2))
    y[1, 1] = "Elapsed time (sec)"
    y[2, 1] = "Type I error rate"
    y[1, 2] = as.character(results$patient_subgroups$`Elapsed time (seconds)`) 
    if(!is.na(results$patient_subgroups$`Type I error rate`)) y[2, 2] = as.character(results$patient_subgroups$`Type I error rate`) else  y[2, 2] = "NA"

    colnames(y) = c("Parameter", "Value")   
    print(y, row.names = FALSE)

}

SaveReport = function(report, report_title) {

  # Create a docx object
  doc = officer::read_docx(system.file(package = "rsides", "template/report_template.docx"))

  dim_doc = officer::docx_dim(doc)

  # Report's title
  doc = officer::set_doc_properties(doc, title = report_title)
  doc = officer::body_add_par(doc, value = report_title, style = "heading 1")

  # Text formatting
  my.text.format = officer::fp_text(font.size = 12, font.family = "Arial")

  # Table formatting
  header.cellProperties = officer::fp_cell(border.left = officer::fp_border(width = 0), border.right = officer::fp_border(width = 0), border.bottom = officer::fp_border(width = 2), border.top = officer::fp_border(width = 2), background.color = "#eeeeee")
  data.cellProperties = officer::fp_cell(border.left = officer::fp_border(width = 0), border.right = officer::fp_border(width = 0), border.bottom = officer::fp_border(width = 0), border.top = officer::fp_border(width = 0))

  header.textProperties = officer::fp_text(font.size = 12, bold = TRUE, font.family = "Arial")
  data.textProperties = officer::fp_text(font.size = 12, font.family = "Arial")

  thick_border = fp_border(color = "black", width = 2)

  leftPar = officer::fp_par(text.align = "left")
  rightPar = officer::fp_par(text.align = "right")
  centerPar = officer::fp_par(text.align = "center")

  # Number of sections in the report (the report's title is not counted)
  n_sections = length(report) 

  # Loop over the sections in the report
  for(section_index in 1:n_sections) {

      # Determine the item's type (text by default)
      type = report[[section_index]]$type

      # Determine the item's label 
      label = report[[section_index]]$label

      # Determine the item's label 
      footnote = report[[section_index]]$footnote

      # Determine the item's value 
      value = report[[section_index]]$value

      # Determine column width 
      column_width = report[[section_index]]$column_width

      # Determine the page break status 
      page_break = report[[section_index]]$page_break
      if (is.null(page_break)) page_break = FALSE

      # Determine the figure's location (for figures only)
      filename = report[[section_index]]$filename

      # Determine the figure's dimensions (for figures only)
      dim = report[[section_index]]$dim

      if (!is.null(type)) {

        # Fully formatted data frame 
        if (type == "table") {

            doc = officer::body_add_par(doc, value = label, style = "heading 2")

            summary_table = flextable::regulartable(data = value)
            summary_table = flextable::style(summary_table, pr_p = leftPar, pr_c = header.cellProperties, pr_t = header.textProperties, part = "header")
            summary_table = flextable::style(summary_table, pr_p = leftPar, pr_c = data.cellProperties, pr_t = data.textProperties, part = "body")

            summary_table = flextable::hline_bottom(summary_table, part = "body", border = thick_border )

            summary_table = flextable::width(summary_table, width = column_width)

            doc = flextable::body_add_flextable(doc, summary_table)

            if (!is.null(footnote)) doc = officer::body_add_par(doc, value = footnote, style = "Normal")

            if (page_break) doc = officer::body_add_break(doc, pos = "after")

        }

        # Enhanced metafile graphics produced by package devEMF 
        if (type == "emf_plot") {

            doc = officer::body_add_par(doc, value = label, style = "heading 2")

            doc = officer::body_add_img(doc, src = filename, width = dim[1], height = dim[2]) 

            if (!is.null(footnote)) doc = officer::body_add_par(doc, value = footnote, style = "Normal")

            if (page_break) doc = officer::body_add_break(doc, pos = "after")

            # Delete the figure
            if (file.exists(filename)) file.remove(filename)   

        }

      }    

  }

  return(doc)       

}
# End SaveReport

CreateTable = function(data_frame, column_names, column_width, title, page_break, footnote = NULL) {

    if (is.null(column_width)) {
      column_width = rep(2, dim(data_frame)[2])
    } 
     
    data_frame = as.data.frame(data_frame)

    data_frame = data.frame(lapply(data_frame, as.character), stringsAsFactors = FALSE)

    colnames(data_frame) = column_names

    item_list = list(label = title, 
                     value = data_frame,
                     column_width = column_width,
                     type = "table",
                     footnote = footnote,
                     page_break = page_break)

    return(item_list)

  }
# End of CreateTable  
  
ExtractSubgroupDescription = function(description, biomarker_names) {

  # List of biomarker labels
  sub_label = list()

  # List of biomarker names
  sub_name = list()

  subs = unlist(strsplit(description, ":"))

  for(i in 1:length(subs)) {

    current_sub = subs[i]

    less = (regexpr('<=', current_sub) > 0)
    greater = (regexpr('>', current_sub) > 0)
    equal = (regexpr('=', current_sub) > 0 & regexpr('<=', current_sub) < 0) 

    if (less) {

      pos = as.numeric(regexpr('<=', current_sub))
      name = substr(current_sub, 1, pos - 1)
      value = as.numeric(substr(current_sub, pos + 2, nchar(current_sub)))
      index = which(biomarker_names == name)
      sub_label[[i]] = paste0(biomarker_names[index], " <= ", value)

    }

    if (greater) {

      pos = as.numeric(regexpr('>', current_sub))
      name = substr(current_sub, 1, pos - 1)
      value = as.numeric(substr(current_sub, pos + 1, nchar(current_sub)))
      index = which(biomarker_names == name)
      sub_label[[i]] = paste0(biomarker_names[index], " > ", value)

    }

    if (equal) {

      pos = as.numeric(regexpr('=', current_sub))
      name = substr(current_sub, 1, pos - 1)
      value = substr(current_sub, pos + 1, nchar(current_sub))
      levels = unlist(strsplit(value, ","))
      index = which(biomarker_names == name)
      sub_label[[i]] = paste0(biomarker_names[index], " = ", paste(levels, collapse=", ")) 

    }

    sub_name[[i]] = name

  }

  result = list(description = paste(sub_label, collapse=" and "))

  return(result) 
    
}  

SubgroupSummary = function(sorted_subgroups, biomarker_names) {

      n_subgroups = dim(sorted_subgroups)[1]

      subgroups = list()

      k = 1

      for (i in 1:n_subgroups) {

        prom_estimate = round(as.numeric(sorted_subgroups[i,5]), 4)
        prom_sderr = round(as.numeric(sorted_subgroups[i,6]), 4)
        prom_sd = round(as.numeric(sorted_subgroups[i,7]), 4)
        pvalue = round(as.numeric(sorted_subgroups[i,8]), 8)
        adj_pvalue = round(as.numeric(sorted_subgroups[i,9]), 4)

        # Overall population
        if (i == 1) {
          description = "Overall population"
          comp_estimate = NA
          comp_sderr = NA
          adj_pvalue = NA
        } else {
          result = ExtractSubgroupDescription(sorted_subgroups[i,1], biomarker_names)         
          description = result$description
        }

        # Save subgroup characteristics
        subgroups[[k]] = list(description = description,
                            size = sorted_subgroups[i,2],
                            size_control = sorted_subgroups[i,3],
                            size_treatment = sorted_subgroups[i,4],
                            prom_estimate = prom_estimate,
                            prom_sderr = prom_sderr,
                            prom_sd = prom_sd,
                            comp_estimate = comp_estimate,
                            comp_sderr= comp_sderr,
                            pvalue = pvalue,
                            adj_pvalue = adj_pvalue)
        k = k + 1

      }

      return(subgroups)

}


GenerateReport = function(results, report_title, report_filename) {

  endpoint_parameters = results$parameters$endpoint_parameters
  subgroup_search_parameters = results$parameters
  patient_subgroups = results$patient_subgroups 
  biomarker_names = results$parameters$data_set_parameters$biomarker_names

  n_biomarkers = length(biomarker_names)

  item_list = list()
  item_index = 1

  ##############################################################################

  # Endpoint list

  column_names = c("Endpoint", "Analysis method")

  if (endpoint_parameters$direction == 1) {
    direction = "Higher value indicates a beneficial effect"
  } else {
    direction = "Lower value indicates a beneficial effect" 
  }
  col1 = c(paste0(endpoint_parameters$label, " (", endpoint_parameters$outcome_variable, ")"), 
           paste0("Type: ", subgroup_search_parameters$data_set_parameters$outcome_variable_type), paste0("Direction: ", direction))

  test_type = endpoint_parameters$analysis_method

  if (!is.null(endpoint_parameters$cont_covariates)) cont_covariates = paste0("Continuous covariates: ", endpoint_parameters$cont_covariates) else cont_covariates = "Continuous covariates: NA"

  if (!is.null(endpoint_parameters$class_covariates)) class_covariates = paste0("Class covariates: ", endpoint_parameters$class_covariates) else class_covariates = "Class covariates: NA"

  col2 = c(test_type, cont_covariates, class_covariates)

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", item_index, ". Endpoint parameters")

  column_width = c(3, 3.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1


  ##############################################################################

  # Biomarker list

  column_names = c("Biomarker ", "Type", "Level")

  col1 = NULL
  col2 = NULL
  col3 = NULL

  for (i in 1:n_biomarkers) {

    col1 = c(col1, paste0("Biomarker ", i, " (", subgroup_search_parameters$data_set_parameters$biomarker_names[i], ")"))
    if (subgroup_search_parameters$data_set_parameters$biomarker_types[i] == "numeric") {
      col2 = c(col2, "Numeric biomarker")
    } else {
      col2 = c(col2, "Nominal biomarker") 
    }
    
    col3 = c(col3, as.character(results$parameters$data_set_parameters$biomarker_levels[i]))

  }

  data_frame = data.frame(col1, col2, col3)
  title = paste0("Table ", item_index, ". List of candidate biomarkers")

  footnote = "A biomarker's level defines the first subgroup search level at which this biomarker is introduced."

  column_width = c(2.5, 2, 2)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1

  ##############################################################################

  # Subgroup search parameters

  column_names = c("Parameter ", "Value")

  n_perms_mult_adjust = subgroup_search_parameters$algorithm_parameters$n_perms_mult_adjust * subgroup_search_parameters$algorithm_parameters$ncores

  col1 = c("Subgroup search algorithm", # 1
           "Minimum subgroup size",     # 2
           "Search depth",              # 3
           "Search width",              # 4
           "Complexity parameters (child-to-parent ratios)", # 5
           "Number of permutations to compute multiplicity-adjusted treatment effect p-values within promising subgroups", # 6
           "Percentile transformation of numeric biomarkers: Maximum number of unique values", # 7
           "Random seed")                                                                      # 8
  col2 = c(subgroup_search_parameters$algorithm_parameters$subgroup_search_algorithm, # 1
           subgroup_search_parameters$algorithm_parameters$min_subgroup_size,         # 2
           subgroup_search_parameters$algorithm_parameters$depth,                     # 3
           subgroup_search_parameters$algorithm_parameters$width,                     # 4
           paste0(subgroup_search_parameters$algorithm_parameters$gamma, collapse=","), # 5
           n_perms_mult_adjust,                                                         # 6
           subgroup_search_parameters$algorithm_parameters$nperc,                       # 7
           subgroup_search_parameters$algorithm_parameters$random_seed)                 # 8

  # Additional parameters
  if (subgroup_search_parameters$algorithm_parameters$subgroup_search_algorithm == "Fixed SIDEScreen procedure") {
    col1 = c(col1, "Number of biomarkers selected for the second stage in Fixed SIDEScreen algorithm") # 9
    col2 = c(col2, subgroup_search_parameters$algorithm_parameters$n_top_biomarkers)                   # 9
  }

  if (subgroup_search_parameters$algorithm_parameters$subgroup_search_algorithm == "Adaptive SIDEScreen procedure") {
    col1 = c(col1, 
             "Multiplier for selecting biomarkers for the second stage in Adaptive SIDEScreen algorithm", # 9
             "Number of permutations for computing the null distribution of the maximum VI score in Adaptive SIDEScreen algorithm" # 10
            )

    # Turn around the empty value
    n_perms_vi_score_prn = subgroup_search_parameters$algorithm_parameters$n_perms_vi_score
    if (is.null(n_perms_vi_score_prn)) { n_perms_vi_score_prn = "(undefined)" }

    col2 = c(col2, 
             subgroup_search_parameters$algorithm_parameters$multiplier, # 9
             n_perms_vi_score_prn # 10
            )
  }

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", item_index, ". Subgroup search parameters")

  column_width = c(5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1

  ##############################################################################

  column_names = c("Parameter", "Value")

  if (!is.na(results$patient_subgroups$`Type I error rate`)) error_rate = as.character(patient_subgroups$`Type I error rate`) else error_rate = "NA"

  col1 = c("Elapsed time (sec)", "Type I error rate")
  col2 = c(as.character(patient_subgroups$`Elapsed time (seconds)`), 
           error_rate)

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", item_index, ". Algorithm's parameters")

  column_width = c(3, 3.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
  item_index = item_index + 1

  ##############################################################################

  # Sort by size
  sorted_subgroups = patient_subgroups$Subgroup[order(-as.numeric(patient_subgroups$Subgroup[,2])), ]

  subgroup_summary = SubgroupSummary(sorted_subgroups, biomarker_names)

  col1 = NULL
  col2 = NULL

  n_subgroups = length(subgroup_summary)

  for (i in 1:n_subgroups) {

    if (i == 1) subgroup_description = "Overall population" else subgroup_description = paste0("Subgroup ", i - 1, ": ", subgroup_summary[[i]]$description)

    if (endpoint_parameters$analysis_method == "T-test") {      

      col1 = c(col1, subgroup_description, rep("", 6))    
      col2 = c(col2, c(paste0("Total number of patients = ", subgroup_summary[[i]]$size),
                       paste0("Number of control patients = ", subgroup_summary[[i]]$size_control),
                       paste0("Number of treatment patients = ", subgroup_summary[[i]]$size_treatment),
                       paste0("Treatment effect = ", subgroup_summary[[i]]$prom_estimate),
                       paste0("Standard deviation = ", subgroup_summary[[i]]$prom_sd),
                       paste0("One-sided unadjusted p-value = ", subgroup_summary[[i]]$pvalue),
                       paste0("One-sided multiplicity adjusted p-value = ", subgroup_summary[[i]]$adj_pvalue)
                       )
              ) 

    }

    if (endpoint_parameters$analysis_method %in% c("Z-test for proportions", "Log-rank test")) {      

      col1 = c(col1, subgroup_description, rep("", 5))    
      col2 = c(col2, c(paste0("Total number of patients = ", subgroup_summary[[i]]$size),
                       paste0("Number of control patients = ", subgroup_summary[[i]]$size_control),
                       paste0("Number of treatment patients = ", subgroup_summary[[i]]$size_treatment),
                       paste0("Treatment effect = ", subgroup_summary[[i]]$prom_estimate),
                       paste0("One-sided unadjusted p-value = ", subgroup_summary[[i]]$pvalue),
                       paste0("One-sided multiplicity adjusted p-value = ", subgroup_summary[[i]]$adj_pvalue)
                       )
              ) 

    }
    if (endpoint_parameters$analysis_method %in% c("ANCOVA", "Logistic regression", "Cox regression")) {      

      col1 = c(col1, subgroup_description, rep("", 6))    
      col2 = c(col2, c(paste0("Total number of patients = ", subgroup_summary[[i]]$size),
                       paste0("Number of control patients = ", subgroup_summary[[i]]$size_control),
                       paste0("Number of treatment patients = ", subgroup_summary[[i]]$size_treatment),
                       paste0("Treatment effect = ", subgroup_summary[[i]]$prom_estimate),
                       paste0("Standard error = ", subgroup_summary[[i]]$prom_sderr),
                       paste0("One-sided unadjusted p-value = ", subgroup_summary[[i]]$pvalue),
                       paste0("One-sided multiplicity adjusted p-value = ", subgroup_summary[[i]]$adj_pvalue)
                       )
              ) 

    }

  }

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", item_index, ". Subgroup summary")

  column_names = c("Subgroup", "Subgroup's characteristics")

  footnote = "The patient subgroups are sorted by the total number of patients in this table. "

  if (endpoint_parameters$type == "continuous") footnote = paste0(footnote, "The treatment effect is defined as the mean difference.")

  if (endpoint_parameters$type == "binary" & endpoint_parameters$analysis_method == "Z-test for proportions") footnote = paste0(footnote, "The treatment effect is defined as the difference in proportions.")

  if (endpoint_parameters$type == "binary" & endpoint_parameters$analysis_method == "Logistic regression") footnote = paste0(footnote, "The treatment effect is defined as the odds ratio.")

  if (endpoint_parameters$type == "survival") footnote = paste0(footnote, "The treatment effect is defined as the hazard ratio.")

  column_width = c(3, 3.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE, footnote)
  item_index = item_index + 1

  vi = patient_subgroups$`Variable Importance`
  vi[, 2] = round(vi[, 2], 3)
  vi = vi[vi[, 2] > 0, ]

  if(dim(vi)[1] > 0) {

    column_names = c("Biomarker", "Variable importance score")

    data_frame = vi

    title = paste0("Table ", item_index, ". Variable importance summary")

    footnote = "Only biomarkers with non-zero variable importance scores are shown in this table."

    if (subgroup_search_parameters$algorithm_parameters$subgroup_search_algorithm == "Adaptive SIDEScreen procedure") {
      vi_threshold = round(results$patient_subgroups$`VI Threshold`, 3)
      footnote = paste0(footnote, " Parameters of the null distribution of the maximum VI score used in Adaptive SIDEScreen algorithm: Mean = ", vi_threshold[1], ", Standard error = ", vi_threshold[2], ", Threshold for selecting biomarkers for the second stage = ", vi_threshold[3], ".")

    }

    column_width = c(3, 3.5)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1

  }

  ##############################################################################

  doc = SaveReport(item_list, report_title)

  # Save the report
  print(doc, target = report_filename)          


}
# End of GenerateReport
