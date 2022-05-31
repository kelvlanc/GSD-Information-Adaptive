library(mvtnorm) # for rmvnorm in iid_centered_mvn()
library(betareg) # for impute_monotone_beta
library(nnet) # for impute_monotone_multinomial
library(boot) # for tmle - bootstrap inference

### compress_range #############################################################
compress_range <- 
  function(
    x,
    min,
    max,
    enforce_range = TRUE
  ) {
    
    if(sum(x > max, na.rm = TRUE) + sum(x < min, na.rm = TRUE) > 0) {
      warning("x contains values outside of specified range")
      if(enforce_range){ 
        warning("values outside of specified range replaced with boundaries")
        x[which(x > max)] <- max
        x[which(x < min)] <- min
      }
    }
    
    (x - min)/(max - min)
  }




### decompress_range ###########################################################
decompress_range <-
  function(
    x,
    min,
    max
  ) {
    x*(max - min) + min
  }




### iid_centered_mvn ###########################################################
iid_centered_mvn <-
  function(
    n,
    n_covariates,
    sigma_x
  ) {
    setNames(
      object =
        data.frame(
          mvtnorm::rmvnorm(
            n = n,
            sigma = sigma_x
          )
        ),
      nm = paste0("x_", 1:n_covariates)
    )
  }




### sim_lme_trial ##############################################################
# Create a data.frame of n simulated participants with p baseline covariates and
#   T continuous outcomes, based on a linear mixed effects (LME) model.
#   Simulated data can incorporate missing data due to study dropout.
#
# n: number of participants
# T: times of outcomes: (T-1) intermediate and one primary outcome
# p: number of baseline covariates
#
# NOTE: First the complete data vector is constructed, and then missingness is
# applied in the last step - missing and complete datasets are returned. If MAR
# mechanism is used, missingness should only depend on previous outcomes and
# baseline covariates.
#
# How to specify outcome and missingness regression models: Suppose the complete
# data vector for each participant consists of 5 baseline covariates and 3
# outcomes. For simplicity, the covariates are named X1, X2, ..., X5.
#
# If outcome_cov is not specified, the covariates will be unrelated to the
# outcome. Elements should be a named vector, with names corresponding to the
# names of baseline covariates or treatment ("tx"). These are used to create
# linear predictors, and the outcome is sampled from a normal distribution given
# the linear predictors and random effects.
#
# Models for dropout are specified similarly. Probabilities are supplied for
# the probability of dropping out at a given visit when all covariates are zero.
# Logistic regression coefficients can be supplied for baseline covariates,
# treatment assignment ("tx"), or outcomes (e.g. "y_1", "y_2", "y_3", etc.),
# however, for the MAR mechanism, missingness should only depend on previously
# measured variables.
# 
# input:
# bl_covariates: (n x p data.frame) of continuous or binary indicators
#   representing covariates measured prior to baseline with no missing data.
# mean_outcomes: (numeric of length T)) - the mean outcome at each assessment
#   for an individual with all covariates and random effects equal to 0.
# outcome_cov: (list of length (T-1)) - regression coefficients
#   for associations between prior data and outcomes, with names corresponding
#   to the columns for baseline covariates, 'treatment', or prior outcomes
#   ('y_1', 'y_2', ...). If NULL (default), outcomes will have no association
#   with either baseline covariates, treatment assignment, or intermediate
#   outcomes.  If a single element in list is NULL, no covariates will have
#   an association with outcomes at that visit
# re_variance: (numeric of length 1-2) - the variance of the random slopes and
#   intercepts.
# re_correlation: correlation between random slopes and random intercepts: only
#   necessary when both random slopes and intercepts are supplied (re_variance
#   has length of 2).
# residual_sd: residual standard deviation
# pr_dropout: (1 x T probability) rate of dropout PRIOR TO each visit for those
#   with the default level of all covariates. If NULL (default), no dropout is
#   introduced.
# dropout_cov: (list of length (T-1)) - logistic regression coefficients
#   for associations between dropout and data (covariates, treatment
#   assignment, or outcomes) at each visit, with names corresponding to the
#   appropriate variables. If NULL (default), dropout will have no association
#   with either baseline covariates, treatment assignment, or intermediate
#   outcomes. If a single element in list is NULL, no covariates will have
#   an association with dropout at that visit.
#
# NOTE: For missingness models, outcomes are CENTERED, but not scaled: If an
#   outcome is included in the model, the intercept is the probability of
#   missingness for individuals with covariates equal to 0 and the population
#   mean value for the outcome.

sim_lme_trial <-
  function(
    bl_covariates,
    visit_times,
    mean_outcomes,
    outcome_cov = NULL,
    re_variance,
    re_correlation = NULL,
    residual_sd,
    pr_dropout = NULL,
    dropout_cov = NULL
  ) {
    
    n_random_effects <-
      nrow(as.matrix(re_variance))
    
    stopifnot(length(visit_times) == length(mean_outcomes))
    
    if(n_random_effects == 1) {
      # Only random intercept is supplied
      re_covariance <- re_variance
    } else if(n_random_effects == 2) {
      # Slope/intercept supplied
      if(is.null(re_correlation)) { # Uncorrelated random effects
        re_correlation <- diag(n_random_effects)
      } else if (nrow(as.matrix(re_correlation)) == 1) { # Scalar correlation
        re_correlation <-
          diag(2) +
          re_correlation*(matrix(data = 1, nrow = 2, ncol = 2) - diag(2))
      } 
      
      re_covariance <- 
        (matrix(sqrt(re_variance), ncol = 1) %*% 
           matrix(sqrt(re_variance), nrow = 1))*re_correlation
      
    } else {
      stop(paste0("Only random intercept, random slope, or random intercept ",
                  "and slope are implemented."))
    }
    
    
    n_outcomes <- length(x = visit_times)
    n_obs <- nrow(bl_covariates)
    
    # 1. Treatment is assigned
    study_data <-
      data.frame(
        id = 1:nrow(bl_covariates),
        bl_covariates,
        tx =
          sample(x = c(rep(0, ceiling(n_obs/2)),
                       rep(1, n_obs - ceiling(n_obs/2))),
                 size = n_obs),
        last_observed = NA,
        dropout_i = NA
      )
    
    # 2. Determine Distribution of outcomes
    
    # 2.1 Determine Fixed Effects
    out_cov_lp <- matrix(0, nrow = n_obs, ncol = n_outcomes)
    
    if(!is.null(outcome_cov)){
      for(i in 1:n_outcomes) {
        if(!is.null(outcome_cov[[i]])) {
          out_cov_lp[, i] <-
            as.matrix(study_data[, names(outcome_cov[[i]])]) %*% 
            outcome_cov[[i]]
        }
      }
    }
    
    # 2.2 Simulate Random effects are simulated
    # Sample random effects
    random_effects <-
      rmvnorm(
        n = n_obs,
        sigma = as.matrix(re_covariance)
      )
    
    # Create (n x T) matrix of disturbances due to random effects
    if(ncol(random_effects) == 1){ # Random Intercept
      random_trajectories <-
        kronecker(random_effects,
                  matrix(1, ncol = length(visit_times)))
      colnames(random_effects) <- c("random_intercept")
    } else { # Centered Random Effects at Mean Visit Time
      random_trajectories <-
        random_effects %*%
        rbind(1, visit_times - mean(visit_times))
      colnames(random_effects) <- c("random_intercept", "random_slope")
    }
    
    # 2.3. Add fixed effects and random effects to residual
    residuals <-
      rmvnorm(n = n_obs,
              mean = matrix(0, ncol = n_outcomes),
              sigma = diag(residual_sd^2, nrow = n_outcomes)) # Residual
    colnames(residuals) <- paste0("residual_", 1:n_outcomes)
    
    study_data[, paste0("y_", 1:n_outcomes)] <-
      study_data[, paste0("y_obs_", 1:n_outcomes)] <-
      # Note: In order to make the dropout model easier to use, the mean
      # outcome is added in later.
      matrix(data = 0, nrow = n_obs, ncol = n_outcomes) +
      out_cov_lp + # Linear Predictor
      random_trajectories + # Random Effects
      residuals
    
    study_data <-
      data.frame(study_data, random_effects, residuals)
    
    # 3. dropout is assigned
    for(i in 1:n_outcomes) {
      
      if(is.null(pr_dropout[[i]])){
      } else if(pr_dropout[[i]] > 0 & pr_dropout[[i]] < 1) {
        
        if(is.null(dropout_cov[[i]])) {
          drop_cov_lp <- qlogis(1 - pr_dropout[[i]])
        } else {
          drop_cov_lp <- qlogis(1 - pr_dropout[[i]]) -
            as.matrix(study_data[names(dropout_cov[[i]])]) %*% 
            dropout_cov[[i]]
        }
        
        study_data$dropout_i <- 
          rlogis(n = n_obs, location = drop_cov_lp, scale = 1) < 0
        
        dropout.rows <-
          with(study_data, which(is.na(last_observed) & dropout_i))
        
        study_data[dropout.rows,
                   paste0("y_obs_", i:n_outcomes)] <- NA
        
        study_data$last_observed[dropout.rows] <- (i - 1)
      }
    }
    
    study_data$last_observed[which(is.na(study_data$last_observed))] <- 
      n_outcomes
    
    study_data$dropout_i <- NULL
    
    # 4. Add in mean outcome
    study_data[, paste0("y_", 1:n_outcomes)] <-
      study_data[, paste0("y_", 1:n_outcomes)] + 
      matrix(data = mean_outcomes, nrow = n_obs, ncol = n_outcomes,
             byrow = TRUE)
    
    study_data[, paste0("y_obs_", 1:n_outcomes)] <-
      study_data[, paste0("y_obs_", 1:n_outcomes)] + 
      matrix(data = mean_outcomes, nrow = n_obs, ncol = n_outcomes,
             byrow = TRUE)
    
    study_data
  }




### impute_covariates_mean_mode ################################################
# Uses mean for numeric covariates, mode for non-numeric covariates.
impute_covariates_mean_mode <-
  function(data, x_columns) {
    
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    any_missing <- function(x){
      any(is.na(x))
    }
    
    impute_mean <- function(x) {
      ifelse(
        test = is.na(x),
        yes = mean(x, na.rm = TRUE),
        no = as.numeric(x)
      )
    }
    
    impute_mode <- function(x) {
      ifelse(
        test = is.na(x),
        yes = Mode(x),
        no = x
      )
    }
    
    numeric_vars <-
      names(which(sapply(X = data[x_columns], FUN = is.numeric)))
    
    non_numeric_vars <-
      setdiff(x_columns, numeric_vars)
    
    for(i in numeric_vars){
      data[, i] <- impute_mean(data[, i])
    }
    
    for(i in non_numeric_vars){
      data[, i] <- impute_mode(data[, i])
    }
    
    return(data)
  }




### non_monotone ###############################################################
# determine which missing values are non-monotone
non_monotone <- function(data, y_columns){
  n_outcomes <- length(y_columns)
  
  non_monotone <-
    matrix(
      data = NA,
      nrow = nrow(data),
      ncol = length(y_columns) - 1
    )
  
  for(i in 1:ncol(non_monotone)){
    
    non_monotone[, i] <-
      is.na(data[, y_columns[i]]) &
      if(i == ncol(non_monotone)) {
        !is.na(data[, y_columns[(i + 1):n_outcomes]])
      } else{
        rowSums(
          !is.na(data[, y_columns[(i + 1):n_outcomes]])
        ) > 0
      }
  }
  
  return(
    setNames(
      object = data.frame(non_monotone),
      nm = head(y_columns, - 1)
    )
  )
}




### absorbing_state_check ######################################################
absorbing_state_check <-
  function(
    data,
    y_columns,
    outcome_type,
    absorbing_state,
    absorbing_outcome
  ) {
    
    if(outcome_type != "multinomial-binomial"){
      stop("Absorbing state only for multinomial-binomial outcomes.")
    } else if(is.null(absorbing_state) | is.null(absorbing_outcome)){
      stop("'absorbing_state' and 'absorbing_outcome' must be supplied.")
    }
    
    final_outcome <- tail(x = y_columns, 1)
    intermediate_outcomes <- setdiff(x = y_columns, y = final_outcome)
    
    if(length(intermediate_outcomes > 0)){
      
      first_absorbing <-
        apply(
          X = data[, head(x = y_columns, -1)],
          FUN = function(x) which(x %in% absorbing_state)[1],
          MARGIN = 1
        )
      
      last_non_absorbing <-
        apply(
          X = data[, head(x = y_columns, -1)],
          FUN = function(x) which(!(x %in% c(NA, absorbing_state)))[1],
          MARGIN = 1
        )
      
      absorbing_conflicts <-
        which(last_non_absorbing > first_absorbing)
      
      if(length(absorbing_conflicts) > 0){
        stop("Transition out of absorbing state in rows: ",
             paste0(absorbing_conflicts, collapse = ", "))
      }
      
      # Carry forward absorbing state
      for(i in 1:(length(y_columns) - 1)){
        rows <- which(first_absorbing == i)
        data[rows, intermediate_outcomes] <- absorbing_state
        data[rows, final_outcome] <- absorbing_outcome
      }
    }
    
    return(data)
  }




### impute_outcomes_to_monotone ################################################
impute_outcomes_to_monotone <-
  function(
    model =
      c("gaussian",
        "beta",
        "binomial",
        "pmm",
        "multinomial-binomial"),
    absorbing_state = NULL,
    ...
  ) {
    arg_list <- as.list(substitute(list(...)))[-1L]
    
    if(model == "gaussian"){
      do.call(
        what = impute_gaussian,
        args = arg_list
      )
    } else if(model == "beta") {
      do.call(
        what = impute_beta,
        args = arg_list
      )
    } else if(model == "binomial") {
      do.call(
        what = impute_binomial,
        args = arg_list
      )
    } else if(model == "pmm") {
      do.call(
        what = impute_pmm,
        args = arg_list
      )
    } else if(model == "multinomial") {
      if(is.null(absorbing_state)){
        do.call(
          what = impute_multinomial,
          args = arg_list
        )
      } else {
        do.call(
          what = impute_multinomial_absorbing,
          args = 
            c(
              arg_list,
              absorbing_state = absorbing_state
            )
        )
      }
    } else{
      stop("`model` not recognized.")
    }
  }




### impute_gaussian ###################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# Gaussian GLM. `impute_columns` and the LHS of `impute_formulas` must
# match temporal ordering of outcomes - This is checked by tmle_precheck().
impute_gaussian <-
  function(
    data = data,
    impute_columns = impute_columns,
    impute_formulas,
    verbose = FALSE,
    alpha = 0.05
  ) {
    
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_to_impute <- colSums(impute_columns)
    
    imputation_models <- list()
    
    for(i in 1:length(impute_formulas)){
      
      # Fit Imputation Model among observed cases
      impute_model <-
        lm(
          formula = impute_formulas[[i]],
          data = data
        )
      
      if(verbose) imputation_models[[i]] <- impute_model
      
      # Sample from Prediction Interval
      pred_interval <-
        with(
          predict(
            object = impute_model,
            newdata = data[which(impute_columns[, i]),],
            interval = "prediction",
            se = TRUE
          ),
          data.frame(fit, se.fit, df)
        )
      
      # Fill in imputed values
      data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
        with(
          pred_interval,
          rnorm(
            n = n_to_impute[i],
            mean = fit,
            sd = (upr - fit)/qt(p = 1 - alpha/2, df = df)
          )
        )
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models
        )
      )
    } else{
      return(data)
    }
  }









### impute_beta #######################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# Beta GLM. `impute_columns` and the LHS of `impute_formulas` must
# match temporal ordering of outcomes - This is checked by tmle_precheck().
impute_beta <-
  function(
    data = data,
    impute_columns = impute_columns,
    impute_formulas,
    verbose = FALSE,
    link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog")[1],
    model_scale = FALSE,
    type = c("BC", "ML", "BR")[1]
  ) {
    
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_to_impute <- colSums(impute_columns)
    
    imputation_models <- list()
    
    for(i in 1:length(impute_formulas)){
      
      ## HERE ###
      
      # If scale is modeled, adjust formula accordingly
      if(model_scale){
        model_rhs <-
          impute_formulas[[i]][[length(impute_formulas[[i]])]]
        
        model_rhs <-
          paste0(Reduce(paste, deparse(model_rhs)), " | ", 
                 Reduce(paste, deparse(model_rhs)))
        
        impute_formulas[[i]] <-
          reformulate(
            termlabels = model_rhs,
            response = impute_formula_lhs[i]
          )
      }
      
      # Fit Imputation Model among observed cases
      impute_model <-
        betareg(
          formula = impute_formulas[[i]],
          link = link,
          data = data,
          type = type
        )
      
      if(verbose) imputation_models[[i]] <- impute_model
      
      beta_params <-
        data.frame(
          mean = 
            predict(
              object = impute_model,
              newdata = data[which(impute_columns[,i]),],
              type = "response"
            ),
          var = 
            predict(
              object = impute_model,
              newdata = data[which(impute_columns[,i]),],
              type = "variance"
            )
        )
      
      beta_params$m <- 
        with(beta_params, (mean*(1 - mean)/var) - 1)
      beta_params$alpha <- 
        with(beta_params, mean*m)
      beta_params$beta <- 
        with(beta_params, (1 - mean)*m)
      
      if(verbose) imputation_parameters[[i]] <- beta_params
      
      # Sample in imputed values
      data[which(impute_columns[,i]), var_to_impute] <-
        rbeta(
          n = n_to_impute[i],
          shape1 = beta_params$alpha,
          shape2 = beta_params$beta,
        )
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models,
          imputation_parameters = imputation_parameters
        )
      )
    } else{
      return(data)
    }
  }




### impute_binomial ###################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# binomial GLM with user-specified link. `impute_columns` and the LHS of
# `impute_formulas` must match temporal ordering of outcomes -
# This is checked by tmle_precheck().
impute_binomial <-
  function(
    data = data,
    impute_columns = impute_columns,
    impute_formulas,
    impute_link = c("logit", "probit", "cauchit", "cloglog")[1],
    verbose = FALSE
  ) {
    
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_to_impute <- colSums(impute_columns)
    
    imputation_models <- list()
    
    for(i in 1:length(impute_formulas)){
      
      # Fit Imputation Model among observed cases
      impute_model <-
        glm(
          formula = impute_formulas[[i]],
          data = data,
          family = binomial(link = impute_link)
        )
      
      if(verbose) imputation_models[[i]] <- impute_model
      
      # Fill in imputed values
      data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
        rbinom(
          n = n_to_impute[i],
          size = 1,
          prob = predict(
            object = impute_model,
            newdata = data[which(impute_columns[, i]), ],
            type = "response"
          )
        )
      
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models
        )
      )
    } else{
      return(data)
    }
  }


### impute_pmm #################################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# binomial GLM with user-specified link. `impute_columns` and the LHS of
# `impute_formulas` must match temporal ordering of outcomes -
# This is checked by tmle_precheck().
impute_pmm <-
  function(
    data = data,
    impute_columns = impute_columns,
    impute_formulas,
    impute_family = c(gaussian, quasipoisson, quasibinomial)[[1]],
    impute_link = c("identity", "log", "inverse", "logit", "probit")[1],
    donors = 10,
    verbose = FALSE
  ){
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_to_impute <- colSums(impute_columns)
    
    imputation_models <- list()
    
    for(i in 1:length(impute_formulas)){
      
      # Fit Imputation Model among observed cases
      impute_model <-
        glm(
          formula = impute_formulas[[i]],
          data = data,
          family = impute_family(link = impute_link)
        )
      
      if(verbose) imputation_models[[i]] <- impute_model
      
      # Get fitted value for missing data
      y_hat <-
        predict(
          object = impute_model,
          newdata = data[which(impute_columns[, i]), ],
          type = "response"
        )
      
      # Get all non-missing values for candidates
      candidates <-
        na.omit(data[, impute_formula_lhs[i]])
      
      # Find distance between fit and candidates 
      candidate_distance <-
        kronecker(
          X = candidates,
          Y = matrix(data = 1, ncol = length(y_hat))
        ) - 
        kronecker(
          X = matrix(data = 1, nrow = length(candidates)),
          Y = matrix(data = y_hat, nrow = 1)
        )
      
      selected_candidates <-
        apply(
          X = candidate_distance^2,
          MARGIN = 2,
          FUN = function(x, donors)
            sample(x = which(rank(x) < donors), size = 1),
          donors = donors
        )
      
      data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
        candidates[selected_candidates]
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models
        )
      )
    } else{
      return(data)
    }
  }




### impute_multinomial #########################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# multinomial GLM with user-specified link. `impute_columns` and the LHS of
# `impute_formulas` must match temporal ordering of outcomes -
# This is checked by tmle_precheck().
# Use `impute_multinomial_absorbing` if the outcome scale has an
# absorbing state (e.g. Death) that should be preserved
impute_multinomial <-
  function(
    data = data,
    impute_columns = impute_columns,
    impute_formulas,
    verbose = FALSE
  ) {
    
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_to_impute <- colSums(impute_columns)
    
    imputation_models <- list()
    
    for(i in 1:length(impute_formulas)){
      
      outcome_levels <- levels(data[, impute_formula_lhs[i]])
      
      # Fit Imputation Model among observed cases
      impute_model <-
        multinom(
          formula = impute_formulas[[i]],
          data = data,
          trace = FALSE
        )
      
      if(verbose) imputation_models[[i]] <- impute_model
      
      # Fill in imputed values
      imputations <-
        apply(
          X = 
            predict(
              object = impute_model,
              newdata = data[which(impute_columns[, i]), ],
              type = "probs"
            ),
          MAR = 1,
          FUN = function(x, outcome_levels) 
            outcome_levels[
              which(rmultinom(n = 1, size = 1, prob = x) == 1)
            ],
          outcome_levels = outcome_levels
        )
      
      data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
        imputations
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models
        )
      )
    } else{
      return(data)
    }
  }




### impute_multinomial_absorbing ################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# multinomial GLM with user-specified link. `impute_columns` and the LHS of
# `impute_formulas` must match temporal ordering of outcomes -
# This is checked by tmle_precheck().
# Use `impute_multinomial_absorbing` if the outcome scale has an
# absorbing state (e.g. Death) that should be preserved
impute_multinomial_absorbing <-
  function(
    data = data,
    impute_columns,
    y_columns,
    impute_formulas,
    absorbing_state,
    verbose = FALSE
  ) {
    
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_outcomes <- ncol(impute_columns)
    n_to_impute <- colSums(impute_columns)
    all_outcome_levels <-
      unique(
        unlist(
          lapply(
            X = data[, names(impute_columns)],
            FUN = levels
          )
        )
      )
    
    imputation_models <- list()
    
    for(i in 1:length(impute_formulas)){
      outcome_levels <- levels(data[, impute_formula_lhs[i]])
      non_absorbing <- setdiff(x = outcome_levels, y = absorbing_state)
      
      outcome_i <- which(y_columns == impute_formula_lhs[i])
      
      # Range of imputations depends on the next observed value:
      # If next observed outcome is absorbing state, use full scale.
      # If next observed outcome is NOT absorbing state, impute from
      # non-absorbing states to avoid inconsistency
      
      next_observed_outcome <-
        apply(
          X = as.matrix(data[, tail(x = y_columns, -outcome_i)]),
          MARGIN = 1,
          FUN = function(x, table = all_outcome_levels)
            min(match(x = x, table = table))
        )
      
      impute_full <-
        (impute_columns[, impute_formula_lhs[i]]) &
        (next_observed_outcome %in% absorbing_state)
      
      impute_non_absorbing <-
        (impute_columns[, impute_formula_lhs[i]]) &
        !(next_observed_outcome %in% absorbing_state)
      
      
      if(sum(impute_full) > 0){
        
        impute_model <-
          multinom(
            formula = impute_formulas[[i]],
            data = data,
            trace = FALSE
          )
        
        if(verbose)
          imputation_models[[i]] <-
            list(full_range = impute_model)
        
        imputations <-
          apply(
            X = predict(
              object = impute_model,
              newdata = data[which(impute_full), ],
              type = "probs"
            ),
            MAR = 1,
            FUN = function(x, outcome_levels) 
              outcome_levels[
                which(rmultinom(n = 1, size = 1, prob = x) == 1)
              ],
            outcome_levels = outcome_levels
          )
        
        
        # Carry absorbing state forward
        if(i < length(impute_formulas)) {
          new_absorbed_rows <-
            which(impute_full)[which(imputations %in% absorbing_state)]
          
          data[new_absorbed_rows, tail(x = y_columns, -outcome_i)] <-
            absorbing_state
          
          impute_columns[new_absorbed_rows, -c(1:outcome_i)] <- TRUE
        }
        
        data[which(impute_full), impute_formula_lhs[i]] <-
          imputations
      }
      
      if(sum(impute_non_absorbing) > 0){
        
        # Subset to those not absorbed in current or next outcome
        non_absorbed <- 
          data[
            which(!data[, y_columns[outcome_i]] %in% absorbing_state), 
          ]
        
        for(j in y_columns[1:outcome_i]){
          non_absorbed[, j] <- droplevels(non_absorbed[, j])
        }
        
        impute_model <-
          multinom(
            formula = impute_formulas[[i]],
            data = non_absorbed,
            trace = FALSE
          )
        
        if(verbose)
          imputation_models[[i]] <-
          list(non_absorbing = impute_model)
        
        imputations <-
          predict(
            object = impute_model,
            newdata = data[which(impute_non_absorbing), ],
            type = "probs"
          )
        
        
        
        
        if(!is.matrix(imputations)){
          # If only 2 classes, multinom produces vector, not matrix
          # Need probability for each class
          if(length(outcome_levels) == 2) {
            imputations <- cbind(1 - imputations, imputations)
          } else {
            imputations <- matrix(data = imputations, nrow = 1)
          }
        }
        
        imputations <-
          apply(
            X = imputations,
            MAR = 1,
            FUN = function(x, outcome_levels) 
              outcome_levels[
                which(rmultinom(n = 1, size = 1, prob = x) == 1)
              ],
            outcome_levels =
              levels(non_absorbed[, impute_formula_lhs[i]])
          )
        
        data[which(impute_non_absorbing), impute_formula_lhs[i]] <-
          imputations
      }
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models
        )
      )
    } else{
      return(data)
    }
  }




### compute_inverse_weights ####################################################
compute_inverse_weights <-
  function(
    data,
    inverse_weight_formulas,
    y_columns,
    verbose = FALSE,
    absorbing_state = NULL
  ) {
    n_outcomes <- length(y_columns)
    
    ipw_models <- list()
    
    missing_outcomes <- 
      colSums(is.na(data[, y_columns]))
    
    # No missing outcomes
    if (all(missing_outcomes == 0)) {
      # No missing data
      ip_weights <-
        matrix(
          data = 1,
          nrow = nrow(data),
          ncol = (n_outcomes + 1)
        )
      
    } else {
      
      # Count non-monotone missingness
      n_non_monotone <-
        sum(colSums(non_monotone(data = data, y_columns = y_columns)))
      
      if(n_non_monotone > 0){
        stop("Missing pattern must be monotone")
      }
      
      ip_weights <-
        matrix(
          data = NA,
          nrow = nrow(data),
          ncol = (n_outcomes + 1)
        )
      
      ip_weights[, 1] <- 1
      
      for(i in 1:n_outcomes){
        
        # If no dropouts or no new dropouts, carry weight forward
        # Find those previously absorbed or censored
        if(i > 1){
          prev_absorbed <-
            data[, y_columns[(i - 1)]] %in% absorbing_state
          
          prev_censored <-
            is.na(data[, y_columns[(i - 1)]])
          
          carry_weights_forward <- 
            missing_outcomes[i] == 0 |
            diff(missing_outcomes[(i - 1):i]) == 0
          
        } else {
          prev_absorbed <- prev_censored <-
            rep(FALSE, nrow(data))
          
          carry_weights_forward <- 
            missing_outcomes[i] == 0
        }
        
        uncensored <- which(!prev_censored)
        
        if(carry_weights_forward) {
          ip_weights[uncensored, (i + 1)] <-
            ip_weights[uncensored, i]
        } else {
          # Newly censored observations
          
          if(is.null(absorbing_state)){
            uncensored_unabsorbed <- uncensored
          } else {
            uncensored_unabsorbed <-
              which(!(prev_absorbed|prev_censored))
            
            # Carry weights forward for those absorbed
            if(sum(prev_absorbed) > 0){
              ip_weights[which(prev_absorbed), (i + 1)] <-
                ip_weights[which(prev_absorbed), i]
            }
          }
          
          ipw_model <-
            glm(
              formula = inverse_weight_formulas[[i]],
              data = data[uncensored_unabsorbed, ],
              family = binomial
            )
          
          if(verbose) ipw_models[[i]] <- ipw_model
          
          pr_cens_i <-
            predict(
              object = ipw_model,
              newdata = data[uncensored_unabsorbed,],
              type = "response"
            )
          
          ip_weights[uncensored_unabsorbed, (i + 1)] <-
            ip_weights[uncensored_unabsorbed, i]*(1/pr_cens_i)
        }
      }
    }
    
    if(verbose){
      return(
        list(
          # Drop leading column
          ipw = ip_weights[, -1],
          ipw_models = ipw_models
        )
      )
    } else{
      return(
        # Drop leading column
        ip_weights[, -1]
      )
    }
  }




### tmle_get_args ##############################################################
# Take arguments from call to TMLE: determine which variables are baseline
# covariates (x_columns), which is the treatment indicator (tx_column),
# and which are outcomes (y_columns) based on the formulas supplied.
# If imputation formulas are supplied, auxilliary covariates (in imputation
# model but not in propensity score, inverse weight, or outcome models) are 
# aslo returned (imputation_x).
tmle_get_args <-
  function(
    data,
    propensity_score_formula,
    inverse_weight_formulas,
    outcome_formulas,
    impute_formulas = NULL
  ) {
    
    # Get names of outcomes from left-hand-side of formulas
    y_columns <- 
      sapply(
        X = outcome_formulas,
        FUN =
          function(x) all.vars(update(x, . ~ 0))
      )
    
    n_outcomes <- length(y_columns)
    
    # Get name of treatment indicator from left-hand-side of formula
    tx_column <-
      all.vars(update(propensity_score_formula, . ~ 0))
    
    
    # Get names of predictors from formulas
    x_columns <-
      unique(
        do.call(
          what = c,
          args = 
            lapply(
              X = c(propensity_score_formula,
                    outcome_formulas,
                    inverse_weight_formulas),
              FUN = function(x, data)
                names(get_all_vars(formula = x, data = data)),
              data = data
            )
        )
      )
    
    # Outcomes may be in inverse weight or outcome formulas: remove them.
    x_columns <-
      setdiff(
        x = x_columns,
        y = c(y_columns, tx_column)
      )
    
    
    if(!is.null(impute_formulas)) {
      imputation_x <-
        do.call(
          what = c,
          args = 
            lapply(
              X = impute_formulas,
              FUN = function(x, data)
                names(get_all_vars(formula = x, data = data)),
              data = data
            )
        )
      
      imputation_x <-
        setdiff(
          x = imputation_x,
          y = c(x_columns, tx_column, y_columns)
        )
    } else {
      imputation_x <- NULL
    }
    
    return(
      list(
        x_columns = x_columns,
        tx_column = tx_column,
        y_columns = y_columns,
        imputation_x = imputation_x
      )
    )
  }




### tmle_precheck ##############################################################
# Check arguments from TMLE for potential problems: missing values in baseline
# covariates or treatment, treatment variable is not binary, missing data that
# is not monotone and imputation arguments are not specified, or the imputation
# arguments are not correctly specified.
# If any issues are identified, the function halts with specific error messages,
# and returns the indices of intermediate outcomes to be imputed if no errors
# are found.
tmle_precheck <-
  function(
    data,
    x_columns,
    y_columns,
    tx_column,
    imputation_x = NULL,
    impute_formulas = NULL,
    impute_model = NULL,
    outcome_type =
      c("gaussian",
        "logistic",
        "binomial")[1],
    outcome_range = NULL,
    absorbing_state = NULL
  ) {
    
    imputation_x <-
      setdiff(
        x = imputation_x,
        y = c(x_columns, tx_column, y_columns)
      )
    
    n_missing_covariates <-
      colSums(is.na(data[, c(x_columns, tx_column, imputation_x)]))
    
    # Check for missing covariates
    if(sum(n_missing_covariates) > 0) {
      n_missing_covariates <-
        n_missing_covariates[n_missing_covariates > 0]
      stop(
        "Missing covariate/treatment data: ",
        paste(
          paste0(names(n_missing_covariates),
                 " (", n_missing_covariates, ")"),
          collapse = ", "
        )
      )
    }
    
    # Check for treatment indicator that is not in {0, 1}
    if(sum(!(data[, tx_column] %in% 0:1)) > 0){
      stop(paste0("Treatment column `", tx_column, "` must be binary."))
    }
    
    
    # Check for logistic data outside [0, 1]
    if(outcome_type == "logistic"){
      if(is.null(outcome_range)){
        stop("Outcome range not specified for logistic outcome model:")
      }
      
      out_of_range <-
        colSums(x = data[, y_columns] > outcome_range[2], na.rm = TRUE) +
        colSums(x = data[, y_columns] < outcome_range[1], na.rm = TRUE)
      
      if(sum(out_of_range) > 0) {
        stop("Outcomes out of range [0, 1]: ",
             paste(names(which(out_of_range > 0)), collapse = ", "))
      }
    }
    
    
    # Check for binary data not in {0, 1}
    if(outcome_type == "binomial"){
      
      out_of_range <-
        apply(
          X = data[, y_columns],
          MARGIN = 2,
          FUN = function(x) sum(!x %in% c(0, 1, NA))
        )
      
      if(sum(out_of_range) > 0) {
        stop("Outcomes out of range {0, 1}: ",
             paste(names(which(out_of_range > 0)), collapse = ", "))
      }
    }
    
    if(outcome_type == "multinomial-binomial"){
      if(
        !all(
          sapply(
            X = data[ , head(y_columns, -1)],
            FUN = class
          ) == "factor"
        )
      ){
        stop("Intermediate multinomial outcomes must be type 'factor'")
      }
      
      
      out_of_range <-
        !(data[, tail(y_columns, 1)] %in% c(0, 1, NA))
      
      if(sum(out_of_range) > 0 ){
        stop("`outcome_type` == 'multinomial-binomial' with final outcome ",
             "values outside {0, 1} in rows ",
             paste(which(out_of_range), collapse = ", "))
      }
    } else{
      if(!is.null(absorbing_state)){
        stop("`absorbing_state` specified without ",
             "`outcome_type` == 'multinomial-binomial'")
      }
    }
    
    # Check for non-monotone missingness
    impute_columns <-
      non_monotone(
        data = data,
        y_columns = y_columns
      )
    
    n_non_monotone <- colSums(impute_columns)
    
    if(sum(n_non_monotone) > 0){
      
      if(is.null(impute_formulas)){
        stop("Missingness is not monotone and no imputation formula is supplied.")
      }
      
      if(!impute_model %in%
         c("gaussian", "binomial", "beta", "pmm", "multinomial")) {
        stop("Invalid specification for `impute_model`:",
             impute_model)
      }
      
      # Find which variables must be imputed: Ensure temporal order
      non_monotone_vars <- names(which(n_non_monotone > 0))
      
      imputation_rhs <- 
        sapply(
          X = impute_formulas,
          FUN =
            function(x) all.vars(update(x, . ~ 0))
        )
      
      imputation_rhs_missing <- 
        setdiff(x = non_monotone_vars, y = imputation_rhs)
      
      if(length(imputation_rhs_missing) > 0) {
        stop("No model specified for outcomes with non-monotone missingness: ",
             paste(imputation_rhs_missing, collapse = ", "))
      }
      
      # Subset formulas to variables with non-monotone missingness
      # Order in temporal sequence if needed
      impute_formulas <-
        impute_formulas[
          match(x = intersect(x = imputation_rhs, y = non_monotone_vars),
                table = y_columns)
        ]
    }
    
    return(
      list(
        impute_columns = impute_columns,
        impute_formulas = impute_formulas
      )
    )
  }




### tmle_get_formulas ##########################################################
# Check specified formulas: missingness formulas can be supplied as a
# right-hand-side formula - use y_columns to add appropriate outcome.
# Since constructing the TMLE algorithm requires constructing new variables,
# new variable names are added to the data.frame to avoid over-writing the
# original values (`..y1`, `..y2`, `..y3` instead of `y1`, `y2`, `y3`). The
# outcome formulas need to be adjusted accordingly.
tmle_get_formulas <-
  function(
    y_columns,
    inverse_weight_formulas,
    outcome_formulas
  ){
    
    for(i in 1:length(inverse_weight_formulas)){
      # Cut off LHS (if supplied) - use y_columns to find response
      inverse_weight_formulas[[i]] <-
        reformulate(
          termlabels =
            Reduce(paste, deparse(inverse_weight_formulas[[i]][[3]])),
          response =
            paste("!is.na(", y_columns[i], ")", collapse = "", sep = "")
        )
    }
    
    for(i in 1:length(outcome_formulas)){
      outcome_formulas[[i]] <-
        reformulate(
          termlabels =
            Reduce(paste, deparse(outcome_formulas[[i]][[3]])),
          response =
            paste0("..", y_columns[i])
        )
    }
    
    return(
      list(
        inverse_weight_formulas = inverse_weight_formulas,
        outcome_formulas = outcome_formulas
      )
    )
  }




### tmle_compute ###############################################################
# NOTE: this is an internal function, meant to be called by tmle(). It has been
# stripped down to speed up bootstrapping at the expense of error handling that
# is done by tmle().
tmle_compute <-
  function(
    data,
    y_columns, tx_column,
    propensity_score_formula,
    inverse_weight_formulas,
    outcome_formulas,
    outcome_type =
      c("gaussian",
        "logistic",
        "binomial")[1],
    outcome_range = NULL,
    absorbing_state = NULL,
    absorbing_outcome = NULL,
    impute = FALSE,
    impute_formulas = NULL,
    impute_model = 
      c("gaussian",
        "binomial",
        "beta",
        "pmm",
        "multinomial")[1],
    imputation_args = NULL,
    verbose = FALSE,
    max_abs_weight = 20
  ) {
    
    if(verbose) data_original <- data
    
    n_outcomes <- length(y_columns)
    
    if(impute){
      impute_columns <-
        non_monotone(data = data, y_columns = y_columns)
      
      imputed_outcomes <- 
        do.call(
          what = impute_outcomes_to_monotone,
          args =
            c(
              list(
                model = impute_model,
                data = data,
                impute_columns = impute_columns,
                impute_formulas = impute_formulas,
                absorbing_state = absorbing_state,
                y_columns = y_columns,
                verbose = verbose
              ),
              imputation_args
            )
        )
      
      if(verbose){
        data <- imputed_outcomes$imputed_data
      } else {
        data <- imputed_outcomes
      }
    } else {
      imputed_outcomes <- NULL
    }
    
    # Fit propensity score model & compute propensity score
    propensity_model <-
      glm(
        formula = propensity_score_formula,
        data = data,
        family = binomial
      )
    
    propensity_score <-
      (data[, tx_column] == 1)*propensity_model$fitted.values +
      (1 - (data[, tx_column] == 1))*(1 - propensity_model$fitted.values)
    
    # Compute inverse probability of censoring weights
    inverse_weights <-
      compute_inverse_weights(
        data = data,
        inverse_weight_formulas = inverse_weight_formulas,
        y_columns = y_columns,
        absorbing_state = absorbing_state,
        verbose = verbose
      )
    
    if(verbose){
      assign(x = "ipw", value = inverse_weights$ipw)
    } else {
      assign(x = "ipw", value = inverse_weights)
    }
    
    ipw <-
      ipw*kronecker(
        X = matrix(data = 1, ncol = n_outcomes),
        Y = 1/propensity_score
      )
    
    # Truncate Weights
    ipw <-
      apply(
        X = ipw,
        MARGIN = 2,
        FUN = function(x) pmin(x, max_abs_weight)
      )
    
    check_ipw <-
      sum(colSums(non_monotone(data = data, y_columns = y_columns)))
    
    if(check_ipw > 0) {
      stop("Error in IPW model.")
    }
    
    # Compute sequence of regression fits:
    regression_sequence <- list()
    
    # Create new columns for fits and IPW:
    data[, "..ipw"] <- NA
    
    fit_y_columns <- paste0("..", c(y_columns))
    data[, fit_y_columns] <- NA
    data[, fit_y_columns[n_outcomes]] <- data[, y_columns[n_outcomes]]
    
    if(outcome_type %in% "gaussian"){
      glm_family <- gaussian
    } else{
      glm_family <- quasibinomial
    }
    
    for(i in n_outcomes:1){
      uncensored <- 
        if(i > 1){
          if(is.null(absorbing_state)){
            which(!is.na(data[, y_columns[(i-1)]]))
          } else {
            which(
              !(is.na(data[, y_columns[(i-1)]]) |
                  data[, y_columns[(i-1)]] %in% absorbing_state)
            )
          }
        } else{
          1:nrow(data)
        }
      
      data$..ipw <- ipw[, i]
      
      outcome_regression <- 
        glm(
          formula = outcome_formulas[[i]],
          data = data[uncensored, ],
          family = glm_family,
          weights = ..ipw
        )
      
      if(verbose) regression_sequence[[i]] <- outcome_regression
      
      if(i > 1) {
        data[uncensored, fit_y_columns[(i - 1)]] <-
          predict(
            object = outcome_regression,
            newdata = data[uncensored, ],
            type = "response"
          )
        
        if(!is.null(absorbing_state)){
          absorbed <- which(data[, y_columns[(i-1)]] %in% absorbing_state)
          if(length(absorbed) > 0){
            data[absorbed, fit_y_columns[(i - 1)]] <- absorbing_outcome
          }
        }
      }
    }
    
    # Compute ATE
    y_tx_1 <-
      predict(
        object = outcome_regression,
        newdata = 
          within(
            data = data,
            expr = {eval(parse(text = paste0(tx_column, " = 1")))}
          ),
        type = "response"
      )
    
    y_tx_0 <-
      predict(
        object = outcome_regression,
        newdata = within(
          data = data,
          expr = {eval(parse(text = paste0(tx_column, " = 0")))}
        ),
        type = "response"
      )
    
    if(outcome_type == "logistic"){
      y_tx_1 <- 
        decompress_range(
          x = y_tx_1,
          min = outcome_range[1],
          max = outcome_range[2]
        )
      
      y_tx_0 <- 
        decompress_range(
          x = y_tx_0,
          min = outcome_range[1],
          max = outcome_range[2]
        )
    }
    
    ate <-
      mean(y_tx_1) - mean(y_tx_0)
    
    if(verbose){
      return(
        list(
          ate = ate,
          y_tx_1 = y_tx_1,
          y_tx_0 = y_tx_0,
          data = data,
          outcome_range = outcome_range,
          outcome_formulas = outcome_formulas,
          outcome_models = regression_sequence,
          ipw = ipw,
          inverse_weights = inverse_weights,
          inverse_weight_formulas = inverse_weight_formulas,
          propensity_model = propensity_model,
          propensity_score = propensity_score,
          imputed_outcomes = imputed_outcomes,
          impute_formulas = impute_formulas,
          imputation_args = imputation_args,
          data_original = data_original
        )
      )
    } else {
      return(ate)
    }
  }


### tmle_boot_wrap #############################################################
# wrapper to be used with boot()
tmle_boot_wrap <-
  function(
    data, indices = NULL, tmle_args
  ) {
    
    if(!is.null(indices)) {
      data <- data[indices, ]
    }
    
    do.call(
      what =
        function(x = data, ...) {
          tmle_compute(data = x, ...)
        },
      args = tmle_args
    )
  }





### tmle #######################################################################
tmle <-
  function(
    data,
    propensity_score_formula,
    inverse_weight_formulas,
    outcome_formulas,
    outcome_type =
      c("gaussian",
        "logistic",
        "binomial",
        "multinomial-binomial")[1],
    outcome_range = NULL,
    absorbing_state = NULL,
    absorbing_outcome = NULL,
    impute_formulas = NULL,
    impute_model = 
      c("gaussian",
        "binomial",
        "beta",
        "pmm",
        "multinomial")[1],
    imputation_args = NULL,
    ci = FALSE,
    verbose = FALSE,
    bootstrap_n = 10000,
    bootstrap_type = c("bca", "norm", "basic", "perc")[1],
    alpha = 0.05,
    ...
  ) {
    
    outcome_type <- tolower(outcome_type)
    impute_model <- tolower(impute_model)
    
    arg_list <- as.list(substitute(list(...)))[-1L]
    
    if(ci) {
      verbose <- FALSE # Override verbose when bootstrapping
      if(!bootstrap_type %in% c("bca", "norm", "basic", "perc")) {
        stop("Unrecognized bootstrap method: ", bootstrap_type)
      }
    }
    
    tmle_args <-
      tmle_get_args(
        data = data,
        propensity_score_formula = propensity_score_formula,
        inverse_weight_formulas = inverse_weight_formulas,
        outcome_formulas = outcome_formulas,
        impute_formulas = impute_formulas
      )
    
    
    # When absorbing state is specified, check for consistency and carry forward
    # absorbing state values
    if(!is.null(absorbing_state)){
      checked_data <-
        absorbing_state_check(
          data = data,
          y_columns = tmle_args$y_columns,
          outcome_type = outcome_type,
          absorbing_state = absorbing_state,
          absorbing_outcome = absorbing_outcome
        )
      
      environment(checked_data) <-
        environment(data)
      
      data <- checked_data
    }
    
    precheck_results <-
      tmle_precheck(
        data = data,
        x_columns = tmle_args$x_columns,
        tx_column = tmle_args$tx_column,
        y_columns = tmle_args$y_columns,
        imputation_x = tmle_args$imputation_x,
        impute_formulas = impute_formulas,
        impute_model = impute_model,
        outcome_type = outcome_type,
        outcome_range = outcome_range,
        absorbing_state = absorbing_state
      )
    
    impute <-
      sum(colSums(non_monotone(data = data, y_columns = tmle_args$y_columns)))
    
    if(outcome_type == "logistic"){
      data[, tmle_args$y_columns] <-
        apply(
          X = data[, tmle_args$y_columns],
          MARGIN = 2,
          FUN = compress_range,
          min = outcome_range[1],
          max = outcome_range[2],
          enforce_range = TRUE
        )
    }
    
    # Construct left-hand-side of formulas
    tmle_formulas <-
      tmle_get_formulas(
        y_columns = tmle_args$y_columns,
        inverse_weight_formulas = inverse_weight_formulas,
        outcome_formulas = outcome_formulas
      )
    
    # Re-set environments for formulas updated within functions
    for (i in 1:length(tmle_formulas$inverse_weight_formulas))
      environment(tmle_formulas$inverse_weight_formulas[[i]]) <-
      environment(inverse_weight_formulas[[i]])
    
    for (i in 1:length(tmle_formulas$outcome_formulas))
      environment(tmle_formulas$outcome_formulas[[i]]) <-
      environment(outcome_formulas[[i]])
    
    for (i in 1:length(precheck_results$impute_formulas))
      environment(precheck_results$impute_formulas[[i]]) <-
      environment(impute_formulas[[i]])
    
    arg_list <-
      c(
        list(
          y_columns = tmle_args$y_columns,
          tx_column = tmle_args$tx_column,
          propensity_score_formula = propensity_score_formula,
          inverse_weight_formulas =
            tmle_formulas$inverse_weight_formulas,
          outcome_formulas =
            tmle_formulas$outcome_formulas,
          outcome_type = outcome_type,
          outcome_range = outcome_range,
          absorbing_state = absorbing_state,
          absorbing_outcome = absorbing_outcome,
          impute = impute,
          impute_formulas = impute_formulas,
          imputation_args = imputation_args,
          impute_model = impute_model,
          verbose = verbose
        ),
        # impute_model, imputation_args, ...
        arg_list
      )
    
    if(ci){
      
      tmle_boot <-
        boot(
          data = data,
          statistic = tmle_boot_wrap,
          R = bootstrap_n,
          tmle_args =
            arg_list
        )
      
      tmle_boot_ci <-
        boot.ci(
          boot.out = tmle_boot,
          conf  = (1 - alpha),
          type = bootstrap_type
        )
      
      lcl_ucl <-
        tail(x = tmle_boot_ci[bootstrap_type][[1]][1,], 2)
      
      return(
        c(
          estimate = tmle_boot$t0,
          se = sd(tmle_boot$t),
          lcl = lcl_ucl[1],
          ucl = lcl_ucl[2],
          alpha = alpha
        )
      )
    } else {
      return(
        tmle_boot_wrap(
          data = data,
          tmle_args = arg_list
        )
      )
    }
  }