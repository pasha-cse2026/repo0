# ==============================================================================
# declaredesign_audit.R
# Protocol for Empirical Validation, Paired Testing, and DeclareDesign Diagnosis
# ==============================================================================

# Required libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(estimatr)
  if (requireNamespace("DeclareDesign", quietly = TRUE)) {
    library(DeclareDesign)
  }
})

# ------------------------------------------------------------------------------
# 1. Load Data and Enforce Integrity Constraints
# ------------------------------------------------------------------------------
csv_file <- "raw_trial_records.csv"

if (!file.exists(csv_file)) {
  message("Note: 'raw_trial_records.csv' not found. Please provide the file to run empirical validation.")
} else {
  dat <- read_csv(csv_file, show_col_types = FALSE)
  
  required_cols <- c(
    "task_id", "model", "interference", "placement",
    "prompt_format", "token_budget", "representation",
    "raw_output", "strict_em", "substring_pass"
  )
  missing <- setdiff(required_cols, names(dat))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Check for NA values in key identification and scoring columns
  critical_fields <- c("task_id", "model", "interference", "prompt_format", "strict_em", "substring_pass")
  if (anyNA(dat[critical_fields])) {
    stop("Critical fields contain NA values; resolve missing data before proceeding.")
  }
  
  # Enforce binary scoring
  if (any(!dat$strict_em %in% 0:1) || any(!dat$substring_pass %in% 0:1)) {
    stop("Scores in 'strict_em' and 'substring_pass' must be binary (0 or 1).")
  }
  
  # Mathematical Invariant: Strict exact match is a strict subset of substring containment
  violations <- which(dat$strict_em == 1 & dat$substring_pass == 0)
  if (length(violations) > 0) {
    stop(sprintf("Invariant violated at rows %s: strict_em=1 but substring_pass=0.", paste(violations, collapse = ", ")))
  }
  
  # Condition column definition
  condition_cols <- c("model", "interference", "placement", "prompt_format", "token_budget", "representation")
  
  # Uniqueness check per task x condition
  duplicates <- dat %>%
    count(task_id, across(all_of(condition_cols))) %>%
    filter(n > 1)
  
  if (nrow(duplicates) > 0) {
    warning(sprintf("Found %d duplicate task-condition keys. Check 'run_id' for multiple runs.", nrow(duplicates)))
  }
  
  # ----------------------------------------------------------------------------
  # 2. Wilson Score Binomial Confidence Intervals
  # ----------------------------------------------------------------------------
  wilson_ci <- function(k, n, conf_level = 0.95) {
    if (n == 0) return(c(rate = NA, lower = NA, upper = NA))
    p <- k / n
    z <- qnorm(1 - (1 - conf_level) / 2)
    den <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / den
    half <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2))) / den
    c(rate = p, lower = max(0, center - half), upper = min(1, center + half))
  }
  
  cat("\n=== 1. Empirical Marginals by Exact Condition Cell ===\n")
  marginal_summary <- dat %>%
    group_by(across(all_of(condition_cols))) %>%
    summarise(
      N = n(),
      Strict_K = sum(strict_em),
      Sub_K = sum(substring_pass),
      Strict_Rate = wilson_ci(sum(strict_em), n())[["rate"]],
      Strict_95CI_Lo = wilson_ci(sum(strict_em), n())[["lower"]],
      Strict_95CI_Hi = wilson_ci(sum(strict_em), n())[["upper"]],
      Sub_Rate = wilson_ci(sum(substring_pass), n())[["rate"]],
      Sub_95CI_Lo = wilson_ci(sum(substring_pass), n())[["lower"]],
      Sub_95CI_Hi = wilson_ci(sum(substring_pass), n())[["upper"]],
      .groups = "drop"
    )
  print(as.data.frame(marginal_summary))
  
  # ----------------------------------------------------------------------------
  # 3. Same-Output Paired Metric Comparison (McNemar & Paired Difference)
  # ----------------------------------------------------------------------------
  cat("\n=== 2. Paired Scoring Gap Analysis (1.5B Decoy, Concat) ===\n")
  metric_arm <- dat %>%
    filter(model == "1.5B", interference == "Decoy", prompt_format == "Concat")
  
  if (nrow(metric_arm) > 0) {
    b <- sum(metric_arm$substring_pass == 1 & metric_arm$strict_em == 0) # Substring only
    c <- sum(metric_arm$substring_pass == 0 & metric_arm$strict_em == 1) # Strict only (should be 0)
    a <- sum(metric_arm$substring_pass == 1 & metric_arm$strict_em == 1) # Both pass
    d <- sum(metric_arm$substring_pass == 0 & metric_arm$strict_em == 0) # Both fail
    
    cat(sprintf("2x2 Table: [Both Pass: %d, Substring-Only: %d, Strict-Only: %d, Both Fail: %d]\n", a, b, c, d))
    cat(sprintf("Paired Gap (Mean Difference): %.4f\n", mean(metric_arm$substring_pass - metric_arm$strict_em)))
    
    if (b + c > 0) {
      mcnemar_exact <- binom.test(b, b + c, p = 0.5)
      cat(sprintf("Exact McNemar Test (Two-Sided): p-value = %g\n", mcnemar_exact$p.value))
    }
    
    # Paired difference-in-proportions with task-clustered robust SEs
    fit_paired <- lm_robust(I(substring_pass - strict_em) ~ 1, data = metric_arm, clusters = task_id)
    print(summary(fit_paired))
  } else {
    message("1.5B Decoy Concat arm not present for paired metric analysis.")
  }
  
  # ----------------------------------------------------------------------------
  # 4. Format-Harmonized Contrasts (Avoiding Confounding)
  # ----------------------------------------------------------------------------
  cat("\n=== 3. Format-Harmonized Scale x Decoy Contrast (Schema Arm) ===\n")
  schema_data <- dat %>% filter(prompt_format == "Schema")
  
  if (n_distinct(schema_data$model) == 2 && n_distinct(schema_data$interference) >= 2) {
    fit_scale_interaction <- lm_robust(
      strict_em ~ model * interference,
      data = schema_data,
      clusters = task_id
    )
    cat("Linear Probability Model on Schema Arm (Unconfounded Scale Interaction):\n")
    print(summary(fit_scale_interaction))
  } else {
    message("Matched Schema data not fully populated for both models.")
  }
  
  cat("\n=== 4. Prompt Remediation Contrast (1.5B Decoy: Schema vs Concat) ===\n")
  remediation_data <- dat %>% filter(model == "1.5B", interference == "Decoy")
  if (n_distinct(remediation_data$prompt_format) == 2) {
    fit_remediation <- lm_robust(strict_em ~ prompt_format, data = remediation_data, clusters = task_id)
    print(summary(fit_remediation))
  }
}

# ------------------------------------------------------------------------------
# 5. DeclareDesign Simulation Model (Design Diagnosis)
# ------------------------------------------------------------------------------
if (requireNamespace("DeclareDesign", quietly = TRUE)) {
  cat("\n=== 5. Formal DeclareDesign Specification & Diagnosis ===\n")
  
  # Step 1: Population of 100 task instances drawn from task grammar
  study_pop <- declare_population(
    task_id = 1:100,
    u_i = rnorm(100, mean = 0, sd = 1.0) # Latent task difficulty
  )
  
  # Step 2: Cartesian Expansion (100 tasks x 4 cells = 400 rows)
  study_expand <- declare_step(
    fabricate,
    data = tidyr::crossing(
      study_pop,
      tibble(
        model = factor(c("1.5B", "3B"), levels = c("1.5B", "3B")),
        interference = factor(c("Clean", "Decoy"), levels = c("Clean", "Decoy"))
      )
    )
  )
  
  # Step 3: Outcomes under greedy decoding (deterministic given task u_i)
  study_outcomes <- declare_step(
    fabricate,
    mu = case_when(
      model == "1.5B" & interference == "Clean" ~  3.18, # ~96%
      model == "1.5B" & interference == "Decoy" ~ -1.90, # ~13%
      model == "3B"   & interference == "Clean" ~  3.89, # ~98%
      model == "3B"   & interference == "Decoy" ~  3.48  # ~97%
    ),
    score = mu + u_i,
    strict_em = as.integer(score > 0),
    sub_shift = if_else(model == "1.5B" & interference == "Decoy", 1.15, 0.20),
    substring_pass = as.integer(score + sub_shift > 0),
    metric_gap = substring_pass - strict_em
  )
  
  # Step 4: Estimands & Inquiries
  inq_scale_int <- declare_inquiry(
    mean((strict_em[model == "3B" & interference == "Decoy"] - 
          strict_em[model == "3B" & interference == "Clean"]) -
         (strict_em[model == "1.5B" & interference == "Decoy"] - 
          strict_em[model == "1.5B" & interference == "Clean"])),
    label = "Risk_Diff_Scale_Interaction"
  )
  
  inq_paired_gap <- declare_inquiry(
    mean(metric_gap[model == "1.5B" & interference == "Decoy"]),
    label = "Paired_Metric_Gap"
  )
  
  # Step 5: Estimators
  est_scale_int <- declare_estimator(
    strict_em ~ model * interference,
    model = lm_robust,
    clusters = task_id,
    term = "model3B:interferenceDecoy",
    inquiry = inq_scale_int,
    label = "LPM_Clustered_Interaction"
  )
  
  est_paired_gap <- declare_estimator(
    metric_gap ~ 1,
    model = lm_robust,
    subset = (model == "1.5B" & interference == "Decoy"),
    inquiry = inq_paired_gap,
    label = "Paired_Difference_Estimator"
  )
  
  declared_design <- study_pop + study_expand + study_outcomes + 
                     inq_scale_int + inq_paired_gap + 
                     est_scale_int + est_paired_gap
  
  cat("DeclareDesign model successfully declared.\n")
  cat("Run 'diagnose_design(declared_design, sims = 100)' to simulate operating characteristics.\n")
}