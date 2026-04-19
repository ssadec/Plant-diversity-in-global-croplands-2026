####STable 20####
rm(list = ls())

# =========================================================
# Stable 20
# Purpose: Fit two-step robust meta-analytic models for the
# five total-response datasets using shared-control covariance
# when valid, with fallback to diagonal sampling variances.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
})

# =========================================================
# 1. Configuration
# =========================================================
cfg <- list(
  project_dir = "/Users/yusha",
  total_names_file = "total_names.txt",
  output_file = "Stable20_TwoStep_Robust_Meta_Results.csv",
  min_rows_required = 5L
)

category_labels <- c(
  "Total_response_of_predator_performance",
  "Total_response_of_parasitoid_performance",
  "Total_response_of_predator_parasitoid_performance",
  "Total_response_of_herbivore_performance",
  "Total_response_of_crop_performance"
)

# =========================================================
# 2. Helper functions
# =========================================================
read_name_vector <- function(file_path) gsub('^"|"$', "", as.character(read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)[, 1]))

calc_shared_block <- function(x) {
  cpart <- x$CONTROL_SD[1]^2 / (x$CONTROL_R[1] * x$CONTROL_M[1]^2)
  v <- matrix(cpart, nrow = nrow(x), ncol = nrow(x))
  diag(v) <- x$vi
  v
}

prepare_dataset <- function(dat) {
  dat <- dat[!is.na(dat$yi), ]
  dat <- dat[is.finite(dat$yi) & is.finite(dat$vi) & dat$vi > 0, ]
  dat
}

build_common_ids <- function(dat) {
  dat %>%
    group_by(S_ID) %>%
    mutate(Common_ID = match(interaction(CONTROL_M, CONTROL_SD, CONTROL_R, drop = TRUE), unique(interaction(CONTROL_M, CONTROL_SD, CONTROL_R, drop = TRUE)))) %>%
    ungroup() %>%
    mutate(Common_ID_global = match(interaction(S_ID, Common_ID, drop = TRUE), unique(interaction(S_ID, Common_ID, drop = TRUE))))
}

can_use_shared_v <- function(dat) {
  check_common <- dat %>%
    group_by(Common_ID_global) %>%
    summarise(n_CONTROL_M = n_distinct(CONTROL_M), n_CONTROL_SD = n_distinct(CONTROL_SD), n_CONTROL_R = n_distinct(CONTROL_R), .groups = "drop")
  all(check_common$n_CONTROL_M == 1) && all(check_common$n_CONTROL_SD == 1) && all(check_common$n_CONTROL_R == 1)
}

build_v_matrix <- function(dat) {
  V <- NULL
  
  if (can_use_shared_v(dat)) {
    V_try <- try(as.matrix(metafor::bldiag(lapply(split(dat, dat$Common_ID_global), calc_shared_block))), silent = TRUE)
    
    if (!inherits(V_try, "try-error")) {
      V_ok <- isSymmetric(V_try) && !any(is.na(V_try)) && !any(is.infinite(V_try))
      
      if (V_ok) {
        eig_min <- try(min(eigen(V_try, symmetric = TRUE, only.values = TRUE)$values), silent = TRUE)
        if (!inherits(eig_min, "try-error") && is.finite(eig_min) && eig_min > 0) V <- V_try
      }
    }
  }
  
  if (is.null(V)) V <- diag(dat$vi)
  V
}

fit_two_step_robust <- function(dat, V) {
  mod_mlfe <- try(rma.mv(yi = yi, V = V, random = ~1 | S_ID, method = "ML", test = "t", dfs = "contain", data = dat), silent = TRUE)
  if (inherits(mod_mlfe, "try-error")) return(NULL)
  
  mod_mlfe_rve <- try(robust(mod_mlfe, cluster = dat$S_ID, adjust = TRUE, clubSandwich = TRUE), silent = TRUE)
  if (inherits(mod_mlfe_rve, "try-error")) return(NULL)
  
  mod_mlfe_rve
}

extract_robust_summary <- function(fit_obj, category_name, dat) {
  if (is.null(fit_obj$b) || length(fit_obj$b) == 0 || is.null(fit_obj$se) || length(fit_obj$se) == 0) return(NULL)
  
  beta_val <- as.numeric(fit_obj$b[1])
  se_val <- as.numeric(fit_obj$se[1])
  tval_val <- if (!is.na(se_val) && se_val != 0) beta_val / se_val else NA
  pval_val <- if (!is.null(fit_obj$pval) && length(fit_obj$pval) > 0) as.numeric(fit_obj$pval[1]) else NA
  cilb_val <- if (!is.null(fit_obj$ci.lb) && length(fit_obj$ci.lb) > 0) as.numeric(fit_obj$ci.lb[1]) else NA
  ciub_val <- if (!is.null(fit_obj$ci.ub) && length(fit_obj$ci.ub) > 0) as.numeric(fit_obj$ci.ub[1]) else NA
  
  data.frame(
    Category = category_name,
    k = fit_obj$k,
    n_cluster = length(unique(dat$S_ID)),
    beta = beta_val,
    SE = se_val,
    tval = tval_val,
    pval = pval_val,
    CIlb = cilb_val,
    CIub = ciub_val,
    stringsAsFactors = FALSE
  )
}

# =========================================================
# 3. Dataset registry
# =========================================================
total_names <- read_name_vector(file.path(cfg$project_dir, cfg$total_names_file))

if (length(total_names) != length(category_labels)) stop("The number of entries in total_names.txt does not match the number of category labels. Please verify the mapping.")

# =========================================================
# 4. Main analysis loop
# =========================================================
res_all <- data.frame()

for (i in seq_along(total_names)) {
  dat <- try(read.table(file.path(cfg$project_dir, paste0(total_names[i], ".txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(dat, "try-error") || nrow(dat) < cfg$min_rows_required) next
  
  dat <- prepare_dataset(dat)
  if (nrow(dat) < cfg$min_rows_required) next
  
  dat <- build_common_ids(dat)
  V <- build_v_matrix(dat)
  fit_obj <- fit_two_step_robust(dat, V)
  if (is.null(fit_obj)) next
  
  tmp <- extract_robust_summary(fit_obj, category_labels[i], dat)
  if (!is.null(tmp)) res_all <- rbind(res_all, tmp)
}

# =========================================================
# 5. Export results
# =========================================================
print(res_all)
write.csv(res_all, file.path(cfg$project_dir, cfg$output_file), row.names = FALSE)

cat("Finished Stable 20.\n")
cat("Results written to: ", file.path(cfg$project_dir, cfg$output_file), "\n")

sessionInfo()
