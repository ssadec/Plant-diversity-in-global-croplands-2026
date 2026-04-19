#### Herbivore-crop paired data cleaning + Stable 12 + SuppTab 13 + Stable 15 + Stable 17 ####
rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(metafor)
  library(Matrix)
  library(corpcor)
  library(nlme)
  library(piecewiseSEM)
  library(openxlsx)
})

#### 1. Paths ####
prefix <- "/Users/yusha/"
path   <- "/Users/yusha/"

#### 2. Shared helpers ####
safe_num <- function(x) suppressWarnings(as.numeric(x))

read_txt_safe <- function(file) {
  read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "", fill = TRUE, blank.lines.skip = FALSE)
}

write_txt_safe <- function(x, file) {
  write.table(x, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
}

round_for_group <- function(x, digits) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA, round(x, digits))
}

make_common_key <- function(m, sd, r) {
  interaction(m, sd, r, drop = TRUE, lex.order = TRUE)
}

fit_rma_mv_try <- function(yi_col, V_input, dat, method = "ML") {
  random_terms <- list(~1 | S_ID, ~1 | ID_global)
  for (opt in c("BFGS", "nlminb", "optim")) {
    fit <- try(
      rma.mv(as.formula(paste0(yi_col, " ~ 1")), V = V_input, random = random_terms, data = dat, test = "t", method = method, control = list(optimizer = opt)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) return(fit)
  }
  "try-error"
}

fit_lme_safe <- function(formula_obj, vi_col, data_obj, nested = TRUE) {
  random_formula <- if (nested) ~1 | S_ID/ID else ~1 | S_ID
  fit <- try(
    lme(fixed = formula_obj, random = random_formula, weights = varFixed(as.formula(paste0("~", vi_col))), data = data_obj, na.action = na.omit, control = lmeControl(sigma = 1, returnObject = TRUE, opt = "optim")),
    silent = TRUE
  )
  fit
}

get_col_safe <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

#### 3. Part I: Data cleaning for herbivore-crop paired datasets ####
filenames <- dir(path = prefix, pattern = "^0402Paired data for herbivore and crop.*\\.xlsx$")
sheets <- c("herbivore", "crop")

dt_total0 <- list()
paired_names_list <- list()

for (i in seq_along(filenames)) {
  file_path <- file.path(prefix, filenames[i])
  pair_id   <- paste0("herbivore_vs_crop_", i)
  dt_new0   <- NULL
  
  for (j in seq_along(sheets)) {
    dt0 <- tryCatch(read_excel(file_path, sheet = sheets[j], range = "A1:AG291"), error = function(e) NULL)
    if (is.null(dt0) || nrow(dt0) == 0) next
    
    colnames(dt0) <- make.names(abbreviate(colnames(dt0), minlength = 4))
    
    dt0 <- dt0 %>%
      filter(!is.na(CONTROL_M), is.finite(CONTROL_M), CONTROL_M > 0,
             !is.na(TREATMENT_M), is.finite(TREATMENT_M), TREATMENT_M > 0,
             !is.na(CONTROL_R), is.finite(CONTROL_R), CONTROL_R > 0,
             !is.na(TREATMENT_R), is.finite(TREATMENT_R), TREATMENT_R > 0,
             !is.na(S_ID))
    
    if (nrow(dt0) == 0) next
    
    dt0$CONTROL_SD <- replace(dt0$CONTROL_SD, !is.na(dt0$CONTROL_SD) & round(dt0$CONTROL_M / dt0$CONTROL_SD, 5) == 10, NA)
    dt0$TREATMENT_SD <- replace(dt0$TREATMENT_SD, is.na(dt0$CONTROL_SD), NA)
    
    ESData <- dt0 %>%
      mutate(cv_Control = na_if(CONTROL_SD / CONTROL_M, Inf),
             cv_Treatment = na_if(TREATMENT_SD / TREATMENT_M, Inf))
    
    ESData <- cv_avg(x = TREATMENT_M, sd = TREATMENT_SD, n = TREATMENT_R, group = S_ID, label = "1", data = ESData)
    ESData <- cv_avg(x = CONTROL_M, sd = CONTROL_SD, n = CONTROL_R, group = S_ID, label = "2", data = ESData)
    
    if (j == 1) {
      ESData <- ESData %>%
        mutate(
          cv2_trea_new1 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2),
          cv2_cont_new1 = if_else(is.na(cv_Control),   b_CV2_2, cv_Control^2),
          yi1 = lnrr_laj(m1 = TREATMENT_M, m2 = CONTROL_M, cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R),
          vi1 = v_lnrr_laj(cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R)
        )
      
      ESData$Trtmnt.s1 <- ESData$TREATMENT_SD
      ESData$Cntrl.s1  <- ESData$CONTROL_SD
      
      idx_t1 <- which(is.na(ESData$Trtmnt.s1))
      if (length(idx_t1) > 0) ESData[idx_t1, "Trtmnt.s1"] <- ESData[idx_t1, "cv2_trea_new1"] * ESData[idx_t1, "TREATMENT_M"]
      
      idx_c1 <- which(is.na(ESData$Cntrl.s1))
      if (length(idx_c1) > 0) ESData[idx_c1, "Cntrl.s1"] <- ESData[idx_c1, "cv2_cont_new1"] * ESData[idx_c1, "CONTROL_M"]
      
      dt_new0 <- ESData
    }
    
    if (j == 2) {
      ESData <- ESData %>%
        mutate(
          cv2_trea_new2 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2),
          cv2_cont_new2 = if_else(is.na(cv_Control),   b_CV2_2, cv_Control^2),
          yi2 = lnrr_laj(m1 = TREATMENT_M, m2 = CONTROL_M, cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R),
          vi2 = v_lnrr_laj(cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R)
        )
      
      ESData$Trtmnt.s2 <- ESData$TREATMENT_SD
      ESData$Cntrl.s2  <- ESData$CONTROL_SD
      
      idx_t2 <- which(is.na(ESData$Trtmnt.s2))
      if (length(idx_t2) > 0) ESData[idx_t2, "Trtmnt.s2"] <- ESData[idx_t2, "cv2_trea_new2"] * ESData[idx_t2, "TREATMENT_M"]
      
      idx_c2 <- which(is.na(ESData$Cntrl.s2))
      if (length(idx_c2) > 0) ESData[idx_c2, "Cntrl.s2"] <- ESData[idx_c2, "cv2_cont_new2"] * ESData[idx_c2, "CONTROL_M"]
      
      if (!is.null(dt_new0)) dt_new0 <- data.frame(dt_new0, ESData[, c("yi2", "vi2", "Trtmnt.s2", "Cntrl.s2", "cv2_trea_new2", "cv2_cont_new2"), drop = FALSE])
    }
  }
  
  if (!is.null(dt_new0)) {
    if (all(c("yi1", "yi2") %in% names(dt_new0))) {
      dt_new0 <- dt_new0 %>% filter(!is.na(yi1), is.finite(yi1), !is.na(yi2), is.finite(yi2))
    }
    if (nrow(dt_new0) == 0) next
    
    sid_tab <- table(dt_new0$S_ID)
    dt_new0$EsID <- unlist(lapply(as.numeric(sid_tab), seq_len), use.names = FALSE)
    
    if ("NCP_N" %in% names(dt_new0)) dt_new0$added_species_number <- log2(dt_new0$NCP_N)
    
    dt_new0$pair_id <- pair_id
    dt_new0$response1 <- "Herbivore"
    dt_new0$response2 <- "Crop"
    dt_new0$trophic_interaction <- "herbivore_vs_crop"
    
    dt_total0[[length(dt_total0) + 1]] <- dt_new0
    paired_names_list[[length(paired_names_list) + 1]] <- data.frame(Pair_ID = pair_id, Trophic_interaction = "herbivore_vs_crop", Response1 = "Herbivore", Response2 = "Crop", Source_file = filenames[i], stringsAsFactors = FALSE)
    
    write_txt_safe(dt_new0, paste0(prefix, pair_id, ".txt"))
  }
}

dt_total0 <- if (length(dt_total0) > 0) bind_rows(dt_total0) else NULL
paired_names <- if (length(paired_names_list) > 0) bind_rows(paired_names_list) else NULL

write_txt_safe(paired_names, paste0(prefix, "paired_names.txt"))
write_txt_safe(dt_total0, paste0(prefix, "herbivore_vs_crop.txt"))

#### 4. Stable 12 ####
paired_names <- read_txt_safe(paste0(prefix, "paired_names.txt"))
round_digits_m <- 4
round_digits_sd <- 4
round_digits_r <- 0
non_pd_list <- character(0)

build_shared_V_st12 <- function(dat, vi_col, sd_col) {
  calc.v <- function(x) {
    cpart <- x[[sd_col]][1]^2 / (x$CONTROL_R[1] * x$CONTROL_M[1]^2)
    v <- matrix(cpart, nrow = nrow(x), ncol = nrow(x))
    diag(v) <- x[[vi_col]]
    v
  }
  
  out <- list(V = NULL, V_type = NA_character_, min_eigen_before_fix = NA_real_, ok = FALSE, msg = NULL)
  
  V_try <- try(as.matrix(bdiag(lapply(split(dat, dat$Common_ID_global), calc.v))), silent = TRUE)
  if (inherits(V_try, "try-error")) { out$msg <- "failed building shared-control V"; return(out) }
  if (!isSymmetric(V_try) || any(is.na(V_try)) || any(is.infinite(V_try))) { out$msg <- "invalid shared-control V"; return(out) }
  
  eig_vals <- try(eigen(V_try, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (inherits(eig_vals, "try-error")) { out$msg <- "failed eigen decomposition"; return(out) }
  
  out$min_eigen_before_fix <- min(eig_vals)
  
  if (!is.finite(out$min_eigen_before_fix) || out$min_eigen_before_fix <= 0) {
    V_fix <- try(make.positive.definite(V_try), silent = TRUE)
    if (inherits(V_fix, "try-error")) { out$msg <- "shared-control V not PD and PD-fix failed"; return(out) }
    out$V <- V_fix
    out$V_type <- "shared_control_pd_fix"
    out$ok <- TRUE
    out$msg <- "shared-control V not PD, fixed before fitting"
  } else {
    out$V <- V_try
    out$V_type <- "shared_control"
    out$ok <- TRUE
  }
  
  out
}

run_one_response_st12 <- function(dat0, yi_col, vi_col, sd_col, category_name, pair_id, trophic_interaction) {
  dat <- dat0 %>%
    filter(!is.na(.data[[yi_col]]), is.finite(.data[[yi_col]]),
           !is.na(.data[[vi_col]]), is.finite(.data[[vi_col]]), .data[[vi_col]] > 0,
           !is.na(.data[[sd_col]]), is.finite(.data[[sd_col]]), .data[[sd_col]] >= 0,
           !is.na(CONTROL_M), is.finite(CONTROL_M), CONTROL_M > 0,
           !is.na(CONTROL_R), is.finite(CONTROL_R), CONTROL_R > 0,
           !is.na(S_ID), !is.na(ID)) %>%
    mutate(CONTROL_M_grp = round_for_group(CONTROL_M, round_digits_m),
           CONTROL_SD_grp = round_for_group(.data[[sd_col]], round_digits_sd),
           CONTROL_R_grp = round_for_group(CONTROL_R, round_digits_r)) %>%
    group_by(S_ID) %>%
    mutate(Common_ID = match(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp), unique(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp)))) %>%
    ungroup() %>%
    mutate(Common_ID_global = match(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE), unique(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE))),
           ID_global = interaction(S_ID, ID, drop = TRUE, lex.order = TRUE))
  
  if (nrow(dat) == 0) {
    return(data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = category_name, Number_of_observations = 0, Number_of_studies = NA, V_type = NA, Min_eigen_before_fix = NA, Effect_size = NA, t_value = NA, P_value = NA, df = NA, CI_lb = NA, CI_ub = NA, stringsAsFactors = FALSE))
  }
  
  check_common <- dat %>% group_by(Common_ID_global) %>% summarise(n_CONTROL_M_grp = n_distinct(CONTROL_M_grp), n_CONTROL_SD_grp = n_distinct(CONTROL_SD_grp), n_CONTROL_R_grp = n_distinct(CONTROL_R_grp), .groups = "drop")
  use_V <- nrow(check_common) > 0 && all(check_common$n_CONTROL_M_grp == 1) && all(check_common$n_CONTROL_SD_grp == 1) && all(check_common$n_CONTROL_R_grp == 1)
  
  dat.fit <- NULL
  V_type <- NA_character_
  min_eigen_before_fix <- NA_real_
  
  if (nrow(dat) == 1) {
    dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "z"), silent = TRUE)
    if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi_single_obs"
  }
  
  if ((is.null(dat.fit) || inherits(dat.fit, "try-error")) && nrow(dat) > 1 && use_V) {
    V_info <- build_shared_V_st12(dat, vi_col = vi_col, sd_col = sd_col)
    min_eigen_before_fix <- V_info$min_eigen_before_fix
    if (!is.null(V_info$msg) && grepl("not PD|invalid|failed", V_info$msg, ignore.case = TRUE)) non_pd_list <<- c(non_pd_list, paste0(pair_id, " | ", category_name, " : ", V_info$msg))
    
    if (isTRUE(V_info$ok)) {
      dat.fit <- fit_rma_mv_try(yi_col, V_info$V, dat, method = "ML")
      if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_try(yi_col, V_info$V, dat, method = "REML")
      if (!inherits(dat.fit, "try-error")) V_type <- V_info$V_type
    }
  }
  
  if (is.null(dat.fit) || inherits(dat.fit, "try-error")) {
    if (length(unique(dat$S_ID)) == 1) {
      dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "t"), silent = TRUE)
      if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi_rma"
    } else {
      dat.fit <- fit_rma_mv_try(yi_col, dat[[vi_col]], dat, method = "ML")
      if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_try(yi_col, dat[[vi_col]], dat, method = "REML")
      if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi"
    }
  }
  
  if (is.null(dat.fit) || inherits(dat.fit, "try-error")) {
    return(data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = category_name, Number_of_observations = nrow(dat), Number_of_studies = length(unique(dat$S_ID)), V_type = NA, Min_eigen_before_fix = min_eigen_before_fix, Effect_size = NA, t_value = NA, P_value = NA, df = NA, CI_lb = NA, CI_ub = NA, stringsAsFactors = FALSE))
  }
  
  tmp_sum <- coef(summary(dat.fit))
  t_value_out <- if ("tval" %in% colnames(tmp_sum)) tmp_sum$tval[1] else if ("zval" %in% colnames(tmp_sum)) tmp_sum$zval[1] else NA
  p_value_out <- if ("pval" %in% colnames(tmp_sum)) tmp_sum$pval[1] else NA
  df_out <- if ("df" %in% colnames(tmp_sum)) tmp_sum$df[1] else NA
  ci_lb_out <- if (!is.null(dat.fit$ci.lb)) dat.fit$ci.lb[1] else NA
  ci_ub_out <- if (!is.null(dat.fit$ci.ub)) dat.fit$ci.ub[1] else NA
  
  data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = category_name, Number_of_observations = nrow(dat), Number_of_studies = length(unique(dat$S_ID)), V_type = V_type, Min_eigen_before_fix = min_eigen_before_fix, Effect_size = as.numeric(dat.fit$b[1]), t_value = t_value_out, P_value = p_value_out, df = df_out, CI_lb = ci_lb_out, CI_ub = ci_ub_out, stringsAsFactors = FALSE)
}

res12 <- list()
for (i in seq_len(nrow(paired_names))) {
  dat0 <- read_txt_safe(paste0(prefix, paired_names$Pair_ID[i], ".txt"))
  res12[[length(res12) + 1]] <- run_one_response_st12(dat0, "yi1", "vi1", "Cntrl.s1", paired_names$Response1[i], paired_names$Pair_ID[i], paired_names$Trophic_interaction[i])
  res12[[length(res12) + 1]] <- run_one_response_st12(dat0, "yi2", "vi2", "Cntrl.s2", paired_names$Response2[i], paired_names$Pair_ID[i], paired_names$Trophic_interaction[i])
}
res12 <- bind_rows(res12)
openxlsx::write.xlsx(res12, paste0(path, "STable12_optimized.xlsx"), rowNames = FALSE)
if (length(non_pd_list) > 0) writeLines(non_pd_list, con = paste0(path, "STable12_non_pd_log_optimized.txt"))

#### 5. SuppTab 13: formal piecewise SEM ####
paired_names <- read_txt_safe(paste0(prefix, "paired_names.txt"))
res13 <- list()

extract_psem_rows <- function(sem_obj, pair_id, n_obs, n_study, sem_type, category_var = NULL, category = NULL) {
  cf <- try(coefs(sem_obj), silent = TRUE)
  if (inherits(cf, "try-error") || is.null(cf) || nrow(cf) == 0) return(NULL)
  
  pred_col <- get_col_safe(cf, c("Predictor", "predictor", "Predictor.name"))
  resp_col <- get_col_safe(cf, c("Response", "response", "Response.name"))
  est_col  <- get_col_safe(cf, c("Estimate", "estimate"))
  se_col   <- get_col_safe(cf, c("Std.Error", "SE", "Std.Error.of.Estimate", "std.error"))
  p_col    <- get_col_safe(cf, c("P.Value", "P-value", "p.value", "P"))
  t_col    <- get_col_safe(cf, c("Crit.Value", "t-value", "z-value", "Crit.Value.", "statistic"))
  df_col   <- get_col_safe(cf, c("DF", "df"))
  if (is.null(pred_col) || is.null(resp_col) || is.null(est_col)) return(NULL)
  
  out <- cf[, c(pred_col, resp_col, est_col), drop = FALSE]
  names(out) <- c("Predictor", "Response", "Estimate")
  out$`Std.Error of Estimate` <- if (!is.null(se_col)) as.numeric(cf[[se_col]]) else NA_real_
  out$`P-value` <- if (!is.null(p_col)) as.numeric(cf[[p_col]]) else NA_real_
  out$`t-value` <- if (!is.null(t_col)) as.numeric(cf[[t_col]]) else NA_real_
  out$df <- if (!is.null(df_col)) as.numeric(cf[[df_col]]) else NA_real_
  out$Estimate <- as.numeric(out$Estimate)
  out$`Lower Bound of Confidence Interval` <- ifelse(is.na(out$df) | is.na(out$`Std.Error of Estimate`), NA_real_, out$Estimate - out$`Std.Error of Estimate` * qt(0.975, out$df))
  out$`Upper Bound of Confidence Interval` <- ifelse(is.na(out$df) | is.na(out$`Std.Error of Estimate`), NA_real_, out$Estimate + out$`Std.Error of Estimate` * qt(0.975, out$df))
  out <- out[!out$Predictor %in% c("(Intercept)", "Intercept"), , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  
  out$`Trophic interaction` <- pair_id
  out$`Number of observations` <- n_obs
  out$`Number of studies` <- n_study
  out$R2 <- NA_real_
  out$SEM_type <- sem_type
  
  if (!is.null(category_var) && !is.null(category)) {
    out$`Category variable` <- category_var
    out$Category <- category
    out <- out[, c("Trophic interaction", "Category variable", "Category", "Predictor", "Response", "Number of observations", "Number of studies", "Estimate", "Std.Error of Estimate", "t-value", "P-value", "df", "Lower Bound of Confidence Interval", "Upper Bound of Confidence Interval", "R2", "SEM_type")]
  } else {
    out <- out[, c("Trophic interaction", "Predictor", "Response", "Number of observations", "Number of studies", "Estimate", "Std.Error of Estimate", "t-value", "P-value", "df", "Lower Bound of Confidence Interval", "Upper Bound of Confidence Interval", "R2", "SEM_type")]
  }
  
  rownames(out) <- NULL
  out
}

for (i in seq_len(nrow(paired_names))) {
  dat <- read_txt_safe(paste0(prefix, paired_names[i, 1], ".txt"))
  num_cols <- c("yi1", "yi2", "vi1", "vi2", "NCP_N", "added_species_number")
  for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- safe_num(dat[[cc]])
  
  dat0 <- dat[!is.na(dat$S_ID) & !is.na(dat$ID) &
                !is.na(dat$yi1) & is.finite(dat$yi1) &
                !is.na(dat$yi2) & is.finite(dat$yi2) &
                !is.na(dat$vi1) & is.finite(dat$vi1) & dat$vi1 > 0 &
                !is.na(dat$vi2) & is.finite(dat$vi2) & dat$vi2 > 0 &
                !is.na(dat$NCP_N) & is.finite(dat$NCP_N) &
                !is.na(dat$added_species_number) & is.finite(dat$added_species_number), , drop = FALSE]
  
  if (nrow(dat0) == 0) next
  dat0$S_ID <- factor(dat0$S_ID)
  dat0$ID <- factor(dat0$ID)
  dat0$plant_diversity <- ifelse(dat0$NCP_N == 1, 0, ifelse(dat0$NCP_N >= 2, 1, NA))
  
  dat_bin <- dat0[!is.na(dat0$plant_diversity), , drop = FALSE]
  if (nrow(dat_bin) > 1 && length(unique(dat_bin$S_ID)) > 1 && length(unique(dat_bin$plant_diversity)) > 1) {
    m1_bin <- fit_lme_safe(yi1 ~ plant_diversity, "vi1", dat_bin, nested = TRUE)
    m2_bin <- fit_lme_safe(yi2 ~ yi1 + plant_diversity, "vi2", dat_bin, nested = TRUE)
    if (!inherits(m1_bin, "try-error") && !inherits(m2_bin, "try-error")) {
      sem_bin <- try(psem(m1_bin, m2_bin), silent = TRUE)
      if (!inherits(sem_bin, "try-error")) res13[[length(res13) + 1]] <- extract_psem_rows(sem_bin, paired_names[i, 1], nrow(dat_bin), length(unique(dat_bin$S_ID)), "binary_piecewiseSEM_nested")
    }
  }
  
  dat_add <- dat0[!is.na(dat0$added_species_number), , drop = FALSE]
  if (nrow(dat_add) > 1 && length(unique(dat_add$S_ID)) > 1 && length(unique(dat_add$added_species_number)) > 1) {
    m1_add <- fit_lme_safe(yi1 ~ added_species_number, "vi1", dat_add, nested = TRUE)
    m2_add <- fit_lme_safe(yi2 ~ yi1 + added_species_number, "vi2", dat_add, nested = TRUE)
    if (!inherits(m1_add, "try-error") && !inherits(m2_add, "try-error")) {
      sem_add <- try(psem(m1_add, m2_add), silent = TRUE)
      if (!inherits(sem_add, "try-error")) res13[[length(res13) + 1]] <- extract_psem_rows(sem_add, paired_names[i, 1], nrow(dat_add), length(unique(dat_add$S_ID)), "continuous_piecewiseSEM_nested")
    }
  }
}

res13 <- bind_rows(res13)
openxlsx::write.xlsx(res13, paste0(path, "SuppTab13_formal_piecewiseSEM_nested_optimized.xlsx"), rowNames = FALSE)

#### 6. Stable 15 ####
paired_names <- read_txt_safe(paste0(prefix, "paired_names.txt"))
category_var <- "Zone"
categories   <- c("Temperate", "Tropic")
res15 <- list()

for (i in seq_len(nrow(paired_names))) {
  dat <- read_txt_safe(paste0(prefix, paired_names[i, 1], ".txt"))
  num_cols <- c("yi1", "yi2", "vi1", "vi2", "NCP_N", "added_species_number")
  for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- safe_num(dat[[cc]])
  
  dat <- dat[!is.na(dat$S_ID) & !is.na(dat$ID) &
               !is.na(dat$yi1) & is.finite(dat$yi1) &
               !is.na(dat$yi2) & is.finite(dat$yi2) &
               !is.na(dat$vi1) & is.finite(dat$vi1) & dat$vi1 > 0 &
               !is.na(dat$vi2) & is.finite(dat$vi2) & dat$vi2 > 0, , drop = FALSE]
  if (nrow(dat) == 0) next
  
  dat$S_ID <- factor(dat$S_ID)
  dat$ID <- factor(dat$ID)
  dat$plant_diversity <- ifelse(dat$NCP_N == 1, 0, ifelse(dat$NCP_N >= 2, 1, NA))
  
  if ("Ex_3" %in% names(dat)) {
    dat$Ex_3 <- trimws(as.character(dat$Ex_3))
    dat$Ex_3 <- ifelse(tolower(dat$Ex_3) == "plot", "Plot", ifelse(tolower(dat$Ex_3) == "pot", "Pot", dat$Ex_3))
  }
  
  if (!(category_var %in% names(dat))) next
  
  for (category in categories) {
    dat_cat <- dat[!is.na(dat[[category_var]]) & dat[[category_var]] == category, , drop = FALSE]
    if (nrow(dat_cat) == 0) next
    
    dat_bin <- dat_cat[!is.na(dat_cat$plant_diversity), , drop = FALSE]
    if (nrow(dat_bin) > 1 && length(unique(dat_bin$S_ID)) > 1 && length(unique(dat_bin$plant_diversity)) > 1) {
      m1_bin <- fit_lme_safe(yi1 ~ plant_diversity, "vi1", dat_bin, nested = TRUE)
      m2_bin <- fit_lme_safe(yi2 ~ yi1 + plant_diversity, "vi2", dat_bin, nested = TRUE)
      if (!inherits(m1_bin, "try-error") && !inherits(m2_bin, "try-error")) {
        sem_bin <- try(psem(m1_bin, m2_bin), silent = TRUE)
        if (!inherits(sem_bin, "try-error")) res14[[length(res14) + 1]] <- extract_psem_rows(sem_bin, paired_names[i, 1], nrow(dat_bin), length(unique(dat_bin$S_ID)), "binary_piecewiseSEM_nested", category_var, category)
      }
    }
    
    dat_add <- dat_cat[!is.na(dat_cat$added_species_number), , drop = FALSE]
    if (nrow(dat_add) > 1 && length(unique(dat_add$S_ID)) > 1 && length(unique(dat_add$added_species_number)) > 1) {
      m1_add <- fit_lme_safe(yi1 ~ added_species_number, "vi1", dat_add, nested = TRUE)
      m2_add <- fit_lme_safe(yi2 ~ yi1 + added_species_number, "vi2", dat_add, nested = TRUE)
      if (!inherits(m1_add, "try-error") && !inherits(m2_add, "try-error")) {
        sem_add <- try(psem(m1_add, m2_add), silent = TRUE)
        if (!inherits(sem_add, "try-error")) res14[[length(res14) + 1]] <- extract_psem_rows(sem_add, paired_names[i, 1], nrow(dat_add), length(unique(dat_add$S_ID)), "continuous_piecewiseSEM_nested", category_var, category)
      }
    }
  }
}

res15 <- bind_rows(res15)
if (nrow(res14) == 0) stop("No valid results were generated.")
openxlsx::write.xlsx(res14, paste0(path, "Stable15_", category_var, "_piecewiseSEM_optimized.xlsx"), rowNames = FALSE)


###### STABLE 17#####
## Bi-trophic subgroup analysis
## (moved later; originally STABLE 16)
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
  library(Matrix)
  library(corpcor)
  library(openxlsx)
})

prefix <- "/Users/yusha/"
path   <- "/Users/yusha/"

#==================================================
# 1. Input
#==================================================
dat_all <- read.table(
  paste0(prefix, "herbivore_vs_crop.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

#==================================================
# 2. Shared-control grouping parameters
#==================================================
round_digits_m  <- 4
round_digits_sd <- 4
round_digits_r  <- 0

#==================================================
# 3. Standardization helpers
#==================================================
standardize_zone <- function(x) {
  x <- trimws(as.character(x))
  x_low <- tolower(x)
  x_out <- x
  x_out[x_low %in% c("tropical", "tropics", "tropic")] <- "Tropic"
  x_out[x_low %in% c("temperate")] <- "Temperate"
  x_out
}

standardize_ex3 <- function(x) {
  x <- trimws(as.character(x))
  x_low <- tolower(x)
  x_out <- x
  x_out[x_low %in% c("plot")] <- "Plot"
  x_out[x_low %in% c("pot")]  <- "Pot"
  x_out
}

#==================================================
# 4. Generic helpers
#==================================================
round_for_group <- function(x, digits = 4) if (is.numeric(x)) round(x, digits = digits) else x
make_common_key <- function(m, sd, r) interaction(m, sd, r, drop = TRUE, lex.order = TRUE)

fit_rma_mv_try <- function(yi_col, V_input, dat, method = "ML") {
  rand_eff <- list(~1 | S_ID, ~1 | ID_global)
  for (opt in c("BFGS", "nlminb", "optim")) {
    fit <- try(
      rma.mv(
        as.formula(paste0(yi_col, " ~ 1")),
        V = V_input,
        random = rand_eff,
        data = dat,
        test = "t",
        control = list(optimizer = opt),
        method = method
      ),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) return(fit)
  }
  "try-error"
}

build_shared_V <- function(dat, vi_col, sd_col) {
  calc.v <- function(x) {
    cpart <- x[[sd_col]][1]^2 / (x$CONTROL_R[1] * x$CONTROL_M[1]^2)
    v <- matrix(cpart, nrow = nrow(x), ncol = nrow(x))
    diag(v) <- x[[vi_col]]
    v
  }
  
  out <- list(V = NULL, V_type = NA_character_, min_eigen_before_fix = NA_real_, ok = FALSE, msg = NULL)
  
  V_try <- try({
    V_tmp <- bldiag(lapply(split(dat, dat$Common_ID_global), calc.v))
    as.matrix(V_tmp)
  }, silent = TRUE)
  
  if (inherits(V_try, "try-error")) {
    out$msg <- "failed building shared-control V"
    return(out)
  }
  
  V_ok <- isSymmetric(V_try) && !any(is.na(V_try)) && !any(is.infinite(V_try))
  if (!V_ok) {
    out$msg <- "invalid shared-control V"
    return(out)
  }
  
  eig_vals <- try(eigen(V_try, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (inherits(eig_vals, "try-error")) {
    out$msg <- "failed eigen decomposition"
    return(out)
  }
  
  out$min_eigen_before_fix <- min(eig_vals)
  
  if (!is.finite(out$min_eigen_before_fix) || out$min_eigen_before_fix <= 0) {
    V_fix <- try(make.positive.definite(V_try), silent = TRUE)
    if (inherits(V_fix, "try-error")) {
      out$msg <- "shared-control V not PD and PD-fix failed"
      return(out)
    }
    out$V <- as.matrix(V_fix)
    out$V_type <- "shared_control_pd_fix"
    out$ok <- TRUE
    out$msg <- "shared-control V not PD, fixed before fitting"
  } else {
    out$V <- V_try
    out$V_type <- "shared_control"
    out$ok <- TRUE
  }
  
  out
}

empty_bi_row <- function(item_name, category_name, n_obs = NA, n_study = NA) {
  data.frame(
    Item = item_name,
    Category = category_name,
    Number_of_observations = n_obs,
    Number_of_studies = n_study,
    Effect_size = NA,
    t_value = NA,
    P_value = NA,
    df = NA,
    CIlb = NA,
    CIub = NA,
    stringsAsFactors = FALSE
  )
}

run_one_response_bi <- function(dat0, yi_col, vi_col, sd_col, item_name, category_name) {
  dat.fit <- NULL
  min_eigen_before_fix <- NA_real_
  
  dat <- dat0 %>%
    filter(
      !is.na(.data[[yi_col]]), is.finite(.data[[yi_col]]),
      !is.na(.data[[vi_col]]), is.finite(.data[[vi_col]]), .data[[vi_col]] > 0,
      !is.na(.data[[sd_col]]), is.finite(.data[[sd_col]]), .data[[sd_col]] >= 0,
      !is.na(CONTROL_M), is.finite(CONTROL_M), CONTROL_M > 0,
      !is.na(CONTROL_R), is.finite(CONTROL_R), CONTROL_R > 0,
      !is.na(S_ID),
      !is.na(ID_global)
    ) %>%
    mutate(
      CONTROL_M_grp  = round_for_group(CONTROL_M, round_digits_m),
      CONTROL_SD_grp = round_for_group(.data[[sd_col]], round_digits_sd),
      CONTROL_R_grp  = round_for_group(CONTROL_R, round_digits_r)
    ) %>%
    group_by(S_ID) %>%
    mutate(
      Common_ID = match(
        make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp),
        unique(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp))
      )
    ) %>%
    ungroup() %>%
    mutate(
      Common_ID_global = match(
        interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE),
        unique(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE))
      )
    )
  
  if (nrow(dat) == 0) return(empty_bi_row(item_name, category_name, n_obs = 0))
  
  check_common <- dat %>%
    group_by(Common_ID_global) %>%
    summarise(
      n_CONTROL_M_grp  = n_distinct(CONTROL_M_grp),
      n_CONTROL_SD_grp = n_distinct(CONTROL_SD_grp),
      n_CONTROL_R_grp  = n_distinct(CONTROL_R_grp),
      .groups = "drop"
    )
  
  use_V <- nrow(check_common) > 0 &&
    all(check_common$n_CONTROL_M_grp == 1) &&
    all(check_common$n_CONTROL_SD_grp == 1) &&
    all(check_common$n_CONTROL_R_grp == 1)
  
  if (nrow(dat) == 1) {
    dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "z"), silent = TRUE)
  }
  
  if ((is.null(dat.fit) || inherits(dat.fit, "try-error")) && nrow(dat) > 1 && use_V) {
    V_info <- build_shared_V(dat, vi_col = vi_col, sd_col = sd_col)
    min_eigen_before_fix <- V_info$min_eigen_before_fix
    
    if (isTRUE(V_info$ok)) {
      dat.fit <- fit_rma_mv_try(yi_col = yi_col, V_input = V_info$V, dat = dat, method = "ML")
      if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_try(yi_col = yi_col, V_input = V_info$V, dat = dat, method = "REML")
    }
  }
  
  if (is.null(dat.fit) || inherits(dat.fit, "try-error")) {
    if (length(unique(dat$S_ID)) == 1) {
      dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "t"), silent = TRUE)
    } else {
      dat.fit <- fit_rma_mv_try(yi_col = yi_col, V_input = dat[[vi_col]], dat = dat, method = "ML")
      if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_try(yi_col = yi_col, V_input = dat[[vi_col]], dat = dat, method = "REML")
    }
  }
  
  if (is.null(dat.fit) || inherits(dat.fit, "try-error")) {
    return(empty_bi_row(item_name, category_name, n_obs = nrow(dat), n_study = length(unique(dat$S_ID))))
  }
  
  tmp_sum <- coef(summary(dat.fit))
  t_value_out <- if ("tval" %in% colnames(tmp_sum)) tmp_sum$tval[1] else if ("zval" %in% colnames(tmp_sum)) tmp_sum$zval[1] else NA
  p_value_out <- if ("pval" %in% colnames(tmp_sum)) tmp_sum$pval[1] else NA
  df_out      <- if ("df" %in% colnames(tmp_sum)) tmp_sum$df[1] else NA
  ci_lb_out   <- if (!is.null(dat.fit$ci.lb)) dat.fit$ci.lb[1] else NA
  ci_ub_out   <- if (!is.null(dat.fit$ci.ub)) dat.fit$ci.ub[1] else NA
  
  data.frame(
    Item = item_name,
    Category = category_name,
    Number_of_observations = nrow(dat),
    Number_of_studies = length(unique(dat$S_ID)),
    Effect_size = as.numeric(dat.fit$b[1]),
    t_value = t_value_out,
    P_value = p_value_out,
    df = df_out,
    CIlb = ci_lb_out,
    CIub = ci_ub_out,
    stringsAsFactors = FALSE
  )
}

#==================================================
# 5. Data preprocessing
#==================================================
dat_all <- dat_all %>%
  mutate(
    Zone = standardize_zone(Zone),
    Ex_3 = standardize_ex3(Ex_3)
  )

if (!"EsID" %in% names(dat_all)) {
  dat_all <- dat_all %>% group_by(S_ID) %>% mutate(EsID = seq_len(n())) %>% ungroup()
}

if (!"ID_global" %in% names(dat_all)) {
  if ("ID" %in% names(dat_all)) {
    dat_all <- dat_all %>% mutate(ID_global = interaction(S_ID, ID, drop = TRUE, lex.order = TRUE))
  } else {
    dat_all <- dat_all %>%
      group_by(S_ID) %>% mutate(ID = seq_len(n())) %>% ungroup() %>%
      mutate(ID_global = interaction(S_ID, ID, drop = TRUE, lex.order = TRUE))
  }
}

#==================================================
# 6. Subgroups
#==================================================
subgroup_list <- list(
  CP_C3 = c("Herbaceous_plant", "Woody_plant"),
  CP_Ct = c("Food_crop", "Cash_crop"),
  Mn__  = c("Generalist", "Specialist"),
  Dv_S = c("Intercropping", "Cover_cropping", "Sown_field_margins"),
  Ex_4  = c("Natural"),
  Ex_3  = c("Plot", "Pot"),
  Zone  = c("Temperate", "Tropic")
)

#==================================================
# 7. Main loop
#==================================================
res <- NULL

for (mod_name in names(subgroup_list)) {
  for (lv in subgroup_list[[mod_name]]) {
    dat_sub <- dat_all %>%
      filter(!is.na(.data[[mod_name]])) %>%
      filter(as.character(.data[[mod_name]]) == lv)
    
    res <- rbind(
      res,
      run_one_response_bi(dat_sub, "yi1", "vi1", "Cntrl.s1", lv, "Herbivore performance"),
      run_one_response_bi(dat_sub, "yi2", "vi2", "Cntrl.s2", lv, "Crop performance")
    )
  }
}

#==================================================
# 8. Output
#==================================================
res <- res %>%
  select(Item, Category, Number_of_observations, Number_of_studies, Effect_size, t_value, P_value, df, CIlb, CIub)

openxlsx::write.xlsx(res, paste0(path, "STABLE17.xlsx"), rowNames = FALSE)

print(res)
cat("Finished STABLE17.\n")