#### data cleaning + Stable 12 + Stable 13 + Stable 14 + Stable 16 ####
rm(list = ls())

suppressPackageStartupMessages({ library(readxl); library(dplyr); library(metafor); library(corpcor); library(openxlsx); library(nlme) })

#### 1. Paths and global settings ####
prefix <- "/Users/yusha/"; path <- "/Users/yusha/"

#### 2. Helper functions ####
safe_num <- function(x) suppressWarnings(as.numeric(x))

read_txt_safe <- function(file_in) read.table(file_in, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "", fill = TRUE, blank.lines.skip = FALSE)

fit_lme_safe <- function(formula_obj, data_obj, vi_col) try(lme(fixed = formula_obj, random = ~1 | S_ID/ID, weights = varFixed(as.formula(paste0("~", vi_col))), data = data_obj, na.action = na.omit, control = lmeControl(sigma = 1, returnObject = TRUE, opt = "optim")), silent = TRUE)

fit_lme_twolevel <- function(fixed_formula, vi_var, data_obj) {
  fit1 <- try(lme(fixed = fixed_formula, random = ~1 | S_ID/ID, weights = varFixed(as.formula(paste0("~", vi_var))), data = data_obj, na.action = na.omit, control = lmeControl(sigma = 1, returnObject = TRUE)), silent = TRUE)
  if (!inherits(fit1, "try-error")) return(list(model = fit1, fit_type = "nested_SID_ID_default"))
  fit2 <- try(lme(fixed = fixed_formula, random = ~1 | S_ID/ID, weights = varFixed(as.formula(paste0("~", vi_var))), data = data_obj, na.action = na.omit, control = lmeControl(sigma = 1, returnObject = TRUE, opt = "optim", msMaxIter = 200, msMaxEval = 400)), silent = TRUE)
  if (!inherits(fit2, "try-error")) return(list(model = fit2, fit_type = "nested_SID_ID_optim"))
  fit3 <- try(lme(fixed = fixed_formula, random = ~1 | S_ID, weights = varFixed(as.formula(paste0("~", vi_var))), data = data_obj, na.action = na.omit, control = lmeControl(sigma = 1, returnObject = TRUE)), silent = TRUE)
  if (!inherits(fit3, "try-error")) return(list(model = fit3, fit_type = "random_SID_default"))
  fit4 <- try(lme(fixed = fixed_formula, random = ~1 | S_ID, weights = varFixed(as.formula(paste0("~", vi_var))), data = data_obj, na.action = na.omit, control = lmeControl(sigma = 1, returnObject = TRUE, opt = "optim", msMaxIter = 200, msMaxEval = 400)), silent = TRUE)
  if (!inherits(fit4, "try-error")) return(list(model = fit4, fit_type = "random_SID_optim"))
  list(model = NULL, fit_type = NA_character_)
}

is_fitworthy <- function(dat_sub, predictor_name) { if (nrow(dat_sub) < 5) return(FALSE); if (length(unique(dat_sub$S_ID)) < 2) return(FALSE); if (!(predictor_name %in% names(dat_sub))) return(FALSE); x <- dat_sub[[predictor_name]]; x <- x[is.finite(x) & !is.na(x)]; if (length(unique(x)) < 2) return(FALSE); TRUE }

make_model_data <- function(dat, vars_needed) { keep <- !is.na(dat$S_ID) & trimws(as.character(dat$S_ID)) != "" & !is.na(dat$ID) & trimws(as.character(dat$ID)) != ""; for (v in vars_needed) { keep <- keep & !is.na(dat[[v]]) & is.finite(dat[[v]]); if (grepl("^vi", v)) keep <- keep & dat[[v]] > 0 }; dat[keep, , drop = FALSE] }

append_result <- function(res_obj, trophic_interaction, category_var, category, predictor, response, n_obs, n_study, fit_obj, coef_name, sem_type) {
  if (is.null(fit_obj$model)) return(res_obj)
  tt <- try(summary(fit_obj$model)$tTable, silent = TRUE)
  if (inherits(tt, "try-error") || is.null(tt) || !(coef_name %in% rownames(tt))) return(res_obj)
  est <- as.numeric(tt[coef_name, "Value"]); se <- as.numeric(tt[coef_name, "Std.Error"]); tv <- as.numeric(tt[coef_name, "t-value"]); pv <- as.numeric(tt[coef_name, "p-value"]); dfv <- as.numeric(tt[coef_name, "DF"]); lb <- est - se * qt(0.975, dfv); ub <- est + se * qt(0.975, dfv)
  out <- data.frame("Trophic interaction" = trophic_interaction, "Category variable" = category_var, "Category" = category, "Predictor" = predictor, "Response" = response, "Number of observations" = n_obs, "Number of studies" = n_study, "Estimate" = est, "Std.Error of Estimate" = se, "t-value" = tv, "P-value" = pv, "df" = dfv, "Lower Bound of Confidence Interval" = lb, "Upper Bound of Confidence Interval" = ub, "R2" = NA_real_, "SEM_type" = sem_type, "Fit_type" = fit_obj$fit_type, check.names = FALSE, stringsAsFactors = FALSE)
  rbind(res_obj, out)
}

get_control_cols <- function(dat, idx) { m_col <- c(paste0("Cntrl.m", idx), "CONTROL_M", "Cntrl.m", "control_m"); s_col <- c(paste0("Cntrl.s", idx), "CONTROL_SD", "Cntrl.s", "control_sd"); r_col <- c(paste0("Cntrl.n", idx), "CONTROL_R", "Cntrl.n", "control_n"); list(m = m_col[m_col %in% names(dat)][1], s = s_col[s_col %in% names(dat)][1], r = r_col[r_col %in% names(dat)][1]) }

fit_rma_mv_safe <- function(yi_col, V_input, dat, method = "ML") { random_terms <- list(~1 | S_ID, ~1 | ID_global); for (opt in c("BFGS", "nlminb", "optim")) { fit <- try(rma.mv(as.formula(paste0(yi_col, " ~ 1")), V = V_input, random = random_terms, data = dat, test = "t", method = method, control = list(optimizer = opt)), silent = TRUE); if (!inherits(fit, "try-error")) return(fit) }; "try-error" }

build_shared_V <- function(dat, vi_col, m_col, s_col, r_col) {
  block_V <- function(x) { cpart <- x[[s_col]][1]^2 / (x[[r_col]][1] * x[[m_col]][1]^2); V <- matrix(cpart, nrow = nrow(x), ncol = nrow(x)); diag(V) <- x[[vi_col]]; V }
  out <- list(V = NULL, V_type = NA_character_, min_eigen = NA_real_, ok = FALSE)
  V_try <- try(as.matrix(Matrix::bdiag(lapply(split(dat, dat$Common_ID_global), block_V))), silent = TRUE)
  if (inherits(V_try, "try-error")) return(out)
  if (!isSymmetric(V_try) || any(is.na(V_try)) || any(is.infinite(V_try))) return(out)
  eig <- try(eigen(V_try, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (inherits(eig, "try-error")) return(out)
  out$min_eigen <- min(eig)
  if (!is.finite(out$min_eigen) || out$min_eigen <= 0) { V_fix <- try(corpcor::make.positive.definite(V_try), silent = TRUE); if (inherits(V_fix, "try-error")) return(out); out$V <- as.matrix(V_fix); out$V_type <- "shared_control_pd_fix" } else { out$V <- V_try; out$V_type <- "shared_control" }
  out$ok <- TRUE; out
}

run_one_response <- function(dat0, yi_col, vi_col, response_name, pair_id, trophic_interaction, idx, round_digits_m = 4, round_digits_sd = 4, round_digits_r = 0) {
  ctrl <- get_control_cols(dat0, idx)
  if (is.na(ctrl$m) || is.na(ctrl$s) || is.na(ctrl$r)) return(data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = response_name, Number_of_observations = NA, Number_of_studies = NA, V_type = NA, Min_eigen_before_fix = NA, Effect_size = NA, t_value = NA, P_value = NA, df = NA, CI_lb = NA, CI_ub = NA, stringsAsFactors = FALSE))
  dat <- dat0; dat[[yi_col]] <- safe_num(dat[[yi_col]]); dat[[vi_col]] <- safe_num(dat[[vi_col]]); dat[[ctrl$m]] <- safe_num(dat[[ctrl$m]]); dat[[ctrl$s]] <- safe_num(dat[[ctrl$s]]); dat[[ctrl$r]] <- safe_num(dat[[ctrl$r]])
  keep <- !is.na(dat[[yi_col]]) & is.finite(dat[[yi_col]]) & !is.na(dat[[vi_col]]) & is.finite(dat[[vi_col]]) & dat[[vi_col]] > 0 & !is.na(dat$S_ID) & trimws(as.character(dat$S_ID)) != "" & !is.na(dat$ID) & trimws(as.character(dat$ID)) != "" & !is.na(dat[[ctrl$m]]) & is.finite(dat[[ctrl$m]]) & dat[[ctrl$m]] > 0 & !is.na(dat[[ctrl$s]]) & is.finite(dat[[ctrl$s]]) & dat[[ctrl$s]] >= 0 & !is.na(dat[[ctrl$r]]) & is.finite(dat[[ctrl$r]]) & dat[[ctrl$r]] > 0
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) == 0) return(data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = response_name, Number_of_observations = 0, Number_of_studies = NA, V_type = NA, Min_eigen_before_fix = NA, Effect_size = NA, t_value = NA, P_value = NA, df = NA, CI_lb = NA, CI_ub = NA, stringsAsFactors = FALSE))
  dat$CONTROL_M_use <- dat[[ctrl$m]]; dat$CONTROL_SD_use <- dat[[ctrl$s]]; dat$CONTROL_R_use <- dat[[ctrl$r]]; dat$CONTROL_M_grp <- round(dat$CONTROL_M_use, round_digits_m); dat$CONTROL_SD_grp <- round(dat$CONTROL_SD_use, round_digits_sd); dat$CONTROL_R_grp <- round(dat$CONTROL_R_use, round_digits_r)
  dat$Common_ID <- ave(seq_len(nrow(dat)), dat$S_ID, FUN = function(ii) match(interaction(dat$CONTROL_M_grp[ii], dat$CONTROL_SD_grp[ii], dat$CONTROL_R_grp[ii], drop = TRUE, lex.order = TRUE), unique(interaction(dat$CONTROL_M_grp[ii], dat$CONTROL_SD_grp[ii], dat$CONTROL_R_grp[ii], drop = TRUE, lex.order = TRUE))))
  dat$Common_ID_global <- match(interaction(dat$S_ID, dat$Common_ID, drop = TRUE, lex.order = TRUE), unique(interaction(dat$S_ID, dat$Common_ID, drop = TRUE, lex.order = TRUE))); dat$ID_global <- interaction(dat$S_ID, dat$ID, drop = TRUE, lex.order = TRUE); dat <- dat[!is.na(dat$ID_global), , drop = FALSE]
  if (nrow(dat) == 0) return(data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = response_name, Number_of_observations = 0, Number_of_studies = NA, V_type = NA, Min_eigen_before_fix = NA, Effect_size = NA, t_value = NA, P_value = NA, df = NA, CI_lb = NA, CI_ub = NA, stringsAsFactors = FALSE))
  check_common <- aggregate(cbind(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp) ~ Common_ID_global, data = dat, FUN = function(x) length(unique(x))); use_V <- all(check_common$CONTROL_M_grp == 1) && all(check_common$CONTROL_SD_grp == 1) && all(check_common$CONTROL_R_grp == 1)
  dat.fit <- NULL; V_type <- NA_character_; min_eigen <- NA_real_
  if (nrow(dat) == 1) { dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "z"), silent = TRUE); if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi_single_obs" } else {
    if (use_V) { V_info <- build_shared_V(dat, vi_col, "CONTROL_M_use", "CONTROL_SD_use", "CONTROL_R_use"); min_eigen <- V_info$min_eigen; if (isTRUE(V_info$ok) && nrow(V_info$V) == nrow(dat)) { dat.fit <- fit_rma_mv_safe(yi_col, V_info$V, dat, method = "ML"); if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_safe(yi_col, V_info$V, dat, method = "REML"); if (!inherits(dat.fit, "try-error")) V_type <- V_info$V_type } }
    if (is.null(dat.fit) || inherits(dat.fit, "try-error")) { if (length(unique(dat$S_ID)) == 1) { dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "t"), silent = TRUE); if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi_rma" } else { dat.fit <- fit_rma_mv_safe(yi_col, dat[[vi_col]], dat, method = "ML"); if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_safe(yi_col, dat[[vi_col]], dat, method = "REML"); if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi" } }
  }
  if (is.null(dat.fit) || inherits(dat.fit, "try-error")) return(data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = response_name, Number_of_observations = nrow(dat), Number_of_studies = length(unique(dat$S_ID)), V_type = NA, Min_eigen_before_fix = min_eigen, Effect_size = NA, t_value = NA, P_value = NA, df = NA, CI_lb = NA, CI_ub = NA, stringsAsFactors = FALSE))
  sm <- coef(summary(dat.fit)); t_value <- if ("tval" %in% colnames(sm)) sm[1, "tval"] else if ("zval" %in% colnames(sm)) sm[1, "zval"] else NA; p_value <- if ("pval" %in% colnames(sm)) sm[1, "pval"] else NA; df_out <- if ("df" %in% colnames(sm)) sm[1, "df"] else NA
  data.frame(Pair_ID = pair_id, Trophic_interaction = trophic_interaction, Category = response_name, Number_of_observations = nrow(dat), Number_of_studies = length(unique(dat$S_ID)), V_type = V_type, Min_eigen_before_fix = min_eigen, Effect_size = as.numeric(dat.fit$b[1]), t_value = as.numeric(t_value), P_value = as.numeric(p_value), df = as.numeric(df_out), CI_lb = if (!is.null(dat.fit$ci.lb)) dat.fit$ci.lb[1] else NA, CI_ub = if (!is.null(dat.fit$ci.ub)) dat.fit$ci.ub[1] else NA, stringsAsFactors = FALSE)
}

#### 3. Data cleaning: enemies + herbivore + crop ####
filenames <- dir(path = prefix, pattern = "^0401Paired data for enemies and herbivore and crop.*\\.xlsx$")
sheets <- c("enemies", "herbivore", "crop")
target_cols_2 <- c("yi2", "vi2", "Trtmnt.s2", "Cntrl.s2", "cv2_trea_new2", "cv2_cont_new2")
target_cols_3 <- c("yi3", "vi3", "Trtmnt.s3", "Cntrl.s3", "cv2_trea_new3", "cv2_cont_new3")
dt_total0 <- NULL; paired_names <- NULL

for (i in seq_along(filenames)) {
  file_path <- file.path(prefix, filenames[i]); pair_id <- paste0("enemies_vs_herbivore_vs_crop_", i); dt_new0 <- NULL
  for (j in 1:3) {
    dt0 <- tryCatch(read_excel(file_path, sheet = sheets[j], range = "A1:AG417"), error = function(e) NULL)
    if (is.null(dt0) || nrow(dt0) == 0) next
    colnames(dt0) <- make.names(abbreviate(colnames(dt0), minlength = 4))
    dt0 <- dt0 %>% filter(!is.na(CONTROL_M), is.finite(CONTROL_M), CONTROL_M > 0, !is.na(TREATMENT_M), is.finite(TREATMENT_M), TREATMENT_M > 0, !is.na(CONTROL_R), is.finite(CONTROL_R), CONTROL_R > 0, !is.na(TREATMENT_R), is.finite(TREATMENT_R), TREATMENT_R > 0, !is.na(S_ID))
    if (nrow(dt0) == 0) next
    dt0$CONTROL_SD <- replace(dt0$CONTROL_SD, !is.na(dt0$CONTROL_SD) & round(dt0$CONTROL_M / dt0$CONTROL_SD, 5) == 10, NA)
    dt0$TREATMENT_SD <- replace(dt0$TREATMENT_SD, is.na(dt0$CONTROL_SD), NA)
    ESData <- dt0 %>% mutate(cv_Control = na_if(CONTROL_SD / CONTROL_M, Inf), cv_Treatment = na_if(TREATMENT_SD / TREATMENT_M, Inf))
    ESData <- cv_avg(x = TREATMENT_M, sd = TREATMENT_SD, n = TREATMENT_R, group = S_ID, label = "1", data = ESData)
    ESData <- cv_avg(x = CONTROL_M, sd = CONTROL_SD, n = CONTROL_R, group = S_ID, label = "2", data = ESData)
    
    if (j == 1) {
      ESData <- ESData %>% mutate(cv2_trea_new1 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2), cv2_cont_new1 = if_else(is.na(cv_Control), b_CV2_2, cv_Control^2), yi1 = lnrr_laj(m1 = TREATMENT_M, m2 = CONTROL_M, cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R), vi1 = v_lnrr_laj(cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R))
      ESData$Trtmnt.s1 <- ESData$TREATMENT_SD; ESData$Cntrl.s1 <- ESData$CONTROL_SD
      idx_t1 <- which(is.na(ESData$Trtmnt.s1)); if (length(idx_t1) > 0) ESData[idx_t1, "Trtmnt.s1"] <- ESData[idx_t1, "cv2_trea_new1"] * ESData[idx_t1, "TREATMENT_M"]
      idx_c1 <- which(is.na(ESData$Cntrl.s1)); if (length(idx_c1) > 0) ESData[idx_c1, "Cntrl.s1"] <- ESData[idx_c1, "cv2_cont_new1"] * ESData[idx_c1, "CONTROL_M"]
      dt_new0 <- ESData
    }
    
    if (j == 2) {
      ESData <- ESData %>% mutate(cv2_trea_new2 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2), cv2_cont_new2 = if_else(is.na(cv_Control), b_CV2_2, cv_Control^2), yi2 = lnrr_laj(m1 = TREATMENT_M, m2 = CONTROL_M, cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R), vi2 = v_lnrr_laj(cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R))
      ESData$Trtmnt.s2 <- ESData$TREATMENT_SD; ESData$Cntrl.s2 <- ESData$CONTROL_SD
      idx_t2 <- which(is.na(ESData$Trtmnt.s2)); if (length(idx_t2) > 0) ESData[idx_t2, "Trtmnt.s2"] <- ESData[idx_t2, "cv2_trea_new2"] * ESData[idx_t2, "TREATMENT_M"]
      idx_c2 <- which(is.na(ESData$Cntrl.s2)); if (length(idx_c2) > 0) ESData[idx_c2, "Cntrl.s2"] <- ESData[idx_c2, "cv2_cont_new2"] * ESData[idx_c2, "CONTROL_M"]
      if (!is.null(dt_new0)) dt_new0 <- data.frame(dt_new0, ESData[, target_cols_2, drop = FALSE])
    }
    
    if (j == 3) {
      ESData <- ESData %>% mutate(cv2_trea_new3 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2), cv2_cont_new3 = if_else(is.na(cv_Control), b_CV2_2, cv_Control^2), yi3 = lnrr_laj(m1 = TREATMENT_M, m2 = CONTROL_M, cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R), vi3 = v_lnrr_laj(cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R))
      ESData$Trtmnt.s3 <- ESData$TREATMENT_SD; ESData$Cntrl.s3 <- ESData$CONTROL_SD
      idx_t3 <- which(is.na(ESData$Trtmnt.s3)); if (length(idx_t3) > 0) ESData[idx_t3, "Trtmnt.s3"] <- ESData[idx_t3, "cv2_trea_new3"] * ESData[idx_t3, "TREATMENT_M"]
      idx_c3 <- which(is.na(ESData$Cntrl.s3)); if (length(idx_c3) > 0) ESData[idx_c3, "Cntrl.s3"] <- ESData[idx_c3, "cv2_cont_new3"] * ESData[idx_c3, "CONTROL_M"]
      if (!is.null(dt_new0)) dt_new0 <- data.frame(dt_new0, ESData[, target_cols_3, drop = FALSE])
    }
  }
  
  if (!is.null(dt_new0)) {
    if (all(c("yi1", "yi2", "yi3") %in% names(dt_new0))) dt_new0 <- dt_new0 %>% filter(!is.na(yi1), is.finite(yi1), !is.na(yi2), is.finite(yi2), !is.na(yi3), is.finite(yi3))
    if (nrow(dt_new0) == 0) next
    S_ID_tab <- table(dt_new0$S_ID); esid <- NULL; for (k in 1:length(S_ID_tab)) esid <- c(esid, seq_len(S_ID_tab[k])); dt_new0$EsID <- esid
    if ("NCP_N" %in% names(dt_new0)) dt_new0$added_species_number <- log2(dt_new0$NCP_N)
    dt_new0$pair_id <- pair_id; dt_new0$response1 <- "Enemies"; dt_new0$response2 <- "Herbivore"; dt_new0$response3 <- "Crop"; dt_new0$trophic_interaction <- "enemies_vs_herbivore_vs_crop"
    dt_total0 <- rbind(dt_total0, dt_new0)
    paired_names <- rbind(paired_names, data.frame(Pair_ID = pair_id, Trophic_interaction = "enemies_vs_herbivore_vs_crop", Response1 = "Enemies", Response2 = "Herbivore", Response3 = "Crop", Source_file = filenames[i], stringsAsFactors = FALSE))
    write.table(dt_new0, file = paste0(prefix, pair_id, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

write.table(paired_names, paste0(prefix, "paired_names.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dt_total0, paste0(prefix, "enemies_vs_herbivore_vs_crop.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

#### 4. Stable 12 ####
paired_names <- read_txt_safe(paste0(prefix, "paired_names.txt"))
round_digits_m <- 4; round_digits_sd <- 4; round_digits_r <- 0
res12 <- data.frame(Pair_ID = character(), Trophic_interaction = character(), Category = character(), Number_of_observations = numeric(), Number_of_studies = numeric(), V_type = character(), Min_eigen_before_fix = numeric(), Effect_size = numeric(), t_value = numeric(), P_value = numeric(), df = numeric(), CI_lb = numeric(), CI_ub = numeric(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(paired_names))) {
  pair_id <- paired_names$Pair_ID[i]; file_in <- paste0(prefix, pair_id, ".txt"); if (!file.exists(file_in)) next
  dat0 <- read_txt_safe(file_in)
  res12 <- rbind(res12, run_one_response(dat0, "yi1", "vi1", paired_names$Response1[i], pair_id, paired_names$Trophic_interaction[i], 1, round_digits_m, round_digits_sd, round_digits_r), run_one_response(dat0, "yi2", "vi2", paired_names$Response2[i], pair_id, paired_names$Trophic_interaction[i], 2, round_digits_m, round_digits_sd, round_digits_r), run_one_response(dat0, "yi3", "vi3", paired_names$Response3[i], pair_id, paired_names$Trophic_interaction[i], 3, round_digits_m, round_digits_sd, round_digits_r))
}

openxlsx::write.xlsx(res12, paste0(path, "Stable12-tri-trophic_clean_NC_style.xlsx"), rowNames = FALSE)

#### 5. Stable 13 ####
paired_names <- read_txt_safe(paste0(prefix, "paired_names.txt"))
res13 <- data.frame(`Trophic interaction` = character(), Predictor = character(), Response = character(), `Number of observations` = numeric(), `Number of studies` = numeric(), Estimate = numeric(), `Std.Error of Estimate` = numeric(), `t-value` = numeric(), `P-value` = numeric(), df = numeric(), `Lower Bound of Confidence Interval` = numeric(), `Upper Bound of Confidence Interval` = numeric(), R2 = numeric(), SEM_type = character(), stringsAsFactors = FALSE)

extract_term_row13 <- function(model, term, pair_id, response, n_obs, n_study, sem_type) {
  tt <- summary(model)$tTable
  data.frame(`Trophic interaction` = pair_id, Predictor = term, Response = response, `Number of observations` = n_obs, `Number of studies` = n_study, Estimate = as.numeric(tt[term, "Value"]), `Std.Error of Estimate` = as.numeric(tt[term, "Std.Error"]), `t-value` = as.numeric(tt[term, "t-value"]), `P-value` = as.numeric(tt[term, "p-value"]), df = as.numeric(tt[term, "DF"]), `Lower Bound of Confidence Interval` = as.numeric(tt[term, "Value"]) - qt(0.975, as.numeric(tt[term, "DF"])) * as.numeric(tt[term, "Std.Error"]), `Upper Bound of Confidence Interval` = as.numeric(tt[term, "Value"]) + qt(0.975, as.numeric(tt[term, "DF"])) * as.numeric(tt[term, "Std.Error"]), R2 = NA_real_, SEM_type = sem_type, stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(paired_names))) {
  pair_id <- paired_names[i, 1]; file_in <- paste0(prefix, pair_id, ".txt"); if (!file.exists(file_in)) next
  dat <- read_txt_safe(file_in)
  num_cols <- c("yi1", "yi2", "yi3", "vi1", "vi2", "vi3", "NCP_N", "added_species_number"); for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- safe_num(dat[[cc]])
  dat$S_ID <- as.character(dat$S_ID); dat$ID <- as.character(dat$ID); dat$plant_diversity <- ifelse(dat$NCP_N == 1, 0, ifelse(dat$NCP_N >= 2, 1, NA))
  
  dat_m1_bin <- dat[!is.na(dat$S_ID) & trimws(dat$S_ID) != "" & !is.na(dat$ID) & trimws(dat$ID) != "" & !is.na(dat$yi1) & is.finite(dat$yi1) & !is.na(dat$vi1) & is.finite(dat$vi1) & dat$vi1 > 0 & !is.na(dat$plant_diversity) & is.finite(dat$plant_diversity), , drop = FALSE]
  if (nrow(dat_m1_bin) > 1) { dat_m1_bin$S_ID <- factor(dat_m1_bin$S_ID); dat_m1_bin$ID <- factor(dat_m1_bin$ID); m1_bin <- fit_lme_safe(yi1 ~ plant_diversity, dat_m1_bin, "vi1"); if (!inherits(m1_bin, "try-error")) res13 <- rbind(res13, extract_term_row13(m1_bin, "plant_diversity", pair_id, "yi1", nrow(dat_m1_bin), length(unique(dat_m1_bin$S_ID)), "binary_psem_nested")) }
  
  dat_m2_bin <- dat[!is.na(dat$S_ID) & trimws(dat$S_ID) != "" & !is.na(dat$ID) & trimws(dat$ID) != "" & !is.na(dat$yi1) & is.finite(dat$yi1) & !is.na(dat$yi2) & is.finite(dat$yi2) & !is.na(dat$vi2) & is.finite(dat$vi2) & dat$vi2 > 0 & !is.na(dat$plant_diversity) & is.finite(dat$plant_diversity), , drop = FALSE]
  if (nrow(dat_m2_bin) > 1) { dat_m2_bin$S_ID <- factor(dat_m2_bin$S_ID); dat_m2_bin$ID <- factor(dat_m2_bin$ID); m2_bin <- fit_lme_safe(yi2 ~ yi1 + plant_diversity, dat_m2_bin, "vi2"); if (!inherits(m2_bin, "try-error")) res13 <- rbind(res13, extract_term_row13(m2_bin, "yi1", pair_id, "yi2", nrow(dat_m2_bin), length(unique(dat_m2_bin$S_ID)), "binary_psem_nested"), extract_term_row13(m2_bin, "plant_diversity", pair_id, "yi2", nrow(dat_m2_bin), length(unique(dat_m2_bin$S_ID)), "binary_psem_nested")) }
  
  dat_m3_bin <- dat[!is.na(dat$S_ID) & trimws(dat$S_ID) != "" & !is.na(dat$ID) & trimws(dat$ID) != "" & !is.na(dat$yi2) & is.finite(dat$yi2) & !is.na(dat$yi3) & is.finite(dat$yi3) & !is.na(dat$vi3) & is.finite(dat$vi3) & dat$vi3 > 0 & !is.na(dat$plant_diversity) & is.finite(dat$plant_diversity), , drop = FALSE]
  if (nrow(dat_m3_bin) > 1) { dat_m3_bin$S_ID <- factor(dat_m3_bin$S_ID); dat_m3_bin$ID <- factor(dat_m3_bin$ID); m3_bin <- fit_lme_safe(yi3 ~ yi2 + plant_diversity, dat_m3_bin, "vi3"); if (!inherits(m3_bin, "try-error")) res13 <- rbind(res13, extract_term_row13(m3_bin, "yi2", pair_id, "yi3", nrow(dat_m3_bin), length(unique(dat_m3_bin$S_ID)), "binary_psem_nested"), extract_term_row13(m3_bin, "plant_diversity", pair_id, "yi3", nrow(dat_m3_bin), length(unique(dat_m3_bin$S_ID)), "binary_psem_nested")) }
  
  dat_m1_add <- dat[!is.na(dat$S_ID) & trimws(dat$S_ID) != "" & !is.na(dat$ID) & trimws(dat$ID) != "" & !is.na(dat$yi1) & is.finite(dat$yi1) & !is.na(dat$vi1) & is.finite(dat$vi1) & dat$vi1 > 0 & !is.na(dat$added_species_number) & is.finite(dat$added_species_number), , drop = FALSE]
  if (nrow(dat_m1_add) > 1) { dat_m1_add$S_ID <- factor(dat_m1_add$S_ID); dat_m1_add$ID <- factor(dat_m1_add$ID); m1_add <- fit_lme_safe(yi1 ~ added_species_number, dat_m1_add, "vi1"); if (!inherits(m1_add, "try-error")) res13 <- rbind(res13, extract_term_row13(m1_add, "added_species_number", pair_id, "yi1", nrow(dat_m1_add), length(unique(dat_m1_add$S_ID)), "continuous_psem_nested")) }
  
  dat_m2_add <- dat[!is.na(dat$S_ID) & trimws(dat$S_ID) != "" & !is.na(dat$ID) & trimws(dat$ID) != "" & !is.na(dat$yi1) & is.finite(dat$yi1) & !is.na(dat$yi2) & is.finite(dat$yi2) & !is.na(dat$vi2) & is.finite(dat$vi2) & dat$vi2 > 0 & !is.na(dat$added_species_number) & is.finite(dat$added_species_number), , drop = FALSE]
  if (nrow(dat_m2_add) > 1) { dat_m2_add$S_ID <- factor(dat_m2_add$S_ID); dat_m2_add$ID <- factor(dat_m2_add$ID); m2_add <- fit_lme_safe(yi2 ~ yi1 + added_species_number, dat_m2_add, "vi2"); if (!inherits(m2_add, "try-error")) res13 <- rbind(res13, extract_term_row13(m2_add, "yi1", pair_id, "yi2", nrow(dat_m2_add), length(unique(dat_m2_add$S_ID)), "continuous_psem_nested"), extract_term_row13(m2_add, "added_species_number", pair_id, "yi2", nrow(dat_m2_add), length(unique(dat_m2_add$S_ID)), "continuous_psem_nested")) }
  
  dat_m3_add <- dat[!is.na(dat$S_ID) & trimws(dat$S_ID) != "" & !is.na(dat$ID) & trimws(dat$ID) != "" & !is.na(dat$yi2) & is.finite(dat$yi2) & !is.na(dat$yi3) & is.finite(dat$yi3) & !is.na(dat$vi3) & is.finite(dat$vi3) & dat$vi3 > 0 & !is.na(dat$added_species_number) & is.finite(dat$added_species_number), , drop = FALSE]
  if (nrow(dat_m3_add) > 1) { dat_m3_add$S_ID <- factor(dat_m3_add$S_ID); dat_m3_add$ID <- factor(dat_m3_add$ID); m3_add <- fit_lme_safe(yi3 ~ yi2 + added_species_number, dat_m3_add, "vi3"); if (!inherits(m3_add, "try-error")) res13 <- rbind(res13, extract_term_row13(m3_add, "yi2", pair_id, "yi3", nrow(dat_m3_add), length(unique(dat_m3_add$S_ID)), "continuous_psem_nested"), extract_term_row13(m3_add, "added_species_number", pair_id, "yi3", nrow(dat_m3_add), length(unique(dat_m3_add$S_ID)), "continuous_psem_nested")) }
}

openxlsx::write.xlsx(res13, file = paste0(path, "Table_S13_tri_trophic_NC_style.xlsx"), rowNames = FALSE)

#### 6. Stable 14 ####
paired_names <- read_txt_safe(paste0(prefix, "paired_names.txt"))
category_var <- "Dv_S"; categories <- c("Intercropping", "Sown_field_margins", "Cover_cropping")
res15 <- data.frame("Trophic interaction" = character(), "Category variable" = character(), "Category" = character(), "Predictor" = character(), "Response" = character(), "Number of observations" = numeric(), "Number of studies" = numeric(), "Estimate" = numeric(), "Std.Error of Estimate" = numeric(), "t-value" = numeric(), "P-value" = numeric(), "df" = numeric(), "Lower Bound of Confidence Interval" = numeric(), "Upper Bound of Confidence Interval" = numeric(), "R2" = numeric(), "SEM_type" = character(), "Fit_type" = character(), check.names = FALSE, stringsAsFactors = FALSE)
log_vec <- character(0)

for (i in seq_len(nrow(paired_names))) {
  pair_id <- paired_names[i, 1]; cat("Processing:", pair_id, "\n"); file_in <- paste0(prefix, pair_id, ".txt")
  if (!file.exists(file_in)) { log_vec <- c(log_vec, paste(pair_id, "| file not found")); next }
  dat <- read_txt_safe(file_in)
  num_cols <- c("yi1", "yi2", "yi3", "vi1", "vi2", "vi3", "NCP_N", "added_species_number"); for (cc in intersect(num_cols, names(dat))) dat[[cc]] <- safe_num(dat[[cc]])
  dat$S_ID <- as.character(dat$S_ID); dat$ID <- as.character(dat$ID); if ("NCP_N" %in% names(dat)) dat$plant_diversity <- ifelse(dat$NCP_N == 1, 0, ifelse(dat$NCP_N >= 2, 1, NA)) else dat$plant_diversity <- NA_real_
  if ("Ex_3" %in% names(dat)) { dat$Ex_3 <- trimws(as.character(dat$Ex_3)); dat$Ex_3 <- ifelse(tolower(dat$Ex_3) == "plot", "Plot", ifelse(tolower(dat$Ex_3) == "pot", "Pot", dat$Ex_3)) }
  if (!(category_var %in% names(dat))) { log_vec <- c(log_vec, paste(pair_id, "| missing category variable:", category_var)); next }
  
  for (category in categories) {
    cat("    Category:", category, "\n"); dat_cat <- dat[!is.na(dat[[category_var]]) & dat[[category_var]] == category, , drop = FALSE]
    if (nrow(dat_cat) == 0) { log_vec <- c(log_vec, paste(pair_id, "|", category, "| no data")); next }
    
    dat_m1_bin <- make_model_data(dat_cat, c("yi1", "vi1", "plant_diversity"))
    if (is_fitworthy(dat_m1_bin, "plant_diversity")) { dat_m1_bin$S_ID <- factor(dat_m1_bin$S_ID); dat_m1_bin$ID <- factor(dat_m1_bin$ID); n_obs_bin <- nrow(dat_m1_bin); n_study_bin <- length(unique(dat_m1_bin$S_ID)); m1_bin <- fit_lme_twolevel(yi1 ~ plant_diversity, "vi1", dat_m1_bin); if (is.null(m1_bin$model)) log_vec <- c(log_vec, paste(pair_id, "|", category, "| m1_bin failed at both levels")) else res15 <- append_result(res15, pair_id, category_var, category, "plant_diversity", "yi1", n_obs_bin, n_study_bin, m1_bin, "plant_diversity", "binary_psem_nested") } else log_vec <- c(log_vec, paste(pair_id, "|", category, "| m1_bin skipped: not fitworthy"))
    
    dat_m2_bin <- make_model_data(dat_cat, c("yi1", "yi2", "vi2", "plant_diversity"))
    if (is_fitworthy(dat_m2_bin, "plant_diversity")) { dat_m2_bin$S_ID <- factor(dat_m2_bin$S_ID); dat_m2_bin$ID <- factor(dat_m2_bin$ID); n_obs_bin <- nrow(dat_m2_bin); n_study_bin <- length(unique(dat_m2_bin$S_ID)); m2_bin <- fit_lme_twolevel(yi2 ~ yi1 + plant_diversity, "vi2", dat_m2_bin); if (is.null(m2_bin$model)) log_vec <- c(log_vec, paste(pair_id, "|", category, "| m2_bin failed at both levels")) else { res15 <- append_result(res15, pair_id, category_var, category, "yi1", "yi2", n_obs_bin, n_study_bin, m2_bin, "yi1", "binary_psem_nested"); res15 <- append_result(res15, pair_id, category_var, category, "plant_diversity", "yi2", n_obs_bin, n_study_bin, m2_bin, "plant_diversity", "binary_psem_nested") } } else log_vec <- c(log_vec, paste(pair_id, "|", category, "| m2_bin skipped: not fitworthy"))
    
    dat_m3_bin <- make_model_data(dat_cat, c("yi2", "yi3", "vi3", "plant_diversity"))
    if (is_fitworthy(dat_m3_bin, "plant_diversity")) { dat_m3_bin$S_ID <- factor(dat_m3_bin$S_ID); dat_m3_bin$ID <- factor(dat_m3_bin$ID); n_obs_bin <- nrow(dat_m3_bin); n_study_bin <- length(unique(dat_m3_bin$S_ID)); m3_bin <- fit_lme_twolevel(yi3 ~ yi2 + plant_diversity, "vi3", dat_m3_bin); if (is.null(m3_bin$model)) log_vec <- c(log_vec, paste(pair_id, "|", category, "| m3_bin failed at both levels")) else { res15 <- append_result(res15, pair_id, category_var, category, "yi2", "yi3", n_obs_bin, n_study_bin, m3_bin, "yi2", "binary_psem_nested"); res15 <- append_result(res15, pair_id, category_var, category, "plant_diversity", "yi3", n_obs_bin, n_study_bin, m3_bin, "plant_diversity", "binary_psem_nested") } } else log_vec <- c(log_vec, paste(pair_id, "|", category, "| m3_bin skipped: not fitworthy"))
    
    if (!("added_species_number" %in% names(dat_cat))) { log_vec <- c(log_vec, paste(pair_id, "|", category, "| missing added_species_number")); next }
    
    dat_m1_add <- make_model_data(dat_cat, c("yi1", "vi1", "added_species_number"))
    if (is_fitworthy(dat_m1_add, "added_species_number")) { dat_m1_add$S_ID <- factor(dat_m1_add$S_ID); dat_m1_add$ID <- factor(dat_m1_add$ID); n_obs_add <- nrow(dat_m1_add); n_study_add <- length(unique(dat_m1_add$S_ID)); m1_add <- fit_lme_twolevel(yi1 ~ added_species_number, "vi1", dat_m1_add); if (is.null(m1_add$model)) log_vec <- c(log_vec, paste(pair_id, "|", category, "| m1_add failed at both levels")) else res15 <- append_result(res15, pair_id, category_var, category, "added_species_number", "yi1", n_obs_add, n_study_add, m1_add, "added_species_number", "continuous_psem_nested") } else log_vec <- c(log_vec, paste(pair_id, "|", category, "| m1_add skipped: not fitworthy"))
    
    dat_m2_add <- make_model_data(dat_cat, c("yi1", "yi2", "vi2", "added_species_number"))
    if (is_fitworthy(dat_m2_add, "added_species_number")) { dat_m2_add$S_ID <- factor(dat_m2_add$S_ID); dat_m2_add$ID <- factor(dat_m2_add$ID); n_obs_add <- nrow(dat_m2_add); n_study_add <- length(unique(dat_m2_add$S_ID)); m2_add <- fit_lme_twolevel(yi2 ~ yi1 + added_species_number, "vi2", dat_m2_add); if (is.null(m2_add$model)) log_vec <- c(log_vec, paste(pair_id, "|", category, "| m2_add failed at both levels")) else { res15 <- append_result(res15, pair_id, category_var, category, "yi1", "yi2", n_obs_add, n_study_add, m2_add, "yi1", "continuous_psem_nested"); res15 <- append_result(res15, pair_id, category_var, category, "added_species_number", "yi2", n_obs_add, n_study_add, m2_add, "added_species_number", "continuous_psem_nested") } } else log_vec <- c(log_vec, paste(pair_id, "|", category, "| m2_add skipped: not fitworthy"))
    
    dat_m3_add <- make_model_data(dat_cat, c("yi2", "yi3", "vi3", "added_species_number"))
    if (is_fitworthy(dat_m3_add, "added_species_number")) { dat_m3_add$S_ID <- factor(dat_m3_add$S_ID); dat_m3_add$ID <- factor(dat_m3_add$ID); n_obs_add <- nrow(dat_m3_add); n_study_add <- length(unique(dat_m3_add$S_ID)); m3_add <- fit_lme_twolevel(yi3 ~ yi2 + added_species_number, "vi3", dat_m3_add); if (is.null(m3_add$model)) log_vec <- c(log_vec, paste(pair_id, "|", category, "| m3_add failed at both levels")) else { res15 <- append_result(res15, pair_id, category_var, category, "yi2", "yi3", n_obs_add, n_study_add, m3_add, "yi2", "continuous_psem_nested"); res15 <- append_result(res15, pair_id, category_var, category, "added_species_number", "yi3", n_obs_add, n_study_add, m3_add, "added_species_number", "continuous_psem_nested") } } else log_vec <- c(log_vec, paste(pair_id, "|", category, "| m3_add skipped: not fitworthy"))
  }
}

if (nrow(res15) == 0) stop("No valid results were generated for this category analysis.")
openxlsx::write.xlsx(res15, paste0(path, "Stable14_NC_style_", category_var, ".xlsx"), rowNames = FALSE)
if (length(log_vec) > 0) writeLines(log_vec, con = paste0(path, "Stable14_NC_style_", category_var, "_log.txt"))


#### 7. Stable 16 ####
############################################################
## STABLE 16
## Tri-trophic subgroup analysis
## (moved forward; originally STABLE 17)
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
paired_names <- read.table(
  paste0(prefix, "paired_names.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE,
  blank.lines.skip = FALSE
)

#==================================================
# 2. Shared-control grouping parameters
#==================================================
round_digits_m  <- 4
round_digits_sd <- 4
round_digits_r  <- 0

#==================================================
# 3. Subgroup definitions
#==================================================
subgroup_list <- list(
  CP_Ct = c("Food_crop", "Cash_crop"),
  CP_C3 = c("Herbaceous_plant", "Woody_plant"),
  Mn__  = c("Generalist", "Specialist"),
  NE_c  = c("Generalist", "Specialist"),
  Dv_S = c("Intercropping", "Cover_cropping", "Sown_field_margins"),
  Ex_4  = c("Natural"),
  Ex_3  = c("Plot", "Pot"),
  Zone  = c("Temperate", "Tropic")
)

#==================================================
# 4. Output container
#==================================================
res <- data.frame(
  Pair_ID = character(),
  Trophic_interaction = character(),
  Moderator = character(),
  Item = character(),
  Category = character(),
  Number_of_observations = numeric(),
  Number_of_studies = numeric(),
  V_type = character(),
  Min_eigen_before_fix = numeric(),
  Effect_size = numeric(),
  t_value = numeric(),
  P_value = numeric(),
  df = numeric(),
  CI_lb = numeric(),
  CI_ub = numeric(),
  stringsAsFactors = FALSE
)

#==================================================
# 5. Standardization helpers
#==================================================
standardize_zone <- function(x) {
  x <- trimws(as.character(x))
  x_low <- tolower(x)
  x_out <- x
  x_out[x_low %in% c("tropical", "tropics", "tropic")] <- "Tropic"
  x_out[x_low %in% c("temperate")] <- "Temperate"
  x_out[!(x_low %in% c("tropical", "tropics", "tropic", "temperate"))] <- x[!(x_low %in% c("tropical", "tropics", "tropic", "temperate"))]
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
# 6. Generic helpers
#==================================================
safe_num <- function(x) suppressWarnings(as.numeric(x))
round_for_group <- function(x, digits = 4) if (is.numeric(x)) round(x, digits = digits) else x
make_common_key <- function(m, sd, r) interaction(m, sd, r, drop = TRUE, lex.order = TRUE)

get_control_cols <- function(dat, idx) {
  m_col <- c(paste0("Cntrl.m", idx), "CONTROL_M", "Cntrl.m", "control_m")
  s_col <- c(paste0("Cntrl.s", idx), "CONTROL_SD", "Cntrl.s", "control_sd")
  r_col <- c(paste0("Cntrl.n", idx), "CONTROL_R", "Cntrl.n", "control_n")
  list(
    m = m_col[m_col %in% names(dat)][1],
    s = s_col[s_col %in% names(dat)][1],
    r = r_col[r_col %in% names(dat)][1]
  )
}

fit_rma_mv_safe <- function(yi_col, V_input, dat, method = "ML") {
  random_terms <- list(~1 | S_ID, ~1 | ID_global)
  for (opt in c("BFGS", "nlminb", "optim")) {
    fit <- try(
      rma.mv(
        as.formula(paste0(yi_col, " ~ 1")),
        V = V_input,
        random = random_terms,
        data = dat,
        test = "t",
        method = method,
        control = list(optimizer = opt)
      ),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) return(fit)
  }
  "try-error"
}

build_shared_V <- function(dat, vi_col, m_col, s_col, r_col) {
  block_V <- function(x) {
    cpart <- x[[s_col]][1]^2 / (x[[r_col]][1] * x[[m_col]][1]^2)
    V <- matrix(cpart, nrow = nrow(x), ncol = nrow(x))
    diag(V) <- x[[vi_col]]
    V
  }
  
  out <- list(V = NULL, V_type = NA_character_, min_eigen = NA_real_, ok = FALSE)
  V_try <- try(as.matrix(Matrix::bdiag(lapply(split(dat, dat$Common_ID_global), block_V))), silent = TRUE)
  if (inherits(V_try, "try-error")) return(out)
  if (!isSymmetric(V_try) || any(is.na(V_try)) || any(is.infinite(V_try))) return(out)
  
  eig <- try(eigen(V_try, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (inherits(eig, "try-error")) return(out)
  
  out$min_eigen <- min(eig)
  
  if (!is.finite(out$min_eigen) || out$min_eigen <= 0) {
    V_fix <- try(corpcor::make.positive.definite(V_try), silent = TRUE)
    if (inherits(V_fix, "try-error")) return(out)
    out$V <- as.matrix(V_fix)
    out$V_type <- "shared_control_pd_fix"
  } else {
    out$V <- V_try
    out$V_type <- "shared_control"
  }
  
  out$ok <- TRUE
  out
}

empty_tri_row <- function(pair_id, trophic_interaction, moderator_name, level_name, response_name, n_obs = NA, n_study = NA, v_type = NA, min_eigen = NA) {
  data.frame(
    Pair_ID = pair_id,
    Trophic_interaction = trophic_interaction,
    Moderator = moderator_name,
    Item = level_name,
    Category = response_name,
    Number_of_observations = n_obs,
    Number_of_studies = n_study,
    V_type = v_type,
    Min_eigen_before_fix = min_eigen,
    Effect_size = NA,
    t_value = NA,
    P_value = NA,
    df = NA,
    CI_lb = NA,
    CI_ub = NA,
    stringsAsFactors = FALSE
  )
}

run_one_response_tri <- function(dat0, yi_col, vi_col, response_name, pair_id, trophic_interaction, idx, moderator_name, level_name) {
  ctrl <- get_control_cols(dat0, idx)
  if (is.na(ctrl$m) || is.na(ctrl$s) || is.na(ctrl$r)) return(empty_tri_row(pair_id, trophic_interaction, moderator_name, level_name, response_name))
  
  dat <- dat0
  dat[[yi_col]]  <- safe_num(dat[[yi_col]])
  dat[[vi_col]]  <- safe_num(dat[[vi_col]])
  dat[[ctrl$m]]  <- safe_num(dat[[ctrl$m]])
  dat[[ctrl$s]]  <- safe_num(dat[[ctrl$s]])
  dat[[ctrl$r]]  <- safe_num(dat[[ctrl$r]])
  
  keep <- !is.na(dat[[yi_col]]) & is.finite(dat[[yi_col]]) &
    !is.na(dat[[vi_col]]) & is.finite(dat[[vi_col]]) & dat[[vi_col]] > 0 &
    !is.na(dat$S_ID) & trimws(as.character(dat$S_ID)) != "" &
    !is.na(dat$ID) & trimws(as.character(dat$ID)) != "" &
    !is.na(dat[[ctrl$m]]) & is.finite(dat[[ctrl$m]]) & dat[[ctrl$m]] > 0 &
    !is.na(dat[[ctrl$s]]) & is.finite(dat[[ctrl$s]]) & dat[[ctrl$s]] >= 0 &
    !is.na(dat[[ctrl$r]]) & is.finite(dat[[ctrl$r]]) & dat[[ctrl$r]] > 0
  
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) == 0) return(empty_tri_row(pair_id, trophic_interaction, moderator_name, level_name, response_name, n_obs = 0))
  
  dat$CONTROL_M_use  <- dat[[ctrl$m]]
  dat$CONTROL_SD_use <- dat[[ctrl$s]]
  dat$CONTROL_R_use  <- dat[[ctrl$r]]
  
  dat$CONTROL_M_grp  <- round(dat$CONTROL_M_use, round_digits_m)
  dat$CONTROL_SD_grp <- round(dat$CONTROL_SD_use, round_digits_sd)
  dat$CONTROL_R_grp  <- round(dat$CONTROL_R_use, round_digits_r)
  
  dat$Common_ID <- ave(
    seq_len(nrow(dat)), dat$S_ID,
    FUN = function(ii) match(
      interaction(dat$CONTROL_M_grp[ii], dat$CONTROL_SD_grp[ii], dat$CONTROL_R_grp[ii], drop = TRUE, lex.order = TRUE),
      unique(interaction(dat$CONTROL_M_grp[ii], dat$CONTROL_SD_grp[ii], dat$CONTROL_R_grp[ii], drop = TRUE, lex.order = TRUE))
    )
  )
  
  dat$Common_ID_global <- match(
    interaction(dat$S_ID, dat$Common_ID, drop = TRUE, lex.order = TRUE),
    unique(interaction(dat$S_ID, dat$Common_ID, drop = TRUE, lex.order = TRUE))
  )
  
  dat$ID_global <- interaction(dat$S_ID, dat$ID, drop = TRUE, lex.order = TRUE)
  dat <- dat[!is.na(dat$ID_global), , drop = FALSE]
  if (nrow(dat) == 0) return(empty_tri_row(pair_id, trophic_interaction, moderator_name, level_name, response_name, n_obs = 0))
  
  check_common <- aggregate(
    cbind(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp) ~ Common_ID_global,
    data = dat,
    FUN = function(x) length(unique(x))
  )
  
  use_V <- all(check_common$CONTROL_M_grp == 1) &&
    all(check_common$CONTROL_SD_grp == 1) &&
    all(check_common$CONTROL_R_grp == 1)
  
  dat.fit <- NULL
  V_type <- NA_character_
  min_eigen <- NA_real_
  
  if (nrow(dat) == 1) {
    dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "z"), silent = TRUE)
    if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi_single_obs"
  } else {
    if (use_V) {
      V_info <- build_shared_V(dat, vi_col, "CONTROL_M_use", "CONTROL_SD_use", "CONTROL_R_use")
      min_eigen <- V_info$min_eigen
      if (isTRUE(V_info$ok) && nrow(V_info$V) == nrow(dat)) {
        dat.fit <- fit_rma_mv_safe(yi_col, V_info$V, dat, method = "ML")
        if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_safe(yi_col, V_info$V, dat, method = "REML")
        if (!inherits(dat.fit, "try-error")) V_type <- V_info$V_type
      }
    }
    
    if (is.null(dat.fit) || inherits(dat.fit, "try-error")) {
      if (length(unique(dat$S_ID)) == 1) {
        dat.fit <- try(rma(yi = dat[[yi_col]], vi = dat[[vi_col]], test = "t"), silent = TRUE)
        if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi_rma"
      } else {
        dat.fit <- fit_rma_mv_safe(yi_col, dat[[vi_col]], dat, method = "ML")
        if (inherits(dat.fit, "try-error")) dat.fit <- fit_rma_mv_safe(yi_col, dat[[vi_col]], dat, method = "REML")
        if (!inherits(dat.fit, "try-error")) V_type <- "diagonal_vi"
      }
    }
  }
  
  if (is.null(dat.fit) || inherits(dat.fit, "try-error")) {
    return(empty_tri_row(pair_id, trophic_interaction, moderator_name, level_name, response_name, n_obs = nrow(dat), n_study = length(unique(dat$S_ID)), min_eigen = min_eigen))
  }
  
  sm <- coef(summary(dat.fit))
  t_value_out <- if ("tval" %in% colnames(sm)) sm[1, "tval"] else if ("zval" %in% colnames(sm)) sm[1, "zval"] else NA
  p_value_out <- if ("pval" %in% colnames(sm)) sm[1, "pval"] else NA
  df_out      <- if ("df"   %in% colnames(sm)) sm[1, "df"] else NA
  
  data.frame(
    Pair_ID = pair_id,
    Trophic_interaction = trophic_interaction,
    Moderator = moderator_name,
    Item = level_name,
    Category = response_name,
    Number_of_observations = nrow(dat),
    Number_of_studies = length(unique(dat$S_ID)),
    V_type = V_type,
    Min_eigen_before_fix = min_eigen,
    Effect_size = as.numeric(dat.fit$b[1]),
    t_value = as.numeric(t_value_out),
    P_value = as.numeric(p_value_out),
    df = as.numeric(df_out),
    CI_lb = if (!is.null(dat.fit$ci.lb)) dat.fit$ci.lb[1] else NA,
    CI_ub = if (!is.null(dat.fit$ci.ub)) dat.fit$ci.ub[1] else NA,
    stringsAsFactors = FALSE
  )
}

#==================================================
# 7. Main loop
#==================================================
for (i in seq_len(nrow(paired_names))) {
  pair_id <- paired_names$Pair_ID[i]
  file_in <- paste0(prefix, pair_id, ".txt")
  if (!file.exists(file_in)) next
  
  dat0 <- read.table(
    file_in,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    blank.lines.skip = FALSE
  )
  
  if ("Ex_3" %in% names(dat0)) dat0$Ex_3 <- standardize_ex3(dat0$Ex_3)
  if ("Zone" %in% names(dat0)) dat0$Zone <- standardize_zone(dat0$Zone)
  
  for (mod_name in names(subgroup_list)) {
    if (!(mod_name %in% names(dat0))) next
    
    for (lv in subgroup_list[[mod_name]]) {
      dat_sub <- dat0[!is.na(dat0[[mod_name]]) & as.character(dat0[[mod_name]]) == lv, , drop = FALSE]
      if (nrow(dat_sub) == 0) next
      
      res <- rbind(
        res,
        run_one_response_tri(dat_sub, "yi1", "vi1", paired_names$Response1[i], pair_id, paired_names$Trophic_interaction[i], 1, mod_name, lv),
        run_one_response_tri(dat_sub, "yi2", "vi2", paired_names$Response2[i], pair_id, paired_names$Trophic_interaction[i], 2, mod_name, lv)
      )
      
      if ("yi3" %in% names(dat_sub) && "vi3" %in% names(dat_sub) && "Response3" %in% names(paired_names)) {
        res <- rbind(
          res,
          run_one_response_tri(dat_sub, "yi3", "vi3", paired_names$Response3[i], pair_id, paired_names$Trophic_interaction[i], 3, mod_name, lv)
        )
      }
    }
  }
}

#==================================================
# 8. Output
#==================================================
openxlsx::write.xlsx(res, paste0(path, "STABLE16_tri_trophic_by_subgroup.xlsx"), rowNames = FALSE)

print(res)
cat("Finished STABLE16.\n")
