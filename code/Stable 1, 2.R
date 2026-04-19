#### STABLE 1 AND STABLE 2####
rm(list = ls())

suppressPackageStartupMessages({ 
  library(dplyr);
  library(metafor);
  library(Matrix); 
  library(orchaRd); 
  library(rotl); 
  library(ape); 
  library(stringr) })

#### 0. Configuration ####
cfg <- list(prefix = "/Users/yusha/", 
            path = "/Users/yusha/", 
            data_file = "/Users/yusha/Alldata.txt", 
            factor_mods = c("CP_C3", "CP_Ct", "Mn__", "NE_c", "Dv_S", "Ex_4", "Ex_3", "Zone"), cont_mod = "added_species_number", 
            round_digits_s1 = c(m = 4, sd = 4, r = 0), round_digits_s2 = c(m = 6, sd = 6, r = 6), 
            phy_round_digits = c(m = 4, sd = 4, r = 0), 
            phy_cache_tax = "/Users/yusha/tax_lookup.rds", 
            phy_cache_tree = "/Users/yusha/Plant_Tree.rds", 
            phy_cache_cor = "/Users/yusha/Plant_Cor.rds")
cfg$mods_to_test <- c(cfg$factor_mods, cfg$cont_mod)

#### 1. Generic helpers ####
standardize_zone <- function(x) { x <- trimws(as.character(x)); x_low <- tolower(x); x_low[x_low %in% c("tropical", "tropics", "tropic")] <- "tropic"; x_low[x_low %in% c("temperate")] <- "temperate"; out <- x; out[x_low == "tropic"] <- "Tropic"; out[x_low == "temperate"] <- "Temperate"; out }
standardize_ex3 <- function(x) { x <- trimws(as.character(x)); x_low <- tolower(x); x_low[x_low %in% c("plot")] <- "plot"; x_low[x_low %in% c("pot")] <- "pot"; out <- x; out[x_low == "plot"] <- "Plot"; out[x_low == "pot"] <- "Pot"; out }
round_for_group <- function(x, digits = 4) { x <- as.numeric(x); ifelse(is.na(x), NA, round(x, digits)) }
make_common_key <- function(m, sd, r) paste(m, sd, r, sep = "_")
clean_species_name <- function(x) { x <- as.character(x); x <- trimws(x); x <- gsub("_", " ", x); x <- gsub("×", " x ", x, fixed = TRUE); x <- gsub("\\bX\\b", " x ", x); x <- stringr::str_squish(x); x <- tolower(x); x[x %in% c("", "na", "n/a", "null")] <- NA; x }
safe_i2 <- function(model) { out <- try(i2_ml(model), silent = TRUE); if (inherits(out, "try-error")) return(NA); out }

extract_anova_row <- function(anv, formula_text, model_label, k_value, i2_value) {
  if (inherits(anv, "try-error") || is.null(anv)) return(data.frame(formula = formula_text, Model = model_label, AIC = NA, logLik = NA, LRT = NA, df = NA, pval = NA, k = k_value, I2 = i2_value, stringsAsFactors = FALSE))
  data.frame(formula = formula_text, Model = model_label, AIC = tryCatch(anv[["fit.stats.f"]][["AIC"]], error = function(e) NA), logLik = tryCatch(anv[["fit.stats.f"]][["ll"]], error = function(e) NA), LRT = tryCatch(anv[["LRT"]], error = function(e) NA), df = tryCatch(anv[["parms.f"]], error = function(e) NA), pval = tryCatch(anv[["pval"]], error = function(e) NA), k = k_value, I2 = i2_value, stringsAsFactors = FALSE)
}

retry_api <- function(fun, ..., max_try = 5, sleep_sec = 8, label = "API call") { last_err <- NULL; for (i in seq_len(max_try)) { cat(label, "- try", i, "of", max_try, "\n"); out <- try(fun(...), silent = TRUE); if (!inherits(out, "try-error")) return(out); last_err <- out; cat("  failed:", as.character(out), "\n"); if (i < max_try) Sys.sleep(sleep_sec) }; stop(paste0(label, " failed after ", max_try, " tries.\n", as.character(last_err))) }

tnrs_match_names_chunked <- function(species_vec, chunk_size = 20, max_try = 5, sleep_sec = 8) { species_vec <- unique(species_vec); chunks <- split(species_vec, ceiling(seq_along(species_vec) / chunk_size)); out_list <- vector("list", length(chunks)); for (i in seq_along(chunks)) { cat("TNRS chunk", i, "of", length(chunks), "\n"); out_list[[i]] <- retry_api(rotl::tnrs_match_names, names = chunks[[i]], max_try = max_try, sleep_sec = sleep_sec, label = paste0("tnrs_match_names chunk ", i)); Sys.sleep(2) }; out <- do.call(rbind, out_list); rownames(out) <- NULL; out }

build_subtree_prune_bad_ott <- function(tax_lookup, max_round = 20) {
  tax_lookup2 <- tax_lookup
  for (round_i in seq_len(max_round)) {
    cat("tol_induced_subtree round", round_i, "\n")
    tr <- try(rotl::tol_induced_subtree(ott_ids = tax_lookup2$ott_id, label_format = "id"), silent = TRUE)
    if (!inherits(tr, "try-error")) return(list(tree = tr, tax_lookup = tax_lookup2))
    err_txt <- as.character(tr); cat("  failed:", err_txt, "\n")
    if (!grepl("pruned_ott_id", err_txt)) stop(tr)
    bad_ott <- regmatches(err_txt, regexpr("ott[0-9]+", err_txt))
    if (length(bad_ott) == 0 || is.na(bad_ott)) stop("Could not parse pruned ott_id from error.")
    bad_ott_num <- suppressWarnings(as.numeric(gsub("^ott", "", bad_ott)))
    if (is.na(bad_ott_num)) stop("Parsed bad ott_id could not be converted to numeric.")
    cat("  removing pruned ott_id:", bad_ott_num, "\n")
    tax_lookup2 <- tax_lookup2[tax_lookup2$ott_id != bad_ott_num, , drop = FALSE]
    if (nrow(tax_lookup2) < 2) stop("Too few taxa remain after removing pruned ott_ids.")
    Sys.sleep(1)
  }
  stop("Failed to build subtree after repeatedly removing pruned ott_ids.")
}

#### 2. Data preparation ####
read_master_data <- function(cfg) read.table(cfg$data_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

prepare_master_data <- function(Data, cfg, include_group_var = c("Trophic_group", "Trophic_group_response", "both"), include_sqrtnum = FALSE, include_phylogeny = FALSE) {
  include_group_var <- match.arg(include_group_var)
  Data <- Data %>% mutate(Zone = standardize_zone(Zone), Ex_3 = standardize_ex3(Ex_3), Trophic_group = trimws(as.character(Trophic_group)), Trophic_group_response = trimws(as.character(Trophic_group_response)), added_species_number = as.numeric(added_species_number))
  if (include_sqrtnum) Data <- Data %>% mutate(SqrtNum = sqrt((CONTROL_R + TREATMENT_R) / (CONTROL_R * TREATMENT_R)))
  if (include_phylogeny) Data <- Data %>% mutate(Plant_Species_raw = clean_species_name(CP_Ltn_nm_f_))
  keep_filter <- c(!is.na(Data$yi) & is.finite(Data$yi), !is.na(Data$vi) & is.finite(Data$vi) & Data$vi > 0, !is.na(Data$S_ID), !is.na(Data$CONTROL_M) & is.finite(Data$CONTROL_M) & Data$CONTROL_M > 0, !is.na(Data$TREATMENT_M) & is.finite(Data$TREATMENT_M) & Data$TREATMENT_M > 0, !is.na(Data$CONTROL_R) & is.finite(Data$CONTROL_R) & Data$CONTROL_R > 0, !is.na(Data$TREATMENT_R) & is.finite(Data$TREATMENT_R) & Data$TREATMENT_R > 0, !is.na(Data$CONTROL_SD) & is.finite(Data$CONTROL_SD) & Data$CONTROL_SD >= 0, !is.na(Data$TREATMENT_SD) & is.finite(Data$TREATMENT_SD) & Data$TREATMENT_SD >= 0)
  if (include_group_var %in% c("Trophic_group", "both")) keep_filter <- c(keep_filter, !is.na(Data$Trophic_group))
  if (include_group_var %in% c("Trophic_group_response", "both")) keep_filter <- c(keep_filter, !is.na(Data$Trophic_group_response))
  if (include_sqrtnum) keep_filter <- c(keep_filter, !is.na(Data$SqrtNum) & is.finite(Data$SqrtNum))
  if (include_phylogeny) keep_filter <- c(keep_filter, !is.na(Data$Plant_Species_raw))
  Data <- Data[Reduce("&", keep_filter), , drop = FALSE] %>% arrange(S_ID) %>% mutate(Obs_ID = seq_len(n()), S_ID = factor(S_ID), Trophic_group = factor(Trophic_group), Trophic_group_response = factor(Trophic_group_response))
  for (v in cfg$factor_mods) if (v %in% names(Data)) Data[[v]] <- factor(Data[[v]])
  Data
}

#### 3. Shared-control V ####
build_shared_V <- function(dat, round_digits = c(m = 4, sd = 4, r = 0)) {
  dat <- dat %>% mutate(CONTROL_M_grp = round_for_group(CONTROL_M, round_digits[["m"]]), CONTROL_SD_grp = round_for_group(CONTROL_SD, round_digits[["sd"]]), CONTROL_R_grp = round_for_group(CONTROL_R, round_digits[["r"]])) %>% group_by(S_ID) %>% mutate(Common_ID = match(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp), unique(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp)))) %>% ungroup() %>% mutate(Common_ID_global = match(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE), unique(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE))))
  calc.v <- function(x) { cpart <- x$CONTROL_SD[1]^2 / (x$CONTROL_R[1] * x$CONTROL_M[1]^2); v <- matrix(cpart, nrow = nrow(x), ncol = nrow(x)); diag(v) <- x$vi; v }
  V <- bldiag(lapply(split(dat, dat$Common_ID_global), calc.v)); list(dat = dat, V = as.matrix(V))
}

#### 4. Model fitting ####
fit_rma_mv <- function(formula_mod, dat, V_mat) try(rma.mv(yi, V = V_mat, mods = formula_mod, random = list(~1 | S_ID, ~1 | ID), data = dat, method = "ML", control = list(optimizer = "BFGS")), silent = TRUE)
fit_rma_mv_phylo <- function(formula_mod, dat, V_mat, phy_cor_sub) try(rma.mv(yi, V = V_mat, mods = formula_mod, random = list(~1 | S_ID, ~1 | Plant_Species, ~1 | ID, ~1 | Plant_Species2), R = list(Plant_Species = phy_cor_sub), data = dat, method = "ML", control = list(optimizer = "BFGS")), silent = TRUE)

#### 5. Phylogeny ####
build_or_load_phylogeny <- function(Whole_data, cfg) {
  species_all <- sort(unique(Whole_data$Plant_Species_raw)); cat("Plant species for phylogeny:", length(species_all), "\n"); if (length(species_all) < 2) stop("Too few valid plant species names after cleaning.")
  if (file.exists(cfg$phy_cache_tax) && file.exists(cfg$phy_cache_tree) && file.exists(cfg$phy_cache_cor)) return(list(tax_lookup = readRDS(cfg$phy_cache_tax), Plant_Tree = readRDS(cfg$phy_cache_tree), Plant_Cor = readRDS(cfg$phy_cache_cor)))
  Taxa_Plant <- tnrs_match_names_chunked(species_all, chunk_size = 20, max_try = 5, sleep_sec = 8)
  if (!("ott_id" %in% colnames(Taxa_Plant))) stop("tnrs_match_names() did not return an ott_id column.")
  Taxa_Plant <- Taxa_Plant[!is.na(Taxa_Plant$ott_id), , drop = FALSE]; if (nrow(Taxa_Plant) < 2) stop("Too few taxa returned valid OTT IDs.")
  matched_name <- Taxa_Plant$unique_name; matched_name[is.na(matched_name) | matched_name == ""] <- Taxa_Plant$search_string; matched_name <- clean_species_name(matched_name)
  tax_lookup <- data.frame(search_string_clean = clean_species_name(Taxa_Plant$search_string), matched_name_clean = matched_name, ott_id = Taxa_Plant$ott_id, stringsAsFactors = FALSE)
  tax_lookup <- tax_lookup[!duplicated(tax_lookup$search_string_clean), , drop = FALSE]; tax_lookup <- tax_lookup[!duplicated(tax_lookup$ott_id), , drop = FALSE]
  subtree_res <- build_subtree_prune_bad_ott(tax_lookup); Plant_Tree <- subtree_res$tree; tax_lookup <- subtree_res$tax_lookup
  tree_ott <- suppressWarnings(as.numeric(gsub("^ott", "", Plant_Tree$tip.label))); tip_match <- match(tree_ott, tax_lookup$ott_id); tip_species <- clean_species_name(tax_lookup$matched_name_clean[tip_match])
  keep_tip <- !is.na(tip_species) & tip_species != ""; if (sum(keep_tip) < 2) stop("After mapping tree tip labels to species names, fewer than 2 tips remain.")
  if (any(!keep_tip)) { Plant_Tree <- ape::drop.tip(Plant_Tree, Plant_Tree$tip.label[!keep_tip]); tip_species <- tip_species[keep_tip] }
  Plant_Tree$tip.label <- tip_species; dup_tip <- duplicated(Plant_Tree$tip.label); if (any(dup_tip)) Plant_Tree <- ape::drop.tip(Plant_Tree, Plant_Tree$tip.label[dup_tip])
  Plant_Cor <- ape::vcv(ape::compute.brlen(Plant_Tree), corr = TRUE); rownames(Plant_Cor) <- clean_species_name(rownames(Plant_Cor)); colnames(Plant_Cor) <- clean_species_name(colnames(Plant_Cor)); Plant_Cor <- Plant_Cor[order(rownames(Plant_Cor)), order(colnames(Plant_Cor)), drop = FALSE]
  saveRDS(tax_lookup, cfg$phy_cache_tax); saveRDS(Plant_Tree, cfg$phy_cache_tree); saveRDS(Plant_Cor, cfg$phy_cache_cor)
  list(tax_lookup = tax_lookup, Plant_Tree = Plant_Tree, Plant_Cor = Plant_Cor)
}

attach_phylogeny_labels <- function(Whole_data, tax_lookup, Plant_Cor, cfg) {
  write.table(data.frame(original_name = unique(Whole_data$CP_Ltn_nm_f_), cleaned_name = clean_species_name(unique(Whole_data$CP_Ltn_nm_f_))), paste0(cfg$path, "Cleaned_plant_names_check.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  Whole_data$Plant_Species_matched <- clean_species_name(tax_lookup$matched_name_clean[match(Whole_data$Plant_Species_raw, tax_lookup$search_string_clean)])
  species_not_in_tree <- setdiff(sort(unique(Whole_data$Plant_Species_raw)), sort(unique(Whole_data$Plant_Species_raw[!is.na(Whole_data$Plant_Species_matched) & Whole_data$Plant_Species_matched %in% rownames(Plant_Cor)])))
  if (length(species_not_in_tree) > 0) { cat("Species not matched into phylogeny:\n"); print(species_not_in_tree) }
  write.table(data.frame(species_not_in_tree = species_not_in_tree), paste0(cfg$path, "Plant_species_not_in_phylogeny.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  Whole_data <- Whole_data %>% filter(!is.na(Plant_Species_matched), Plant_Species_matched %in% rownames(Plant_Cor))
  cat("Rows retained after phylogeny matching:", nrow(Whole_data), "\n"); cat("Unique matched plant species:", length(unique(Whole_data$Plant_Species_matched)), "\n"); if (nrow(Whole_data) == 0) stop("No rows retained after phylogeny matching.")
  Whole_data
}

align_phylo_subset <- function(dat, phy_cor_full) { sp_keep <- sort(unique(dat$Plant_Species_matched)); sp_keep <- sp_keep[sp_keep %in% rownames(phy_cor_full)]; dat <- dat %>% filter(Plant_Species_matched %in% sp_keep); phy_sub <- phy_cor_full[sp_keep, sp_keep, drop = FALSE]; dat$Plant_Species <- factor(dat$Plant_Species_matched, levels = sp_keep); dat$Plant_Species2 <- dat$Plant_Species; list(dat = dat, phy_cor = phy_sub) }

#### 6. Stable 1 engine: moderator model comparison ####
run_stable1_core <- function(Whole_data, group_var, mods_to_test, factor_mods, round_digits, output_file, phy_cor_full = NULL) {
  use_phy <- !is.null(phy_cor_full)
  base_data <- Whole_data %>% filter(!is.na(.data[[group_var]])) %>% droplevels(); base_data[[group_var]] <- factor(base_data[[group_var]])
  if (use_phy) { phy0 <- align_phylo_subset(base_data, phy_cor_full); base_data <- phy0$dat; phy_cor0 <- phy0$phy_cor }
  tmp0 <- build_shared_V(base_data, round_digits = round_digits); base_data0 <- tmp0$dat; V0 <- tmp0$V
  if (use_phy) { phy0b <- align_phylo_subset(base_data0, phy_cor0); base_data0 <- phy0b$dat; phy_cor0 <- phy0b$phy_cor }
  null_formula <- as.formula(paste0("~ ", group_var, " - 1"))
  model.null <- if (use_phy) fit_rma_mv_phylo(null_formula, dat = base_data0, V_mat = V0, phy_cor_sub = phy_cor0) else fit_rma_mv(null_formula, dat = base_data0, V_mat = V0)
  Results <- data.frame(formula = group_var, Model = NA, AIC = if (!inherits(model.null, "try-error")) AIC(model.null) else NA, logLik = if (!inherits(model.null, "try-error")) as.numeric(logLik(model.null)) else NA, LRT = NA, df = NA, pval = NA, k = if (!inherits(model.null, "try-error")) summary(model.null)$k else NA, I2 = if (!inherits(model.null, "try-error")) safe_i2(model.null) else NA, stringsAsFactors = FALSE)
  
  for (i in seq_along(mods_to_test)) {
    mod <- mods_to_test[i]; cat("Processing:", group_var, "|", mod, "\n")
    Tem_data <- Whole_data %>% filter(!is.na(.data[[group_var]]), !is.na(.data[[mod]])) %>% droplevels()
    if (mod %in% factor_mods) { Tem_data[[mod]] <- factor(Tem_data[[mod]]); invalid_mod <- nlevels(Tem_data[[mod]]) < 2 } else { Tem_data[[mod]] <- as.numeric(Tem_data[[mod]]); invalid_mod <- all(is.na(Tem_data[[mod]])) || length(unique(Tem_data[[mod]])) < 2 }
    if (invalid_mod) {
      Results <- rbind(Results, data.frame(formula = paste0(group_var, " + ", mod, " - 1"), Model = "1", AIC = NA, logLik = NA, LRT = NA, df = NA, pval = NA, k = nrow(Tem_data), I2 = NA, stringsAsFactors = FALSE), data.frame(formula = paste0(group_var, " + ", mod, " + ", group_var, ":", mod, " - 1"), Model = paste0("1.", i, "A"), AIC = NA, logLik = NA, LRT = NA, df = NA, pval = NA, k = nrow(Tem_data), I2 = NA, stringsAsFactors = FALSE))
      next
    }
    Tem_data[[group_var]] <- factor(Tem_data[[group_var]]); Tem_data$S_ID <- factor(Tem_data$S_ID)
    if (use_phy) { phy_sub <- align_phylo_subset(Tem_data, phy_cor_full); Tem_data <- phy_sub$dat; phy_cor_sub <- phy_sub$phy_cor }
    tmp <- build_shared_V(Tem_data, round_digits = round_digits); Tem_data2 <- tmp$dat; V_mod <- tmp$V
    if (use_phy) { phy_sub2 <- align_phylo_subset(Tem_data2, phy_cor_sub); Tem_data2 <- phy_sub2$dat; phy_cor_sub <- phy_sub2$phy_cor }
    formula.A <- as.formula(paste0("~ ", group_var, " + ", mod, " - 1")); formula.B <- as.formula(paste0("~ ", group_var, " * ", mod, " - 1")); formula_A_text <- paste0(group_var, " + ", mod, " - 1"); formula_B_text <- paste0(group_var, " + ", mod, " + ", group_var, ":", mod, " - 1")
    model.null.mod <- if (use_phy) fit_rma_mv_phylo(as.formula(paste0("~ ", group_var, " - 1")), dat = Tem_data2, V_mat = V_mod, phy_cor_sub = phy_cor_sub) else fit_rma_mv(as.formula(paste0("~ ", group_var, " - 1")), dat = Tem_data2, V_mat = V_mod)
    model.A <- if (use_phy) fit_rma_mv_phylo(formula.A, dat = Tem_data2, V_mat = V_mod, phy_cor_sub = phy_cor_sub) else fit_rma_mv(formula.A, dat = Tem_data2, V_mat = V_mod)
    model.B <- if (use_phy) fit_rma_mv_phylo(formula.B, dat = Tem_data2, V_mat = V_mod, phy_cor_sub = phy_cor_sub) else fit_rma_mv(formula.B, dat = Tem_data2, V_mat = V_mod)
    I2.A <- if (!inherits(model.A, "try-error")) safe_i2(model.A) else NA; I2.B <- if (!inherits(model.B, "try-error")) safe_i2(model.B) else NA
    A_res <- extract_anova_row(try(anova(model.A, model.null.mod), silent = TRUE), formula_A_text, "1", if (!inherits(model.A, "try-error")) summary(model.A)$k else nrow(Tem_data2), I2.A)
    BA_res <- if (!inherits(model.A, "try-error") && !inherits(model.B, "try-error") && length(model.A$beta) != length(model.B$beta)) extract_anova_row(try(anova(model.B, model.A), silent = TRUE), formula_B_text, paste0("1.", i, "A"), summary(model.B)$k, I2.B) else data.frame(formula = formula_B_text, Model = paste0("1.", i, "A"), AIC = NA, logLik = NA, LRT = NA, df = NA, pval = NA, k = if (!inherits(model.B, "try-error")) summary(model.B)$k else nrow(Tem_data2), I2 = I2.B, stringsAsFactors = FALSE)
    Results <- rbind(Results, A_res, BA_res)
  }
  
  Results <- Results[, c("formula", "Model", "AIC", "logLik", "LRT", "df", "pval", "k", "I2")]
  write.csv(Results, output_file, row.names = FALSE)
  Results
}

#### 7. Stable 2: publication-bias assessment ####
egger_from_model <- function(model, predictor_label, group_type) {
  if (inherits(model, "try-error") || is.null(model)) return(data.frame(Group_type = group_type, Predictor_variables = predictor_label, N = NA, test_value = NA, P_value = NA, stringsAsFactors = FALSE))
  rr <- try(residuals(model), silent = TRUE); vv <- try(vcov(model, type = "resid"), silent = TRUE)
  if (inherits(rr, "try-error") || inherits(vv, "try-error")) return(data.frame(Group_type = group_type, Predictor_variables = predictor_label, N = summary(model)$k, test_value = NA, P_value = NA, stringsAsFactors = FALSE))
  s_vec <- try(sign(rr) * sqrt(diag(vv)), silent = TRUE)
  if (inherits(s_vec, "try-error")) return(data.frame(Group_type = group_type, Predictor_variables = predictor_label, N = summary(model)$k, test_value = NA, P_value = NA, stringsAsFactors = FALSE))
  egger_dt <- data.frame(r = rr, s = s_vec); egger_dt <- egger_dt[is.finite(egger_dt$r) & is.finite(egger_dt$s), , drop = FALSE]
  if (nrow(egger_dt) < 3 || length(unique(egger_dt$s)) < 2) return(data.frame(Group_type = group_type, Predictor_variables = predictor_label, N = nrow(egger_dt), test_value = NA, P_value = NA, stringsAsFactors = FALSE))
  egg_test <- try(lm(r ~ s, data = egger_dt), silent = TRUE)
  if (inherits(egg_test, "try-error")) return(data.frame(Group_type = group_type, Predictor_variables = predictor_label, N = nrow(egger_dt), test_value = NA, P_value = NA, stringsAsFactors = FALSE))
  sm <- summary(egg_test); data.frame(Group_type = group_type, Predictor_variables = predictor_label, N = nrow(egg_test$model), test_value = coef(sm)[1, 3], P_value = coef(sm)[1, 4], stringsAsFactors = FALSE)
}

run_pb_assessment <- function(data_all, group_var, mods_to_test, factor_mods, round_digits, phy_cor_full = NULL) {
  use_phy <- !is.null(phy_cor_full); res <- NULL
  dat0 <- data_all %>% filter(!is.na(.data[[group_var]])) %>% droplevels(); dat0[[group_var]] <- factor(dat0[[group_var]])
  if (use_phy) { phy0 <- align_phylo_subset(dat0, phy_cor_full); dat0 <- phy0$dat; phy0_cor <- phy0$phy_cor }
  tmp0 <- build_shared_V(dat0, round_digits = round_digits); dat_base <- tmp0$dat; V_base <- tmp0$V
  if (use_phy) { phy0b <- align_phylo_subset(dat_base, phy0_cor); dat_base <- phy0b$dat; phy0_cor <- phy0b$phy_cor }
  formula_null <- as.formula(paste0("~ ", group_var, " + SqrtNum - 1"))
  model.null <- if (use_phy) fit_rma_mv_phylo(formula_null, dat = dat_base, V_mat = V_base, phy_cor_sub = phy0_cor) else fit_rma_mv(formula_null, dat = dat_base, V_mat = V_base)
  res <- rbind(res, egger_from_model(model.null, predictor_label = group_var, group_type = group_var))
  
  for (i in seq_along(mods_to_test)) {
    mod <- mods_to_test[i]; cat("Publication bias:", group_var, "|", mod, "\n")
    Tem_data <- data_all %>% filter(!is.na(.data[[group_var]]), !is.na(.data[[mod]])) %>% droplevels(); Tem_data[[group_var]] <- factor(Tem_data[[group_var]])
    if (mod %in% factor_mods) { Tem_data[[mod]] <- factor(Tem_data[[mod]]); invalid_mod <- nlevels(Tem_data[[mod]]) < 2 } else { Tem_data[[mod]] <- as.numeric(Tem_data[[mod]]); invalid_mod <- all(is.na(Tem_data[[mod]])) || length(unique(Tem_data[[mod]])) < 2 }
    if (invalid_mod) {
      res <- rbind(res, data.frame(Group_type = group_var, Predictor_variables = paste0(group_var, " + ", mod, " + SqrtNum - 1"), N = nrow(Tem_data), test_value = NA, P_value = NA, stringsAsFactors = FALSE), data.frame(Group_type = group_var, Predictor_variables = paste0(group_var, " + ", mod, " + ", group_var, ":", mod, " + SqrtNum - 1"), N = nrow(Tem_data), test_value = NA, P_value = NA, stringsAsFactors = FALSE))
      next
    }
    if (use_phy) { phy_sub <- align_phylo_subset(Tem_data, phy_cor_full); Tem_data <- phy_sub$dat; phy_cor_sub <- phy_sub$phy_cor }
    tmp <- build_shared_V(Tem_data, round_digits = round_digits); dat_mod <- tmp$dat; V_mod <- tmp$V
    if (use_phy) { phy_sub2 <- align_phylo_subset(dat_mod, phy_cor_sub); dat_mod <- phy_sub2$dat; phy_cor_sub <- phy_sub2$phy_cor }
    formula.A <- as.formula(paste0("~ ", group_var, " + ", mod, " + SqrtNum - 1")); formula.B <- as.formula(paste0("~ ", group_var, " * ", mod, " + SqrtNum - 1"))
    model.A <- if (use_phy) fit_rma_mv_phylo(formula.A, dat = dat_mod, V_mat = V_mod, phy_cor_sub = phy_cor_sub) else fit_rma_mv(formula.A, dat = dat_mod, V_mat = V_mod)
    model.B <- if (use_phy) fit_rma_mv_phylo(formula.B, dat = dat_mod, V_mat = V_mod, phy_cor_sub = phy_cor_sub) else fit_rma_mv(formula.B, dat = dat_mod, V_mat = V_mod)
    res <- rbind(res, egger_from_model(model.A, predictor_label = paste0(group_var, " + ", mod, " + SqrtNum - 1"), group_type = group_var), egger_from_model(model.B, predictor_label = paste0(group_var, " + ", mod, " + ", group_var, ":", mod, " + SqrtNum - 1"), group_type = group_var))
  }
  rownames(res) <- NULL; res
}

#### 8. STABLE 1.1: moderator tests without phylogeny ####
master_data_s1 <- prepare_master_data(read_master_data(cfg), cfg, include_group_var = "both", include_sqrtnum = FALSE, include_phylogeny = FALSE)
res_s11_tg <- run_stable1_core(master_data_s1, group_var = "Trophic_group", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$round_digits_s1, output_file = paste0(cfg$path, "SupplementaryTable_Trophic_group_moderator_tests.csv"))
res_s11_tgr <- run_stable1_core(master_data_s1, group_var = "Trophic_group_response", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$round_digits_s1, output_file = paste0(cfg$path, "SupplementaryTable_Trophic_group_response_moderator_tests.csv"))
print(res_s11_tg); print(res_s11_tgr)

#### 9. STABLE 1.2: moderator tests with phylogeny ####
master_data_s12 <- prepare_master_data(read_master_data(cfg), cfg, include_group_var = "both", include_sqrtnum = FALSE, include_phylogeny = TRUE)
phy_obj <- build_or_load_phylogeny(master_data_s12, cfg)
master_data_s12 <- attach_phylogeny_labels(master_data_s12, phy_obj$tax_lookup, phy_obj$Plant_Cor, cfg)
res_s12_tg <- run_stable1_core(master_data_s12, group_var = "Trophic_group", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$phy_round_digits, output_file = paste0(cfg$path, "SupplementaryTable_Trophic_group_phylogeny_moderator_tests2.csv"), phy_cor_full = phy_obj$Plant_Cor)
res_s12_tgr <- run_stable1_core(master_data_s12, group_var = "Trophic_group_response", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$phy_round_digits, output_file = paste0(cfg$path, "SupplementaryTable_Trophic_group_response_phylogeny_moderator_tests2.csv"), phy_cor_full = phy_obj$Plant_Cor)
print(res_s12_tg); print(res_s12_tgr)

#### 10. STABLE 2.1: publication-bias assessment without phylogeny ####
master_data_s2 <- prepare_master_data(read_master_data(cfg), cfg, include_group_var = "both", include_sqrtnum = TRUE, include_phylogeny = FALSE)
res_s21_tg <- run_pb_assessment(master_data_s2, group_var = "Trophic_group", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$round_digits_s2)
res_s21_tgr <- run_pb_assessment(master_data_s2, group_var = "Trophic_group_response", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$round_digits_s2)
res_s21_all <- rbind(res_s21_tg, res_s21_tgr); res_s21_all <- res_s21_all[!duplicated(res_s21_all[, c("Group_type", "Predictor_variables")]), ]
write.csv(res_s21_tg, paste0(cfg$path, "SupplementTable_publication_bias_Trophic_group.csv"), row.names = FALSE)
write.csv(res_s21_tgr, paste0(cfg$path, "SupplementTable_publication_bias_Trophic_group_response.csv"), row.names = FALSE)
write.csv(res_s21_all, paste0(cfg$path, "SupplementTable_publication_bias_combined.csv"), row.names = FALSE)
print(res_s21_tg); print(res_s21_tgr); print(res_s21_all)

#### 11. STABLE 2.2: publication-bias assessment with phylogeny ####
master_data_s22 <- prepare_master_data(read_master_data(cfg), cfg, include_group_var = "both", include_sqrtnum = TRUE, include_phylogeny = TRUE)
phy_obj2 <- build_or_load_phylogeny(master_data_s22, cfg)
master_data_s22 <- attach_phylogeny_labels(master_data_s22, phy_obj2$tax_lookup, phy_obj2$Plant_Cor, cfg)
res_s22_tg <- run_pb_assessment(master_data_s22, group_var = "Trophic_group", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$phy_round_digits, phy_cor_full = phy_obj2$Plant_Cor)
res_s22_tgr <- run_pb_assessment(master_data_s22, group_var = "Trophic_group_response", mods_to_test = cfg$mods_to_test, factor_mods = cfg$factor_mods, round_digits = cfg$phy_round_digits, phy_cor_full = phy_obj2$Plant_Cor)
res_s22_all <- rbind(res_s22_tg, res_s22_tgr); res_s22_all <- res_s22_all[!duplicated(res_s22_all[, c("Group_type", "Predictor_variables")]), ]
write.csv(res_s22_tg, paste0(cfg$path, "SupplementTable_publication_bias_Trophic_group_phylogeny.csv"), row.names = FALSE)
write.csv(res_s22_tgr, paste0(cfg$path, "SupplementTable_publication_bias_Trophic_group_response_phylogeny.csv"), row.names = FALSE)
write.csv(res_s22_all, paste0(cfg$path, "SupplementTable_publication_bias_combined_phylogeny.csv"), row.names = FALSE)
print(res_s22_tg); print(res_s22_tgr); print(res_s22_all)

#### 12. Final message ####
cat("Finished: STABLE 1.1, 1.2, 2.1, and 2.2 completed.\n")