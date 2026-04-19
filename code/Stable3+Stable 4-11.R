rm(list = ls())

# =========================================================
# Script: 01_prepare_effect_sizes.R
# Purpose: Read raw Excel sheets, calculate effect sizes,
#          apply Geary filtering, and export cleaned datasets.
# =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(metafor)
  library(nlme)
  library(corpcor)
  library(Matrix)
})

# =========================================================
# 1. Configuration
# =========================================================
cfg <- list(
  project_dir = "/Users/yusha",
  input_excel = "data of all trophic groups_2026.xlsx",
  n_max = 1200L,
  round_digits_m = 4L,
  round_digits_sd = 4L,
  round_digits_r = 0L,
  geary_threshold = 3
)

input_file <- file.path(cfg$project_dir, cfg$input_excel)
output_dir <- cfg$project_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Helper functions
# =========================================================
geary_eq13 <- function(mean, sd, n) {
  ifelse(is.na(mean) | is.na(sd) | is.na(n) | sd == 0, NA_real_, (mean / sd) * (4 * n^(3 / 2) / (1 + 4 * n)))
}

infer_trophic_group <- function(dataset_name) {
  if (grepl("^Predator_parasitoid_", dataset_name)) return("Predator_parasitoid")
  if (grepl("^Predator_", dataset_name)) return("Predator")
  if (grepl("^Parasitoid_", dataset_name)) return("Parasitoid")
  if (grepl("^Herbivore_", dataset_name)) return("Herbivore")
  if (grepl("^Crop_", dataset_name)) return("Crop")
  return(NA_character_)
}

safe_write_table <- function(x, file) {
  write.table(x, file = file, sep = "\t", row.names = FALSE, quote = TRUE)
}

read_sheet_subset <- function(file, sheet, n_max, n_cols = 28L) {
  read_xlsx(file, sheet = sheet, n_max = n_max)[, seq_len(n_cols)]
}

# =========================================================
# 3. Sheet registry
# =========================================================
sheet_registry <- tibble::tribble(
  ~sheet_id, ~dataset_name,
  1L,  "Predator_reproduction_response",
  2L,  "Predator_diversity_response",
  3L,  "Predator_predation_response",
  4L,  "Parasitoid_growth_response",
  5L,  "Parasitoid_reproduction_response",
  6L,  "Parasitoid_diversity_response",
  7L,  "Parasitoid_parasitism_response",
  8L,  "Predator_parasitoid_reproduction_response",
  9L,  "Predator_parasitoid_diversity_response",
  10L, "Herbivore_growth_response",
  11L, "Herbivore_reproduction_response",
  12L, "Herbivore_damage_response",
  13L, "Crop_growth_response",
  14L, "Crop_reproduction_response",
  15L, "Crop_quality_response"
)

# =========================================================
# 4. Read all sheets into a named list
# =========================================================
raw_data_list <- lapply(sheet_registry$sheet_id, function(i) read_sheet_subset(input_file, sheet = i, n_max = cfg$n_max, n_cols = 38L))
names(raw_data_list) <- sheet_registry$dataset_name

colnames_abbr <- abbreviate(colnames(raw_data_list[["Crop_growth_response"]]), minlength = 4)

# =========================================================
# 5. Initialize merged outputs
# =========================================================
all_data_geary <- NULL
herbivore_data_geary <- NULL
crop_data_geary <- NULL

# =========================================================
# 6. Main loop
# =========================================================
for (dataset_name in names(raw_data_list)) {
  message("Processing: ", dataset_name)
  
  dat <- raw_data_list[[dataset_name]]
  colnames(dat) <- make.names(colnames_abbr)
  
  dat <- subset(
    dat,
    !is.na(CONTROL_M) & !is.na(TREATMENT_M) &
      !is.na(CONTROL_R) & !is.na(TREATMENT_R) &
      CONTROL_M > 0 & TREATMENT_M > 0 &
      CONTROL_R > 0 & TREATMENT_R > 0
  )
  
  if (nrow(dat) == 0) {
    safe_write_table(data.frame(), file.path(output_dir, paste0(dataset_name, ".txt")))
    next
  }
  
  es_data <- dat %>%
    mutate(
      cv_Control = na_if(CONTROL_SD / CONTROL_M, Inf),
      cv_Treatment = na_if(TREATMENT_SD / TREATMENT_M, Inf)
    )
  
  es_data <- cv_avg(x = TREATMENT_M, sd = TREATMENT_SD, n = TREATMENT_R, group = S_ID, label = "1", data = es_data)
  es_data <- cv_avg(x = CONTROL_M, sd = CONTROL_SD, n = CONTROL_R, group = S_ID, label = "2", data = es_data)
  
  es_data <- es_data %>%
    mutate(
      cv2_trea_new = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2),
      cv2_cont_new = if_else(is.na(cv_Control), b_CV2_2, cv_Control^2),
      yi = lnrr_laj(m1 = TREATMENT_M, m2 = CONTROL_M, cv1_2 = b_CV2_1, cv2_2 = b_CV2_2, n1 = TREATMENT_R, n2 = CONTROL_R),
      vi = v_lnrr_laj(cv1_2 = b_CV2_1, n1 = TREATMENT_R, cv2_2 = b_CV2_2, n2 = CONTROL_R)
    )
  
  idx_t <- which(is.na(es_data$TREATMENT_SD))
  if (length(idx_t) > 0) es_data[idx_t, "TREATMENT_SD"] <- es_data[idx_t, "cv2_trea_new"] * es_data[idx_t, "TREATMENT_M"]
  
  idx_c <- which(is.na(es_data$CONTROL_SD))
  if (length(idx_c) > 0) es_data[idx_c, "CONTROL_SD"] <- es_data[idx_c, "cv2_cont_new"] * es_data[idx_c, "CONTROL_M"]
  
  dat <- subset(es_data, !is.na(yi) & !is.na(vi) & is.finite(yi) & is.finite(vi) & vi > 0)
  
  if (nrow(dat) == 0) {
    safe_write_table(data.frame(), file.path(output_dir, paste0(dataset_name, ".txt")))
    next
  }
  
  dat <- dat %>%
    group_by(S_ID) %>%
    mutate(EsID = row_number()) %>%
    ungroup()
  
  dat$added_species_number <- ifelse(!is.na(dat$NCP_N) & dat$NCP_N > 0, log2(dat$NCP_N), NA_real_)
  dat$Trophic_group <- infer_trophic_group(dataset_name)
  dat$Trophic_group_response <- dataset_name
  
  dat$geary_control <- geary_eq13(dat$CONTROL_M, dat$CONTROL_SD, dat$CONTROL_R)
  dat$geary_treatment <- geary_eq13(dat$TREATMENT_M, dat$TREATMENT_SD, dat$TREATMENT_R)
  
  dat$pass_control <- !is.na(dat$geary_control) & dat$geary_control >= cfg$geary_threshold
  dat$pass_treatment <- !is.na(dat$geary_treatment) & dat$geary_treatment >= cfg$geary_threshold
  dat$pass_geary <- dat$pass_control & dat$pass_treatment
  
  dat_geary <- subset(dat, pass_geary)
  
  safe_write_table(dat_geary, file.path(output_dir, paste0(dataset_name, ".txt")))
  
  all_data_geary <- rbind(all_data_geary, dat_geary)
  
  if (identical(dat$Trophic_group[1], "Herbivore")) herbivore_data_geary <- rbind(herbivore_data_geary, dat_geary)
  if (identical(dat$Trophic_group[1], "Crop")) crop_data_geary <- rbind(crop_data_geary, dat_geary)
}

# =========================================================
# 7. Export merged datasets and registries
# =========================================================
safe_write_table(all_data_geary, file.path(output_dir, "Alldata.txt"))
safe_write_table(herbivore_data_geary, file.path(output_dir, "Antadata.txt"))
safe_write_table(crop_data_geary, file.path(output_dir, "Cropdata.txt"))

safe_write_table(data.frame(raw_name = names(colnames_abbr), abbr_name = colnames_abbr), file.path(output_dir, "name_index.txt"))
safe_write_table(data.frame(data_names = names(raw_data_list)), file.path(output_dir, "data_names.txt"))

total_names <- data.frame(
  total_name = c(
    "Total_response_of_predator",
    "Total_response_of_parasitoid",
    "Total_response_of_predator_parasitoid",
    "Total_response_of_herbivore",
    "Total_response_of_crop"
  ),
  stringsAsFactors = FALSE
)

safe_write_table(total_names, file.path(output_dir, "total_names.txt"))

# =========================================================
# 8. Session info
# =========================================================
sessionInfo()



rm(list = ls())

# =========================================================
##### Stable 3####
# Purpose: Fit overall meta-analytic models for total and non-total
# response datasets, construct shared-control variance matrices,
# generate orchard plots, and export summary tables.
# =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(metafor)
  library(orchaRd)
  library(gridExtra)
  library(ggplot2)
  library(Matrix)
  library(corpcor)
  library(openxlsx)
})

# =========================================================
# 1. Configuration
# =========================================================
cfg <- list(
  project_dir = "/Users/yusha",
  total_names_file = "total_names.txt",
  data_names_file = "data_names.txt",
  suppfig_total_file = "SuppFig10_final.png",
  suppfig_data_file = "SuppFig11_final.png",
  supptab_file = "SuppTab3_final.xlsx",
  round_digits_m = 4L,
  round_digits_sd = 4L,
  round_digits_r = 0L,
  page_width_total_cm = 80,
  page_height_total_cm = 30,
  page_width_data_cm = 80,
  page_height_data_cm = 60,
  dpi = 300
)

# =========================================================
# 2. Helper functions
# =========================================================
read_name_vector <- function(file_path) {
  x <- read.table(file = file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
  gsub('^"|"$', "", as.character(x[[1]]))
}

round_for_group <- function(x, digits = 4L) if (is.numeric(x)) round(x, digits = digits) else x

capitalize_words <- function(x) tools::toTitleCase(x)

make_common_key <- function(m, sd, r) interaction(m, sd, r, drop = TRUE, lex.order = TRUE)

safe_i2 <- function(fit) {
  out <- try(i2_ml(fit), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_)
  if (length(out) >= 2L) return(out[2])
  NA_real_
}

safe_read_dataset <- function(file_path) {
  read.table(file = file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = "", fill = TRUE)
}

empty_result_row <- function(category, n_obs = NA, n_studies = NA, effect = NA, t_value = NA, p_value = NA, df = NA, ci_lb = NA, ci_ub = NA, i2 = NA, v_type = NA, min_eigen = NA) {
  data.frame(
    Category = category,
    Number_of_observations = n_obs,
    Number_of_studies = n_studies,
    Effect_size = effect,
    t_value = t_value,
    P_value = p_value,
    df = df,
    CI_lb = ci_lb,
    CI_ub = ci_ub,
    I2 = i2,
    V_type = v_type,
    Min_eigen_before_fix = min_eigen,
    stringsAsFactors = FALSE
  )
}

save_plot_pages <- function(plot_list, file_name, ncol = 4, nrow = 2, width_cm = 80, height_cm = 30, dpi = 300) {
  if (length(plot_list) == 0L) return(invisible(NULL))
  plots_per_page <- ncol * nrow
  n_pages <- ceiling(length(plot_list) / plots_per_page)
  page_grobs <- vector("list", n_pages)
  
  for (pg in seq_len(n_pages)) {
    start_idx <- (pg - 1L) * plots_per_page + 1L
    end_idx <- min(pg * plots_per_page, length(plot_list))
    page_plots <- plot_list[start_idx:end_idx]
    layout_mat <- matrix(seq_len(nrow * ncol), nrow = nrow, ncol = ncol, byrow = TRUE)
    layout_mat[layout_mat > length(page_plots)] <- NA
    page_grobs[[pg]] <- arrangeGrob(grobs = page_plots, layout_matrix = layout_mat, top = NULL)
  }
  
  ggsave(filename = file_name, plot = marrangeGrob(grobs = page_grobs, nrow = 1, ncol = 1, top = NULL), dpi = dpi, width = width_cm, height = height_cm, units = "cm")
}

calc_shared_block <- function(x) {
  cpart <- x$CONTROL_SD[1]^2 / (x$CONTROL_R[1] * x$CONTROL_M[1]^2)
  v <- matrix(cpart, nrow = nrow(x), ncol = nrow(x))
  diag(v) <- x$vi
  v
}

build_common_ids <- function(dat, round_digits_m, round_digits_sd, round_digits_r) {
  dat %>%
    mutate(
      CONTROL_M_grp = round_for_group(CONTROL_M, round_digits_m),
      CONTROL_SD_grp = round_for_group(CONTROL_SD, round_digits_sd),
      CONTROL_R_grp = round_for_group(CONTROL_R, round_digits_r)
    ) %>%
    group_by(S_ID) %>%
    mutate(Common_ID = match(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp), unique(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp)))) %>%
    ungroup() %>%
    mutate(Common_ID_global = match(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE), unique(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE))))
}

can_use_shared_v <- function(dat) {
  check_common <- dat %>%
    group_by(Common_ID_global) %>%
    summarise(
      n_CONTROL_M_grp = n_distinct(CONTROL_M_grp),
      n_CONTROL_SD_grp = n_distinct(CONTROL_SD_grp),
      n_CONTROL_R_grp = n_distinct(CONTROL_R_grp),
      .groups = "drop"
    )
  all(check_common$n_CONTROL_M_grp == 1L) && all(check_common$n_CONTROL_SD_grp == 1L) && all(check_common$n_CONTROL_R_grp == 1L)
}

fit_rma_mv_safe <- function(dat, V_input, method = "ML") {
  fit <- try(rma.mv(yi, V = V_input, random = list(~1 | S_ID, ~1 | ID), data = dat, test = "t", control = list(optimizer = "BFGS"), method = method), silent = TRUE)
  if (!inherits(fit, "try-error")) return(fit)
  fit <- try(rma.mv(yi, V = V_input, random = list(~1 | S_ID, ~1 | ID), data = dat, test = "t", control = list(optimizer = "nlminb"), method = method), silent = TRUE)
  fit
}

extract_model_summary <- function(fit, category, n_studies, v_type, min_eigen_before_fix) {
  sm <- coef(summary(fit))
  t_value_out <- if ("tval" %in% colnames(sm)) sm$tval[1] else if ("zval" %in% colnames(sm)) sm$zval[1] else NA
  p_value_out <- if ("pval" %in% colnames(sm)) sm$pval[1] else NA
  df_out <- if ("df" %in% colnames(sm)) sm$df[1] else NA
  
  empty_result_row(
    category = category,
    n_obs = fit$k,
    n_studies = n_studies,
    effect = sm$estimate[1],
    t_value = t_value_out,
    p_value = p_value_out,
    df = df_out,
    ci_lb = sm$ci.lb[1],
    ci_ub = sm$ci.ub[1],
    i2 = safe_i2(fit),
    v_type = v_type,
    min_eigen = min_eigen_before_fix
  )
}

make_orchard_plot <- function(fit, dataset_name, tag_letter, is_total = TRUE) {
  x_label_text <- if (is_total) capitalize_words(gsub("Total_response_of_", "", gsub("_", " ", dataset_name))) else capitalize_words(gsub("_response", "", gsub("_", " ", dataset_name)))
  orchard_plot(fit, group = "S_ID", xlab = "Effect size (lnRR)", transfm = "none", twig.size = 2, trunk.size = 2) +
    theme(
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      plot.tag = element_text(size = 30)
    ) +
    labs(tag = tag_letter) +
    scale_x_discrete("", labels = x_label_text)
}

# =========================================================
# 3. Input files and dataset registry
# =========================================================
total_names_vec <- read_name_vector(file.path(cfg$project_dir, cfg$total_names_file))
data_names_vec <- read_name_vector(file.path(cfg$project_dir, cfg$data_names_file))
all_names <- c(total_names_vec, data_names_vec)

# =========================================================
# 4. Initialize outputs
# =========================================================
results_table <- data.frame(
  Category = character(),
  Number_of_observations = numeric(),
  Number_of_studies = numeric(),
  Effect_size = numeric(),
  t_value = numeric(),
  P_value = numeric(),
  df = numeric(),
  CI_lb = numeric(),
  CI_ub = numeric(),
  I2 = numeric(),
  V_type = character(),
  Min_eigen_before_fix = numeric(),
  stringsAsFactors = FALSE
)

total_plots <- list()
non_total_plots <- list()
total_plot_idx <- 1L
non_total_plot_idx <- 1L

# =========================================================
# 5. Main analysis loop
# =========================================================
for (dataset_name in all_names) {
  file_path <- file.path(cfg$project_dir, paste0(dataset_name, ".txt"))
  message("Processing: ", file_path)
  
  if (!file.exists(file_path)) {
    warning("File does not exist: ", file_path)
    results_table <- rbind(results_table, empty_result_row(category = dataset_name))
    next
  }
  
  dat <- safe_read_dataset(file_path)
  if (!is.data.frame(dat)) stop("Input is not a data.frame: ", file_path)
  
  required_cols <- c("S_ID", "ID", "yi", "vi", "CONTROL_M", "CONTROL_SD", "CONTROL_R")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    warning("Missing required columns in ", dataset_name, ": ", paste(missing_cols, collapse = ", "))
    results_table <- rbind(results_table, empty_result_row(category = dataset_name, n_obs = nrow(dat), n_studies = if ("S_ID" %in% names(dat)) length(unique(dat$S_ID)) else NA, v_type = "missing_required_columns"))
    next
  }
  
  dat <- dat %>%
    filter(
      !is.na(yi), is.finite(yi),
      !is.na(vi), is.finite(vi), vi > 0,
      !is.na(S_ID), !is.na(ID),
      !is.na(CONTROL_M), is.finite(CONTROL_M), CONTROL_M > 0,
      !is.na(CONTROL_SD), is.finite(CONTROL_SD), CONTROL_SD >= 0,
      !is.na(CONTROL_R), is.finite(CONTROL_R), CONTROL_R > 0
    )
  
  if (nrow(dat) == 0L) {
    results_table <- rbind(results_table, empty_result_row(category = dataset_name, n_obs = 0))
    next
  }
  
  dat <- build_common_ids(dat, round_digits_m = cfg$round_digits_m, round_digits_sd = cfg$round_digits_sd, round_digits_r = cfg$round_digits_r)
  
  dat_fit <- NULL
  v_type <- NA_character_
  min_eigen_before_fix <- NA_real_
  
  if (nrow(dat) == 1L) {
    dat_fit <- try(rma(yi, vi, data = dat, test = "z"), silent = TRUE)
    if (!inherits(dat_fit, "try-error")) v_type <- "diagonal_vi_single_obs"
  } else {
    if (can_use_shared_v(dat)) {
      V_try <- try(as.matrix(bdiag(lapply(split(dat, dat$Common_ID_global), calc_shared_block))), silent = TRUE)
      
      if (!inherits(V_try, "try-error")) {
        V_ok <- isSymmetric(V_try) && !any(is.na(V_try)) && !any(is.infinite(V_try))
        if (V_ok) {
          eig_vals <- try(eigen(V_try, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
          if (!inherits(eig_vals, "try-error")) {
            min_eigen_before_fix <- min(eig_vals)
            if (!is.finite(min_eigen_before_fix) || min_eigen_before_fix <= 0) {
              V_try <- as.matrix(make.positive.definite(V_try))
              v_type <- "shared_control_pd_fix"
            } else {
              v_type <- "shared_control"
            }
            dat_fit <- fit_rma_mv_safe(dat, V_try, method = "ML")
          }
        }
      }
    }
    
    if (is.null(dat_fit) || inherits(dat_fit, "try-error")) {
      if (length(unique(dat$S_ID)) == 1L) {
        dat_fit <- try(rma(yi, vi, data = dat, test = "t"), silent = TRUE)
        if (!inherits(dat_fit, "try-error")) v_type <- "diagonal_vi_rma"
      } else {
        dat_fit <- fit_rma_mv_safe(dat, dat$vi, method = "ML")
        if (!inherits(dat_fit, "try-error")) v_type <- "diagonal_vi"
      }
    }
  }
  
  if (inherits(dat_fit, "try-error") || is.null(dat_fit)) {
    warning("Model fitting failed for: ", dataset_name)
    results_table <- rbind(results_table, empty_result_row(category = dataset_name, n_obs = nrow(dat), n_studies = length(unique(dat$S_ID)), v_type = v_type, min_eigen = min_eigen_before_fix))
    next
  }
  
  results_table <- rbind(results_table, extract_model_summary(dat_fit, category = dataset_name, n_studies = length(unique(dat$S_ID)), v_type = v_type, min_eigen_before_fix = min_eigen_before_fix))
  
  if (dataset_name %in% total_names_vec) {
    total_plots[[total_plot_idx]] <- make_orchard_plot(dat_fit, dataset_name, tag_letter = letters[total_plot_idx], is_total = TRUE)
    total_plot_idx <- total_plot_idx + 1L
  } else {
    non_total_plots[[non_total_plot_idx]] <- make_orchard_plot(dat_fit, dataset_name, tag_letter = letters[non_total_plot_idx], is_total = FALSE)
    non_total_plot_idx <- non_total_plot_idx + 1L
  }
}

# =========================================================
# 6. Export figures
# =========================================================
if (length(total_plots) > 0L) {
  save_plot_pages(
    plot_list = total_plots,
    file_name = file.path(cfg$project_dir, cfg$suppfig_total_file),
    ncol = 4,
    nrow = 2,
    width_cm = cfg$page_width_total_cm,
    height_cm = cfg$page_height_total_cm,
    dpi = cfg$dpi
  )
}

if (length(non_total_plots) > 0L) {
  save_plot_pages(
    plot_list = non_total_plots,
    file_name = file.path(cfg$project_dir, cfg$suppfig_data_file),
    ncol = 4,
    nrow = 4,
    width_cm = cfg$page_width_data_cm,
    height_cm = cfg$page_height_data_cm,
    dpi = cfg$dpi
  )
}

# =========================================================
# 7. Export summary table
# =========================================================
results_table[] <- lapply(results_table, function(x) if (is.list(x)) as.character(x) else x)
names(results_table) <- make.unique(trimws(names(results_table)))
names(results_table)[names(results_table) == ""] <- paste0("V", which(names(results_table) == ""))

openxlsx::write.xlsx(results_table, file = file.path(cfg$project_dir, cfg$supptab_file), rowNames = FALSE)

cat("Finished Stable 3.\n")
sessionInfo()





rm(list = ls())

# =========================================================
##### Stable 4-11####
# Purpose: Run subgroup meta-analytic models across predefined
# moderator levels, generate orchard plots, and export summary tables.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
  library(orchaRd)
  library(gridExtra)
  library(ggplot2)
  library(Matrix)
  library(corpcor)
  library(openxlsx)
})

# =========================================================
# 1. Configuration
# =========================================================
cfg <- list(
  project_dir = "/Users/yusha",
  total_names_file = "total_names.txt",
  data_names_file = "data_names.txt",
  output_root = "Stable4_11_compare_outputs_final",
  round_digits_m = 4L,
  round_digits_sd = 4L,
  round_digits_r = 0L,
  per_plot_width_cm = 20,
  per_plot_height_cm = 15,
  dpi = 300,
  max_per_page = 8L,
  max_ncol = 4L
)

output_root <- file.path(cfg$project_dir, cfg$output_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Helper functions
# =========================================================
read_name_vector <- function(file_path) {
  x <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = "", fill = TRUE)
  gsub('^"|"$', "", as.character(x[[1]]))
}

round_for_group <- function(x, digits = 4L) if (is.numeric(x)) round(x, digits = digits) else x

make_common_key <- function(m, sd, r) interaction(m, sd, r, drop = TRUE, lex.order = TRUE)

capitalize_words <- function(x) tools::toTitleCase(x)

safe_i2 <- function(fit) {
  out <- try(i2_ml(fit), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_)
  if (length(out) >= 2L) return(out[2])
  NA_real_
}

standardize_zone <- function(x) {
  x <- trimws(as.character(x))
  x_low <- tolower(x)
  x_low[x_low %in% c("tropical", "tropics", "tropic")] <- "tropic"
  x_low[x_low %in% c("temperate")] <- "temperate"
  x_out <- x
  x_out[x_low == "tropic"] <- "Tropic"
  x_out[x_low == "temperate"] <- "Temperate"
  x_out[!(x_low %in% c("tropic", "temperate"))] <- NA_character_
  x_out
}

standardize_ex3 <- function(x) {
  x <- trimws(as.character(x))
  x_low <- tolower(x)
  x_low[x_low %in% c("plot")] <- "plot"
  x_low[x_low %in% c("pot")] <- "pot"
  x_out <- x
  x_out[x_low == "plot"] <- "Plot"
  x_out[x_low == "pot"] <- "Pot"
  x_out[!(x_low %in% c("plot", "pot"))] <- NA_character_
  x_out
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)

safe_read_dataset <- function(file_path, need_zone_standardize = FALSE, need_ex3_standardize = FALSE) {
  dat <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = "", fill = TRUE)
  if (!is.data.frame(dat)) return(NULL)
  names(dat) <- gsub('^"|"$', "", names(dat))
  if (need_zone_standardize && "Zone" %in% names(dat)) dat$Zone <- standardize_zone(dat$Zone)
  if (need_ex3_standardize && "Ex_3" %in% names(dat)) dat$Ex_3 <- standardize_ex3(dat$Ex_3)
  dat
}

save_plots_dynamic <- function(plot_list, out_dir, prefix_name, per_plot_width = 20, per_plot_height = 15, dpi = 300, max_per_page = 8L, max_ncol = 4L) {
  if (length(plot_list) == 0L) return(invisible(NULL))
  label_index <- 1L
  
  for (i in seq(1L, length(plot_list), by = max_per_page)) {
    plots_this_page <- plot_list[i:min(i + max_per_page - 1L, length(plot_list))]
    n_plots <- length(plots_this_page)
    
    for (j in seq_along(plots_this_page)) {
      plots_this_page[[j]] <- plots_this_page[[j]] + labs(tag = letters[(label_index - 1L) %% 26L + 1L])
      label_index <- label_index + 1L
    }
    
    ncol_layout <- min(max_ncol, n_plots)
    nrow_layout <- ceiling(n_plots / ncol_layout)
    layout_mat <- matrix(seq_len(nrow_layout * ncol_layout), nrow = nrow_layout, ncol = ncol_layout, byrow = TRUE)
    layout_mat[layout_mat > n_plots] <- NA
    
    g <- arrangeGrob(grobs = plots_this_page, layout_matrix = layout_mat, top = NULL)
    
    ggsave(
      filename = file.path(out_dir, paste0(prefix_name, "_page", ceiling(i / max_per_page), ".png")),
      plot = g,
      dpi = dpi,
      width = per_plot_width * ncol_layout,
      height = per_plot_height * nrow_layout,
      units = "cm",
      limitsize = FALSE
    )
  }
}

build_common_ids <- function(dat, round_digits_m, round_digits_sd, round_digits_r) {
  dat %>%
    mutate(
      CONTROL_M_grp = round_for_group(CONTROL_M, round_digits_m),
      CONTROL_SD_grp = round_for_group(CONTROL_SD, round_digits_sd),
      CONTROL_R_grp = round_for_group(CONTROL_R, round_digits_r)
    ) %>%
    group_by(S_ID) %>%
    mutate(Common_ID = match(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp), unique(make_common_key(CONTROL_M_grp, CONTROL_SD_grp, CONTROL_R_grp)))) %>%
    ungroup() %>%
    mutate(Common_ID_global = match(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE), unique(interaction(S_ID, Common_ID, drop = TRUE, lex.order = TRUE))))
}

build_shared_v <- function(dat) {
  calc_v_block <- function(x) {
    cpart <- x$CONTROL_SD[1]^2 / (x$CONTROL_R[1] * x$CONTROL_M[1]^2)
    v <- matrix(cpart, nrow = nrow(x), ncol = nrow(x))
    diag(v) <- x$vi
    v
  }
  
  out <- list(V = NULL, V_type = NA_character_, min_eigen_before_fix = NA_real_, ok = FALSE, msg = NULL)
  
  V_try <- try(as.matrix(bdiag(lapply(split(dat, dat$Common_ID_global), calc_v_block))), silent = TRUE)
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
    out$msg <- "shared-control V fixed by make.positive.definite"
  } else {
    out$V <- V_try
    out$V_type <- "shared_control"
    out$ok <- TRUE
  }
  
  out
}

fit_rma_mv_safe <- function(V_input, dat, method = "ML") {
  random_terms <- list(~1 | S_ID, ~1 | ID)
  for (opt in c("BFGS", "nlminb", "optim")) {
    fit <- try(rma.mv(yi, V = V_input, random = random_terms, data = dat, test = "t", method = method, control = list(optimizer = opt)), silent = TRUE)
    if (!inherits(fit, "try-error")) return(fit)
  }
  "try-error"
}

empty_subset_result <- function(subgroup_level, category_name, v_type = NA, min_eigen = NA, n_obs = NA, n_studies = NA) {
  data.frame(
    subgroup = subgroup_level,
    Category = category_name,
    Num_obs = n_obs,
    Num_studies = n_studies,
    Effect_size = NA,
    t_value = NA,
    P_value = NA,
    df = NA,
    CI_lb = NA,
    CI_ub = NA,
    I2 = NA,
    V_type = v_type,
    Min_eigen_before_fix = min_eigen,
    stringsAsFactors = FALSE
  )
}

extract_subset_summary <- function(fit, subgroup_level, category_name, v_type, min_eigen_before_fix, dat) {
  sm <- coef(summary(fit))
  est <- if ("estimate" %in% colnames(sm)) sm[1, "estimate"] else NA
  tval <- if ("tval" %in% colnames(sm)) sm[1, "tval"] else if ("zval" %in% colnames(sm)) sm[1, "zval"] else NA
  pval <- if ("pval" %in% colnames(sm)) sm[1, "pval"] else NA
  cilb <- if ("ci.lb" %in% colnames(sm)) sm[1, "ci.lb"] else NA
  ciub <- if ("ci.ub" %in% colnames(sm)) sm[1, "ci.ub"] else NA
  df_val <- if ("df" %in% colnames(sm)) sm[1, "df"] else nrow(dat) - 1L
  
  data.frame(
    subgroup = subgroup_level,
    Category = category_name,
    Num_obs = fit$k,
    Num_studies = length(unique(dat$S_ID)),
    Effect_size = est,
    t_value = tval,
    P_value = pval,
    df = df_val,
    CI_lb = cilb,
    CI_ub = ciub,
    I2 = safe_i2(fit),
    V_type = v_type,
    Min_eigen_before_fix = min_eigen_before_fix,
    stringsAsFactors = FALSE
  )
}

run_one_subset <- function(dat0, category_name, subgroup_var, subgroup_level, round_digits_m, round_digits_sd, round_digits_r) {
  if (is.null(dat0) || !is.data.frame(dat0) || nrow(dat0) == 0L) return(list(result = empty_subset_result(subgroup_level, category_name, n_obs = 0), fit = NULL, note = "empty_input"))
  if (!(subgroup_var %in% names(dat0))) return(list(result = empty_subset_result(subgroup_level, category_name, v_type = "missing_subgroup_column"), fit = NULL, note = "missing_subgroup_column"))
  
  dat <- subset(dat0, dat0[[subgroup_var]] == subgroup_level)
  
  required_cols <- c("S_ID", "ID", "yi", "vi", "CONTROL_M", "CONTROL_SD", "CONTROL_R")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) return(list(result = empty_subset_result(subgroup_level, category_name, v_type = "missing_required_columns", n_obs = nrow(dat), n_studies = if ("S_ID" %in% names(dat)) length(unique(dat$S_ID)) else NA), fit = NULL, note = "missing_required_columns"))
  
  dat <- dat %>%
    filter(
      !is.na(yi), is.finite(yi),
      !is.na(vi), is.finite(vi), vi > 0,
      !is.na(S_ID), trimws(as.character(S_ID)) != "",
      !is.na(ID), trimws(as.character(ID)) != "",
      !is.na(CONTROL_M), is.finite(CONTROL_M), CONTROL_M > 0,
      !is.na(CONTROL_SD), is.finite(CONTROL_SD), CONTROL_SD >= 0,
      !is.na(CONTROL_R), is.finite(CONTROL_R), CONTROL_R > 0
    )
  
  if (nrow(dat) == 0L) return(list(result = empty_subset_result(subgroup_level, category_name, n_obs = 0, n_studies = 0), fit = NULL, note = "filtered_to_zero"))
  
  dat <- build_common_ids(dat, round_digits_m = round_digits_m, round_digits_sd = round_digits_sd, round_digits_r = round_digits_r)
  
  dat_fit <- NULL
  v_type <- NA_character_
  min_eigen_before_fix <- NA_real_
  note <- NULL
  
  if (nrow(dat) == 1L) {
    dat_fit <- try(rma(yi, vi, data = dat, test = "z"), silent = TRUE)
    if (!inherits(dat_fit, "try-error")) v_type <- "diagonal_vi_single_obs"
  } else {
    check_common <- dat %>%
      group_by(Common_ID_global) %>%
      summarise(
        n_CONTROL_M_grp = n_distinct(CONTROL_M_grp),
        n_CONTROL_SD_grp = n_distinct(CONTROL_SD_grp),
        n_CONTROL_R_grp = n_distinct(CONTROL_R_grp),
        .groups = "drop"
      )
    
    use_V <- all(check_common$n_CONTROL_M_grp == 1L) && all(check_common$n_CONTROL_SD_grp == 1L) && all(check_common$n_CONTROL_R_grp == 1L)
    
    if (use_V) {
      V_info <- build_shared_v(dat)
      min_eigen_before_fix <- V_info$min_eigen_before_fix
      note <- V_info$msg
      
      if (isTRUE(V_info$ok)) {
        dat_fit <- fit_rma_mv_safe(V_info$V, dat, method = "ML")
        if (inherits(dat_fit, "try-error")) dat_fit <- fit_rma_mv_safe(V_info$V, dat, method = "REML")
        if (!inherits(dat_fit, "try-error")) v_type <- V_info$V_type
      }
    }
    
    if (is.null(dat_fit) || inherits(dat_fit, "try-error")) {
      if (length(unique(dat$S_ID)) == 1L) {
        dat_fit <- try(rma(yi, vi, data = dat, test = "t"), silent = TRUE)
        if (!inherits(dat_fit, "try-error")) v_type <- "diagonal_vi_rma"
      } else {
        dat_fit <- fit_rma_mv_safe(dat$vi, dat, method = "ML")
        if (inherits(dat_fit, "try-error")) dat_fit <- fit_rma_mv_safe(dat$vi, dat, method = "REML")
        if (!inherits(dat_fit, "try-error")) v_type <- "diagonal_vi"
      }
    }
  }
  
  if (inherits(dat_fit, "try-error") || is.null(dat_fit)) return(list(result = empty_subset_result(subgroup_level, category_name, v_type = v_type, min_eigen = min_eigen_before_fix, n_obs = nrow(dat), n_studies = length(unique(dat$S_ID))), fit = NULL, note = note))
  
  list(result = extract_subset_summary(dat_fit, subgroup_level, category_name, v_type, min_eigen_before_fix, dat), fit = dat_fit, note = note)
}

make_orchard_plot <- function(fit_obj, subgroup_label, dataset_name, is_total = TRUE) {
  x_label_text <- if (is_total) capitalize_words(gsub("Total_response_of_", "", gsub("_", " ", dataset_name))) else capitalize_words(gsub("_response", "", gsub("_", " ", dataset_name)))
  orchard_plot(fit_obj, group = "S_ID", xlab = "Effect size (lnRR)", transfm = "none", twig.size = 2, trunk.size = 2) +
    theme(
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      plot.tag = element_text(size = 30),
      plot.title = element_text(size = 18, hjust = 0.5)
    ) +
    ggtitle(subgroup_label) +
    scale_x_discrete("", labels = x_label_text)
}

run_stable_analysis <- function(stable_name, subgroup_var, subgroup_levels, start_fig_num, first_col_name, output_csv_name, datasets, total_names, output_root, project_dir, round_digits_m, round_digits_sd, round_digits_r, need_zone_standardize = FALSE, need_ex3_standardize = FALSE, per_plot_width = 20, per_plot_height = 15, dpi = 300, max_per_page = 8L, max_ncol = 4L) {
  message("========================================")
  message("Running ", stable_name)
  message("Subgroup variable: ", subgroup_var)
  message("========================================")
  
  stable_dir <- file.path(output_root, stable_name)
  figure_dir <- file.path(stable_dir, "figures")
  table_dir <- file.path(stable_dir, "tables")
  dir.create(stable_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  fig_index <- start_fig_num
  non_pd_log <- character()
  
  res <- data.frame(
    subgroup = character(),
    Category = character(),
    Num_obs = numeric(),
    Num_studies = numeric(),
    Effect_size = numeric(),
    t_value = numeric(),
    P_value = numeric(),
    df = numeric(),
    CI_lb = numeric(),
    CI_ub = numeric(),
    I2 = numeric(),
    V_type = character(),
    Min_eigen_before_fix = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (lev in subgroup_levels) {
    message("------------------------------")
    message("Processing ", subgroup_var, " = ", lev)
    message("------------------------------")
    
    total_plot_list <- list()
    data_plot_list <- list()
    total_idx <- 1L
    data_idx <- 1L
    
    for (dataset_name in datasets) {
      file_path <- file.path(project_dir, paste0(dataset_name, ".txt"))
      message("  Dataset: ", dataset_name)
      
      if (!file.exists(file_path)) {
        res <- rbind(res, data.frame(subgroup = lev, Category = dataset_name, Num_obs = NA, Num_studies = NA, Effect_size = NA, t_value = NA, P_value = NA, df = NA, CI_lb = NA, CI_ub = NA, I2 = NA, V_type = "file_not_found", Min_eigen_before_fix = NA, stringsAsFactors = FALSE))
        next
      }
      
      dat0 <- safe_read_dataset(file_path, need_zone_standardize = need_zone_standardize, need_ex3_standardize = need_ex3_standardize)
      
      fit_out <- run_one_subset(
        dat0 = dat0,
        category_name = dataset_name,
        subgroup_var = subgroup_var,
        subgroup_level = lev,
        round_digits_m = round_digits_m,
        round_digits_sd = round_digits_sd,
        round_digits_r = round_digits_r
      )
      
      res <- rbind(res, fit_out$result)
      
      if (!is.null(fit_out$note) && grepl("fixed|invalid|failed", fit_out$note, ignore.case = TRUE)) non_pd_log <- c(non_pd_log, paste0(lev, " | ", dataset_name, " : ", fit_out$note))
      
      if (!is.null(fit_out$fit)) {
        plot_try <- try(make_orchard_plot(fit_out$fit, subgroup_label = lev, dataset_name = dataset_name, is_total = dataset_name %in% total_names), silent = TRUE)
        
        if (!inherits(plot_try, "try-error")) {
          if (dataset_name %in% total_names) {
            total_plot_list[[total_idx]] <- plot_try
            total_idx <- total_idx + 1L
          } else {
            data_plot_list[[data_idx]] <- plot_try
            data_idx <- data_idx + 1L
          }
        } else {
          message("Plot failed for: ", lev, " | ", dataset_name)
        }
      }
    }
    
    if (length(total_plot_list) > 0L) {
      save_plots_dynamic(
        plot_list = total_plot_list,
        out_dir = figure_dir,
        prefix_name = paste0("SuppFig", fig_index, "_", sanitize_name(lev), "_total"),
        per_plot_width = per_plot_width,
        per_plot_height = per_plot_height,
        dpi = dpi,
        max_per_page = max_per_page,
        max_ncol = max_ncol
      )
    }
    fig_index <- fig_index + 1L
    
    if (length(data_plot_list) > 0L) {
      save_plots_dynamic(
        plot_list = data_plot_list,
        out_dir = figure_dir,
        prefix_name = paste0("SuppFig", fig_index, "_", sanitize_name(lev), "_data"),
        per_plot_width = per_plot_width,
        per_plot_height = per_plot_height,
        dpi = dpi,
        max_per_page = max_per_page,
        max_ncol = max_ncol
      )
    }
    fig_index <- fig_index + 1L
  }
  
  colnames(res) <- c(first_col_name, "Category", "Number of observations", "Number of studies", "Effect size", "t-value", "P-value", "df", "Lower Bound", "Upper Bound", "I2", "V_type", "Min_eigen_before_fix")
  write.csv(res, file.path(table_dir, output_csv_name), row.names = FALSE)
  
  if (length(non_pd_log) > 0L) writeLines(unique(non_pd_log), con = file.path(table_dir, paste0(stable_name, "_non_pd_log.txt")))
  
  res
}

# =========================================================
# 3. Dataset registry
# =========================================================
total_names <- read_name_vector(file.path(cfg$project_dir, cfg$total_names_file))
data_names <- read_name_vector(file.path(cfg$project_dir, cfg$data_names_file))
datasets <- c(total_names, data_names)

# =========================================================
# 4. Stable 4-11 analyses
# =========================================================
res_stable4 <- run_stable_analysis(
  stable_name = "Stable4",
  subgroup_var = "CP_C3",
  subgroup_levels = c("Herbaceous_plant", "Woody_plant"),
  start_fig_num = 11,
  first_col_name = "Crop form",
  output_csv_name = "SuppTab4.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

res_stable5 <- run_stable_analysis(
  stable_name = "Stable5",
  subgroup_var = "CP_Ct",
  subgroup_levels = c("Food_crop", "Cash_crop"),
  start_fig_num = 17,
  first_col_name = "Crop type",
  output_csv_name = "SuppTab5.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

res_stable6 <- run_stable_analysis(
  stable_name = "Stable6",
  subgroup_var = "Mn__",
  subgroup_levels = c("Generalist", "Specialist"),
  start_fig_num = 21,
  first_col_name = "Diet breadth",
  output_csv_name = "SuppTab6.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

res_stable7 <- run_stable_analysis(
  stable_name = "Stable7",
  subgroup_var = "NE_c",
  subgroup_levels = c("Generalist", "Specialist"),
  start_fig_num = 25,
  first_col_name = "Natural enemy host range",
  output_csv_name = "SuppTab7.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

res_stable8 <- run_stable_analysis(
  stable_name = "Stable8",
  subgroup_var = "Dv_S",
  subgroup_levels = c("Intercropping", "Cover_cropping", "Sown_field_margins"),
  start_fig_num = 29,
  first_col_name = "Diversification strategy",
  output_csv_name = "SuppTab8.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

res_stable9 <- run_stable_analysis(
  stable_name = "Stable9",
  subgroup_var = "Ex_4",
  subgroup_levels = c("Natural", "Semi_natural", "Controlled"),
  start_fig_num = 35,
  first_col_name = "Experimental context",
  output_csv_name = "SuppTab9.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

res_stable10 <- run_stable_analysis(
  stable_name = "Stable10",
  subgroup_var = "Ex_3",
  subgroup_levels = c("Plot", "Pot"),
  start_fig_num = 41,
  first_col_name = "Experimental setting",
  output_csv_name = "SuppTab10.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  need_ex3_standardize = TRUE,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

res_stable11 <- run_stable_analysis(
  stable_name = "Stable11",
  subgroup_var = "Zone",
  subgroup_levels = c("Temperate", "Tropic"),
  start_fig_num = 45,
  first_col_name = "Climatic zone",
  output_csv_name = "SuppTab11.csv",
  datasets = datasets,
  total_names = total_names,
  output_root = output_root,
  project_dir = cfg$project_dir,
  round_digits_m = cfg$round_digits_m,
  round_digits_sd = cfg$round_digits_sd,
  round_digits_r = cfg$round_digits_r,
  need_zone_standardize = TRUE,
  per_plot_width = cfg$per_plot_width_cm,
  per_plot_height = cfg$per_plot_height_cm,
  dpi = cfg$dpi,
  max_per_page = cfg$max_per_page,
  max_ncol = cfg$max_ncol
)

cat("========================================\n")
cat("All Stable 4-11 analyses finished.\n")
cat("Outputs saved in:\n")
cat(output_root, "\n")
cat("========================================\n")

sessionInfo()