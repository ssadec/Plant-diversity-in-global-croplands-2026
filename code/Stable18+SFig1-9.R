##### Stable 18####
# Purpose: Fit weighted regression or mixed-effects models
# for added_species_number across trophic-group datasets,
# export summary tables, generate individual plots, and
# assemble combined multi-panel figures.
# =========================================================

suppressPackageStartupMessages({
  library(nlme)
  library(effects)
  library(scales)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(rsq)
  library(gridExtra)
})

# =========================================================
# 1. Configuration
# =========================================================
cfg <- list(
  project_dir = "/Users/yusha",
  total_names_file = "total_names.txt",
  all_data_file = "Alldata.txt",
  output_dir = "Stable18_output_final",
  plots_dir = "plots",
  combined_dir = "combined_plots",
  result_csv = "Stable18.csv",
  figure_prefix = "Supplementary Data Fig.",
  x_limits = c(0, 6),
  y_limits = c(-8, 8),
  individual_plot_width = 6.5,
  individual_plot_height = 5.2,
  individual_plot_dpi = 300,
  combined_png_width_cm = 8,
  combined_png_height_cm = 9,
  combined_png_res = 300
)

out_dir <- file.path(cfg$project_dir, cfg$output_dir)
plot_dir <- file.path(out_dir, cfg$plots_dir)
combined_dir <- file.path(out_dir, cfg$combined_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. Dataset registry and labels
# =========================================================
total_names_df <- read_tsv(file.path(cfg$project_dir, cfg$total_names_file), col_names = TRUE, show_col_types = FALSE)
total_names <- if ("total_name" %in% names(total_names_df)) total_names_df$total_name else total_names_df[[1]]
total_names <- gsub('^"|"$', "", total_names)[1:5]

epdt <- read.table(file.path(cfg$project_dir, cfg$all_data_file), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = "", fill = TRUE)
if (!("S_ID" %in% names(epdt))) epdt$S_ID <- paste0("Athr", seq_len(nrow(epdt)))

category_labels <- c("Predator performance", "Parasitoid performance", "Predator-parasitoid performance", "Herbivore performance", "Crop performance")
category_short <- c("Predator", "Parasitoid", "Predator-parasitoid", "Herbivore", "Crop")
category_colors <- c("#BC3A24", "#EFB882", "#06798F", "#7A9A01", "#384D73")
names(category_colors) <- category_short

strata <- c(
  "Global",
  "Herbaceous_plant", "Woody_plant",
  "Food_crop", "Cash_crop",
  "Mn__Generalist", "Mn__Specialist",
  "NE_c_Generalist", "NE_c_Specialist",
  "Intercropping", "Cover_cropping", "Sown_field_margins",
  "Natural", "Semi_natural", "Controlled",
  "Plot", "Pot",
  "Temperate", "Tropic"
)

fig_ids <- paste(rep(seq_along(strata), each = length(total_names)), rep(letters[seq_along(total_names)], length(strata)), sep = "")
fig_sources <- paste(cfg$figure_prefix, fig_ids)

combined_plot_groups <- list(
  Global = c("Global"),
  CP_C3 = c("Herbaceous_plant", "Woody_plant"),
  CP_Ct = c("Food_crop", "Cash_crop"),
  Mn__ = c("Mn__Generalist", "Mn__Specialist"),
  NE_c = c("NE_c_Generalist", "NE_c_Specialist"),
  Dv_S = c("Intercropping", "Cover_cropping", "Sown_field_margins"),
  Ex_4 = c("Natural", "Semi_natural", "Controlled"),
  Ex_3 = c("Plot", "Pot"),
  Zone = c("Temperate", "Tropic")
)

# =========================================================
# 3. Helper functions
# =========================================================
standardize_ex3 <- function(x) {
  x <- trimws(as.character(x)); x_low <- tolower(x); x_low[x_low %in% "plot"] <- "plot"; x_low[x_low %in% "pot"] <- "pot"
  x_out <- x; x_out[x_low == "plot"] <- "Plot"; x_out[x_low == "pot"] <- "Pot"; x_out[!(x_low %in% c("plot", "pot"))] <- NA
  x_out
}

standardize_zone <- function(x) {
  x <- trimws(as.character(x)); x_low <- tolower(x); x_low[x_low %in% c("tropical", "tropics", "tropic")] <- "tropic"; x_low[x_low %in% "temperate"] <- "temperate"
  x_out <- x; x_out[x_low == "tropic"] <- "Tropic"; x_out[x_low == "temperate"] <- "Temperate"; x_out[!(x_low %in% c("tropic", "temperate"))] <- NA
  x_out
}

safe_point_size <- function(vi, default_size = 3, min_size = 1.5, max_size = 6.5) {
  vi_num <- suppressWarnings(as.numeric(vi)); wi <- 1 / sqrt(vi_num); ok <- is.finite(wi) & !is.na(wi) & wi > 0
  if (sum(ok) == 0) return(rep(default_size, length(vi)))
  if (sum(ok) == 1) { out <- rep(default_size, length(vi)); out[!ok] <- default_size; return(out) }
  wi_ok <- wi[ok]; wmin <- min(wi_ok, na.rm = TRUE); wmax <- max(wi_ok, na.rm = TRUE)
  if (!is.finite(wmin) || !is.finite(wmax) || wmax == wmin) { out <- rep(default_size, length(vi)); out[!ok] <- default_size; return(out) }
  out <- rep(default_size, length(vi)); out[ok] <- min_size + (max_size - min_size) * (wi_ok - wmin) / (wmax - wmin); out[!ok] <- default_size
  out
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)

make_eq <- function(intercept, slope) paste0("Y = ", round(intercept, 4), ifelse(slope >= 0, " + ", " - "), abs(round(slope, 4)), "X")

format_p_label <- function(p) {
  if (is.null(p) || length(p) != 1 || is.na(p) || !is.finite(p)) return(NA_character_)
  if (p < 0.001) "P < 0.001" else paste0("P = ", sprintf("%.3f", p))
}

format_r2_label <- function(r2) {
  if (is.null(r2) || length(r2) != 1 || is.na(r2) || !is.finite(r2)) return(NA_character_)
  paste0("R² = ", sprintf("%.4f", r2))
}

empty_fit_result <- function() list(model_type = NA, fit = NULL, intercept = NA, slope = NA, lower = NA, upper = NA, se = NA, df = NA, t_val = NA, p_val = NA, R2 = NA, regression_eq = NA)

prepare_data <- function(dat, epdt) {
  colnames(dat) <- gsub("\\.", "_", trimws(colnames(dat))); names(dat) <- gsub('^"|"$', "", names(dat))
  dat$Epdt <- if ("S_ID" %in% names(dat) && "S_ID" %in% names(epdt)) epdt[match(dat$S_ID, epdt$S_ID), 2] else NA
  if ("yi" %in% names(dat)) dat$yi <- suppressWarnings(as.numeric(dat$yi))
  if ("vi" %in% names(dat)) dat$vi <- suppressWarnings(as.numeric(dat$vi))
  dat$added_species_number <- if ("added_species_number" %in% names(dat)) suppressWarnings(as.numeric(dat$added_species_number)) else NA
  dat$Ex_3_std <- if ("Ex_3" %in% names(dat)) standardize_ex3(dat$Ex_3) else NA
  dat$Zone_std <- if ("Zone" %in% names(dat)) standardize_zone(dat$Zone) else NA
  
  dat %>% filter(!is.na(yi), !is.na(vi), !is.na(added_species_number), is.finite(yi), is.finite(vi), is.finite(added_species_number), vi > 0)
}

subset_by_strat <- function(dat, strat_name) {
  if (strat_name == "Global") return(dat)
  if (strat_name %in% c("Herbaceous_plant", "Woody_plant")) return(if ("CP_C3" %in% names(dat)) subset(dat, CP_C3 == strat_name) else dat[0, , drop = FALSE])
  if (strat_name %in% c("Food_crop", "Cash_crop")) return(if ("CP_Ct" %in% names(dat)) subset(dat, CP_Ct == strat_name) else dat[0, , drop = FALSE])
  if (strat_name %in% c("Mn__Generalist", "Mn__Specialist")) return(if ("Mn__" %in% names(dat)) subset(dat, Mn__ == sub("Mn__", "", strat_name)) else dat[0, , drop = FALSE])
  if (strat_name %in% c("NE_c_Generalist", "NE_c_Specialist")) return(if ("NE_c" %in% names(dat)) subset(dat, NE_c == sub("NE_c_", "", strat_name)) else dat[0, , drop = FALSE])
  if (strat_name %in% c("Intercropping", "Cover_cropping", "Sown_field_margins")) return(if ("Dv_S" %in% names(dat)) subset(dat, Dv_S == strat_name) else dat[0, , drop = FALSE])
  if (strat_name %in% c("Natural", "Semi_natural", "Controlled")) return(if ("Ex_4" %in% names(dat)) subset(dat, Ex_4 == strat_name) else dat[0, , drop = FALSE])
  if (strat_name %in% c("Plot", "Pot")) return(if ("Ex_3_std" %in% names(dat)) subset(dat, Ex_3_std == strat_name) else dat[0, , drop = FALSE])
  if (strat_name %in% c("Temperate", "Tropic")) return(if ("Zone_std" %in% names(dat)) subset(dat, Zone_std == strat_name) else dat[0, , drop = FALSE])
  dat[0, , drop = FALSE]
}

fit_one_model <- function(dat_sub) {
  out <- empty_fit_result()
  if (nrow(dat_sub) <= 1) return(out)
  
  x <- dat_sub$added_species_number
  x <- x[is.finite(x)]
  if (length(unique(x)) <= 1) return(out)
  
  n_study <- length(unique(na.omit(dat_sub$S_ID)))
  
  if (n_study <= 1) {
    fit <- try(lm(yi ~ added_species_number, data = dat_sub, weights = 1 / vi), silent = TRUE)
    if (inherits(fit, "try-error")) return(out)
    coefs <- try(coef(summary(fit)), silent = TRUE)
    if (inherits(coefs, "try-error") || is.null(coefs) || !("(Intercept)" %in% rownames(coefs)) || !("added_species_number" %in% rownames(coefs))) return(out)
    
    intercept <- coefs["(Intercept)", "Estimate"]; slope <- coefs["added_species_number", "Estimate"]; se <- coefs["added_species_number", "Std. Error"]
    df_val <- fit$df.residual; t_val <- coefs["added_species_number", "t value"]; p_val <- coefs["added_species_number", "Pr(>|t|)"]
    lower <- slope - qt(0.975, df = df_val) * se; upper <- slope + qt(0.975, df = df_val) * se
    
    out$model_type <- "lm"; out$fit <- fit; out$intercept <- intercept; out$slope <- slope; out$lower <- lower; out$upper <- upper
    out$se <- se; out$df <- df_val; out$t_val <- t_val; out$p_val <- p_val; out$R2 <- tryCatch(summary(fit)$r.squared, error = function(e) NA); out$regression_eq <- make_eq(intercept, slope)
    return(out)
  }
  
  fit <- try(lme(yi ~ added_species_number, random = ~1 | S_ID, data = dat_sub, weights = varFixed(~vi), na.action = na.omit, control = lmeControl(opt = "optim", sigma = 1, returnObject = TRUE)), silent = TRUE)
  if (inherits(fit, "try-error")) return(out)
  coefs <- try(summary(fit)$tTable, silent = TRUE)
  if (inherits(coefs, "try-error") || is.null(coefs) || !("(Intercept)" %in% rownames(coefs)) || !("added_species_number" %in% rownames(coefs))) return(out)
  
  intercept <- coefs["(Intercept)", "Value"]; slope <- coefs["added_species_number", "Value"]; se <- coefs["added_species_number", "Std.Error"]
  df_val <- coefs["added_species_number", "DF"]; t_val <- coefs["added_species_number", "t-value"]; p_val <- coefs["added_species_number", "p-value"]
  lower <- slope - qt(0.975, df = df_val) * se; upper <- slope + qt(0.975, df = df_val) * se
  
  out$model_type <- "lme"; out$fit <- fit; out$intercept <- intercept; out$slope <- slope; out$lower <- lower; out$upper <- upper
  out$se <- se; out$df <- df_val; out$t_val <- t_val; out$p_val <- p_val; out$R2 <- tryCatch(rsq.lmm(fit)$model, error = function(e) NA); out$regression_eq <- make_eq(intercept, slope)
  out
}

get_effect_prediction <- function(fit_obj, xvar = "added_species_number", xseq = NULL, dat_sub = NULL) {
  if (is.null(fit_obj)) return(NULL)
  if (is.null(xseq)) {
    if (is.null(dat_sub)) return(NULL)
    xseq <- seq(min(dat_sub[[xvar]], na.rm = TRUE), max(dat_sub[[xvar]], na.rm = TRUE), length.out = 200)
  }
  pred <- try(Effect(xvar, fit_obj, xlevels = setNames(list(xseq), xvar)), silent = TRUE)
  if (inherits(pred, "try-error")) return(NULL)
  data.frame(x = pred$x[, 1], fit = pred$fit[, 1], lower = pred$lower[, 1], upper = pred$upper[, 1])
}

save_reg_plot <- function(dat_sub, fit_res, fig_label, category_name, strat_name, out_file, point_color, obs_n, study_n) {
  if (obs_n <= 0 || study_n <= 0) return(invisible(NULL))
  
  show_reg <- !is.na(fit_res$p_val) && is.finite(fit_res$p_val) && fit_res$p_val < 0.05 && !is.na(fit_res$regression_eq)
  dat_sub$point_size_raw <- 1 / dat_sub$vi
  
  p <- ggplot() +
    theme_bw(base_size = 12) +
    labs(title = paste0(fig_label, " | ", category_name, " | ", strat_name), x = expression(log[2] ~ Number ~ of ~ added ~ plant ~ species), y = "Effect size (yi)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank()) +
    annotate("text", x = Inf, y = Inf, label = paste0("n = ", obs_n, "\nStudy = ", study_n), hjust = 1.05, vjust = 1.2, size = 3.8) +
    geom_point(data = dat_sub, aes(x = added_species_number, y = yi, size = point_size_raw), shape = 21, fill = point_color, colour = "black", alpha = 0.8) +
    scale_size_continuous(guide = "none")
  
  if (show_reg) {
    pred_df <- get_effect_prediction(fit_res$fit, dat_sub = dat_sub)
    if (!is.null(pred_df) && nrow(pred_df) > 1) {
      p <- p +
        geom_ribbon(data = pred_df, aes(x = x, ymin = lower, ymax = upper), inherit.aes = FALSE, fill = "grey70", alpha = 0.5) +
        geom_line(data = pred_df, aes(x = x, y = fit), inherit.aes = FALSE, linewidth = 0.9, colour = "black") +
        annotate("text", x = Inf, y = Inf, label = paste(fit_res$regression_eq, format_p_label(fit_res$p_val), format_r2_label(fit_res$R2), sep = "\n"), hjust = 1.05, vjust = 4.5, size = 3.8)
    }
  } else {
    reason_txt <- if (nrow(dat_sub) <= 1) "No estimable regression" else if (is.na(fit_res$regression_eq)) "" else NULL
    if (!is.null(reason_txt)) p <- p + annotate("text", x = Inf, y = Inf, label = reason_txt, hjust = 1.05, vjust = 3.0, size = 3.8)
  }
  
  ggsave(filename = out_file, plot = p, width = cfg$individual_plot_width, height = cfg$individual_plot_height, dpi = cfg$individual_plot_dpi)
}

build_combined_plot_objects <- function(strat_values, total_names, category_labels, res_table, project_dir, epdt, category_short) {
  all_plots <- list()
  plot_labels <- letters
  label_index <- 1L
  
  for (strat_name in strat_values) {
    for (i in seq_along(category_labels)) {
      dat_path <- file.path(project_dir, paste0(total_names[i], ".txt"))
      if (!file.exists(dat_path)) next
      
      dat <- read.table(dat_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = "", fill = TRUE)
      dat <- prepare_data(dat, epdt)
      dat_sub <- subset_by_strat(dat, strat_name)
      
      obs <- nrow(dat_sub)
      study_n <- if (nrow(dat_sub) > 0) length(unique(na.omit(dat_sub$S_ID))) else 0
      if (obs <= 0 || study_n <= 0) next
      
      dat_sub$size <- safe_point_size(dat_sub$vi)
      dat_sub$Category_short <- factor(category_labels[i], levels = category_labels, labels = category_short)
      
      res_row <- subset(res_table, Strat == strat_name & Category == category_labels[i])
      eq <- if (nrow(res_row) == 1) as.character(res_row$Regression_equation) else NA
      pval <- if (nrow(res_row) == 1) suppressWarnings(as.numeric(as.character(res_row$P_value))) else NA
      r2 <- if (nrow(res_row) == 1) suppressWarnings(as.numeric(as.character(res_row$R2))) else NA
      
      all_plots[[length(all_plots) + 1L]] <- list(data = dat_sub, eq = eq, pval = pval, R2 = r2, obs_label = paste0("(", obs, "/", study_n, ")"), label = plot_labels[label_index], ylab = category_labels[i], strat = strat_name)
      label_index <- label_index + 1L
    }
  }
  
  all_plots
}

draw_combined_plot <- function(all_plots, outfile, category_colors) {
  if (length(all_plots) == 0L) return(NULL)
  
  ncol_plot <- 5L
  nrow_plot <- ceiling(length(all_plots) / ncol_plot)
  
  png(outfile, units = "cm", width = cfg$combined_png_width_cm * ncol_plot, height = cfg$combined_png_height_cm * nrow_plot, res = cfg$combined_png_res)
  par(mai = c(1, 0.8, 0.3, 0.3), mfrow = c(nrow_plot, ncol_plot), mgp = c(3, 0.6, 0), family = "Arial")
  
  for (i in seq_along(all_plots)) {
    p <- all_plots[[i]]
    dat_sub <- p$data
    
    plot(dat_sub$added_species_number, dat_sub$yi, type = "n", xlab = "", ylab = p$ylab, ylim = cfg$y_limits, xlim = cfg$x_limits, cex.lab = 1.5, font.lab = 1, cex.axis = 1.2)
    mtext(expression(atop("Log"[2] * " Number of added", "plant species")), side = 1, line = 5, cex = 1, family = "Arial")
    
    show_reg <- !is.na(p$pval) && is.finite(p$pval) && p$pval < 0.05 && !is.na(p$eq) && trimws(p$eq) != ""
    
    if (show_reg) {
      fit_res <- tryCatch(fit_one_model(dat_sub), error = function(e) empty_fit_result())
      if (!is.null(fit_res$fit) && !is.na(fit_res$p_val) && is.finite(fit_res$p_val) && fit_res$p_val < 0.05 && !is.na(fit_res$regression_eq)) {
        pred_df <- get_effect_prediction(fit_res$fit, dat_sub = dat_sub)
        if (!is.null(pred_df) && nrow(pred_df) > 1) {
          polygon(c(rev(pred_df$x), pred_df$x), c(rev(pred_df$upper), pred_df$lower), col = alpha("grey", 0.5), border = NA)
          lines(pred_df$x, pred_df$fit, col = "black", lwd = 3)
          text(5.5, 6.8, labels = fit_res$regression_eq, cex = 1.0, adj = 1)
          text(5.5, 6.0, labels = format_p_label(fit_res$p_val), cex = 1.0, adj = 1)
          r2_lab <- format_r2_label(fit_res$R2)
          if (!is.na(r2_lab)) text(5.5, 5.2, labels = r2_lab, cex = 1.0, adj = 1)
        }
      }
    } else {
      reason_txt <- if (nrow(dat_sub) <= 1) "" else if (is.na(p$eq) || trimws(p$eq) == "" || tolower(trimws(p$eq)) == "na") "" else NULL
      if (!is.null(reason_txt)) text(5.5, 6.0, labels = reason_txt, cex = 1.0, adj = 1)
    }
    
    points(dat_sub$added_species_number, dat_sub$yi, bg = category_colors[as.character(dat_sub$Category_short)], col = "black", pch = 21, cex = dat_sub$size)
    text(0.4, 7, labels = p$label, cex = 1.5, font = 2)
    text(5.5, 4.4, labels = p$obs_label, cex = 1.0, adj = 1)
  }
  
  dev.off()
}

# =========================================================
# 4. Main loop: fit models and export individual plots
# =========================================================
res <- NULL
k <- 1L

for (j in strata) {
  for (i in seq_along(total_names)) {
    dat_path <- file.path(cfg$project_dir, paste0(total_names[i], ".txt"))
    
    if (!file.exists(dat_path)) {
      warning("File does not exist: ", dat_path)
      res <- rbind(res, data.frame(Figsource = fig_sources[k], Dataset = total_names[i], Category = category_labels[i], Strat = j, n = 0, n_study = 0, Regression_equation = NA, Intercept = NA, Slope = NA, Lower = NA, Upper = NA, SE = NA, DF = NA, T_value = NA, P_value = NA, R2 = NA, stringsAsFactors = FALSE))
      k <- k + 1L
      next
    }
    
    dat <- read.table(dat_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", comment.char = "", fill = TRUE)
    dat <- prepare_data(dat, epdt)
    dat_sub <- subset_by_strat(dat, j)
    
    obs_n <- nrow(dat_sub)
    study_n <- length(unique(na.omit(dat_sub$S_ID)))
    
    message("Checking: Strat = ", j, " | Dataset = ", total_names[i], " | Category = ", category_labels[i], " | n = ", obs_n, " | n_study = ", study_n, " | unique_x = ", length(unique(dat_sub$added_species_number)))
    
    fit_res <- tryCatch(fit_one_model(dat_sub), error = function(e) { message("Model failed: Strat = ", j, " | Dataset = ", total_names[i], " | Category = ", category_labels[i], " | n = ", obs_n, " | n_study = ", study_n, " | error = ", conditionMessage(e)); empty_fit_result() })
    
    res <- rbind(res, data.frame(Figsource = fig_sources[k], Dataset = total_names[i], Category = category_labels[i], Strat = j, n = obs_n, n_study = study_n, Regression_equation = fit_res$regression_eq, Intercept = fit_res$intercept, Slope = fit_res$slope, Lower = fit_res$lower, Upper = fit_res$upper, SE = fit_res$se, DF = fit_res$df, T_value = fit_res$t_val, P_value = fit_res$p_val, R2 = fit_res$R2, stringsAsFactors = FALSE))
    
    if (obs_n > 0 && study_n > 0) {
      fig_file <- file.path(plot_dir, paste0("Stable18_", sanitize_name(fig_sources[k]), "_", sanitize_name(total_names[i]), "_", sanitize_name(j), ".png"))
      save_reg_plot(dat_sub = dat_sub, fit_res = fit_res, fig_label = fig_sources[k], category_name = category_labels[i], strat_name = j, out_file = fig_file, point_color = category_colors[category_short[i]], obs_n = obs_n, study_n = study_n)
    }
    
    k <- k + 1L
  }
}

# =========================================================
# 5. Export summary table
# =========================================================
write.csv(res, file.path(out_dir, cfg$result_csv), row.names = FALSE)

# =========================================================
# 6. Export combined plots
# =========================================================
for (grp_name in names(combined_plot_groups)) {
  strat_values <- combined_plot_groups[[grp_name]]
  plot_objs <- build_combined_plot_objects(strat_values = strat_values, total_names = total_names, category_labels = category_labels, res_table = res, project_dir = cfg$project_dir, epdt = epdt, category_short = category_short)
  if (length(plot_objs) == 0L) next
  draw_combined_plot(all_plots = plot_objs, outfile = file.path(combined_dir, paste0("Stable18_", grp_name, "_combined.png")), category_colors = category_colors)
}

# =========================================================
# 7. Completion message
# =========================================================
cat("Done!\n")
cat("Output folder: ", out_dir, "\n")
cat("Result table: ", file.path(out_dir, cfg$result_csv), "\n")
cat("Individual plots folder: ", plot_dir, "\n")
cat("Combined plots folder: ", combined_dir, "\n")

sessionInfo()
