#### SuppTab19a-1 ####
rm(list = ls())

library(nlme)
library(MuMIn)

prefix <- '/Users/yusha/'
paired_names <- read.table(paste0(prefix, "paired_names.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#==================================================
# 1. Safe helpers
#==================================================

safe_AICc <- function(fit) tryCatch(as.numeric(MuMIn::AICc(fit)), error = function(e) NA_real_)
safe_logLik <- function(fit) tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)

safe_r2_glmm <- function(fit) {
  out <- tryCatch(MuMIn::r.squaredGLMM(fit), error = function(e) NULL)
  if (is.null(out)) return(c(R2_marginal = NA_real_, R2_conditional = NA_real_))
  out <- as.data.frame(out)
  if (all(c("R2m", "R2c") %in% colnames(out))) return(c(R2_marginal = as.numeric(out$R2m[1]), R2_conditional = as.numeric(out$R2c[1])))
  if (ncol(out) >= 2) return(c(R2_marginal = as.numeric(out[1, 1]), R2_conditional = as.numeric(out[1, 2])))
  c(R2_marginal = NA_real_, R2_conditional = NA_real_)
}

#==================================================
# 2. paired_names section
#==================================================

res <- NULL

for (i in 1:nrow(paired_names)) {
  dat <- read.table(paste0(prefix, paired_names[i, 1], ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  if (!("ID" %in% names(dat))) dat$ID <- seq_len(nrow(dat))
  
  dat.fit1 <- lme(yi1 ~ 1, dat, random = ~1 | S_ID/ID, weights = varFixed(~vi1), control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2 ~ 1, dat, random = ~1 | S_ID/ID, weights = varFixed(~vi2), control = lmeControl(sigma = 1))
  dat.fit3 <- lme(yi3 ~ 1, dat, random = ~1 | S_ID/ID, weights = varFixed(~vi3), control = lmeControl(sigma = 1))
  
  dat.residual <- data.frame(r1 = residuals(dat.fit1, type = "pearson"), r2 = residuals(dat.fit2, type = "pearson"), r3 = residuals(dat.fit3, type = "pearson"), S_ID = dat$S_ID, ID = dat$ID, stringsAsFactors = FALSE)
  
  yi1_on_yi2 <- lme(r2 ~ -1 + r1, dat.residual, random = ~ -1 + r1 | S_ID/ID)
  yi2_on_yi3 <- lme(r3 ~ -1 + r2, dat.residual, random = ~ -1 + r2 | S_ID/ID)
  
  r2_1 <- safe_r2_glmm(yi1_on_yi2)
  r2_2 <- safe_r2_glmm(yi2_on_yi3)
  
  res <- rbind(
    res,
    data.frame(
      Item = paired_names[i, 1],
      AIC_1 = AIC(yi1_on_yi2),
      AIC_2 = AIC(yi2_on_yi3),
      AICc_1 = safe_AICc(yi1_on_yi2),
      AICc_2 = safe_AICc(yi2_on_yi3),
      BIC_1 = BIC(yi1_on_yi2),
      BIC_2 = BIC(yi2_on_yi3),
      logLik_1 = safe_logLik(yi1_on_yi2),
      logLik_2 = safe_logLik(yi2_on_yi3),
      R2_marginal_1 = r2_1["R2_marginal"],
      R2_conditional_1 = r2_1["R2_conditional"],
      R2_marginal_2 = r2_2["R2_marginal"],
      R2_conditional_2 = r2_2["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      fix.empty.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

#==================================================
# 3. Strat section
#==================================================

Strat <- c("Global", "Herbaceous_plant", "Food_crop", "Cash_crop", "Mn__-Generalist", "NE_c-Generalist", "Intercropping", "Cover_cropping", "Natural", "Plot", "Temperate", "Tropic")

dat <- read.table(paste0(prefix, "enemies_vs_herbivore_vs_crop.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!("ID" %in% names(dat))) dat$ID <- seq_len(nrow(dat))

for (i in Strat) {
  if (i == "Global") {
    dat_sub <- dat
  } else if (grepl("Mn__-", i)) {
    level <- sub("Mn__-", "", i)
    dat_sub <- subset(dat, Mn__ == level)
  } else if (grepl("NE_c-", i)) {
    level <- sub("NE_c-", "", i)
    dat_sub <- subset(dat, NE_c == level)
  } else {
    dat_sub <- subset(dat, CP_C3 == i | CP_Ct == i | Mn__ == i | NE_c == i | Dv_s == i | Ex_4 == i | Ex_3 == i | Zone == i)
  }
  
  if (nrow(dat_sub) < 2) next
  if (!("ID" %in% names(dat_sub))) dat_sub$ID <- seq_len(nrow(dat_sub))
  
  dat.fit1 <- lme(yi1 ~ 1, dat_sub, random = ~1 | S_ID/ID, weights = varFixed(~vi1), control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2 ~ 1, dat_sub, random = ~1 | S_ID/ID, weights = varFixed(~vi2), control = lmeControl(sigma = 1))
  dat.fit3 <- lme(yi3 ~ 1, dat_sub, random = ~1 | S_ID/ID, weights = varFixed(~vi3), control = lmeControl(sigma = 1))
  
  dat.residual <- data.frame(r1 = residuals(dat.fit1, type = "pearson"), r2 = residuals(dat.fit2, type = "pearson"), r3 = residuals(dat.fit3, type = "pearson"), S_ID = dat_sub$S_ID, ID = dat_sub$ID, stringsAsFactors = FALSE)
  
  yi1_on_yi2 <- lme(r2 ~ -1 + r1, dat.residual, random = ~ -1 + r1 | S_ID/ID)
  yi2_on_yi3 <- lme(r3 ~ -1 + r2, dat.residual, random = ~ -1 + r2 | S_ID/ID)
  
  r2_1 <- safe_r2_glmm(yi1_on_yi2)
  r2_2 <- safe_r2_glmm(yi2_on_yi3)
  
  res <- rbind(
    res,
    data.frame(
      Item = i,
      AIC_1 = AIC(yi1_on_yi2),
      AIC_2 = AIC(yi2_on_yi3),
      AICc_1 = safe_AICc(yi1_on_yi2),
      AICc_2 = safe_AICc(yi2_on_yi3),
      BIC_1 = BIC(yi1_on_yi2),
      BIC_2 = BIC(yi2_on_yi3),
      logLik_1 = safe_logLik(yi1_on_yi2),
      logLik_2 = safe_logLik(yi2_on_yi3),
      R2_marginal_1 = r2_1["R2_marginal"],
      R2_conditional_1 = r2_1["R2_conditional"],
      R2_marginal_2 = r2_2["R2_marginal"],
      R2_conditional_2 = r2_2["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      fix.empty.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

View(res)
write.csv(res, '/Users/yusha/SuppTab19a-1也许派.csv', row.names = FALSE)



#### SuppTab19a-2 ####
rm(list = ls())

library(nlme)
library(MuMIn)

prefix <- "/Users/yusha/"

paired_names <- read.table(paste0(prefix, "paired_names.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

safe_AICc <- function(fit) tryCatch(as.numeric(MuMIn::AICc(fit)), error = function(e) NA_real_)
safe_logLik <- function(fit) tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)

safe_r2 <- function(fit) {
  out <- tryCatch(MuMIn::r.squaredGLMM(fit), error = function(e) NULL)
  if (is.null(out)) return(c(R2_marginal = NA_real_, R2_conditional = NA_real_))
  out <- as.data.frame(out)
  if (all(c("R2m", "R2c") %in% colnames(out))) return(c(R2_marginal = as.numeric(out$R2m[1]), R2_conditional = as.numeric(out$R2c[1])))
  if (ncol(out) >= 2) return(c(R2_marginal = as.numeric(out[1, 1]), R2_conditional = as.numeric(out[1, 2])))
  c(R2_marginal = NA_real_, R2_conditional = NA_real_)
}

fit_lme_safe <- function(formula_obj, data_obj, vi_col) {
  tryCatch(lme(fixed = formula_obj, data = data_obj, random = ~1 | S_ID/ID, weights = varFixed(as.formula(paste0("~", vi_col))), control = lmeControl(sigma = 1)), error = function(e) NULL)
}

fit_resid_safe <- function(formula_obj, data_obj, rand_term) {
  tryCatch(lme(fixed = formula_obj, data = data_obj, random = rand_term), error = function(e) NULL)
}

res <- NULL

for (ii in seq_len(nrow(paired_names))) {
  dat <- read.table(paste0(prefix, paired_names[ii, 1], ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  if (!("ID" %in% names(dat))) dat$ID <- seq_len(nrow(dat))
  
  fit1 <- fit_lme_safe(yi1 ~ 1 + added_species_number, dat, "vi1")
  fit2 <- fit_lme_safe(yi2 ~ 1 + added_species_number, dat, "vi2")
  fit3 <- fit_lme_safe(yi3 ~ 1 + added_species_number, dat, "vi3")
  
  if (is.null(fit1) || is.null(fit2) || is.null(fit3)) next
  
  dat.residual <- data.frame(r1 = residuals(fit1, type = "pearson"), r2 = residuals(fit2, type = "pearson"), r3 = residuals(fit3, type = "pearson"), S_ID = dat$S_ID, ID = dat$ID, stringsAsFactors = FALSE)
  
  yi1_on_yi2 <- fit_resid_safe(r2 ~ -1 + r1, dat.residual, ~ -1 + r1 | S_ID/ID)
  yi2_on_yi3 <- fit_resid_safe(r3 ~ -1 + r2, dat.residual, ~ -1 + r2 | S_ID/ID)
  
  if (is.null(yi1_on_yi2) || is.null(yi2_on_yi3)) next
  
  r2_1 <- safe_r2(yi1_on_yi2)
  r2_2 <- safe_r2(yi2_on_yi3)
  
  res <- rbind(
    res,
    data.frame(
      Item = paired_names[ii, 1],
      AIC_1 = AIC(yi1_on_yi2),
      AICc_1 = safe_AICc(yi1_on_yi2),
      BIC_1 = BIC(yi1_on_yi2),
      logLik_1 = safe_logLik(yi1_on_yi2),
      R2_marginal_1 = r2_1["R2_marginal"],
      R2_conditional_1 = r2_1["R2_conditional"],
      AIC_2 = AIC(yi2_on_yi3),
      AICc_2 = safe_AICc(yi2_on_yi3),
      BIC_2 = BIC(yi2_on_yi3),
      logLik_2 = safe_logLik(yi2_on_yi3),
      R2_marginal_2 = r2_2["R2_marginal"],
      R2_conditional_2 = r2_2["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      stringsAsFactors = FALSE
    )
  )
}

Strat <- c("Global", "Herbaceous_plant", "Food_crop", "Cash_crop", "Mn__-Generalist", "NE_c-Generalist", "Intercropping", "Cover_cropping", "Natural", "Plot", "Temperate", "Tropic")

dat_all <- read.table(paste0(prefix, "enemies_vs_herbivore_vs_crop.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!("ID" %in% names(dat_all))) dat_all$ID <- seq_len(nrow(dat_all))

for (ii in Strat) {
  if (ii == "Global") {
    dat_sub <- dat_all
  } else if (grepl("Mn__-", ii)) {
    lev <- sub("Mn__-", "", ii)
    dat_sub <- subset(dat_all, Mn__ == lev)
  } else if (grepl("NE_c-", ii)) {
    lev <- sub("NE_c-", "", ii)
    dat_sub <- subset(dat_all, NE_c == lev)
  } else {
    dat_sub <- subset(dat_all, CP_C3 == ii | CP_Ct == ii | Mn__ == ii | NE_c == ii | Dv_s == ii | Ex_4 == ii | Ex_3 == ii | Zone == ii)
  }
  
  if (nrow(dat_sub) < 2) next
  if (!("ID" %in% names(dat_sub))) dat_sub$ID <- seq_len(nrow(dat_sub))
  
  fit1 <- fit_lme_safe(yi1 ~ 1 + added_species_number, dat_sub, "vi1")
  fit2 <- fit_lme_safe(yi2 ~ 1 + added_species_number, dat_sub, "vi2")
  fit3 <- fit_lme_safe(yi3 ~ 1 + added_species_number, dat_sub, "vi3")
  
  if (is.null(fit1) || is.null(fit2) || is.null(fit3)) next
  
  dat.residual <- data.frame(r1 = residuals(fit1, type = "pearson"), r2 = residuals(fit2, type = "pearson"), r3 = residuals(fit3, type = "pearson"), S_ID = dat_sub$S_ID, ID = dat_sub$ID, stringsAsFactors = FALSE)
  
  yi1_on_yi2 <- fit_resid_safe(r2 ~ -1 + r1, dat.residual, ~ -1 + r1 | S_ID/ID)
  yi2_on_yi3 <- fit_resid_safe(r3 ~ -1 + r2, dat.residual, ~ -1 + r2 | S_ID/ID)
  
  if (is.null(yi1_on_yi2) || is.null(yi2_on_yi3)) next
  
  r2_1 <- safe_r2(yi1_on_yi2)
  r2_2 <- safe_r2(yi2_on_yi3)
  
  res <- rbind(
    res,
    data.frame(
      Item = ii,
      AIC_1 = AIC(yi1_on_yi2),
      AICc_1 = safe_AICc(yi1_on_yi2),
      BIC_1 = BIC(yi1_on_yi2),
      logLik_1 = safe_logLik(yi1_on_yi2),
      R2_marginal_1 = r2_1["R2_marginal"],
      R2_conditional_1 = r2_1["R2_conditional"],
      AIC_2 = AIC(yi2_on_yi3),
      AICc_2 = safe_AICc(yi2_on_yi3),
      BIC_2 = BIC(yi2_on_yi3),
      logLik_2 = safe_logLik(yi2_on_yi3),
      R2_marginal_2 = r2_2["R2_marginal"],
      R2_conditional_2 = r2_2["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      stringsAsFactors = FALSE
    )
  )
}

View(res)
write.csv(res, "/Users/yusha/SuppTab19a-2也许派.csv", row.names = FALSE)



#### STable 19b-1 ####
rm(list = ls())

library(nlme)
library(MuMIn)

prefix <- "/Users/yusha/"

paired_names <- read.table(paste0(prefix, "paired_names.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

safe_AICc <- function(fit) {
  tryCatch(as.numeric(MuMIn::AICc(fit)), error = function(e) NA_real_)
}

safe_logLik <- function(fit) {
  tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
}

safe_r2 <- function(fit) {
  out <- tryCatch(MuMIn::r.squaredGLMM(fit), error = function(e) NULL)
  
  if (is.null(out)) return(c(R2_marginal = NA_real_, R2_conditional = NA_real_))
  
  out <- as.data.frame(out)
  
  if (all(c("R2m", "R2c") %in% colnames(out))) {
    return(c(R2_marginal = as.numeric(out$R2m[1]), R2_conditional = as.numeric(out$R2c[1])))
  }
  
  if (ncol(out) >= 2) {
    return(c(R2_marginal = as.numeric(out[1, 1]), R2_conditional = as.numeric(out[1, 2])))
  }
  
  c(R2_marginal = NA_real_, R2_conditional = NA_real_)
}

fit_lme_safe <- function(formula_obj, data_obj, vi_col) {
  tryCatch(
    lme(fixed = formula_obj, data = data_obj, random = ~1 | S_ID/ID, weights = varFixed(as.formula(paste0("~", vi_col))), control = lmeControl(sigma = 1)),
    error = function(e) NULL
  )
}

fit_resid_safe <- function(formula_obj, data_obj) {
  tryCatch(
    lme(fixed = formula_obj, data = data_obj, random = ~ -1 + r1 | S_ID/ID),
    error = function(e) NULL
  )
}

res <- NULL

## paired_names part
for (ii in seq_len(nrow(paired_names))) {
  dat <- read.table(paste0(prefix, paired_names[ii, 1], ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  if (!("ID" %in% names(dat))) dat$ID <- seq_len(nrow(dat))
  
  fit1 <- fit_lme_safe(yi1 ~ 1, dat, "vi1")
  fit2 <- fit_lme_safe(yi2 ~ 1, dat, "vi2")
  
  if (is.null(fit1) || is.null(fit2)) next
  
  dat.residual <- data.frame(r1 = residuals(fit1, type = "pearson"), r2 = residuals(fit2, type = "pearson"), vi1 = dat$vi1, vi2 = dat$vi2, S_ID = dat$S_ID, ID = dat$ID, stringsAsFactors = FALSE)
  
  yi1_on_yi2 <- fit_resid_safe(r2 ~ -1 + r1, dat.residual)
  if (is.null(yi1_on_yi2)) next
  
  r2_out <- safe_r2(yi1_on_yi2)
  
  res <- rbind(
    res,
    data.frame(
      Item = paired_names[ii, 1],
      AIC = AIC(yi1_on_yi2),
      AICc = safe_AICc(yi1_on_yi2),
      BIC = BIC(yi1_on_yi2),
      logLik = safe_logLik(yi1_on_yi2),
      R2_marginal = r2_out["R2_marginal"],
      R2_conditional = r2_out["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      stringsAsFactors = FALSE
    )
  )
}

## Strat part
Strat <- c("Global", "Herbaceous_plant", "Woody_plant", "Food_crop", "Cash_crop", "Generalist", "Specialist", "Intercropping", "Cover_cropping", "Sown_field_margins", "Natural", "Plot", "Temperate", "Tropic")

for (ii in Strat) {
  dat <- read.table(paste0(prefix, "herbivore_vs_crop.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  if (!("ID" %in% names(dat))) dat$ID <- seq_len(nrow(dat))
  
  if (ii != "Global") {
    dat <- subset(dat, (CP_C3 == ii) | (CP_Ct == ii) | (Mn__ == ii) | (Dv_s == ii) | (Ex_4 == ii) | (Ex_3 == ii) | (Zone == ii))
  }
  
  fit1 <- fit_lme_safe(yi1 ~ 1, dat, "vi1")
  fit2 <- fit_lme_safe(yi2 ~ 1, dat, "vi2")
  
  if (is.null(fit1) || is.null(fit2)) next
  
  dat.residual <- data.frame(r1 = residuals(fit1, type = "pearson"), r2 = residuals(fit2, type = "pearson"), vi1 = dat$vi1, vi2 = dat$vi2, S_ID = dat$S_ID, ID = dat$ID, stringsAsFactors = FALSE)
  
  yi1_on_yi2 <- fit_resid_safe(r2 ~ -1 + r1, dat.residual)
  if (is.null(yi1_on_yi2)) next
  
  r2_out <- safe_r2(yi1_on_yi2)
  
  res <- rbind(
    res,
    data.frame(
      Item = ii,
      AIC = AIC(yi1_on_yi2),
      AICc = safe_AICc(yi1_on_yi2),
      BIC = BIC(yi1_on_yi2),
      logLik = safe_logLik(yi1_on_yi2),
      R2_marginal = r2_out["R2_marginal"],
      R2_conditional = r2_out["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      stringsAsFactors = FALSE
    )
  )
}

View(res)
write.csv(res, "/Users/yusha/STable19b-1也许派.csv", row.names = FALSE)



#### STable 19b-2 ####
rm(list = ls())

library(nlme)
library(MuMIn)

prefix <- "/Users/yusha/"

paired_names <- read.table(paste0(prefix, "paired_names.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

safe_AICc <- function(fit) {
  tryCatch(as.numeric(AICc(fit)), error = function(e) NA_real_)
}

safe_logLik <- function(fit) {
  tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
}

safe_r2 <- function(fit) {
  out <- tryCatch(MuMIn::r.squaredGLMM(fit), error = function(e) NULL)
  
  if (is.null(out)) return(c(R2_marginal = NA_real_, R2_conditional = NA_real_))
  
  out <- as.data.frame(out)
  
  if (all(c("R2m", "R2c") %in% colnames(out))) {
    return(c(R2_marginal = as.numeric(out$R2m[1]), R2_conditional = as.numeric(out$R2c[1])))
  }
  
  if (ncol(out) >= 2) {
    return(c(R2_marginal = as.numeric(out[1, 1]), R2_conditional = as.numeric(out[1, 2])))
  }
  
  c(R2_marginal = NA_real_, R2_conditional = NA_real_)
}

fix_id_sid <- function(dat) {
  if (!("ID" %in% names(dat))) {
    dat$ID <- seq_len(nrow(dat))
  } else {
    bad_id <- is.na(dat$ID) | trimws(as.character(dat$ID)) == ""
    if (any(bad_id)) dat$ID[bad_id] <- paste0("tmpID_", seq_len(sum(bad_id)))
  }
  
  dat$S_ID <- as.character(dat$S_ID)
  dat$ID   <- as.character(dat$ID)
  
  dat
}

fit_main_safe <- function(dat) {
  dat_use <- dat[
    !is.na(dat$S_ID) & trimws(dat$S_ID) != "" &
      !is.na(dat$ID) & trimws(dat$ID) != "" &
      !is.na(dat$yi1) & is.finite(dat$yi1) &
      !is.na(dat$yi2) & is.finite(dat$yi2) &
      !is.na(dat$vi1) & is.finite(dat$vi1) & dat$vi1 > 0 &
      !is.na(dat$vi2) & is.finite(dat$vi2) & dat$vi2 > 0 &
      !is.na(dat$added_species_number) & is.finite(dat$added_species_number),
    , drop = FALSE
  ]
  
  if (nrow(dat_use) < 2) return(NULL)
  
  dat_use$S_ID <- factor(dat_use$S_ID)
  dat_use$ID   <- factor(dat_use$ID)
  
  fit1 <- tryCatch(
    lme(yi1 ~ 1 + added_species_number, dat_use, random = ~1 | S_ID/ID, weights = varFixed(~vi1), control = lmeControl(sigma = 1)),
    error = function(e) NULL
  )
  
  fit2 <- tryCatch(
    lme(yi2 ~ 1 + added_species_number, dat_use, random = ~1 | S_ID/ID, weights = varFixed(~vi2), control = lmeControl(sigma = 1)),
    error = function(e) NULL
  )
  
  if (is.null(fit1) || is.null(fit2)) return(NULL)
  
  dat.residual <- data.frame(r1 = residuals(fit1, type = "pearson"), r2 = residuals(fit2, type = "pearson"), S_ID = dat_use$S_ID, ID = dat_use$ID, stringsAsFactors = FALSE)
  
  fit_resid <- tryCatch(
    lme(r2 ~ -1 + r1, dat.residual, random = ~1 | S_ID/ID),
    error = function(e) NULL
  )
  
  if (is.null(fit_resid)) return(NULL)
  
  list(fit = fit_resid)
}

res <- NULL

## paired_names
for (i in seq_len(nrow(paired_names))) {
  dat <- read.table(paste0(prefix, paired_names[i, 1], ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  dat <- fix_id_sid(dat)
  out <- fit_main_safe(dat)
  if (is.null(out)) next
  
  r2_out <- safe_r2(out$fit)
  
  res <- rbind(
    res,
    data.frame(
      Item = paired_names[i, 1],
      AIC = AIC(out$fit),
      AICc = safe_AICc(out$fit),
      BIC = BIC(out$fit),
      logLik = safe_logLik(out$fit),
      R2_marginal = r2_out["R2_marginal"],
      R2_conditional = r2_out["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      stringsAsFactors = FALSE
    )
  )
}

## stratified
Strat <- c("Global", "Herbaceous_plant", "Woody_plant", "Food_crop", "Cash_crop", "Generalist", "Specialist", "Intercropping", "Cover_cropping", "Sown_field_margins", "Natural", "Plot", "Pot", "Temperate", "Tropic")

dat_all <- read.table(paste0(prefix, "herbivore_vs_crop.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dat_all <- fix_id_sid(dat_all)

for (i in Strat) {
  dat_sub <- dat_all
  
  if (i != "Global") {
    dat_sub <- subset(dat_sub, (CP_C3 == i) | (CP_Ct == i) | (Mn__ == i) | (Dv_s == i) | (Ex_4 == i) | (Ex_3 == i) | (Zone == i))
  }
  
  out <- fit_main_safe(dat_sub)
  if (is.null(out)) next
  
  r2_out <- safe_r2(out$fit)
  
  res <- rbind(
    res,
    data.frame(
      Item = i,
      AIC = AIC(out$fit),
      AICc = safe_AICc(out$fit),
      BIC = BIC(out$fit),
      logLik = safe_logLik(out$fit),
      R2_marginal = r2_out["R2_marginal"],
      R2_conditional = r2_out["R2_conditional"],
      R2_method = "MuMIn::r.squaredGLMM",
      stringsAsFactors = FALSE
    )
  )
}

View(res)
write.csv(res, '/Users/yusha/STable19b-2也许派.csv', row.names = FALSE)