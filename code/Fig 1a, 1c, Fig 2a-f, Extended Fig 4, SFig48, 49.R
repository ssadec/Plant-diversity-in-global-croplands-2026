####Fig1####
rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(rnaturalearth)
  library(stringi)
})

# =========================================================
##### Figure 1a####
# Purpose: Plot the global distribution of study locations.
# =========================================================

cfg <- list(
  input_file = "/Users/yusha/Fig1a.csv",
  output_file = "/Users/yusha/SFig1a_map.png",
  point_size = 3,
  point_alpha = 0.8,
  jitter_radius = 50000,
  seed = 1234,
  width = 14,
  height = 8,
  dpi = 600
)

cat_colors <- c("A" = "#BC3A24", "B" = "#EFB882", "C" = "#06798F", "D" = "#7A9A01", "E" = "#384D73")
cat_labels <- c("A" = "Invertebrate predators", "B" = "Invertebrate parasitoids", "C" = "Invertebrate predator-parasitoids", "D" = "Invertebrate herbivores", "E" = "Crops")

data_clean <- read.csv(cfg$input_file, stringsAsFactors = FALSE) %>%
  mutate(Category = trimws(Category), Category = gsub("[[:space:]]+", "_", Category), Category = stri_replace_all_regex(Category, "[\\p{C}]", ""), Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>%
  filter(!is.na(Category), Category != "", !grepl("^[[:space:]]*$", Category), !is.na(Latitude), !is.na(Longitude), Category %in% names(cat_colors)) %>%
  mutate(Category = factor(Category, levels = names(cat_colors)))

world_sf <- ne_countries(scale = "medium", returnclass = "sf")
set.seed(cfg$seed)

points_sf <- st_as_sf(data_clean, coords = c("Longitude", "Latitude"), crs = 4326)
robinson_crs <- "+proj=robin +datum=WGS84"
points_proj <- st_transform(points_sf, crs = robinson_crs)
world_proj <- st_transform(world_sf, crs = robinson_crs)

data_coords <- cbind(data_clean, st_coordinates(points_proj)) %>%
  group_by(X, Y) %>%
  mutate(n_dup = n(), dup_id = row_number()) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(angle = ifelse(n_dup > 1, 2 * pi * (dup_id - 1) / n_dup, 0), radius = ifelse(n_dup > 1, cfg$jitter_radius * sqrt(n_dup), 0), X_jit = X + radius * cos(angle), Y_jit = Y + radius * sin(angle)) %>%
  ungroup()

p1 <- ggplot() +
  geom_sf(data = world_proj, fill = "gainsboro", color = "white", linewidth = 0.5) +
  geom_point(data = data_coords, aes(x = X_jit, y = Y_jit, color = Category), size = cfg$point_size, alpha = cfg$point_alpha) +
  scale_color_manual(values = cat_colors, labels = cat_labels, name = "Trophic groups:") +
  coord_sf(crs = robinson_crs, xlim = st_bbox(world_proj)[c("xmin", "xmax")] * 1.02, ylim = st_bbox(world_proj)[c("ymin", "ymax")] * 1.02, expand = FALSE) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#87CEFA"),
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
    legend.background = element_rect(fill = "gainsboro", color = "white", linewidth = 0.6),
    legend.margin = margin(5, 5, 5, 5), legend.box.margin = margin(5, 5, 5, 5),
    legend.text = element_text(size = 22), legend.title = element_text(face = "bold", size = 22),
    legend.spacing.y = unit(0.5, "cm"), legend.key.height = unit(1, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4), keyheight = unit(1.5, "cm")))

ggsave(filename = cfg$output_file, plot = p1, width = cfg$width, height = cfg$height, dpi = cfg$dpi, bg = "white")


##### Figure 1c####
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gghalves)
  library(ggsci)
  library(stringr)
  library(scales)
  library(tibble)
})

# Define input and output paths in one place for easier reuse.
cfg_fig1c <- list(
  input_file = "/Users/yusha/Fig.1c.csv",
  output_file = "/Users/yusha/Fig1c.png",
  width = 14,
  height = 8,
  dpi = 600
)

# Read source data.
df <- read.csv(cfg_fig1c$input_file, stringsAsFactors = FALSE)

# Standardize grouping variables used in the figure.
df_clean <- df %>%
  mutate(
    Zone = str_to_lower(trimws(Zone)),
    Zone = case_when(
      Zone %in% c("tropic", "tropical") ~ "Tropical zone",
      Zone == "temperate" ~ "Temperate zone",
      TRUE ~ NA_character_
    ),
    CP_Ct = case_when(
      CP_Ct == "Food_crop" ~ "Food crop",
      CP_Ct == "Cash_crop" ~ "Cash crop",
      TRUE ~ NA_character_
    ),
    CP_C3 = case_when(
      CP_C3 == "Herbaceous_plant" ~ "Herbaceous crop",
      CP_C3 == "Woody_plant" ~ "Woody crop",
      TRUE ~ NA_character_
    )
  )

# Build the four grouping layers shown on the x-axis.
df_all  <- df_clean %>% mutate(Group = "All studies")
df_zone <- df_clean %>% filter(!is.na(Zone)) %>% transmute(S_ID, NCP_N, Group = Zone)
df_ct   <- df_clean %>% filter(!is.na(CP_Ct)) %>% transmute(S_ID, NCP_N, Group = CP_Ct)
df_c3   <- df_clean %>% filter(!is.na(CP_C3)) %>% transmute(S_ID, NCP_N, Group = CP_C3)

df_plot <- bind_rows(df_all %>% select(S_ID, NCP_N, Group), df_zone, df_ct, df_c3)

# Count studies rather than observations for the x-axis labels.
n_table <- bind_rows(
  tibble(Group = "All studies", N = n_distinct(df_all$S_ID)),
  df_zone %>% group_by(Group) %>% summarise(N = n_distinct(S_ID), .groups = "drop"),
  df_ct   %>% group_by(Group) %>% summarise(N = n_distinct(S_ID), .groups = "drop"),
  df_c3   %>% group_by(Group) %>% summarise(N = n_distinct(S_ID), .groups = "drop")
)

# Fix the x-axis order explicitly so the panel is stable across reruns.
group_order <- c(
  "All studies",
  "Herbaceous crop", "Woody crop",
  "Food crop", "Cash crop",
  "Temperate zone", "Tropical zone"
)

x_levels <- paste0(group_order, "\n(N=", n_table$N[match(group_order, n_table$Group)], ")")

df_plot <- df_plot %>%
  left_join(n_table, by = "Group") %>%
  mutate(x_label = factor(paste0(Group, "\n(N=", N, ")"), levels = x_levels))

# Keep one point per study within each x-axis group for the left-side jitter layer.
df_points <- df_plot %>% distinct(S_ID, x_label, NCP_N, Group)

# Summarize mean values for the yellow diamond and numeric label.
summary_df <- df_plot %>%
  group_by(Group, x_label) %>%
  summarise(mean_val = mean(NCP_N, na.rm = TRUE), .groups = "drop") %>%
  mutate(mean_lab = sprintf("%.2f", mean_val))

print(summary_df)

# Draw the figure with a right half-violin, boxplot, left jittered points, and mean labels.
p2 <- ggplot(df_plot, aes(x = x_label, y = NCP_N, fill = Group, color = Group)) +
  geom_half_violin(
    color = NA, side = "r", trim = FALSE, adjust = 1.2, width = 0.8,
    position = position_nudge(x = 0.12)
  ) +
  geom_boxplot(
    color = "black", width = 0.15, size = 0.6, outlier.size = 1,
    position = position_nudge(x = -0.05)
  ) +
  geom_point(
    data = df_points,
    aes(x = as.numeric(x_label) - 0.15, y = NCP_N, color = Group),
    inherit.aes = FALSE, size = 2, alpha = 0.9,
    position = position_jitter(width = 0.05)
  ) +
  geom_point(
    data = summary_df,
    aes(x = x_label, y = mean_val),
    inherit.aes = FALSE, shape = 23, size = 3.2,
    fill = "yellow", color = "black",
    position = position_nudge(x = -0.05)
  ) +
  geom_text(
    data = summary_df,
    aes(x = x_label, y = mean_val, label = mean_lab),
    inherit.aes = FALSE, family = "Arial", size = 6,
    fontface = "bold", color = "black", vjust = -0.9,
    position = position_nudge(x = -0.05)
  ) +
  scale_fill_jama() +
  scale_color_jama() +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = c(0, 2, 5, 10, 20, 40),
    limits = c(0, 40),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(x = NULL, y = "Number of added plant species") +
  theme_test(base_family = "Arial") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, family = "Arial", color = "black"),
    axis.text.x  = element_text(size = 24, hjust = 0.5, vjust = 0.5, family = "Arial", color = "black"),
    axis.text.y  = element_text(size = 24, family = "Arial", color = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 10, b = 50, l = 10),
    aspect.ratio = 1 / 3
  ) +
  coord_cartesian(clip = "off")

print(p2)

ggsave(
  filename = cfg_fig1c$output_file,
  plot = p2,
  width = cfg_fig1c$width,
  height = cfg_fig1c$height,
  dpi = cfg_fig1c$dpi,
  bg = "white"
)




#### SFig.49 ####
rm(list = ls())

cfg_49 <- list(input_file = "/Users/yusha/Fig.1c.csv", out_a = "/Users/yusha/SFig49a.png", out_b = "/Users/yusha/SFig49b.png", width = 14, height = 4.8, dpi = 600, y_max = 40)

group_colors_49 <- c(
  "Generalist herbivore" = "#1b9e77", "Specialist herbivore" = "#d95f02",
  "Generalist natural enemy" = "#7570b3", "Specialist natural enemy" = "#e7298a",
  "Intercropping" = "#66a61e", "Cover cropping" = "#e6ab02", "Sown field margins" = "#a6761d",
  "Natural experiment" = "#1f78b4", "Semi-natural experiment" = "#b2df8a", "Controlled experiment" = "#fb9a99",
  "Plot experiment" = "#6a3d9a", "Pot experiment" = "#ff7f00"
)

df <- read.csv(cfg_49$input_file, stringsAsFactors = FALSE)
if (!("NE_c" %in% names(df)) && ("NE_C" %in% names(df))) df$NE_c <- df$NE_C

df_clean <- df %>%
  mutate(
    S_ID = as.character(S_ID), NCP_N = as.numeric(NCP_N),
    Mn__ = case_when(trimws(as.character(Mn__)) %in% c("Generalist", "generalist") ~ "Generalist", trimws(as.character(Mn__)) %in% c("Specialist", "specialist") ~ "Specialist", TRUE ~ NA_character_),
    NE_c = case_when(trimws(as.character(NE_c)) %in% c("Generalist", "generalist") ~ "Generalist", trimws(as.character(NE_c)) %in% c("Specialist", "specialist") ~ "Specialist", TRUE ~ NA_character_),
    Ex_4 = case_when(trimws(as.character(Ex_4)) %in% c("Intercropping", "intercropping") ~ "Intercropping", trimws(as.character(Ex_4)) %in% c("Cover_cropping", "cover_cropping", "Cover cropping", "cover cropping") ~ "Cover cropping", trimws(as.character(Ex_4)) %in% c("Sown_field_margins", "sown_field_margins", "Sown field margins", "sown field margins") ~ "Sown field margins", TRUE ~ NA_character_),
    Ex_4 = case_when(tolower(trimws(as.character(Ex_4))) %in% c("natural") ~ "Natural experiment", tolower(trimws(as.character(Ex_4))) %in% c("semi_natural", "semi-natural", "seminatural") ~ "Semi-natural experiment", tolower(trimws(as.character(Ex_4))) %in% c("controlled") ~ "Controlled experiment", TRUE ~ NA_character_),
    Ex_3 = case_when(tolower(trimws(as.character(Ex_3))) == "plot" ~ "Plot experiment", tolower(trimws(as.character(Ex_3))) == "pot" ~ "Pot experiment", TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(NCP_N), is.finite(NCP_N), NCP_N >= 0)

make_combined_half_violin_plot <- function(df_long, group_levels, outfile, y_max = 40) {
  df_long <- df_long %>% filter(Group %in% group_levels) %>% mutate(Group = factor(Group, levels = group_levels))
  if (nrow(df_long) == 0) stop("No valid data after combining groups.")
  
  n_table <- df_long %>% group_by(Group) %>% summarise(N = n_distinct(S_ID), .groups = "drop")
  df_long <- df_long %>% left_join(n_table, by = "Group") %>% mutate(x_label = paste0(as.character(Group), "\n(N=", N, ")"))
  x_levels <- paste0(group_levels, "\n(N=", n_table$N[match(group_levels, n_table$Group)], ")")
  df_long$x_label <- factor(df_long$x_label, levels = x_levels)
  
  df_points <- df_long %>% distinct(S_ID, x_label, NCP_N, Group)
  summary_df <- df_long %>% group_by(Group, x_label) %>% summarise(mean_val = mean(NCP_N, na.rm = TRUE), .groups = "drop") %>% mutate(mean_lab = sprintf("%.2f", mean_val))
  df_violin <- df_long %>% group_by(x_label) %>% filter(n() >= 2, dplyr::n_distinct(NCP_N) >= 2) %>% ungroup()
  
  p_base <- ggplot(df_long, aes(x = x_label, y = NCP_N, fill = Group, color = Group))
  p <- tryCatch(
    p_base + geom_half_violin(data = df_violin, aes(group = x_label), color = NA, side = "r", trim = FALSE, adjust = 1.2, width = 0.8, position = position_nudge(x = 0.12), na.rm = TRUE),
    error = function(e) p_base + geom_violin(data = df_violin, aes(group = x_label), color = NA, trim = FALSE, adjust = 1.2, width = 0.8, na.rm = TRUE)
  )
  
  p <- p +
    geom_boxplot(color = "black", width = 0.15, size = 0.5, outlier.size = 0.8, position = position_nudge(x = -0.05), na.rm = TRUE) +
    geom_point(data = df_points, aes(x = as.numeric(x_label) - 0.15, y = NCP_N, color = Group), inherit.aes = FALSE, size = 1.8, alpha = 0.9, position = position_jitter(width = 0.05), na.rm = TRUE) +
    scale_fill_manual(values = group_colors_49, limits = names(group_colors_49), drop = FALSE) +
    scale_color_manual(values = group_colors_49, limits = names(group_colors_49), drop = FALSE) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10), breaks = c(0, 2, 5, 10, 20, 40), limits = c(0, y_max), expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Number of added plant species") +
    theme_test(base_family = "Arial") +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 18, family = "Arial", color = "black"), axis.text.x = element_text(size = 18, angle = 20, hjust = 1, vjust = 1, family = "Arial", color = "black"), axis.text.y = element_text(size = 14, family = "Arial", color = "black"), legend.position = "none", panel.background = element_rect(fill = "#FFFFFF"), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(t = 15, r = 10, b = 40, l = 10), aspect.ratio = 1 / 2.8) +
    coord_cartesian(clip = "off") +
    geom_point(data = summary_df, aes(x = x_label, y = mean_val), inherit.aes = FALSE, shape = 23, size = 2.8, fill = "yellow", color = "black", position = position_nudge(x = -0.05)) +
    geom_text(data = summary_df, aes(x = x_label, y = mean_val, label = mean_lab), inherit.aes = FALSE, family = "Arial", size = 3.4, fontface = "bold", color = "black", vjust = -0.9, position = position_nudge(x = -0.05))
  
  ggsave(filename = outfile, plot = p, width = cfg_49$width, height = cfg_49$height, dpi = cfg_49$dpi, bg = "white")
  p
}

df_49a_long <- bind_rows(
  df_clean %>% filter(!is.na(Mn__)) %>% mutate(Group = case_when(Mn__ == "Generalist" ~ "Generalist herbivore", Mn__ == "Specialist" ~ "Specialist herbivore", TRUE ~ NA_character_)) %>% filter(!is.na(Group)) %>% select(S_ID, NCP_N, Group),
  df_clean %>% filter(!is.na(NE_c)) %>% mutate(Group = case_when(NE_c == "Generalist" ~ "Generalist natural enemy", NE_c == "Specialist" ~ "Specialist natural enemy", TRUE ~ NA_character_)) %>% filter(!is.na(Group)) %>% select(S_ID, NCP_N, Group),
  df_clean %>% filter(!is.na(Ex_4)) %>% mutate(Group = Ex_4) %>% filter(!is.na(Group)) %>% select(S_ID, NCP_N, Group)
)

group_levels_49a <- c("Generalist herbivore", "Specialist herbivore", "Generalist natural enemy", "Specialist natural enemy", "Intercropping", "Cover cropping", "Sown field margins")
p49a <- make_combined_half_violin_plot(df_long = df_49a_long, group_levels = group_levels_49a, outfile = cfg_49$out_a, y_max = cfg_49$y_max)

df_49b_long <- bind_rows(
  df_clean %>% filter(!is.na(Ex_4)) %>% mutate(Group = Ex_4) %>% filter(!is.na(Group)) %>% select(S_ID, NCP_N, Group),
  df_clean %>% filter(!is.na(Ex_3)) %>% mutate(Group = Ex_3) %>% filter(!is.na(Group)) %>% select(S_ID, NCP_N, Group)
)

group_levels_49b <- c("Natural experiment", "Semi-natural experiment", "Controlled experiment", "Plot experiment", "Pot experiment")
p49b <- make_combined_half_violin_plot(df_long = df_49b_long, group_levels = group_levels_49b, outfile = cfg_49$out_b, y_max = cfg_49$y_max)


#### Fig2a-Fig2f ####
rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

cfg_fig2 <- list(
  base_dir = "/Users/yusha",
  colors = c("#BC3A24", "#EFB882", "#06798F", "#7A9A01", "#384D73"),
  family = "Arial"
)

read_fig2_data <- function(file) {
  dat <- read.csv(file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  dat$V2 <- as.numeric(dat$V2); dat$V3 <- as.numeric(dat$V3); dat$V4 <- as.numeric(dat$V4); dat$V5 <- as.factor(dat$V5)
  dat
}

make_fig2_plot <- function(dat, y_limits_vec, x_limits, x_breaks, ratio, point_size = 3.5, error_width = 0.3, vline_width = 2, output_file, unique_y = FALSE) {
  if (unique_y) {
    dat$V1_unique <- make.unique(as.character(dat$V1))
    y_order <- rev(dat$V1_unique)
    y_labels <- rev(dat$V1)
    p <- ggplot(dat, aes(x = V2, y = V1_unique, color = V5)) +
      geom_point(size = point_size) +
      geom_errorbar(aes(xmin = V3, xmax = V4), linewidth = 1, width = error_width) +
      scale_y_discrete(limits = y_order, labels = y_labels)
  } else {
    p <- ggplot(dat, aes(x = V2, y = as.factor(V1), color = V5)) +
      geom_point(size = point_size) +
      geom_errorbar(aes(xmin = V3, xmax = V4), linewidth = 1, width = error_width) +
      scale_y_discrete(limits = y_limits_vec)
  }
  
  p <- p +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    xlab("") + ylab("") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = vline_width) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = unit(c(1, 1, 1, 0), "cm"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 3),
      panel.spacing = unit(0.1, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(size = 48),
      axis.title.x = element_text(hjust = 0.5),
      text = element_text(color = "black", family = cfg_fig2$family, size = 24),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black", margin = margin(r = 1, unit = "pt"), lineheight = 50, size = 18),
      axis.ticks = element_line(color = "black")
    ) +
    scale_color_manual(values = cfg_fig2$colors) +
    coord_fixed(ratio = ratio)
  
  ggsave(output_file, p, width = 18, height = 18, units = "cm", dpi = 300)
  p
}

y_fig2a <- c(
  "Crop quality (26/4)", "Crop reproduction (271/75)", "Crop growth (95/35)", "",
  "Herbivore damage (239/48)", "Herbivore reproduction (445/109)", "Herbivore growth (1/1)", "",
  "Predator-parasitoid diversity (6/3)", "Predator-parasitoid reproduction (13/5)", "",
  "Parasitoid parasitism (36/19)", "Parasitoid diversity (14/2)", "Parasitoid reproduction (52/21)", "Parasitoid growth (3/1)", "",
  "Predator predation (39/12)", "Predator diversity (62/16)", "Predator reproduction (221/66)", "",
  "Total crop response (392/105)", "Total herbivore response (686/132)", "Total predator-parasitoid response (19/7)", "Total parasitoid response (105/33)", "Total predator response (322/75)"
)
fig2a <- make_fig2_plot(dat = read_fig2_data(file.path(cfg_fig2$base_dir, "Fig.2a.csv")), y_limits_vec = y_fig2a, x_limits = c(-2.5, 2.5), x_breaks = seq(-2, 2, by = 1), ratio = 0.265, point_size = 3, error_width = 0.4, output_file = file.path(cfg_fig2$base_dir, "fig2a.png"))

y_fig2b <- c(
  "Total crop response (82/26)", "Total herbivore response (171/137)", "Total predator-parasitoid response (3/2)", "Total parasitoid response (26/5)", "Total predator response (105/20)", "",
  "Total crop response (310/81)", "Total herbivore response (515/97)", "Total predator-parasitoid response (11/5)", "Total parasitoid response (73/26)", "Total predator response (215/55)", "",
  "Total crop response (25/6)", "Total herbivore response (107/24)", "Total predator-parasitoid response (7/2)", "Total parasitoid response (26/7)", "Total predator response (108/17)", "",
  "Total crop response (367/99)", "Total herbivore response (579/108)", "Total predator-parasitoid response (12/5)", "Total parasitoid response (79/26)", "Total predator response (214/58)"
)
fig2b <- make_fig2_plot(dat = read_fig2_data(file.path(cfg_fig2$base_dir, "Fig.2b.csv")), y_limits_vec = y_fig2b, x_limits = c(-4, 4), x_breaks = seq(-4, 4, by = 2), ratio = 0.475, point_size = 3.5, error_width = 0.3, output_file = file.path(cfg_fig2$base_dir, "fig2b.crop.png"))

fig2c <- make_fig2_plot(dat = read_fig2_data(file.path(cfg_fig2$base_dir, "Fig.2c.csv")), y_limits_vec = NULL, x_limits = c(-2.5, 2.5), x_breaks = seq(-2, 2, by = 1), ratio = 0.295, point_size = 4, error_width = 0.4, vline_width = 1, output_file = file.path(cfg_fig2$base_dir, "fig2c.png"), unique_y = TRUE)

y_fig2d <- c(
  "Total crop response (86/22)", "Total herbivore response (84/24)", "Total predator-parasitoid response (10/3)", "Total parasitoid response (32/11)", "Total predator response (76/17)", "",
  "Total crop response (76/16)", "Total herbivore response (146/30)", "Total predator-parasitoid response (0/0)", "Total parasitoid response (14/1)", "Total predator response (120/21)", "",
  "Total crop response (230/69)", "Total herbivore response (456/79)", "Total predator-parasitoid response (9/4)", "Total parasitoid response (59/22)", "Total predator response (126/39)"
)
fig2d <- make_fig2_plot(dat = read_fig2_data(file.path(cfg_fig2$base_dir, "Fig.2d.csv")), y_limits_vec = y_fig2d, x_limits = c(-2.5, 2.5), x_breaks = seq(-2, 2, by = 1), ratio = 0.4, point_size = 3.5, error_width = 0.3, output_file = file.path(cfg_fig2$base_dir, "fig2d.png"))

y_fig2e <- c(
  "Total crop response (0/0)", "Total herbivore response (53/5)", "Total predator-parasitoid response (0/0)", "Total parasitoid response (12/2)", "Total predator response (15/3)", "",
  "Total crop response (1/1)", "Total herbivore response (14/4)", "Total predator-parasitoid response (0/0)", "Total parasitoid response (0/0)", "Total predator response (1/1)", "",
  "Total crop response (391/104)", "Total herbivore response (617/124)", "Total predator-parasitoid response (19/7)", "Total parasitoid response (93/32)", "Total predator response (306/71)"
)
fig2e <- make_fig2_plot(dat = read_fig2_data(file.path(cfg_fig2$base_dir, "Fig.2e.csv")), y_limits_vec = y_fig2e, x_limits = c(-2.5, 2.5), x_breaks = seq(-2, 2, by = 1), ratio = 0.415, point_size = 3.5, error_width = 0.2, output_file = file.path(cfg_fig2$base_dir, "fig2e.png"))

y_fig2f <- c(
  "Total crop response (149/35)", "Total herbivore response (256/35)", "Total predator-parasitoid response (0/0)", "Total parasitoid response (15/6)", "Total predator response (15/7)", "",
  "Total crop response (241/68)", "Total herbivore response (355/89)", "Total predator-parasitoid response (19/7)", "Total parasitoid response (78/26)", "Total predator response (290/63)"
)
fig2f <- make_fig2_plot(dat = read_fig2_data(file.path(cfg_fig2$base_dir, "Fig.2f.csv")), y_limits_vec = y_fig2f, x_limits = c(-2.5, 2.5), x_breaks = seq(-2, 2, by = 1), ratio = 0.61, point_size = 3.5, error_width = 0.2, output_file = file.path(cfg_fig2$base_dir, "fig2f.png"))



#### Extended data Fig. 4 ####
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(gghalves)
  library(ggsci)
})

cfg_ed4 <- list(input_file = "Alldata.txt", y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12)

data <- read.table(cfg_ed4$input_file, header = TRUE, stringsAsFactors = FALSE)

prepare_plot_data <- function(dat) {
  dat %>%
    mutate(
      rel.per = ((TREATMENT_M - CONTROL_M) / CONTROL_M) * 100,
      Zone = case_when(str_to_lower(trimws(Zone)) %in% c("tropic", "tropical") ~ "Tropical zone", str_to_lower(trimws(Zone)) == "temperate" ~ "Temperate zone", TRUE ~ NA_character_),
      CP_Ct = case_when(CP_Ct == "Food_crop" ~ "Food crop", CP_Ct == "Cash_crop" ~ "Cash crop", TRUE ~ NA_character_),
      CP_C3 = case_when(CP_C3 == "Herbaceous_plant" ~ "Herbaceous crop", CP_C3 == "Woody_plant" ~ "Woody crop", TRUE ~ NA_character_)
    )
}

build_df_plot <- function(dat) {
  df_all <- dat %>% mutate(Group = "All studies")
  df_zone <- dat %>% filter(!is.na(Zone)) %>% mutate(Group = Zone)
  df_ct <- dat %>% filter(!is.na(CP_Ct)) %>% mutate(Group = CP_Ct)
  df_c3 <- dat %>% filter(!is.na(CP_C3)) %>% mutate(Group = CP_C3)
  df_plot <- bind_rows(df_all, df_zone, df_ct, df_c3)
  
  n_table <- bind_rows(
    tibble(Group = "All studies", N = nrow(df_all)),
    df_c3 %>% group_by(Group) %>% summarise(N = n(), .groups = "drop"),
    df_ct %>% group_by(Group) %>% summarise(N = n(), .groups = "drop"),
    df_zone %>% group_by(Group) %>% summarise(N = n(), .groups = "drop")
  ) %>% distinct(Group, .keep_all = TRUE)
  
  group_order <- c("All studies", "Herbaceous crop", "Woody crop", "Food crop", "Cash crop", "Temperate zone", "Tropical zone")
  level_labels <- paste0(group_order, "\n(N=", n_table$N[match(group_order, n_table$Group)], ")")
  
  df_plot %>% left_join(n_table, by = "Group") %>% mutate(x_label = factor(paste0(Group, "\n(N=", N, ")"), levels = level_labels))
}

build_summary_values <- function(df_plot) df_plot %>% group_by(x_label) %>% summarise(mean_rel_per = round(mean(rel.per, na.rm = TRUE), 2), median_rel_per = round(median(rel.per, na.rm = TRUE), 2), lower_25 = round(quantile(rel.per, 0.25, na.rm = TRUE), 2), upper_75 = round(quantile(rel.per, 0.75, na.rm = TRUE), 2), .groups = "drop")

build_mean_label_positions <- function(df_plot, text_offset = 12) df_plot %>% group_by(x_label) %>% summarise(mean_plot = mean(rel.per, na.rm = TRUE), upper_75 = quantile(rel.per, 0.75, na.rm = TRUE), label_y = upper_75 + text_offset, .groups = "drop")

make_extended_fig4_plot <- function(df_plot, y_title, y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12) {
  label_df <- build_mean_label_positions(df_plot, text_offset)
  
  ggplot(df_plot, aes(x = x_label, y = rel.per, fill = Group, color = Group)) +
    geom_half_violin(color = NA, side = "r", trim = FALSE, adjust = 1.1, width = 0.5, alpha = 0.4, position = position_nudge(x = 0.1)) +
    geom_boxplot(color = "black", width = 0.15, size = 0.6, outlier.size = 1, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "yellow", color = "black") +
    geom_point(data = df_plot %>% distinct(S_ID, x_label, rel.per, Group), aes(x = as.numeric(x_label) - 0.2, y = rel.per, fill = Group, color = Group), inherit.aes = FALSE, size = 2, alpha = 0.6, position = position_jitter(width = 0.08)) +
    scale_fill_jama() +
    scale_color_jama() +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = y_limits, breaks = y_breaks) +
    labs(y = y_title) +
    theme_test(base_family = "Arial") +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 24, family = "Arial", color = "black"), axis.text.x = element_text(size = 24, angle = 0, hjust = 0.5, family = "Arial", color = "black"), axis.text.y = element_text(size = 24, family = "Arial", color = "black"), legend.position = "none", panel.background = element_rect(fill = "#FFFFFF"), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(t = 25, r = 10, b = 50, l = 10)) +
    coord_cartesian(clip = "off") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(data = label_df, aes(x = x_label, y = label_y, label = sprintf("%.2f", mean_plot)), inherit.aes = FALSE, color = "black", size = 6, fontface = "bold")
}

run_extended_fig4 <- function(data, subset_expr, y_title, png_file, csv_file, y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12) {
  dat <- data %>% filter({{ subset_expr }}) %>% prepare_plot_data()
  df_plot <- build_df_plot(dat)
  summary_values <- build_summary_values(df_plot)
  p <- make_extended_fig4_plot(df_plot = df_plot, y_title = y_title, y_limits = y_limits, y_breaks = y_breaks, text_offset = text_offset)
  ggsave(png_file, p, width = 12, height = 5.5, dpi = 300)
  write.csv(summary_values, csv_file, row.names = FALSE)
  list(data = dat, df_plot = df_plot, summary_values = summary_values, plot = p)
}

res_crop <- run_extended_fig4(data = data, subset_expr = Trophic_group == "Crop", y_title = "Increased percentage of crop performance", png_file = "ExtendedFig.4_CropPerformance_gghalves_final.png", csv_file = "ExtendedFig.4_CropPerformance_summary_values_2decimals.csv", y_limits = cfg_ed4$y_limits, y_breaks = cfg_ed4$y_breaks, text_offset = cfg_ed4$text_offset)
res_enemy <- run_extended_fig4(data = data, subset_expr = Trophic_group %in% c("Predator", "Parasitoid", "Predator_parasitoid"), y_title = "Increased percentage of natural enemy performance", png_file = "ExtendedFig.4_EnemyPerformance_gghalves_final.png", csv_file = "ExtendedFig.4_EnemyPerformance_summary_values_2decimals.csv", y_limits = cfg_ed4$y_limits, y_breaks = cfg_ed4$y_breaks, text_offset = cfg_ed4$text_offset)
res_herbivore <- run_extended_fig4(data = data, subset_expr = Trophic_group == "Herbivore", y_title = "Decreased percentage of herbivore performance", png_file = "ExtendedFig.4_HerbivorePerformance_gghalves_final.png", csv_file = "ExtendedFig.4_HerbivorePerformance_summary_values_2decimals.csv", y_limits = cfg_ed4$y_limits, y_breaks = cfg_ed4$y_breaks, text_offset = cfg_ed4$text_offset)


#### SFig 48 ####
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(gghalves)
})

data <- read.table("Alldata.txt", header = TRUE, stringsAsFactors = FALSE)

group_order <- c(
  "Generalist herbivore", "Specialist herbivore",
  "Generalist natural enemy", "Specialist natural enemy",
  "Intercropping", "Cover cropping", "Sown field margins",
  "Natural experiment", "Semi-natural experiment", "Controlled experiment",
  "Plot experiment", "Pot experiment"
)

group_colors <- c(
  "Generalist herbivore" = "#1b9e77", "Specialist herbivore" = "#d95f02",
  "Generalist natural enemy" = "#7570b3", "Specialist natural enemy" = "#e7298a",
  "Intercropping" = "#66a61e", "Cover cropping" = "#e6ab02", "Sown field margins" = "#a6761d",
  "Natural experiment" = "#1f78b4", "Semi-natural experiment" = "#b2df8a", "Controlled experiment" = "#fb9a99",
  "Plot experiment" = "#6a3d9a", "Pot experiment" = "#ff7f00"
)

prepare_plot_data <- function(dat) {
  dat %>%
    mutate(
      rel.per = ((TREATMENT_M - CONTROL_M) / CONTROL_M) * 100,
      Mn__ = case_when(str_trim(Mn__) == "Generalist" ~ "Generalist herbivore", str_trim(Mn__) == "Specialist" ~ "Specialist herbivore", TRUE ~ NA_character_),
      NE_c = case_when(str_trim(NE_c) == "Generalist" ~ "Generalist natural enemy", str_trim(NE_c) == "Specialist" ~ "Specialist natural enemy", TRUE ~ NA_character_),
      Ex_4 = case_when(str_trim(Ex_4) == "Intercropping" ~ "Intercropping", str_trim(Ex_4) == "Cover_cropping" ~ "Cover cropping", str_trim(Ex_4) == "Sown_field_margins" ~ "Sown field margins", TRUE ~ NA_character_),
      Ex_4 = case_when(str_to_lower(str_trim(Ex_4)) == "natural" ~ "Natural experiment", str_to_lower(str_trim(Ex_4)) %in% c("semi_natural", "semi-natural", "seminatural") ~ "Semi-natural experiment", str_to_lower(str_trim(Ex_4)) == "controlled" ~ "Controlled experiment", TRUE ~ NA_character_),
      Ex_3 = case_when(str_to_lower(str_trim(Ex_3)) == "plot" ~ "Plot experiment", str_to_lower(str_trim(Ex_3)) == "pot" ~ "Pot experiment", TRUE ~ NA_character_)
    )
}

build_df_plot_12groups <- function(dat) {
  df_plot <- bind_rows(
    dat %>% filter(!is.na(Mn__)) %>% mutate(Group = Mn__),
    dat %>% filter(!is.na(NE_c)) %>% mutate(Group = NE_c),
    dat %>% filter(!is.na(Ex_4)) %>% mutate(Group = Ex_4),
    dat %>% filter(!is.na(Ex_4)) %>% mutate(Group = Ex_4),
    dat %>% filter(!is.na(Ex_3)) %>% mutate(Group = Ex_3)
  ) %>% filter(!is.na(rel.per), is.finite(rel.per))
  
  n_table <- df_plot %>% group_by(Group) %>% summarise(N = n(), .groups = "drop") %>% filter(Group %in% group_order)
  present_groups <- group_order[group_order %in% n_table$Group]
  
  label_table <- tibble(Group = present_groups, x_id = seq_along(present_groups)) %>% left_join(n_table, by = "Group") %>% mutate(x_label = paste0(Group, "\n(N=", N, ")"))
  
  list(
    df_plot = df_plot %>% filter(Group %in% present_groups) %>% left_join(label_table, by = "Group") %>% mutate(Group = factor(Group, levels = present_groups)),
    label_table = label_table
  )
}

build_summary_values <- function(df_plot, label_table) df_plot %>% group_by(Group, x_id) %>% summarise(mean_rel_per = round(mean(rel.per, na.rm = TRUE), 2), median_rel_per = round(median(rel.per, na.rm = TRUE), 2), lower_25 = round(quantile(rel.per, 0.25, na.rm = TRUE), 2), upper_75 = round(quantile(rel.per, 0.75, na.rm = TRUE), 2), .groups = "drop") %>% left_join(label_table, by = c("Group", "x_id"))

build_mean_label_positions <- function(df_plot, text_offset = 12) df_plot %>% group_by(Group, x_id) %>% summarise(mean_plot = mean(rel.per, na.rm = TRUE), upper_75 = quantile(rel.per, 0.75, na.rm = TRUE), label_y = upper_75 + text_offset, .groups = "drop")

make_sfig48_plot <- function(df_plot, label_table, y_title, y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12) {
  df_violin <- df_plot %>% group_by(x_id) %>% filter(n() >= 2, dplyr::n_distinct(rel.per) >= 2) %>% ungroup()
  label_df <- build_mean_label_positions(df_plot, text_offset)
  p_base <- ggplot(df_plot, aes(x = x_id, y = rel.per, fill = Group, color = Group))
  
  p <- tryCatch(
    p_base + geom_half_violin(data = df_violin, aes(group = x_id), color = NA, side = "r", trim = FALSE, adjust = 1.1, width = 0.5, alpha = 0.4, position = position_nudge(x = 0.1), na.rm = TRUE),
    error = function(e) p_base + geom_violin(data = df_violin, aes(group = x_id), color = NA, trim = FALSE, adjust = 1.1, width = 0.5, alpha = 0.4, na.rm = TRUE)
  )
  
  p +
    geom_boxplot(aes(group = x_id), color = "black", width = 0.15, size = 0.7, outlier.size = 1.2, alpha = 0.5, na.rm = TRUE) +
    stat_summary(aes(group = x_id), fun = mean, geom = "point", shape = 23, size = 3.5, fill = "yellow", color = "black", na.rm = TRUE) +
    geom_point(data = df_plot %>% distinct(S_ID, Group, x_id, rel.per), aes(x = x_id - 0.2, y = rel.per, fill = Group, color = Group), inherit.aes = FALSE, size = 2.4, alpha = 0.65, position = position_jitter(width = 0.08, height = 0), na.rm = TRUE) +
    scale_fill_manual(values = group_colors, limits = group_order, drop = FALSE) +
    scale_color_manual(values = group_colors, limits = group_order, drop = FALSE) +
    scale_x_continuous(breaks = label_table$x_id, labels = label_table$x_label) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = y_limits, breaks = y_breaks) +
    labs(y = y_title, x = NULL) +
    theme_test(base_family = "Arial") +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 22, family = "Arial", color = "black"), axis.text.x = element_text(size = 17, angle = 30, hjust = 1, vjust = 1, family = "Arial", color = "black", lineheight = 0.95), axis.text.y = element_text(size = 24, family = "Arial", color = "black"), legend.position = "none", panel.background = element_rect(fill = "#FFFFFF"), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(t = 30, r = 15, b = 70, l = 15)) +
    coord_cartesian(clip = "off") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(data = label_df, aes(x = x_id, y = label_y, label = sprintf("%.2f", mean_plot)), inherit.aes = FALSE, color = "black", size = 6.2, fontface = "bold")
}

run_sfig48 <- function(data, subset_expr, y_title, png_file, csv_file, y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12) {
  dat <- data %>% filter({{ subset_expr }}) %>% prepare_plot_data()
  built <- build_df_plot_12groups(dat)
  df_plot <- built$df_plot
  label_table <- built$label_table
  if (nrow(df_plot) == 0) stop("No valid rows available for plotting after cleaning.")
  summary_values <- build_summary_values(df_plot, label_table)
  p <- make_sfig48_plot(df_plot = df_plot, label_table = label_table, y_title = y_title, y_limits = y_limits, y_breaks = y_breaks, text_offset = text_offset)
  ggsave(filename = png_file, plot = p, width = 22, height = 9.5, dpi = 300)
  write.csv(summary_values, csv_file, row.names = FALSE)
  list(data = dat, df_plot = df_plot, label_table = label_table, summary_values = summary_values, plot = p)
}

res_enemy <- run_sfig48(data = data, subset_expr = Trophic_group %in% c("Predator", "Parasitoid", "Predator_parasitoid"), y_title = "Relative percentage change in natural enemy performance", png_file = "SFig48_NaturalEnemy.png", csv_file = "SFig48_NaturalEnemy_summary.csv", y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12)
res_herbivore <- run_sfig48(data = data, subset_expr = Trophic_group == "Herbivore", y_title = "Relative percentage change in herbivore performance", png_file = "SFig48_Herbivore.png", csv_file = "SFig48_Herbivore_summary.csv", y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12)
res_crop <- run_sfig48(data = data, subset_expr = Trophic_group == "Crop", y_title = "Relative percentage change in crop performance", png_file = "SFig48_Crop.png", csv_file = "SFig48_Crop_summary.csv", y_limits = c(-100, 500), y_breaks = seq(-100, 500, 100), text_offset = 12)