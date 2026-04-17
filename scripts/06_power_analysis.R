# ============================================================
# 05_power_analysis.R
# Simulation-based power analysis for paired PERMANOVA design.
#
# Approach: within-pair signal injection.
# For each simulation, a synthetic CLR-space shift is added to
# the E+ member of each pair only, scaled as multiples of the
# empirical within-pair SD. This matches what the marginal
# endophyte term in adonis2 is actually estimating -- the E+/E-
# difference after pair identity is accounted for.
#
# INPUTS (from 01_clr_processing.R):
#   gh_clr, gh_meta
#
# OUTPUTS:
#   power_results  - tibble of effect size vs. power
#   type1_error    - empirical type I error rate at effect = 0
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "01_clr_processing.R"))

library(tidyverse)
library(vegan)

set.seed(4721)

# ============================================================
# PARAMETERS
# ============================================================

n_sim        <- 200
effect_sizes <- seq(0, 0.3, by = 0.01)

# ============================================================
# PREP: align CLR matrix with metadata
# ============================================================

clr_mat <- as.matrix(gh_clr)   # rows = samples, cols = compounds

# Confirm row order matches gh_meta
stopifnot(identical(rownames(clr_mat), rownames(gh_meta)))

# Identify E+ and E- indices within each pair
# Each pair should have exactly one E+ and one E- member
pair_ids   <- unique(gh_meta$pair)
ep_idx     <- which(gh_meta$swa == "1")   # E+ row indices
em_idx     <- which(gh_meta$swa == "0")   # E- row indices


# Identify complete pairs (exactly one E+ and one E-)
complete_pairs <- gh_meta %>%
  group_by(pair) %>%
  summarise(
    n_ep = sum(swa == "1"),
    n_em = sum(swa == "0"),
    .groups = "drop"
  ) %>%
  filter(n_ep == 1, n_em == 1) %>%
  pull(pair)

cat("Complete pairs:", length(complete_pairs), "\n")
cat("Dropped pairs:", sum(!unique(gh_meta$pair) %in% complete_pairs), "\n")

# Subset metadata and CLR matrix to complete pairs only
keep_rows <- gh_meta$pair %in% complete_pairs
gh_meta_pw <- gh_meta[keep_rows, ]
clr_mat_pw <- clr_mat[keep_rows, ]

# Rebuild pair lookup on complete pairs only
pair_lookup <- gh_meta_pw %>%
  mutate(row_idx = row_number()) %>%
  group_by(pair) %>%
  summarise(
    ep_row = row_idx[swa == "1"],
    em_row = row_idx[swa == "0"],
    .groups = "drop"
  )

# Sanity check: every pair has exactly one E+ and one E-
stopifnot(nrow(pair_lookup) == length(complete_pairs))

# ============================================================
# EMPIRICAL WITHIN-PAIR SD
# ============================================================
# For each compound, compute the SD of (E+ - E-) differences
# across pairs. This is the scale unit for signal injection.

within_pair_diff <- clr_mat_pw[pair_lookup$ep_row, ] -
  clr_mat_pw[pair_lookup$em_row, ]

compound_sd <- apply(within_pair_diff, 2, sd)

cat("Mean within-pair SD across compounds:",
    round(mean(compound_sd), 4), "\n")

# ============================================================
# SIMULATION FUNCTION
# ============================================================
# Injects a uniform shift (effect * compound_sd) into all E+
# rows, then runs the same adonis2 model as 02_permanova.R.
# Returns the p-value for the swa term.

run_sim <- function(effect_size) {
  
  sim_mat <- clr_mat_pw
  
  signal <- matrix(
    rep(effect_size * compound_sd, each = nrow(pair_lookup)),
    nrow = nrow(pair_lookup),
    ncol = ncol(clr_mat_pw)
  )
  sim_mat[pair_lookup$ep_row, ] <- sim_mat[pair_lookup$ep_row, ] + signal
  
  sim_dist <- vegdist(sim_mat, method = "euclidean")
  
  fit <- adonis2(
    as.dist(sim_dist) ~ swa + pair,
    data         = gh_meta_pw,
    permutations = 499,
    strata       = gh_meta_pw$year,
    by           = "margin"
  )
  
  as.data.frame(fit)["swa", "Pr(>F)"]
}

# ============================================================
# RUN SIMULATIONS
# ============================================================

power_results <- map_dfr(effect_sizes, function(es) {
  
  cat("Running effect size:", es, "\n")
  
  p_vals <- replicate(n_sim, run_sim(es))
  
  tibble(
    effect_size = es,
    power       = mean(p_vals < 0.05, na.rm = TRUE),
    n_na        = sum(is.na(p_vals))
  )
})

# Type I error rate is power at effect_size == 0
type1_error <- power_results %>%
  filter(effect_size == 0) %>%
  pull(power)

cat("\nEmpirical type I error rate (effect = 0):",
    round(type1_error, 3), "\n")
cat("Nominal alpha: 0.05\n")

print(power_results, n = 100)

# ============================================================
# PLOT
# ============================================================

ggplot(power_results, aes(x = effect_size, y = power)) +
  geom_hline(yintercept = 0.05,  linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 0.80,  linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = type1_error, linetype = "dotted",
             color = "firebrick", linewidth = 0.5) +
  geom_line(linewidth = 1, color = "steelblue4") +
  geom_point(size = 2,   color = "steelblue4") +
  annotate("text", x = max(effect_sizes) * 0.95, y = 0.08,
           label = "α = 0.05", hjust = 1, size = 3.2,
           color = "gray40") +
  annotate("text", x = max(effect_sizes) * 0.95, y = 0.84,
           label = "Power = 0.80", hjust = 1, size = 3.2,
           color = "gray40") +
  annotate("text", x = max(effect_sizes) * 0.95,
           y = type1_error + 0.04,
           label = paste0("Empirical type I error = ",
                          round(type1_error, 2)),
           hjust = 1, size = 3.2, color = "firebrick") +
  scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = "Effect size (multiples of within-pair SD)",
    y     = "Power (proportion p < 0.05)",
    title = "Simulation-based power analysis",
    subtitle = paste0("Paired PERMANOVA, within-pair signal injection\n",
                      "n = ", n_sim, " simulations per effect size")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3)
  )



n_sim        <- 200
effect_sizes <- seq(0.13, 0.17, by = 0.005)