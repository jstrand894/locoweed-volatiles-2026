# ============================================================
# fig3_simper.R
# Figure 3: SIMPER compound contributions to VOC dissimilarity
# between endophyte-absent and endophyte-present plants.
#
# Two panels: GH (combined) and Field living tissue (combined).
# Bars ranked by average contribution to Bray-Curtis dissimilarity.
# Fill color indicates direction: which group has higher abundance.
#
# INPUTS (from 04_simper.R):
#   simper_gh, simper_live
# ============================================================

source(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "04_simper.R"))

summary(simper_gh)[[1]] %>% rownames()
# ============================================================
# PREPARE DATA
# ============================================================

name_lookup <- c(
  "Caryophyllene"      = "\u03b2-Caryophyllene",
  "Bpinene"            = "\u03b2-Pinene",
  "Bocimene"           = "\u03b2-Ocimene",
  "Dlimonene"          = "(R)-Limonene",
  "Methylsalicylate"   = "Methyl salicylate",
  "Hexenylisovalerate" = "(Z)-3-Hexenyl isovalerate",
  "Cisthreehexiso"     = "(Z)-3-Hexenyl isobutyrate",
  "Z3hex"              = "(Z)-3-Hexenyl acetate",
  "Sixmethyl"          = "6-Methyl-5-hepten-2-one",
  "Octanal"            = "Octanal",
  "Nonanal"            = "Nonanal",
  "Decanal"            = "Decanal",
  "Heptanal"           = "Heptanal",
  "Manoyl oxide"       = "Manoyl oxide"
)


# Pull all compounds from each simper object, add direction and dataset label.
prep_simper_plot <- function(sim, label, n = 10) {
  summary(sim)[[1]] |>
    rownames_to_column("compound") |>
    arrange(desc(average)) |>
    slice_head(n = n) |>
    mutate(
      dataset   = label,
      direction = if_else(ava > avb, "Higher in Absent", "Higher in Present"),
      compound  = str_to_sentence(compound),
      compound  = recode(compound, !!!name_lookup),
      compound  = fct_reorder(compound, average)
    )
}

simper_plot_data <- bind_rows(
  prep_simper_plot(simper_gh,   "Greenhouse"),
  prep_simper_plot(simper_live, "Field (living tissue)")
) %>%
  mutate(dataset = factor(dataset, levels = c("Greenhouse", "Field (living tissue)")))


# ============================================================
# PLOT
# ============================================================

simper_plot_data <- bind_rows(gh_data, live_data) %>%
  mutate(dataset = factor(dataset, levels = c("Greenhouse", "Field (living tissue)")))

# Build decoupled levels per panel
gh_levels   <- gh_data   %>% arrange(average) %>% mutate(cp = paste("Greenhouse", compound, sep = "__")) %>% pull(cp)
live_levels <- live_data %>% arrange(average) %>% mutate(cp = paste("Field (living tissue)", compound, sep = "__")) %>% pull(cp)

simper_plot_data <- simper_plot_data %>%
  mutate(
    compound_panel = paste(dataset, compound, sep = "__"),
    compound_panel = factor(compound_panel, levels = c(live_levels, gh_levels))
  )

SIMPER_Plot <-
  ggplot(simper_plot_data,
         aes(x = average, y = compound_panel, fill = direction)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.7) +
  geom_errorbar(aes(xmin = average - sd, xmax = average + sd),
                orientation = "y",
                width = 0.25, linewidth = 0.4) +
  scale_fill_manual(
    values = c("Higher in Absent"  = "#D55E00",
               "Higher in Present" = "#009E73")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_discrete(labels = \(x) sub(".*__", "", x)) +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(
    x    = "Average contribution to dissimilarity (Bray-Curtis)",
    y    = NULL,
    fill = "Direction"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = NA),
    legend.background  = element_rect(color = "black", fill = "white", linewidth = 0.3),
    axis.text.y        = element_text(size = 9)
  )

SIMPER_Plot

ggsave("figures/Fig3_SIMPER.png", SIMPER_Plot,
       width = 9, height = 5, dpi = 300)


